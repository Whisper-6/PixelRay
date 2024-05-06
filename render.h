#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <algorithm>
#include <vector>
#include "threadpool.h"
#include "geometry.h"
#include "graphics.h"

struct Voxel {
	Point pos; // 空间位置
	Vector normal; // 法向量
	Color col; // 颜色
};

unsigned timeStamp;

// 三角面，记录三个顶点的位置、法向、颜色，中间的点的数据由线性插值而得.
struct Triangle {
	Voxel v1, v2, v3;
	unsigned stamp;
};

// 管理三角面的数据结构（目前是 std::vector）
class ChunkSet {
private:
public:
	std::vector<Triangle> tris;
	void appendTriangle(Triangle& tri) {
		tris.push_back(tri);
	}
};

#undef min
#undef max

using std::min;
using std::max;

// 深度图信息，记录每个 Pixel 被哪个 Triangle 遮挡，以及遮挡处的 depth
struct DepthData {
	Triangle* tri;
	float depth;
	void clear() {
		tri = nullptr;
		depth = 1e9;
	}
};

// 正交相机
struct Camera {
	Vector Veci, Vecj, Vecd;
	float deti, detj;
	// 设定相机朝向. alpha 为地面角, beta 为仰角
	void set(float alpha, float beta, Point Lpos) {
		Veci = make_vec(-sin(beta) * sin(alpha), sin(beta) * cos(alpha), -cos(beta));
		Vecj = make_vec(cos(alpha), sin(alpha), 0);
		Vecd = vec3::cross(Veci, Vecj);
		deti = -(Lpos * Veci);
		detj = -(Lpos * Vecj);
	}
	// 令取景框在屏幕上移动 (_deti, _detj)
	void move(float _deti, float _detj) {
		deti -= _deti;
		detj -= _detj;
	}
	// 得到点在面板上的 i 坐标
	inline float i(Point pos) const {
		return pos * Veci + deti;
	}
	// 得到点在面板上的 j 坐标
	inline float j(Point pos) const {
		return pos * Vecj + detj;
	}
	// 得到点的 depth
	inline float d(Point pos) const {
		return pos * Vecd;
	}
	// 还原点的世界坐标
	inline Point pos(float i, float j, float d) const {
		return (i - deti) * Veci + (j - detj) * Vecj + d * Vecd;
	}
};

template<int N, int M>
class DepthMap {
public:
	DepthData map[N][M];
	Camera camera;
	void project(ChunkSet& chunks);
	inline bool shadowTest(Point& pos);
	float accurateShadowTest(Point& pos);
	void getScreen(Color col[N][M], Point cameraPosMap[N][M], Vector normalMap[N][M]);
};

inline float intersection(float i1, float j1, float i2, float j2, float i) {
	return j1 + (j2 - j1) * (i - i1) / (i2 - i1);
}

// 将 chunks 中所有 Triangle 投影到 DepthMap 上
template<int N, int M>
void DepthMap<N, M>::project(ChunkSet& chunks) {
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			map[i][j].clear();
	for (Triangle& tri : chunks.tris) {
		if (camera.Vecd * tri.v1.normal > 0 &&
			camera.Vecd * tri.v2.normal > 0 &&
			camera.Vecd * tri.v3.normal > 0) continue;
		float
			i1 = camera.i(tri.v1.pos), j1 = camera.j(tri.v1.pos),
			i2 = camera.i(tri.v2.pos), j2 = camera.j(tri.v2.pos),
			i3 = camera.i(tri.v3.pos), j3 = camera.j(tri.v3.pos);
		float cross = (i1 - i2) * (j3 - j2) - (i3 - i2) * (j1 - j2);
		if (fabsf(cross) < 1e-4) continue;
		float
			icoef1 = (j3 - j2) / cross,
			jcoef1 = -(i3 - i2) / cross,
			ccoef1 = 1 - i1 * icoef1 - j1 * jcoef1,
			icoef2 = (j1 - j3) / cross,
			jcoef2 = -(i1 - i3) / cross,
			ccoef2 = 1 - i2 * icoef2 - j2 * jcoef2;
		if (i3 < i2) { std::swap(i3, i2); std::swap(j3, j2); }
		if (i2 < i1) { std::swap(i2, i1); std::swap(j2, j1); }
		if (i3 < i2) { std::swap(i3, i2); std::swap(j3, j2); }
		if (i2 - i1 < 1e-4f) i2 += 1e-3f;
		if (i3 - i2 < 1e-4f) i3 += 1e-3f;
		int iL = i1 - 1e-4f, iR = i3;
		iR = min(iR, N - 1);
		float d1 = camera.d(tri.v1.pos), d2 = camera.d(tri.v2.pos), d3 = camera.d(tri.v3.pos);
		for (int i = max(iL + 1, 0); i <= iR; i++) {
			float
				_jL = intersection(i1, j1, i3, j3, i),
				_jR = i < i2 ? intersection(i1, j1, i2, j2, i) : intersection(i2, j2, i3, j3, i);
			if (_jL > _jR)std::swap(_jL, _jR);
			int jL = _jL - 1e-4f, jR = _jR + 1e-4f;
			float
				dk = jcoef1 * (d1 - d3) + jcoef2 * (d2 - d3),
				db = (icoef1 * i + ccoef1) * (d1 - d3) + (icoef2 * i + ccoef2) * (d2 - d3) + d3;
			jR = min(jR, M - 1);
			for (int j = max(jL + 1, 0); j <= jR; j++) {
				float depth = dk * j + db;
				if (depth < map[i][j].depth) {
					map[i][j].depth = depth;
					map[i][j].tri = &tri;
				}
			}
		}
	}
}


// 计算每个 Pixel 的 Color 和 normal
template<int N, int M>
void DepthMap<N, M>::getScreen(Color colMap[N][M], Point posMap[N][M], Vector normalMap[N][M]) {
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
			if (map[i][j].tri == nullptr) {
				colMap[i][j] = make_col(0, 0, 0, 1);
				continue;
			}
			Triangle& tri = *map[i][j].tri;
			float
				i1 = camera.i(tri.v1.pos), j1 = camera.j(tri.v1.pos),
				i2 = camera.i(tri.v2.pos), j2 = camera.j(tri.v2.pos),
				i3 = camera.i(tri.v3.pos), j3 = camera.j(tri.v3.pos);
			float
				cross = (i1 - i2) * (j3 - j2) - (i3 - i2) * (j1 - j2),
				k1 = ((i - i2) * (j3 - j2) - (j - j2) * (i3 - i2)) / cross,
				k2 = ((i - i3) * (j1 - j3) - (j - j3) * (i1 - i3)) / cross,
				k3 = 1 - k1 - k2;
			colMap[i][j] = k1 * tri.v1.col + k2 * tri.v2.col + k3 * tri.v3.col;
			posMap[i][j] = k1 * tri.v1.pos + k2 * tri.v2.pos + k3 * tri.v3.pos;
			normalMap[i][j] = k1 * tri.v1.normal + k2 * tri.v2.normal + k3 * tri.v3.normal;
		}
}

// 测试 pos 是否在阴影中.
template<int N, int M>
inline bool DepthMap<N, M>::shadowTest(Point& pos) {
	float i = camera.i(pos), j = camera.j(pos), depth = camera.d(pos);
	int i0 = i + 0.5f, j0 = j + 0.5f;
	if (i0 < 0 || i0 >= N || j0 < 0 || j0 >= M)
		return false;
	return depth > map[i0][j0].depth;
}

// 计算三角形 tri 对点 (i,j,depth) 的覆盖度.
// 若点在 tri 内, 覆盖度为 1.0; 否则计算点到 tri 的距离, 覆盖度为 1.0 - dis/BlurDistance
float cover(Triangle* _tri, float i, float j, float depth, Camera& camera) {
	if (_tri == nullptr) return .0f;
	Triangle& tri = *_tri;
	float
		i1 = camera.i(tri.v1.pos), j1 = camera.j(tri.v1.pos),
		i2 = camera.i(tri.v2.pos), j2 = camera.j(tri.v2.pos),
		i3 = camera.i(tri.v3.pos), j3 = camera.j(tri.v3.pos);
	float
		cross = (i1 - i2) * (j3 - j2) - (i3 - i2) * (j1 - j2),
		k1 = ((i - i2) * (j3 - j2) - (j - j2) * (i3 - i2)) / cross,
		k2 = ((i - i3) * (j1 - j3) - (j - j3) * (i1 - i3)) / cross,
		k3 = 1 - k1 - k2;
	float depth2 =
		k1 * camera.d(tri.v1.pos) +
		k2 * camera.d(tri.v2.pos) +
		k3 * camera.d(tri.v3.pos);
	if (depth < depth2) return .0;
	static constexpr float eps = 1e-4f;
	char mask = (k1 < -eps) | ((k2 < -eps) << 1) | ((k3 < -eps) << 2);
	float BlurDistance = min(1.0f, (depth - depth2) * 0.04f); // 软阴影宽度
	return
		(mask == 0) ? 1.0 :
		(mask == 4) ? max(1.0f - vec2::dis(i1, j1, i2, j2, i, j) / BlurDistance, .0f) :
		(mask == 2) ? max(1.0f - vec2::dis(i1, j1, i3, j3, i, j) / BlurDistance, .0f) :
		(mask == 1) ? max(1.0f - vec2::dis(i2, j2, i3, j3, i, j) / BlurDistance, .0f) :
		(mask == 6) ? max(1.0f - vec2::dis(i1, j1, i, j) / BlurDistance, .0f) :
		(mask == 5) ? max(1.0f - vec2::dis(i2, j2, i, j) / BlurDistance, .0f) :
		(mask == 3) ? max(1.0f - vec2::dis(i3, j3, i, j) / BlurDistance, .0f) : .0f;
}

// 枚举附近的三角面, 精确地测试 pos 是否在阴影中.
template<int N, int M>
float DepthMap<N, M>::accurateShadowTest(Point& pos) {
	float i = camera.i(pos), j = camera.j(pos), depth = camera.d(pos);
	int i0 = i + 0.5f, j0 = j + 0.5f;
	static constexpr int DetectionRadius = 1; // 枚举三角面的范围
	if (i < DetectionRadius || i + DetectionRadius >= N || j < DetectionRadius || j + DetectionRadius >= M)
		return .0f;
	if (map[i0][j0].depth < depth &&
		map[i0 + 1][j0].depth < depth &&
		map[i0 - 1][j0].depth < depth &&
		map[i0][j0 + 1].depth < depth &&
		map[i0][j0 - 1].depth < depth)
		return 1.0f;
	timeStamp++;
	float ret = .0f;
	for (int i2 = i0 - DetectionRadius; i2 <= i0 + DetectionRadius; i2++)
		for (int j2 = j0 - DetectionRadius; j2 <= j0 + DetectionRadius; j2++)
			if (map[i2][j2].depth < depth) {
				if (map[i2][j2].tri->stamp == timeStamp)
					continue;
				map[i2][j2].tri->stamp = timeStamp;
				ret = max(ret, cover(map[i2][j2].tri, i, j, depth, camera));
			}
	return ret;
}

// 枚举附近的三角面, 计算点 Pixel (i,j) 的 AO (环境光遮蔽).
template<int N, int M>
float rangeAOTest(Point& pos, DepthData depthMap[N][M], Vector normalMap[N][M], Camera& camera) {
	static constexpr struct {
		int di, dj;
		float weight;
	} Filter[] = {
		{-3,-1,1.0000},{-3,1,1.0000},{-2,-2,1.7321},{-2,0,2.6458},{-2,2,1.7321},{-1,-3,1.0000},{-1,-1,3.0000},
		{-1,1,3.0000},{-1,3,1.0000},{0,-2,2.6458},{0,0,3.3166},{0,2,2.6458},{1,-3,1.0000},{1,-1,3.0000},
		{1,1,3.0000},{1,3,1.0000},{2,-2,1.7321},{2,0,2.6458},{2,2,1.7321},{3,-1,1.0000},{3,1,1.0000},
	}; // 棋盘采样
	float i = camera.i(pos), j = camera.j(pos), depth = camera.d(pos);
	int i0 = i, j0 = j;
	if (i0 < 3 || i0 + 3 >= N || j0 < 3 || j0 + 3 >= M)
		return .0f;
	static constexpr float
		DistanceCutoff = 8.0f,	// 深度差达到 CutoffDistance 时，AO 阴影淡出
		AngleCutoff = 0.2f;		// 角度接近垂直时，AO 阴影淡出
	float light = .0f, total = .0f;
	for (auto& sample : Filter) {
		int i2 = i0 + sample.di,
			j2 = j0 + sample.dj;
		float
			weight = sample.weight,
			depthDet = depth - depthMap[i2][j2].depth,
			dot = normalMap[i2][j2] * camera.Vecd;
		if (depthDet > DistanceCutoff)
			weight *= DistanceCutoff / depthDet;
		if (dot > -AngleCutoff)
			weight *= dot / (-AngleCutoff);
		light += max(min(- depthDet, weight), -weight);
		total += weight;
	}
	return (light + total) / (2.0f * total + 1e-4f);
}

constexpr int TileSize = 4; // 分类 Tile 的大小.

template<int N, int M>
inline char checkNeighbors(int i, int j, char tile[N / TileSize][M / TileSize]) {
	char mask = 0;
	for (int i2 = max(i - 1, 0); i2 <= min(i + 1, N / TileSize - 1); i2++)
		for (int j2 = max(j - 1, 0); j2 <= min(j + 1, M / TileSize - 1); j2++)
			mask |= tile[i2][j2];
	return mask;
}

template<int N, int M>
void getAOTiles(DepthData depthMap[N][M], Vector normalMap[N][M], char tile[N / TileSize][M / TileSize]) {
	static constexpr float AngularCutoff = 0.8f, DepthCutoff = 8.0f; // 判断两个 Pixel 是否发生转折
	memset(tile, 0, sizeof(char) * (N / TileSize) * (M / TileSize));
	for (int i = 1; i < N; i++)
		for (int j = 1; j < M; j++)
			tile[i / TileSize][j / TileSize] |=
			fabsf(depthMap[i][j].depth - depthMap[i - 1][j].depth) > DepthCutoff ||
			normalMap[i][j] * normalMap[i - 1][j] < AngularCutoff ||
			fabsf(depthMap[i][j].depth - depthMap[i][j - 1].depth) > DepthCutoff ||
			normalMap[i][j] * normalMap[i][j - 1] < AngularCutoff;
}

#define forTile(i,j,i0,j0) \
	for (int ti = (i0) * TileSize, tj = (j0) * TileSize, i = ti; i < ti + TileSize; i++) \
		for (int j = tj; j < tj + TileSize; j++)

int SSAO_TIME, PCSS_TIME;

template<int N, int M>
void SSAO(DepthMap<N, M>& depthMap, Point posMap[N][M], Vector normalMap[N][M], float lightMap[N][M]) {
	SSAO_TIME -= clock();
	static char tile[N / TileSize][M / TileSize];
	static constexpr float SampleRadius = 3.5f; // 采样半径，即法向偏移量
	getAOTiles<N, M>(depthMap.map, normalMap, tile);
	for (int i0 = 0; i0 < N / TileSize; i0++)
		for (int j0 = 0; j0 < M / TileSize; j0++) {
			int ti = i0 * TileSize, tj = j0 * TileSize;
			if (!checkNeighbors<N, M>(i0, j0, tile)) {
				forTile(i, j, i0, j0)
					lightMap[i][j] = 1.0f;
				continue;
			}
			forTile(i, j, i0, j0) {
				if (depthMap.map[i][j].depth > 1e8f)
					continue;
				Point pos = posMap[i][j] + (SampleRadius + 0.5f) * normalMap[i][j];
				lightMap[i][j] = rangeAOTest<N, M>(pos, depthMap.map, normalMap, depthMap.camera);
			}
		}
	SSAO_TIME += clock();
}

template<int N, int M>
void Shadow(DepthMap<N, M>& cameraDepthMap, Vector cameraNormalMap[N][M], DepthMap<N, M>& sunDepthMap, float lightMap[N][M]) {
	PCSS_TIME -= clock();
	static float NdotL[N][M];
	static char tile[N / TileSize][M / TileSize];
	static constexpr float BIAS = 2.0f; // 法向偏移量
	Camera& camera = cameraDepthMap.camera, & sun = sunDepthMap.camera;
	memset(tile, 0, sizeof(char) * (N / TileSize) * (M / TileSize));
	// 本影采样
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			NdotL[i][j] = -(cameraNormalMap[i][j] * sun.Vecd);
	for (int i = 0; i < N; i++)
		for (int j = i & 1; j < M; j += 2) {
			if (NdotL[i][j] > 0) {
				Point pos =
					camera.pos(i, j, cameraDepthMap.map[i][j].depth)
					+ BIAS * cameraNormalMap[i][j];
				tile[i / TileSize][j / TileSize] |= (1 << sunDepthMap.shadowTest(pos));
			}
		}
	for (int i0 = 0; i0 < N / TileSize; i0++)
		for (int j0 = 0; j0 < M / TileSize; j0++) {
			char mask = checkNeighbors<N, M>(i0, j0, tile);
			if (mask == 3) { // mask == 3 表示该区域可能有半影，进行昂贵的半影采样
				forTile(i, j, i0, j0) {
					if (NdotL[i][j] < 0)
						lightMap[i][j] = .0f;
					else {
						Point pos =
							camera.pos(i, j, cameraDepthMap.map[i][j].depth)
							+ BIAS * cameraNormalMap[i][j];
						lightMap[i][j] = NdotL[i][j] * (1.0f - sunDepthMap.accurateShadowTest(pos));
					}
				}
			}
			else { // 该区域无半影，按照 mask 信息直接赋值
				if (mask == 1)
					forTile(i, j, i0, j0)
					lightMap[i][j] = max(NdotL[i][j], .0f);
				else
					forTile(i, j, i0, j0)
					lightMap[i][j] = .0f;
			}
		}
	PCSS_TIME += clock();
}

// 3x3 高斯模糊
template<int N, int M>
inline void blur(float map[N][M]) {
	static float sav[N][M];
	for (int i = 1; i + 1 < N; i++)
		for (int j = 1; j + 1 < M; j++)
			sav[i][j] = (map[i - 1][j] + 2.0f * map[i][j] + map[i + 1][j]) / 4.0f;
	for (int i = 1; i + 1 < N; i++)
		for (int j = 1; j + 1 < M; j++)
			map[i][j] = (sav[i][j - 1] + 2.0f * sav[i][j] + sav[i][j + 1]) / 4.0f;
}

template<int N, int M>
class RayBoard {
private:
	DepthMap<N, M> cameraDepthMap;
	DepthMap<N, M> sunDepthMap;
	Point cameraPosMap[N][M];
	Vector cameraNormalMap[N][M]; // cameraMap 的法向
	float
		environmentLight[N][M], // 环境光强度
		sunLight[N][M]; // 阳光强度
	ThreadPool<2> threadPool;
	Camera camera, sun;
public:
	RayBoard() {
		threadPool.init();
	}
	void setCamera(float alpha, float beta, Point Mpos) {
		camera.set(alpha, beta, Mpos);
		camera.move(-0.5f * N, -0.5f * M);
	}
	void setSun(float alpha, float beta, Point Mpos) {
		sun.set(alpha, beta, Mpos);
		sun.move(-0.5f * N, -0.5f * M);
	}
	void render(ChunkSet& chunks, Picture<N, M>& pic);
};

int cameraRender_time, sunRender_time;

template<int N, int M>
void cameraRender(DepthMap<N, M>& cameraDepthMap, Point cameraPosMap[N][M], Vector cameraNormalMap[N][M], ChunkSet& chunks, Picture<N, M>& pic, Camera& camera) {
	cameraRender_time -= clock();
	cameraDepthMap.camera = camera;
	cameraDepthMap.project(chunks);
	cameraDepthMap.getScreen(pic.col, cameraPosMap, cameraNormalMap);
	cameraRender_time += clock();
}

template<int N, int M, typename Args>
void run_cameraRender(void* _args) {
	Args* args = (Args*)_args;
	cameraRender<N, M>(*args->cameraDepthMap, args->cameraPosMap, args->cameraNormalMap, *args->chunks, *args->pic, *args->camera);
}

template<int N, int M>
void sunRender(DepthMap<N, M>& sunDepthMap, ChunkSet& chunks, Camera& sun) {
	sunRender_time -= clock();
	sunDepthMap.camera = sun;
	sunDepthMap.project(chunks);
	sunRender_time += clock();
}

template<int N, int M, typename Args>
void run_sunRender(void* _args) {
	Args* args = (Args*)_args;
	sunRender<N, M>(*args->sunDepthMap, *args->chunks, *args->sun);
}

template<int N, int M, typename Args>
void run_Shadow(void* _args) {
	Args* args = (Args*)_args;
	Shadow<N, M>(*args->cameraDepthMap, args->cameraNormalMap, *args->sunDepthMap, args->lightMap);
}

template<int N, int M, typename Args>
void run_SSAO(void* _args) {
	Args* args = (Args*)_args;
	SSAO<N, M>(*args->depthMap, args->posMap, args->normalMap, args->lightMap);
	blur<N, M>(args->lightMap);
}

// 将 chunks 的渲染结果输出到 pic 中
template<int N, int M>
void RayBoard<N, M>::render(ChunkSet& chunks, Picture<N, M>& pic) {

	struct CameraRenderArgs {
		DepthMap<N, M>* cameraDepthMap;
		Point(*cameraPosMap)[M];
		Vector(*cameraNormalMap)[M];
		ChunkSet* chunks;
		Picture<N, M>* pic;
		Camera* camera;
	} camera_args = { &cameraDepthMap, cameraPosMap, cameraNormalMap, &chunks, &pic, &camera };
	threadPool.setTask(0, run_cameraRender<N, M, CameraRenderArgs>, (void*)&camera_args);

	struct SunRenderArgs {
		DepthMap<N, M>* sunDepthMap;
		ChunkSet* chunks;
		Camera* sun;
	} sun_args = { &sunDepthMap, &chunks, &sun };
	threadPool.setTask(1, run_sunRender<N, M, SunRenderArgs>, (void*)&sun_args);

	threadPool.work();

	struct ShadowArgs {
		DepthMap<N, M>* cameraDepthMap;
		Vector(*cameraNormalMap)[M];
		DepthMap<N, M>* sunDepthMap;
		float(*lightMap)[M];
	} Shadow_args = { &cameraDepthMap, cameraNormalMap, &sunDepthMap, sunLight };
	threadPool.setTask(0, run_Shadow<N, M, ShadowArgs>, (void*)&Shadow_args);

	struct SSAOArgs {
		DepthMap<N, M>* depthMap;
		Point(*posMap)[M];
		Vector(*normalMap)[M];
		float(*lightMap)[M];
	} SSAO_args = { &cameraDepthMap, cameraPosMap, cameraNormalMap, environmentLight };
	threadPool.setTask(1, run_SSAO<N, M, SSAOArgs>, (void*)&SSAO_args);

	threadPool.work();

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			pic.col[i][j] = (0.01 + 0.24 * environmentLight[i][j] + 0.7 * sunLight[i][j]) * pic.col[i][j];
}

#endif