#ifndef _DISPLAY_H_
#define _DISPLAY_H_

#include <immintrin.h>
#include <graphics.h>
#include <algorithm>
#include <math.h>
#include <cstdio>
#include "geometry.h"

#undef min
#undef max

using std::min;
using std::max;

// 颜色表示为向量
using Color = Vector;

Color make_col(float r, float g, float b, float alpha) {
	return { _mm_setr_ps(r, g, b, alpha) };
}

// 按权重归一化，用于计算平均颜色
inline Color normalize(const Color &col) {
	float buf[4];
	_mm_store_ps(buf, col.vec);
	if (buf[3] == 0.0)
		return { _mm_setzero_ps() };
	return { _mm_setr_ps(
		buf[0] / buf[3],
		buf[1] / buf[3],
		buf[2] / buf[3],
		1.0
	) };
}

// 用 Color 记录的图片
template<int N, int M>
struct Picture {
	Color col[N][M];
	void print(IMAGE& img, int PixelScale);
};

const int ColorDivision = 32768;
int gammaTable[ColorDivision];

void initGamma() {
	for (int i = 0; i < ColorDivision; i++) {
		float x = 1.0f / ColorDivision * i;
		gammaTable[i] = 256 * pow(x, 1.0f / 2.2f);
	}
}

// 将 [0,1] 内的浮点 rgb 转化为 0~255 的整数.
inline int _channel(float x) {
	int id = x * ColorDivision;
	return gammaTable[min(id,32767)];
}

// 将 Color 转化为 DWORD
inline DWORD getColorid(Color& col) {
	float buf[4];
	_mm_store_ps(buf, col.vec);
	return
		(_channel(buf[0]) << 16) |
		(_channel(buf[1]) << 8) |
		_channel(buf[2]);
}

// 从 Picture 输出到图形库的 IMAGE 类，像素颗粒大小为 PixelScale
template<int N, int M>
void Picture<N, M>::print(IMAGE& img, int PixelScale) {
	DWORD* buf = GetImageBuffer(&img);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
			DWORD nowCol = getColorid(col[i][j]);
			for (int i2 = 0; i2 < PixelScale; i2++)
				for (int j2 = 0; j2 < PixelScale; j2++)
					buf[(i * PixelScale + i2) * M * PixelScale + (j * PixelScale + j2)] = nowCol;
		}
}

template<int N, int M>
void add(Picture<N, M>& A, Picture<N, M>& B) {
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			A.col[i][j] = A.col[i][j] + B.col[i][j];
}

#endif
