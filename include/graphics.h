#ifndef _GRAPHICS_H_
#define _GRAPHICS_H_

#include <SDL2/SDL.h>
#include <algorithm>
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
	void print(SDL_Surface* img, int PixelScale);
};

constexpr int ColorDivision = 32768;
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

// 将 Color 转化为 Uint32
inline Uint32 getColorid(Color& col) {
	float buf[4];
	_mm_store_ps(buf, col.vec);
	return
		(_channel(buf[0]) << 16) |
		(_channel(buf[1]) << 8) |
		_channel(buf[2]);
}

// 从 Picture 输出到图形库的 SDL_Surface 类，像素颗粒大小为 PixelScale
template<int N, int M>
void Picture<N, M>::print(SDL_Surface* img, int PixelScale) {
	auto pixels = (Uint32*)img->pixels;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
			Uint32 nowCol = getColorid(col[i][j]);
			for (int i2 = 0; i2 < PixelScale; i2++)
				for (int j2 = 0; j2 < PixelScale; j2++)
					pixels[(i * PixelScale + i2) * M * PixelScale + (j * PixelScale + j2)] = nowCol;
		}
}

struct Initializer {
	Initializer() {
		if (SDL_Init(SDL_INIT_EVERYTHING) == -1) {
			SDL_Log("SDL initialization failed!%s\n", SDL_GetError());
			exit(-1);
		}
		initGamma();
	}
};

// SDL 库创建的图形窗口
class Window {
private:
	SDL_Window* window;
public:
	SDL_Surface* surface;
	// surface 是 Uint32 数组, 代表屏幕

	Window(int width, int height, const char* title) {
		window = SDL_CreateWindow(
			title,
			SDL_WINDOWPOS_CENTERED,
			SDL_WINDOWPOS_CENTERED,
			width,
			height,
			SDL_WINDOW_SHOWN
		);
		if (!window) {
			SDL_Log("create window failed!%s\n", SDL_GetError());
			exit(-1);
		}
		surface = SDL_GetWindowSurface(window);
	}
	~Window() {
		SDL_DestroyWindow(window);
	}
	void flush() {
		SDL_UpdateWindowSurface(window);
	}
};

#endif
