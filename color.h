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

// ��ɫ��ʾΪ����
using Color = Vector;

Color make_col(float r, float g, float b, float alpha) {
	return { _mm_setr_ps(r, g, b, alpha) };
}

// ��Ȩ�ع�һ�������ڼ���ƽ����ɫ
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

// �� Color ��¼��ͼƬ
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

// �� [0,1] �ڵĸ��� rgb ת��Ϊ 0~255 ������.
inline int _channel(float x) {
	int id = x * ColorDivision;
	return gammaTable[min(id,32767)];
}

// �� Color ת��Ϊ DWORD
inline DWORD getColorid(Color& col) {
	float buf[4];
	_mm_store_ps(buf, col.vec);
	return
		(_channel(buf[0]) << 16) |
		(_channel(buf[1]) << 8) |
		_channel(buf[2]);
}

// �� Picture �����ͼ�ο�� IMAGE �࣬���ؿ�����СΪ PixelScale
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
