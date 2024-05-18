#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <immintrin.h>
#include <cstdio>
#include <cmath>

const float pi = acos(-1);

struct Vector {
	__m128 vec;
};

// SIMD 加速向量的基本运算
inline Vector operator + (const Vector& ColA, const Vector& ColB) {
	return { _mm_add_ps(ColA.vec, ColB.vec) };
}

inline Vector operator - (const Vector& VecA, const Vector& VecB) {
	return { _mm_sub_ps(VecA.vec, VecB.vec) };
}

inline Vector operator * (const float& scale, const Vector& Vec) {
	return { _mm_mul_ps(_mm_set1_ps(scale), Vec.vec) };
}

void print(Vector p) {
	float a[4];
	_mm_store_ps(a, p.vec);
	printf("(%.3f,%.3f,%.3f,%.3f)\n", a[0], a[1], a[2], a[3]);
}

// 三维向量，第四位总是 0
Vector make_vec(float x, float y, float z) {
	return { _mm_setr_ps(x, y, z, .0f) };
}

// 向量点积
inline float operator * (const Vector& VecA, const Vector& VecB) {
	return _mm_cvtss_f32(_mm_dp_ps(VecA.vec, VecB.vec, 0xF1));
}

using Point = Vector;

inline float InvSqrt(float x) {
	float xhalf = 0.5f * x;
	int i = *(int*)&x; // get bits for floating value
	i = 0x5f3759df - (i >> 1); // gives initial guess y0
	x = *(float*)&i; // convert bits back to float
	x = x * (1.5f - xhalf * x * x); // Newton step, repeating increases accuracy
	return x;
}

// 二维计算几何
namespace vec2 {
	float norm(float x, float y) {
		return sqrtf(x * x + y * y);
	}

	float dis(float x1, float y1, float x2, float y2) {
		return norm(x1 - x2, y1 - y2);
	}

	// 点 (x3,y3) 到线段 (x1,y1)-(x2,y2) 的距离
	inline float dis(float x1, float y1, float x2, float y2, float x3, float y3) {
		x2 -= x1; y2 -= y1;
		x3 -= x1; y3 -= y1;
		float dot = x2 * x3 + y2 * y3;
		if (dot < 0) return norm(x3, y3);
		if (dot > x2 * x2 + y2 * y2) return dis(x2, y2, x3, y3);
		return fabsf(x2 * y3 - x3 * y2) / norm(x2, y2);
	}
};

// 三维计算几何
namespace vec3 {
	inline float norm(Point p) {
		return sqrtf(p * p);
	}

	inline float dis(Point p1, Point p2) {
		return norm(p1 - p2);
	}

	// 向量叉积
	inline Vector cross(const Vector& VecA, const Vector& VecB) {
		float A[4], B[4];
		_mm_store_ps(A, VecA.vec);
		_mm_store_ps(B, VecB.vec);
		return make_vec(
			A[1] * B[2] - B[1] * A[2],
			B[0] * A[2] - A[0] * B[2],
			A[0] * B[1] - B[0] * A[1]
		);
	}

	// 点 p3 到线段 p1-p2 的距离
	inline float dis(Point p1, Point p2, Point p3) {
		p2 = p2 - p1;
		p3 = p3 - p1;
		float dot = p2 * p3;
		if (dot < 0) return norm(p3);
		if (dot > p2 * p2) return dis(p2, p3);
		return norm(cross(p3, p2)) / norm(p2);
	}
};

#endif