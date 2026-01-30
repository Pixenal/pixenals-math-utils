/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#pragma once

#include <math.h>

#include "../../pixenals-types/include/pixenals_types.h"

#define PIX_MATH_EPSILON .0000001f
#ifdef PIX_MATH_USE_EPSILON
#define F32_EQL_INTERN(a, b) (fabsf((a) - (b)) <= PIX_MATH_EPSILON)
#define F64_EQL_INTERN(a, b) (fabs((a) - (b)) <= PIX_MATH_EPSILON)
#define F32_NOTEQL_INTERN(a, b) (fabsf((a) - (b)) > PIX_MATH_EPSILON)
#define F64_NOTEQL_INTERN(a, b) (fabs((a) - (b)) > PIX_MATH_EPSILON)
#define F32_GREAT_INTERN(a, b) ((a) - (b) > PIX_MATH_EPSILON)
#define F64_GREAT_INTERN(a, b) F32_GREAT_INTERN(a, b)
#define F32_LESS_INTERN(a, b) ((a) - (b) < -PIX_MATH_EPSILON)
#define F64_LESS_INTERN(a, b) F32_LESS_INTERN(a, b)
#define F32_GREATEQL_INTERN(a, b) (F32_GREAT_INTERN(a, b) || F32_EQL_INTERN(a, b))
#define F64_GREATEQL_INTERN(a, b) (F64_GREAT_INTERN(a, b) || F64_EQL_INTERN(a, b))
#define F32_LESSEQL_INTERN(a, b) (F32_LESS_INTERN(a, b) || F32_EQL_INTERN(a, b))
#define F64_LESSEQL_INTERN(a, b) (F64_LESS_INTERN(a, b) || F64_EQL_INTERN(a, b))
#else
#define pixmF32_EQL_INTERN(a, b) ((a) == (b))
#define pixmF64_EQL_INTERN(a, b) ((a) == (b))
#define pixmF32_NOTEQL_INTERN(a, b) ((a) != (b))
#define pixmF64_NOTEQL_INTERN(a, b) ((a) != (b))
#define pixmF32_GREAT_INTERN(a, b) ((a) > (b))
#define pixmF64_GREAT_INTERN(a, b) ((a) > (b))
#define pixmF32_LESS_INTERN(a, b) ((a) < (b))
#define pixmF64_LESS_INTERN(a, b) ((a) < (b))
#define pixmF32_GREATEQL_INTERN(a, b) ((a) >= (b))
#define pixmF64_GREATEQL_INTERN(a, b) ((a) >= (b))
#define pixmF32_LESSEQL_INTERN(a, b) ((a) <= (b))
#define pixmF64_LESSEQL_INTERN(a, b) ((a) <= (b))
#endif
#define F32_EQL ,F32_EQL_INTERN,
#define F64_EQL ,F64_EQL_INTERN,
#define F32_NOTEQL ,F32_NOTEQL_INTERN,
#define F64_NOTEQL ,F64_NOTEQL_INTERN,
#define F32_GREAT ,F32_GREAT_INTERN,
#define F64_GREAT ,F64_GREAT_INTERN,
#define F32_LESS ,F32_LESS_INTERN,
#define F64_LESS ,F64_LESS_INTERN,
#define F32_GREATEQL ,F32_GREATEQL_INTERN,
#define F64_GREATEQL ,F64_GREATEQL_INTERN,
#define F32_LESSEQL ,F32_LESSEQL_INTERN,
#define F64_LESSEQL ,F64_LESSEQL_INTERN,


#define PIX_MATH_IDENT_MAT4X4 (PixtyM4x4) {\
	1.0, .0, .0, .0,\
	.0, 1.0, .0, .0,\
	.0, .0, 1.0, .0,\
	.0, .0, .0, 1.0\
}

#define PIXM_MIN(a, b) (a < b ? a : b)
#define PIXM_MAX(a, b) (a > b ? a : b)

static inline
float pixmF32Lerp(float a, float b, float alpha) {
	return b * alpha + (1.0f - alpha) * a;
}

static inline
double pixmF64Lerp(double a, double b, double alpha) {
	return b * alpha + (1.0 - alpha) * a;
}

PixtyM3x3 pixmM3x3Adjugate(const PixtyM3x3 *pA);
PixtyV3_F32 pixmCartesianToBarycentric(
	const PixtyV3_F32 *pTri32,
	const PixtyV3_F32 *pPoint32,
	const PixtyV3_F32 *pNormal
);
PixtyV3_F32 pixmBarycentricToCartesian(const PixtyV3_F32 *pTri, PixtyV3_F32 point);

static inline
PixtyV4_F32 pixmV4F32MultiplyM4x4(PixtyV4_F32 a, const PixtyM4x4 *pB) {
	PixtyV4_F32 c = {0};
	c.d[0] = 
		a.d[0] * pB->d[0][0] +
		a.d[1] * pB->d[1][0] +
		a.d[2] * pB->d[2][0] +
		a.d[3] * pB->d[3][0];
	c.d[1] =
		a.d[0] * pB->d[0][1] +
		a.d[1] * pB->d[1][1] +
		a.d[2] * pB->d[2][1] +
		a.d[3] * pB->d[3][1];
	c.d[2] =
		a.d[0] * pB->d[0][2] +
		a.d[1] * pB->d[1][2] +
		a.d[2] * pB->d[2][2] +
		a.d[3] * pB->d[3][2];
	return c;
}

static inline
void pixmV4F32MultiplyEqualPixtyM4x4(PixtyV4_F32 *pA, const PixtyM4x4 *pB) {
	*pA = pixmV4F32MultiplyM4x4(*pA, pB);
}
#define V4MULEQLM4X4 ,V4F32MultiplyEqualPixtyM4x4,

static inline
PixtyV3_F32 pixmV3F32DivideScalar(PixtyV3_F32 a, float b) {
	return (PixtyV3_F32) {a.d[0] / b, a.d[1] / b, a.d[2] / b};
}
#define V3DIVS ,V3F32DivideScalar,

static inline
PixtyV3_F32 pixmDivideByW(const PixtyV4_F32 *pA) {
	return pixmV3F32DivideScalar((*(PixtyV3_F32 *)pA), (pA->d[3]));
}

static inline
PixtyV3_F32 pixmV3F32MultiplyScalar(PixtyV3_F32 a, float b) {
	return (PixtyV3_F32) {a.d[0] * b, a.d[1] * b, a.d[2] * b};
}
#define V3MULS ,V3F32MultiplyScalar,

static inline
PixtyV3_F64 pixmV3F64MultiplyScalar(PixtyV3_F64 a, double b) {
	return (PixtyV3_F64) {a.d[0] * b, a.d[1] * b, a.d[2] * b};
}

static inline
void pixmV3F32DivideEqualScalar(PixtyV3_F32 *pA, float b) {
	pA->d[0] /= b;
	pA->d[1] /= b;
	pA->d[2] /= b;
}
#define V3DIVEQLS ,V3F32DivideEqualScalar,

static inline
PixtyV3_F32 pixmV3F32SubtractScalar(PixtyV3_F32 a, float b) {
	return (PixtyV3_F32) {a.d[0] - b, a.d[1] - b, a.d[2] - b};
}
#define V3SUBS ,V3F32SubtractScalar,

static inline
PixtyV3_F32 pixmV3F32Subtract(PixtyV3_F32 a, PixtyV3_F32 b) {
	return (PixtyV3_F32) {a.d[0] - b.d[0], a.d[1] - b.d[1], a.d[2] - b.d[2]};
}
#define V3SUB ,V3F32Subtract,

static inline
PixtyV3_F64 pixmV3F64Subtract(PixtyV3_F64 a, PixtyV3_F64 b) {
	return (PixtyV3_F64) {a.d[0] - b.d[0], a.d[1] - b.d[1], a.d[2] - b.d[2]};
}

static inline
PixtyV3_F32 pixmV3F32AddScalar(PixtyV3_F32 a, float b) {
	return (PixtyV3_F32) {a.d[0] + b, a.d[1] + b, a.d[2] + b};
}
#define V3ADDS ,V3F32AddScalar,

static inline
PixtyV3_F32 pixmV3F32Add(PixtyV3_F32 a, PixtyV3_F32 b) {
	return (PixtyV3_F32) {a.d[0] + b.d[0], a.d[1] + b.d[1], a.d[2] + b.d[2]};
}
#define V3ADD ,V3F32Add,

static inline
PixtyV3_F64 pixmV3F64Add(PixtyV3_F64 a, PixtyV3_F64 b) {
	return (PixtyV3_F64) {a.d[0] + b.d[0], a.d[1] + b.d[1], a.d[2] + b.d[2]};
}

static inline
void pixmV3F32AddEqual(PixtyV3_F32 *pA, PixtyV3_F32 b) {
	pA->d[0] += b.d[0];
	pA->d[1] += b.d[1];
	pA->d[2] += b.d[2];
}
#define V3ADDEQL ,V3F32AddEqual,

static inline
bool pixmV3F32Equal(PixtyV3_F32 a, PixtyV3_F32 b) {
	return
		pixmF32_EQL_INTERN((a.d[0]), (b.d[0])) && pixmF32_EQL_INTERN((a.d[1]), (b.d[1])) && pixmF32_EQL_INTERN((a.d[2]), (b.d[2]));
}
#define V3EQL ,V3F32Equal,

static inline
bool pixmV3F64Equal(PixtyV3_F64 a, PixtyV3_F64 b) {
	return
		pixmF64_EQL_INTERN((a.d[0]), (b.d[0])) && pixmF64_EQL_INTERN((a.d[1]), (b.d[1])) && pixmF64_EQL_INTERN((a.d[2]), (b.d[2]));
}
#define V364EQL ,V3F64Equal,

static inline
bool pixmV3F32NotEqual(PixtyV3_F32 a, PixtyV3_F32 b) {
	return
		pixmF32_NOTEQL_INTERN((a.d[0]), (b.d[0])) ||
		pixmF32_NOTEQL_INTERN((a.d[1]), (b.d[1])) ||
		pixmF32_NOTEQL_INTERN((a.d[2]), (b.d[2]));
}
#define V3NOTEQL ,V3F32NotEqual,

static inline
bool pixmV3F32GreaterThan(PixtyV3_F32 a, PixtyV3_F32 b) {
	return
		pixmF32_GREAT_INTERN((a.d[0]), (b.d[0])) &&
		pixmF32_GREAT_INTERN((a.d[1]), (b.d[1])) &&
		pixmF32_GREAT_INTERN((a.d[2]), (b.d[2]));
}
#define V3GREAT ,V3F32GreaterThan,

static inline
bool pixmV3F32LessThan(PixtyV3_F32 a, PixtyV3_F32 b) {
	return
		pixmF32_LESS_INTERN((a.d[0]), (b.d[0])) &&
		pixmF32_LESS_INTERN((a.d[1]), (b.d[1])) &&
		pixmF32_LESS_INTERN((a.d[2]), (b.d[2]));
}
#define V3LESS ,V3F32LessThan,

static inline
PixtyV3_F32 pixmV3F32Lerp(PixtyV3_F32 a, PixtyV3_F32 b, float alpha) {
	float alphaInverse = 1.0f - alpha;
	return (PixtyV3_F32) {.d = {
		a.d[0] * alphaInverse + b.d[0] * alpha,
		a.d[1] * alphaInverse + b.d[1] * alpha,
		a.d[2] * alphaInverse + b.d[2] * alpha
	}};
}

static inline
PixtyV3_F32 pixmV3F32UnitFromPoints(PixtyV3_F32 a, PixtyV3_F32 b) {
	PixtyV3_F32 dir = pixmV3F32Subtract((b), (a)); //direction
	float magnitude = sqrtf(dir.d[0] * dir.d[0] + dir.d[1] * dir.d[1]);
	return pixmV3F32DivideScalar((dir), (magnitude));
}

static inline
PixtyV3_F32 pixmV3F32Cross(PixtyV3_F32 a, PixtyV3_F32 b) {
	return (PixtyV3_F32) {
		.d[0] = a.d[1] * b.d[2] - a.d[2] * b.d[1],
		.d[1] = a.d[2] * b.d[0] - a.d[0] * b.d[2],
		.d[2] = a.d[0] * b.d[1] - a.d[1] * b.d[0]
	};
}
#define V3CROSS ,V3F32Cross,

static inline
PixtyV3_F64 pixmV3F64Cross(PixtyV3_F64 a, PixtyV3_F64 b) {
	return (PixtyV3_F64) {
		.d[0] = a.d[1] * b.d[2] - a.d[2] * b.d[1],
		.d[1] = a.d[2] * b.d[0] - a.d[0] * b.d[2],
		.d[2] = a.d[0] * b.d[1] - a.d[1] * b.d[0]
	};
}

static inline
PixtyV3_F32 pixmV3F32MultiplyM3x3(PixtyV3_F32 a, const PixtyM3x3 *pB) {
	return (PixtyV3_F32) {.d = {
		a.d[0] * pB->d[0][0] + a.d[1] * pB->d[1][0] + a.d[2] * pB->d[2][0],
		a.d[0] * pB->d[0][1] + a.d[1] * pB->d[1][1] + a.d[2] * pB->d[2][1],
		a.d[0] * pB->d[0][2] + a.d[1] * pB->d[1][2] + a.d[2] * pB->d[2][2]
	}};
}
#define V3MULM3X3 ,V3F32MultiplyM3x3,

static inline
void pixmV3F32MultiplyEqualPixtyM3x3(PixtyV3_F32 *pA, const PixtyM3x3 *pB) {
	*pA = pixmV3F32MultiplyM3x3(*pA, pB);
}
#define V3MULEQLM3X3 ,V3F32MultiplyEqualPixtyM3x3,

static inline
PixtyV3_F32 pixmV3F32Normalize(PixtyV3_F32 a) {
	float magnitude = sqrtf(a.d[0] * a.d[0] + a.d[1] * a.d[1] + a.d[2] * a.d[2]);
	return pixmV3F32DivideScalar((a), (magnitude));
}

static inline
bool pixmV3F32IsFinite(PixtyV3_F32 a) {
	return isfinite(a.d[0]) && isfinite(a.d[1]) && isfinite(a.d[2]);
}

static inline
float pixmV3F32SquareLen(PixtyV3_F32 a) {
	return a.d[0] * a.d[0] + a.d[1] * a.d[1] + a.d[2] * a.d[2];
}

static inline
float pixmV3F32Len(PixtyV3_F32 a) {
	return sqrtf(pixmV3F32SquareLen(a));
}

static inline
float pixmV3F32TriArea(PixtyV3_F32 a, PixtyV3_F32 b, PixtyV3_F32 c) {
	PixtyV3_F32 ba = pixmV3F32Subtract((a), (b));
	PixtyV3_F32 bc = pixmV3F32Subtract((c), (b));
	PixtyV3_F32 cross = pixmV3F32Cross((ba), (bc));
	return pixmV3F32Len(cross) / 2.0f;
}

static inline
bool pixmV3F32DegenerateTri(PixtyV3_F32 a, PixtyV3_F32 b, PixtyV3_F32 c, float threshold) {
	PixtyV3_F32 ac = pixmV3F32Subtract((a), (c));
	PixtyV3_F32 bc = pixmV3F32Subtract((b), (c));
	PixtyV3_F32 cross = pixmV3F32Cross((ac), (bc));
	float len = pixmV3F32Len(cross);
	return pixmF32_LESSEQL_INTERN((len), (threshold)) && pixmF32_GREATEQL_INTERN((len), (-threshold));
}

static inline
float pixmV3F32TriHeight(PixtyV3_F32 a, PixtyV3_F32 b, PixtyV3_F32 c) {
	PixtyV3_F32 ac = pixmV3F32Subtract((a), (c));
	PixtyV3_F32 bc = pixmV3F32Subtract((b), (c));
	PixtyV3_F32 cross = pixmV3F32Cross((ac), (bc));
	return pixmV3F32Len(cross);
}

static inline
float pixmV3F32Dot(PixtyV3_F32 a, PixtyV3_F32 b) {
	return a.d[0] * b.d[0] + a.d[1] * b.d[1] + a.d[2] * b.d[2];
}
#define V3DOT ,V3F32Dot,

static inline
double pixmV3F64Dot(PixtyV3_F64 a, PixtyV3_F64 b) {
	return a.d[0] * b.d[0] + a.d[1] * b.d[1] + a.d[2] * b.d[2];
}

static inline
PixtyV2_F32 pixmV2F32Abs(PixtyV2_F32 a) {
	return (PixtyV2_F32) {fabsf(a.d[0]), fabsf(a.d[1])};
}

static inline
PixtyV2_F32 pixmV2F32Multiply(PixtyV2_F32 a, PixtyV2_F32 b) {
	return (PixtyV2_F32) {.d = {a.d[0] * b.d[0], a.d[1] * b.d[1]}};
}
#define V2MUL ,V2F32Multiply,

static inline
void pixmV2F32MultiplyEqual(PixtyV2_F32 *pA, PixtyV2_F32 b) {
	pA->d[0] *= b.d[0];
	pA->d[1] *= b.d[1];
}
#define V2MULEQL ,V2F32MultiplyEqual,

static inline
PixtyV2_F32 pixmV2F32DivideScalar(PixtyV2_F32 a, float b) {
	return (PixtyV2_F32) {.d = {a.d[0] / b, a.d[1] / b}};
}
#define V2DIVS ,V2F32DivideScalar,

static inline
void pixmV2F32DivideEqualScalar(PixtyV2_F32 *pA, float b) {
	pA->d[0] /= b;
	pA->d[1] /= b;
}
#define V2DIVSEQL ,V2F32DivideEqualScalar,

static inline
PixtyV2_F32 pixmV2F32Subtract(PixtyV2_F32 a, PixtyV2_F32 b) {
	return (PixtyV2_F32) {a.d[0] - b.d[0], a.d[1] - b.d[1]};
}
#define V2SUB ,V2F32Subtract,

static inline
PixtyV2_F64 pixmV2F64Subtract(PixtyV2_F64 a, PixtyV2_F64 b) {
	return (PixtyV2_F64) {a.d[0] - b.d[0], a.d[1] - b.d[1]};
}

static inline
void pixmV2F32SubtractEqual(PixtyV2_F32 *pA, PixtyV2_F32 b) {
	pA->d[0] -= b.d[0];
	pA->d[1] -= b.d[1];
}
#define V2SUBEQL ,V2F32SubtractEqual,

static inline
PixtyV2_F32 pixmV2F32SubtractScalar(PixtyV2_F32 a, float b) {
	return (PixtyV2_F32) {a.d[0] - b, a.d[1] - b};
}
#define V2SUBS ,V2F32SubtractScalar,

static inline
PixtyV2_F32 pixmV2F32Add(PixtyV2_F32 a, PixtyV2_F32 b) {
	return (PixtyV2_F32) {a.d[0] + b.d[0], a.d[1] + b.d[1]};
}
#define V2ADD ,V2F32Add,

static inline
PixtyV2_F32 pixmV2F32AddScalar(PixtyV2_F32 a, float b) {
	return (PixtyV2_F32) {a.d[0] + b, a.d[1] + b};
}
#define V2ADDS ,V2F32AddScalar,

static inline
void pixmV2F32AddEqual(PixtyV2_F32 *pA, PixtyV2_F32 b) {
	pA->d[0] += b.d[0];
	pA->d[1] += b.d[1];
}
#define V2ADDEQL ,V2F32AddEqual,

static inline
void pixmV2F32AddEqualScalar(PixtyV2_F32 *pA, float b) {
	pA->d[0] += b;
	pA->d[1] += b;
}
#define V2ADDEQLS ,V2F32AddEqualScalar,

static inline
void pixmV2F32MultiplyEqualScalar(PixtyV2_F32 *pA, float b) {
	pA->d[0] *= b;
	pA->d[1] *= b;
}
#define V2MULSEQL ,V2F32MultiplyEqualScalar,

static inline
PixtyV2_F32 pixmV2F32MultiplyScalar(PixtyV2_F32 a, float b) {
	return (PixtyV2_F32) {a.d[0] * b, a.d[1] * b};
}
#define V2MULS ,V2F32MultiplyScalar,

static inline
float pixmV2F32Dot(PixtyV2_F32 a, PixtyV2_F32 b) {
	return a.d[0] * b.d[0] + a.d[1] * b.d[1];
}
#define V2DOT ,V2F32Dot,

static inline
double pixmV2F64Dot(PixtyV2_F64 a, PixtyV2_F64 b) {
	return a.d[0] * b.d[0] + a.d[1] * b.d[1];
}

static inline
float pixmV2F32Cross(PixtyV2_F32 a, PixtyV2_F32 b) {
	return a.d[0] * b.d[1] - a.d[1] * b.d[0];
}
#define V2CROSS ,V2F32Cross,

static inline
double pixmV2F64Cross(PixtyV2_F64 a, PixtyV2_F64 b) {
	return a.d[0] * b.d[1] - a.d[1] * b.d[0];
}

static inline
PixtyV2_F32 pixmV2F32LineNormal(PixtyV2_F32 a) {
	PixtyV2_F32 b = {a.d[1], -a.d[0]};
	return b;
}

static inline
PixtyV2_F64 pixmV2F64LineNormal(PixtyV2_F64 a) {
	PixtyV2_F64 b = {a.d[1], -a.d[0]};
	return b;
}

static inline
PixtyV2_F32 pixmV2F32ModScalar(PixtyV2_F32 a, float b) {
	return (PixtyV2_F32) {fmodf(a.d[0], b), fmodf(a.d[1], b)};
}
#define V2MODS ,V2F32ModScalar,

static inline
float pixmV2F32SquareLen(PixtyV2_F32 a) {
	return a.d[0] * a.d[0] + a.d[1] * a.d[1];
}

static inline
double pixmV2F64SquareLen(PixtyV2_F64 a) {
	return a.d[0] * a.d[0] + a.d[1] * a.d[1];
}

static inline
float pixmV2F32Len(PixtyV2_F32 a) {
	return sqrtf(pixmV2F32SquareLen(a));
}

static inline
double pixmV2F64Len(PixtyV2_F64 a) {
	return sqrt(pixmV2F64SquareLen(a));
}

static inline
float pixmV2F32TriArea(PixtyV2_F32 a, PixtyV2_F32 b, PixtyV2_F32 c) {
	PixtyV2_F32 ba = pixmV2F32Subtract((a), (b));
	PixtyV2_F32 bc = pixmV2F32Subtract((c), (b));
	PixtyV3_F32 ba3 = {ba.d[0], ba.d[1], .0f};
	PixtyV3_F32 bc3 = {bc.d[0], bc.d[1], .0f};
	PixtyV3_F32 cross = pixmV3F32Cross((ba3), (bc3));
	return fabsf(cross.d[2]) / 2.0f;
}

static inline
float pixmV2F32Determinate(PixtyV2_F32 a, PixtyV2_F32 b) {
	return a.d[0] * b.d[1] - a.d[1] * b.d[0];
}
#define V2DET ,V2F32Determinate,

static inline
double pixmV2F64Determinate(PixtyV2_F64 a, PixtyV2_F64 b) {
	return a.d[0] * b.d[1] - a.d[1] * b.d[0];
}

static inline
void pixmV2F32ModEqualScalar(PixtyV2_F32 *a, float b) {
	a->d[0] = fmodf(a->d[0], b);
	a->d[1] = fmodf(a->d[1], b);
}
#define V2MODEQLS ,V2F32ModEqualScalar,

static inline
bool pixmV2F32GreaterThanEqualTo(PixtyV2_F32 a, PixtyV2_F32 b) {
	return pixmF32_GREATEQL_INTERN((a.d[0]), (b.d[0])) && pixmF32_GREATEQL_INTERN((a.d[1]), (b.d[1]));
}
#define V2GREATEQL ,V2F32GreaterThanEqualTo,

static inline
bool pixmV2F32GreaterThan(PixtyV2_F32 a, PixtyV2_F32 b) {
	return pixmF32_GREAT_INTERN((a.d[0]), (b.d[0])) && pixmF32_GREAT_INTERN((a.d[1]), (b.d[1]));
}
#define V2GREAT ,V2F32GreaterThan,

static inline
bool pixmV2F32GreaterThanScalar(PixtyV2_F32 a, float b) {
	return pixmF32_GREAT_INTERN((a.d[0]), (b)) && pixmF32_GREAT_INTERN((a.d[1]), (b));
}
#define V2GREATS ,V2F32GreaterThanScalar,

static inline
bool pixmV2F32LessThan(PixtyV2_F32 a, PixtyV2_F32 b) {
	return pixmF32_LESS_INTERN((a.d[0]), (b.d[0])) && pixmF32_LESS_INTERN((a.d[1]), (b.d[1]));
}
#define V2LESS ,V2F32LessThan,

static inline
bool pixmV2F32LessThanScalar(PixtyV2_F32 a, float b) {
	return pixmF32_LESS_INTERN((a.d[0]), (b)) && pixmF32_LESS_INTERN((a.d[1]), (b));
}
#define V2LESSS ,V2F32LessThanScalar,

static inline
bool pixmV2F32LessThanEqualTo(PixtyV2_F32 a, PixtyV2_F32 b) {
	return pixmF32_LESSEQL_INTERN((a.d[0]), (b.d[0])) && pixmF32_LESSEQL_INTERN((a.d[1]), (b.d[1]));
}
#define V2LESSEQL ,V2F32LessThanEqualTo,

static inline
bool pixmV2F32NotEqual(PixtyV2_F32 a, PixtyV2_F32 b) {
	return pixmF32_NOTEQL_INTERN((a.d[0]), (b.d[0])) || pixmF32_NOTEQL_INTERN((a.d[1]), (b.d[1]));
}
#define V2NOTEQL ,V2F32NotEqual,

//TODO replace return with bool in comparison funcs like this
static inline
bool pixmV2F32Equal(PixtyV2_F32 a, PixtyV2_F32 b) {
	return pixmF32_EQL_INTERN((a.d[0]), (b.d[0])) && pixmF32_EQL_INTERN((a.d[1]), (b.d[1]));
}
#define V2EQL ,V2F32Equal,

static inline
bool pixmV2F64Equal(PixtyV2_F64 a, PixtyV2_F64 b) {
	return pixmF64_EQL_INTERN((a.d[0]), (b.d[0])) && pixmF64_EQL_INTERN((a.d[1]), (b.d[1]));
}
#define V264EQL ,V2F64Equal,

static inline
bool pixmV2F32AproxEqual(PixtyV2_F32 a, PixtyV2_F32 b) {
	PixtyV2_F32 bLow = pixmV2F32SubtractScalar((b), (.0000001f));
	PixtyV2_F32 bHigh = pixmV2F32AddScalar((b), (.0000001f));
	return pixmV2F32LessThanEqualTo((a), (bHigh)) && pixmV2F32GreaterThanEqualTo((a), (bLow));
}
#define V2APROXEQL ,V2F32AproxEqual,

static inline
bool pixmV2F32AproxEqualThres(PixtyV2_F32 a, PixtyV2_F32 b, float threshold) {
	PixtyV2_F32 bLow = pixmV2F32SubtractScalar((b), (threshold));
	PixtyV2_F32 bHigh = pixmV2F32AddScalar((b), (threshold));
	return pixmV2F32LessThanEqualTo((a), (bHigh)) && pixmV2F32GreaterThanEqualTo((a), (bLow));
}

static inline
bool pixmV2F32DegenerateTri(PixtyV2_F32 a, PixtyV2_F32 b, PixtyV2_F32 c, float threshold) {
	PixtyV2_F32 ac = pixmV2F32Subtract((a), (c));
	PixtyV2_F32 bc = pixmV2F32Subtract((b), (c));
	float cross = ac.d[0] * bc.d[1] - bc.d[0] * ac.d[1];
	return pixmF32_LESSEQL_INTERN((cross), (threshold)) && pixmF32_GREATEQL_INTERN((cross), (-threshold));
}

static inline
float pixmV2F32TriHeight(PixtyV2_F32 a, PixtyV2_F32 b, PixtyV2_F32 c) {
	PixtyV2_F32 ac = pixmV2F32Subtract((a), (c));
	PixtyV2_F32 bc = pixmV2F32Subtract((b), (c));
	return ac.d[0] * bc.d[1] - bc.d[0] * ac.d[1];
}

static inline
bool pixmV2F32IsFinite(PixtyV2_F32 a) {
	return isfinite(a.d[0]) && isfinite(a.d[1]);
}

static inline
bool pixmV2I8Equal(PixtyV2_I8 a, PixtyV2_I8 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1];
}
#define V2I8EQL ,V2I8Equal,

static inline
bool pixmV2I16Equal(PixtyV2_I16 a, PixtyV2_I16 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1];
}
#define V2I16EQL ,V2I16Equal,

static inline
bool pixmV2I32Equal(PixtyV2_I32 a, PixtyV2_I32 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1];
}
#define V2IEQL ,V2I32Equal,

static inline
bool pixmV2I64Equal(PixtyV2_I64 a, PixtyV2_I64 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1];
}
#define V2I64EQL ,V2I64Equal,

static inline
bool pixmV3I8Equal(PixtyV3_I8 a, PixtyV3_I8 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1] && a.d[2] == b.d[2];
}
#define V3I8EQL ,V3I8Equal,

static inline
bool pixmV3I16Equal(PixtyV3_I16 a, PixtyV3_I16 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1] && a.d[2] == b.d[2];
}
#define V3I16EQL ,V3I16Equal,

static inline
bool pixmV3I32Equal(PixtyV3_I32 a, PixtyV3_I32 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1] && a.d[2] == b.d[2];
}
#define V3IEQL ,V3I32Equal,

static inline
bool pixmV3I64Equal(PixtyV3_I64 a, PixtyV3_I64 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1] && a.d[2] == b.d[2];
}
#define V3I64EQL ,V3I64Equal,

static inline
bool pixmV4I8Equal(PixtyV4_I8 a, PixtyV4_I8 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1] && a.d[2] == b.d[2] && a.d[3] == b.d[3];
}
#define V4I8EQL ,V4I8Equal,

static inline
bool pixmV4I16Equal(PixtyV4_I16 a, PixtyV4_I16 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1] && a.d[2] == b.d[2] && a.d[3] == b.d[3];
}
#define V4I16EQL ,V4I16Equal,

static inline
bool pixmV4I32Equal(PixtyV4_I32 a, PixtyV4_I32 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1] && a.d[2] == b.d[2] && a.d[3] == b.d[3];
}
#define V4IEQL ,V4I32Equal,

static inline
bool pixmV4I64Equal(PixtyV4_I64 a, PixtyV4_I64 b) {
	return a.d[0] == b.d[0] && a.d[1] == b.d[1] && a.d[2] == b.d[2] && a.d[3] == b.d[3];
}
#define V4I64EQL ,V4I64Equal,

static inline
PixtyM2x2 pixmM2x2Adjugate(PixtyM2x2 a) {
	PixtyM2x2 c = {0};
	c.d[0][0] = a.d[1][1];
	c.d[0][1] = -a.d[0][1];
	c.d[1][0] = -a.d[1][0];
	c.d[1][1] = a.d[0][0];
	return c;
}

static inline
float pixmM2x2Determinate(PixtyM2x2 a) {
	return a.d[0][0] * a.d[1][1] - a.d[0][1] * a.d[1][0];
}

static inline
void pixmM2x2MultiplyEqualScalar(PixtyM2x2 *pA, float b) {
	pA->d[0][0] *= b;
	pA->d[0][1] *= b;
	pA->d[1][0] *= b;
	pA->d[1][1] *= b;
}

static inline
PixtyM2x2 pixmM2x2Invert(PixtyM2x2 a) {
	float determinate = pixmM2x2Determinate(a);
	PixtyM2x2 inverse = pixmM2x2Adjugate(a);
	pixmM2x2MultiplyEqualScalar(&inverse, 1.0f / determinate);
	return inverse;
}

static inline
bool pixmM2x2IsFinite(const PixtyM2x2 *pA) {
	return isfinite(pA->d[0][0]) && isfinite(pA->d[0][1]) &&
	       isfinite(pA->d[1][0]) && isfinite(pA->d[1][1]);
}

static inline
float pixmM3x3Determinate(const PixtyM3x3 *pA) {
	float aDet = pA->d[1][1] * pA->d[2][2] - pA->d[2][1] * pA->d[1][2];
	float bDet = pA->d[0][1] * pA->d[2][2] - pA->d[2][1] * pA->d[0][2];
	float cDet = pA->d[0][1] * pA->d[1][2] - pA->d[1][1] * pA->d[0][2];
	return pA->d[0][0] * aDet - pA->d[1][0] * bDet + pA->d[2][0] * cDet;
}

static inline
void pixmM3x3MultiplyEqualScalar(PixtyM3x3 *pA, float b) {
	for (int32_t i = 0; i < 3; ++i) {
		for (int32_t j = 0; j < 3; ++j) {
			pA->d[i][j] *= b;
		}
	}
}

static inline
PixtyM3x3 pixmM3x3FromV3_F32(PixtyV3_F32 a, PixtyV3_F32 b, PixtyV3_F32 c) {
	PixtyM3x3 mat = {0};
	mat.d[0][0] = a.d[0];
	mat.d[0][1] = a.d[1];
	mat.d[0][2] = a.d[2];
	mat.d[1][0] = b.d[0];
	mat.d[1][1] = b.d[1];
	mat.d[1][2] = b.d[2];
	mat.d[2][0] = c.d[0];
	mat.d[2][1] = c.d[1];
	mat.d[2][2] = c.d[2];
	return mat;
}

static inline
PixtyM3x3 pixmM3x3FromM4x4(const PixtyM4x4 *pA) {
	PixtyM3x3 b = {
		pA->d[0][0], pA->d[0][1], pA->d[0][2],
		pA->d[1][0], pA->d[1][1], pA->d[1][2],
		pA->d[2][0], pA->d[2][1], pA->d[2][2]
	};
	return b;
}

static inline
PixtyM3x3 pixmM3x3Invert(const PixtyM3x3 *pA) {
	float determinate = pixmM3x3Determinate(pA);
	PixtyM3x3 inverse = pixmM3x3Adjugate(pA);
	pixmM3x3MultiplyEqualScalar(&inverse, 1.0f / determinate);
	return inverse;
}

static inline
bool pixmM3x3IsFinite(const PixtyM3x3 *pA) {
	return 
		isfinite(pA->d[0][0]) && isfinite(pA->d[0][1]) && isfinite(pA->d[0][2]) &&
		isfinite(pA->d[1][0]) && isfinite(pA->d[1][1]) && isfinite(pA->d[1][2]) &&
		isfinite(pA->d[2][0]) && isfinite(pA->d[2][1]) && isfinite(pA->d[2][2]);
}

static inline
PixtyM2x3 pixmM2x2MultiplyM2x3(PixtyM2x2 a, PixtyM2x3 b) {
	PixtyM2x3 c = {0};
	c.d[0][0] = a.d[0][0] * b.d[0][0] + a.d[0][1] * b.d[1][0];
	c.d[0][1] = a.d[0][0] * b.d[0][1] + a.d[0][1] * b.d[1][1];
	c.d[0][2] = a.d[0][0] * b.d[0][2] + a.d[0][1] * b.d[1][2];
	c.d[1][0] = a.d[1][0] * b.d[0][0] + a.d[1][1] * b.d[1][0];
	c.d[1][1] = a.d[1][0] * b.d[0][1] + a.d[1][1] * b.d[1][1];
	c.d[1][2] = a.d[1][0] * b.d[0][2] + a.d[1][1] * b.d[1][2];
	return c;
}

static inline
float pixmFloor(float a) {
	int32_t aTrunc = (int32_t)a;
	aTrunc -= pixmF32_NOTEQL_INTERN(((float)aTrunc), (a)) && pixmF32_LESS_INTERN((a), (.0f));
	return (float)aTrunc;
}

static inline
PixtyV2_I32 pixmV2F32FloorAssign(PixtyV2_F32 *pA) {
	PixtyV2_I32 c = {0};
	pA->d[0] = pixmFloor((float)pA->d[0]);
	pA->d[1] = pixmFloor((float)pA->d[1]);
	c.d[0] = (int32_t)pA->d[0];
	c.d[1] = (int32_t)pA->d[1];
	return c;
}

static inline
bool pixmV4F32Equal(PixtyV4_F32 a, PixtyV4_F32 b) {
	return
		pixmF32_EQL_INTERN((a.d[0]), (b.d[0])) &&
		pixmF32_EQL_INTERN((a.d[1]), (b.d[1])) &&
		pixmF32_EQL_INTERN((a.d[2]), (b.d[2])) &&
		pixmF32_EQL_INTERN((a.d[3]), (b.d[3]));
}
#define V4EQL ,V4F32Equal,

static inline
bool pixmV4F64Equal(PixtyV4_F64 a, PixtyV4_F64 b) {
	return
		pixmF64_EQL_INTERN((a.d[0]), (b.d[0])) &&
		pixmF64_EQL_INTERN((a.d[1]), (b.d[1])) &&
		pixmF64_EQL_INTERN((a.d[2]), (b.d[2])) &&
		pixmF64_EQL_INTERN((a.d[3]), (b.d[3]));
}
#define V464EQL ,V4F64Equal,

static inline
void pixmV4F32DivideEqualScalar(PixtyV4_F32 *pA, float b) {
	pA->d[0] /= b;
	pA->d[1] /= b;
	pA->d[2] /= b;
	pA->d[3] /= b;
}
#define V4DIVSEQL ,V4F32DivideEqualScalar,

static inline
PixtyV4_F32 pixmV4F32DivideScalar(PixtyV4_F32 a, float b) {
	return (PixtyV4_F32){
		.d = {a.d[0] / b, a.d[1] / b, a.d[2] / b, a.d[3] / b}
	};
}
#define V4DIVS ,V4F32DivideScalar,

#define PIXM_INFIX(a,o,b) pixm##o((a),(b))
#define _(a) PIXM_INFIX(a)
