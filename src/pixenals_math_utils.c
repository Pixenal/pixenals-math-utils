/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#include <pixenals_math_utils.h>

typedef PixtyV2_I8 V2_I8;
typedef PixtyV2_I16 V2_I16;
typedef PixtyV2_I32 V2_I32;
typedef PixtyV2_I64 V2_I64;
typedef PixtyV2_F32 V2_F32;
typedef PixtyV2_F64 V2_F64;
typedef PixtyV3_I8 V3_I8;
typedef PixtyV3_I16 V3_I16;
typedef PixtyV3_I32 V3_I32;
typedef PixtyV3_I64 V3_I64;
typedef PixtyV3_F32 V3_F32;
typedef PixtyV3_F64 V3_F64;
typedef PixtyV4_I8 V4_I8;
typedef PixtyV4_I16 V4_I16;
typedef PixtyV4_I32 V4_I32;
typedef PixtyV4_I64 V4_I64;
typedef PixtyV4_F32 V4_F32;
typedef PixtyV4_F64 V4_F64;
typedef PixtyM4x4 M4x4;
typedef PixtyM3x3 M3x3;
typedef PixtyM2x2 M2x2;
typedef PixtyM2x3 M2x3;

typedef int8_t I8;
typedef int16_t I16;
typedef int32_t I32;
typedef int64_t I64;

typedef uint8_t U8;
typedef uint16_t U16;
typedef uint32_t U32;
typedef uint64_t U64;

typedef float F32;
typedef double F64;

M3x3 pixmM3x3Adjugate(const M3x3 *pA) {
	M3x3 c = {0};
	c.d[0][0] = pA->d[1][1] * pA->d[2][2] - pA->d[2][1] * pA->d[1][2];
	c.d[0][1] = pA->d[0][1] * pA->d[2][2] - pA->d[2][1] * pA->d[0][2];
	c.d[0][2] = pA->d[0][1] * pA->d[1][2] - pA->d[1][1] * pA->d[0][2];
	c.d[1][0] = pA->d[1][0] * pA->d[2][2] - pA->d[2][0] * pA->d[1][2];
	c.d[1][1] = pA->d[0][0] * pA->d[2][2] - pA->d[2][0] * pA->d[0][2];
	c.d[1][2] = pA->d[0][0] * pA->d[1][2] - pA->d[1][0] * pA->d[0][2];
	c.d[2][0] = pA->d[1][0] * pA->d[2][1] - pA->d[2][0] * pA->d[1][1];
	c.d[2][1] = pA->d[0][0] * pA->d[2][1] - pA->d[2][0] * pA->d[0][1];
	c.d[2][2] = pA->d[0][0] * pA->d[1][1] - pA->d[1][0] * pA->d[0][1];
	c.d[1][0] *= -1.0f;
	c.d[0][1] *= -1.0f;
	c.d[2][1] *= -1.0f;
	c.d[1][2] *= -1.0f;
	return c;
}

V3_F32 pixmBarycentricToCartesian(const V3_F32 *pTri, V3_F32 point) {
	V3_F32 pointCartesian = {0};
	pointCartesian.d[0] =
		(point.d[0] * pTri[0].d[0]) +
		(point.d[1] * pTri[1].d[0]) +
		(point.d[2] * pTri[2].d[0]);
	pointCartesian.d[1] =
		(point.d[0] * pTri[0].d[1]) +
		(point.d[1] * pTri[1].d[1]) +
		(point.d[2] * pTri[2].d[1]);
	pointCartesian.d[2] =
		(point.d[0] * pTri[0].d[2]) +
		(point.d[1] * pTri[1].d[2]) +
		(point.d[2] * pTri[2].d[2]);
	return pointCartesian;
}

V3_F32 pixmCartesianToBarycentric(
	const V3_F32 *pTri32,
	const V3_F32 *pPoint32,
	const V3_F32 *pNormal
) {
	F64 derta = .0;
	F64 dertau = .0;
	F64 dertav = .0;

	//Convert to F64
	V2_F64 point = {.d = {pPoint32->d[0], pPoint32->d[1]}};
	V2_F64 tri[3] = {0};
	for (I32 i = 0; i < 3; ++i) {
		if ((pNormal->d[2] == 1.0) || (pNormal->d[2] == -1.0)) {
			tri[i] = (V2_F64){.d = {(F64)pTri32[i].d[0], (F64)pTri32[i].d[1]}};
			point = (V2_F64){.d = {(F64)pPoint32->d[0], (F64)pPoint32->d[1]}};
		}
		else if ((pNormal->d[1] == 1.0) || (pNormal->d[1] == -1.0) || (pNormal->d[0] == .0f)) {
			tri[i] = (V2_F64){.d = {(F64)pTri32[i].d[0], (F64)pTri32[i].d[2]}};
			point = (V2_F64){.d = {(F64)pPoint32->d[0], (F64)pPoint32->d[2]}};
		}
		else {
			tri[i] = (V2_F64){.d = {(F64)pTri32[i].d[1], (F64)pTri32[i].d[2]}};
			point = (V2_F64){.d = {(F64)pPoint32->d[1], (F64)pPoint32->d[2]}};
		}
	}
	//Perform cramers rule
	derta = (tri[0].d[0] * tri[1].d[1]) - (tri[0].d[0] * tri[2].d[1]) -
	        (tri[1].d[0] * tri[0].d[1]) + (tri[1].d[0] * tri[2].d[1]) +
	        (tri[2].d[0] * tri[0].d[1]) - (tri[2].d[0] * tri[1].d[1]);
	//Get determinate of Au
	dertau = (point.d[0] * tri[1].d[1]) - (point.d[0] * tri[2].d[1]) -
	         (tri[1].d[0] * point.d[1]) + (tri[1].d[0] * tri[2].d[1]) +
	         (tri[2].d[0] * point.d[1]) - (tri[2].d[0] * tri[1].d[1]);
	//Get determinate of Av
	dertav = (tri[0].d[0] * point.d[1]) - (tri[0].d[0] * tri[2].d[1]) -
	         (point.d[0] * tri[0].d[1]) + (point.d[0] * tri[2].d[1]) +
	         (tri[2].d[0] * tri[0].d[1]) - (tri[2].d[0] * point.d[1]);
	
	V3_F32 pointBc = {0};
	//u = dert(Au) / dert(A)
	pointBc.d[0] = (F32)(dertau / derta);
	//u = dert(Av) / dert(A)
	pointBc.d[1] = (F32)(dertav / derta);
	//w can be derived from u and v
	pointBc.d[2] = 1.0f - pointBc.d[0] - pointBc.d[1];

	return pointBc;
}
