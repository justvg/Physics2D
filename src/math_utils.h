#pragma once

#include <math.h>
#include <stdint.h>

#define PI 3.14159265358979323846f

#define Epsilon (1.19e-7f)
#define I32_MIN (-2147483647 - 1)
#define I32_MAX 2147483647
#define uint32_t_MAX 0xFFFFFFFF
#define FLT_MAX 3.402823466e+38F

#define Max(A, B) (((A) > (B)) ? (A) : (B))
#define Min(A, B) (((A) < (B)) ? (A) : (B))

union vec2
{
	struct
	{
		float x, y;
	};

	struct
	{
		float u, v;
	};

	float E[2];

	vec2() { x = y = 0.0f; }
	vec2(float X, float Y) { x = X; y = Y; }
};

union vec3
{
	struct
	{
		float x, y, z;
	};

	struct
	{
		float u, v, w;
	};

	struct
	{
		float r, g, b;
	};

	struct
	{
		vec2 xy;
		float Ignored_0;
	};
	struct
	{
		float Ignored_1;
		vec2 yz;
	};
	struct
	{
		vec2 uv;
		float Ignored_2;
	};
	struct
	{
		float Ignored_3;
		vec2 vw;
	};

	float E[3];

	vec3() {}
	vec3(float X, float Y, float Z) { x = X; y = Y; z = Z; }
	vec3(vec2 XY, float Z) { x = XY.x; y = XY.y; z = Z; }
};

union vec4
{
	struct
	{
		float x, y, z, w;
	};

	struct
	{
		float r, g, b, a;
	};

	struct
	{
		vec3 xyz;
		float Ignored_0;
	};
	struct
	{
		vec2 xy;
		float Ignored_1;
		float Ignored_2;
	};
	struct
	{
		float Ignored_3;
		vec2 yz;
		float Ignored_4;
	};
	struct
	{
		float Ignored_5;
		float Ignored_6;
		vec2 zw;
	};

	struct
	{
		vec3 rgb;
		float Ignored_7;
	};

	float E[4];

	vec4() {}
	vec4(float X, float Y, float Z, float W) { x = X; y = Y; z = Z; w = W; }
};

inline vec2
vec2i(int32_t X, int32_t Y)
{
	vec2 Result = vec2((float)X, (float)Y);

	return(Result);
}

inline vec2
vec2i(uint32_t X, uint32_t Y)
{
	vec2 Result = vec2((float)X, (float)Y);

	return(Result);
}

inline vec3
vec3i(int32_t X, int32_t Y, int32_t Z)
{
	vec3 Result = vec3((float)X, (float)Y, (float)Z);

	return(Result);
}

inline vec3
vec3i(uint32_t X, uint32_t Y, uint32_t Z)
{
	vec3 Result = vec3((float)X, (float)Y, (float)Z);

	return(Result);
}

inline vec4
vec4i(int32_t X, int32_t Y, int32_t Z, int32_t W)
{
	vec4 Result = vec4((float)X, (float)Y, (float)Z, (float)W);

	return(Result);
}

inline vec4
vec4i(uint32_t X, uint32_t Y, uint32_t Z, uint32_t W)
{
	vec4 Result = vec4((float)X, (float)Y, (float)Z, (float)W);

	return(Result);
}

union mat2
{
	float E[4];
	struct
	{
		float a11, a21;
		float a12, a22;
	};
};

union mat3
{
	float E[9];
	struct
	{
		float a11, a21, a31;
		float a12, a22, a32;
		float a13, a23, a33;
	};
};

union mat4
{
	float E[16];

	struct
	{
		float a11, a21, a31, a41;
		float a12, a22, a32, a42;
		float a13, a23, a33, a43;
		float a14, a24, a34, a44;
	};
};

inline float
Square(float Value)
{
	Value = Value * Value;
	return(Value);
}

inline float
Absolute(float Value)
{
	if(Value < 0.0f) Value = -Value;
	return(Value);
}

inline float
Radians(float Angle)
{
	float Result = (Angle / 180.0f) * PI;;

	return(Result);
}

inline float
Degrees(float Rad)
{
	float Result = (Rad / PI) * 180.0f;

	return(Result);
}

inline float
Lerp(float A, float B, float t)
{
	float Result = A + (B - A) * t;

	return(Result);
}

inline float
Clamp(float Value, float Min, float Max)
{
	if(Min > Max)
	{
		float Temp = Min;
		Min = Max;
		Max = Temp;
	}

	if(Value < Min) Value = Min;
	else if(Value > Max) Value = Max;

	return(Value);
}

// 
// NOTE(georgy): vec2
// 

inline vec2
operator+ (vec2 A, vec2 B)
{
	vec2 Result;

	Result.x = A.x + B.x;
	Result.y = A.y + B.y;

	return(Result);
}

inline vec2
operator- (vec2 A, vec2 B)
{
	vec2 Result;

	Result.x = A.x - B.x;
	Result.y = A.y - B.y;

	return(Result);
}

inline vec2
Hadamard(vec2 A, vec2 B)
{
	vec2 Result;

	Result.x = A.x * B.x;
	Result.y = A.y * B.y;

	return(Result);
}

inline vec2
operator* (vec2 A, float B)
{
	vec2 Result;

	Result.x = A.x * B;
	Result.y = A.y * B;

	return(Result);
}

inline vec2
operator* (float B, vec2 A)
{
	vec2 Result = A * B;

	return(Result);
}

inline vec2&
operator+= (vec2& A, vec2 B)
{
	A = A + B;

	return(A);
}

inline vec2&
operator-= (vec2& A, vec2 B)
{
	A = A - B;

	return(A);
}

inline vec2&
operator*= (vec2& A, float B)
{
	A = A * B;

	return(A);
}

inline vec2
operator- (vec2 A)
{
	vec2 Result;

	Result.x = -A.x;
	Result.y = -A.y;

	return(Result);
}

inline float
Dot(vec2 A, vec2 B)
{
	float Result = A.x * B.x + A.y * B.y;

	return(Result);
}

inline float
LengthSq(vec2 A)
{
	float Result = Dot(A, A);

	return(Result);
}

inline float
Length(vec2 A)
{
	float Result = sqrtf(LengthSq(A));

	return(Result);
}

inline vec2
Normalize(vec2 A)
{
	vec2 Result = A * (1.0f / Length(A));

	return(Result);
}

inline vec2
NOZ(vec2 A)
{
	vec2 Result;

	float Len = Length(A);
	if (Len <= Epsilon)
	{
		Result = vec2(0.0f, 0.0f);
	}
	else
	{
		Result = A * (1.0f / Length(A));
	}

	return(Result);
}

inline vec2
Perp(vec2 A)
{
	vec2 Result = vec2(-A.y, A.x);

	return(Result);
}

inline float
Cross2D(vec2 A, vec2 B)
{
	float Result = Dot(Perp(A), B);

	return(Result);
}

inline vec2
Lerp(vec2 A, vec2 B, float t)
{
	vec2 Result = A + (B - A) * t;

	return(Result);
}

inline bool
VectorsAreEqual(vec2 A, vec2 B)
{
	bool XIsEqual = (fabs(A.x - B.x) <= Epsilon);
	bool YIsEqual = (fabs(A.y - B.y) <= Epsilon);

	return(XIsEqual && YIsEqual);
}

// 
// NOTE(georgy): vec3
// 

inline vec3
operator+ (vec3 A, vec3 B)
{
	vec3 Result;

	Result.x = A.x + B.x;
	Result.y = A.y + B.y;
	Result.z = A.z + B.z;

	return(Result);
}

inline vec3
operator- (vec3 A, vec3 B)
{
	vec3 Result;

	Result.x = A.x - B.x;
	Result.y = A.y - B.y;
	Result.z = A.z - B.z;

	return(Result);
}

inline vec3
Hadamard(vec3 A, vec3 B)
{
	vec3 Result;

	Result.x = A.x * B.x;
	Result.y = A.y * B.y;
	Result.z = A.z * B.z;

	return(Result);
}

inline vec3
operator* (vec3 A, float B)
{
	vec3 Result;

	Result.x = A.x * B;
	Result.y = A.y * B;
	Result.z = A.z * B;

	return(Result);
}

inline vec3
operator* (float B, vec3 A)
{
	vec3 Result = A * B;

	return(Result);
}

inline vec3&
operator+= (vec3& A, vec3 B)
{
	A = A + B;

	return(A);
}

inline vec3&
operator-= (vec3& A, vec3 B)
{
	A = A - B;

	return(A);
}

inline vec3&
operator*= (vec3& A, float B)
{
	A = A * B;

	return(A);
}

inline vec3
operator- (vec3 A)
{
	vec3 Result;

	Result.x = -A.x;
	Result.y = -A.y;
	Result.z = -A.z;

	return(Result);
}

inline float
Dot(vec3 A, vec3 B)
{
	float Result = A.x * B.x + A.y * B.y + A.z * B.z;

	return(Result);
}

inline float
LengthSq(vec3 A)
{
	float Result = Dot(A, A);

	return(Result);
}

inline float
Length(vec3 A)
{
	float Result = sqrtf(LengthSq(A));

	return(Result);
}

inline vec3
Normalize(vec3 A)
{
	vec3 Result = A * (1.0f / Length(A));

	return(Result);
}

inline vec3
NOZ(vec3 A)
{
	vec3 Result;

	float Len = Length(A);
	if (Len <= Epsilon)
	{
		Result = vec3(0.0f, 0.0f, 0.0f);
	}
	else
	{
		Result = A * (1.0f / Length(A));
	}

	return(Result);
}

inline vec3
Cross(vec3 A, vec3 B)
{
	vec3 Result;

	Result.x = A.y * B.z - B.y * A.z;
	Result.y = A.z * B.x - B.z * A.x;
	Result.z = A.x * B.y - B.x * A.y;

	return(Result);
}

inline vec3
Lerp(vec3 A, vec3 B, float t)
{
	vec3 Result = A + (B - A) * t;

	return(Result);
}

static bool
PointInTriangle(vec3 P, vec3 A, vec3 B, vec3 C)
{
	A -= P; B -= P; C -= P;
	vec3 U = Cross(B, C);
	vec3 V = Cross(C, A);
	vec3 W = Cross(A, B);
	bool Result = (Dot(U, V) >= 0.0f) && (Dot(U, W) >= 0.0f);

	return(Result);
}

static vec3
ClosestPointInSegment(vec3 P, vec3 A, vec3 B)
{
	float t = Dot(P - A, (B - A) * (1.0f / Length(B - A))) / Length(B - A);
	vec3 Result = A + t * (B - A);

	return(Result);
}

enum voronoi_region
{
	VoronoiRegion_Vertex,
	VoronoiRegion_Edge,
	VoronoiRegion_Triangle,
};

struct closest_voronoi_region_triangle_point
{
	voronoi_region VR;
	vec3 P;

	// NOTE(georgy): These are valid if VR = VoronoiRegion_Edge
	vec3 EdgeP0;
	vec3 EdgeP1;

	// NOTE(georgy): These are valid if VR = VoronoiRegion_Vertex or VR = VoronoiRegion_Edge
	uint32_t VertIndex0;
	uint32_t VertIndex1;
};

static closest_voronoi_region_triangle_point
ClosestPointInTriangleVR(vec3 P, vec3 A, vec3 B, vec3 C)
{
	closest_voronoi_region_triangle_point Result;

	vec3 AB = B - A;
	vec3 AC = C - A;
	vec3 BC = C - B;

	// NOTE(georgy): Parametric position s for projection P' of P on AB,
	// P' = A + s*AB, s = sNom / (sNom + sDenom)
	float sNom = Dot(P - A, AB);
	float sDenom = Dot(P - B, A - B);

	// NOTE(georgy): Parametric position t for projection P' of P on AC,
	// P' = A + t*AC, s = tNom / (tNom + tDenom)
	float tNom = Dot(P - A, AC);
	float tDenom = Dot(P - C, A - C);

	if ((sNom <= 0.0f) && (tNom <= 0.0f))
	{
		Result.VR = VoronoiRegion_Vertex;
		Result.P = A;
		Result.VertIndex0 = 0;
		return(Result);
	}

	// NOTE(georgy): Parametric position u for projection P' of P on BC,
	// P' = B + u*BC, u = uNom / (uNom + uDenom)
	float uNom = Dot(P - B, BC);
	float uDenom = Dot(P - C, B - C);

	if ((sDenom <= 0.0f) && (uNom <= 0.0f))
	{
		Result.VR = VoronoiRegion_Vertex;
		Result.P = B;
		Result.VertIndex0 = 1;
		return(Result);
	}

	if ((tDenom <= 0.0f) && (uDenom <= 0.0f))
	{
		Result.VR = VoronoiRegion_Vertex;
		Result.P = C;
		Result.VertIndex0 = 2;
		return(Result);
	}

	vec3 N = Cross(B - A, C - A);
	float RAB = Dot(N, Cross(A - P, B - P));
	if ((RAB <= 0.0f) && (sNom >= 0.0f) && (sDenom >= 0.0f))
	{
		Result.VR = VoronoiRegion_Edge;
		Result.P = A + AB * (sNom / (sNom + sDenom));
		Result.EdgeP0 = A;
		Result.EdgeP1 = B;
		Result.VertIndex0 = 0;
		Result.VertIndex1 = 1;
		return(Result);
	}

	float RBC = Dot(N, Cross(B - P, C - P));
	if ((RBC <= 0.0f) && (uNom >= 0.0f) && (uDenom >= 0.0f))
	{
		Result.VR = VoronoiRegion_Edge;
		Result.P = B + BC * (uNom / (uNom + uDenom));
		Result.EdgeP0 = B;
		Result.EdgeP1 = C;
		Result.VertIndex0 = 1;
		Result.VertIndex1 = 2;
		return(Result);
	}

	float RCA = Dot(N, Cross(C - P, A - P));
	if ((RCA <= 0.0f) && (tNom >= 0.0f) && (tDenom >= 0.0f))
	{
		Result.VR = VoronoiRegion_Edge;
		Result.P = A + AC * (tNom / (tNom + tDenom));
		Result.EdgeP0 = C;
		Result.EdgeP1 = A;
		Result.VertIndex0 = 2;
		Result.VertIndex1 = 0;
		return(Result);
	}

	float u = RBC / (RAB + RBC + RCA);
	float v = RCA / (RAB + RBC + RCA);
	float w = 1.0f - u - v;
	Result.VR = VoronoiRegion_Triangle;
	Result.P = u * A + v * B + w * C;
	return(Result);

#if 0
	// NOTE(georgy): Slow but easy to understand way to do this

	b32 NegativeHalfspaceOfABPlane = (Dot(P - A, B - A) <= 0.0f);
	b32 NegativeHalfspaceOfACPlane = (Dot(P - A, C - A) <= 0.0f);
	if (NegativeHalfspaceOfABPlane && NegativeHalfspaceOfACPlane)
	{
		Result.VR = VoronoiRegion_Vertex;
		Result.P = A;
	}
	else
	{
		b32 NegativeHalfspaceOfBAPlane = (Dot(P - B, A - B) <= 0.0f);
		b32 NegativeHalfspaceOfBCPlane = (Dot(P - B, C - B) <= 0.0f);
		if (NegativeHalfspaceOfBAPlane && NegativeHalfspaceOfBCPlane)
		{
			Result.VR = VoronoiRegion_Vertex;
			Result.P = B;
		}
		else
		{
			b32 NegativeHalfspaceOfCAPlane = (Dot(P - C, A - C) <= 0.0f);
			b32 NegativeHalfspaceOfCBPlane = (Dot(P - C, B - C) <= 0.0f);
			if (NegativeHalfspaceOfCAPlane && NegativeHalfspaceOfCBPlane)
			{
				Result.VR = VoronoiRegion_Vertex;
				Result.P = C;
			}
			else
			{
				vec3 N = Cross(B - A, C - A);

				float RAB = Dot(N, Cross(A - P, B - P));
				if ((RAB <= 0.0f) && !NegativeHalfspaceOfABPlane && !NegativeHalfspaceOfBAPlane)
				{
					Result.VR = VoronoiRegion_Edge;
					Result.P = ClosestPointInSegment(P, A, B);
					Result.EdgeP0 = A;
					Result.EdgeP1 = B;
				}
				else
				{
					float RBC = Dot(N, Cross(B - P, C - P));
					if ((RBC <= 0.0f) && !NegativeHalfspaceOfBCPlane && !NegativeHalfspaceOfCBPlane)
					{
						Result.VR = VoronoiRegion_Edge;
						Result.P = ClosestPointInSegment(P, B, C);
						Result.EdgeP0 = B;
						Result.EdgeP1 = C;
					}
					else
					{
						float RCA = Dot(N, Cross(C - P, A - P));
						if ((RCA <= 0.0f) && !NegativeHalfspaceOfCAPlane && !NegativeHalfspaceOfACPlane)
						{
							Result.VR = VoronoiRegion_Edge;
							Result.P = ClosestPointInSegment(P, C, A);
							Result.EdgeP0 = C;
							Result.EdgeP1 = A;
						}
						else
						{
							Result.VR = VoronoiRegion_Triangle;
							Result.P = P - N * (Dot(P - A, N) / Dot(N, N));
						}
					}
				}
			}
		}
	}
#endif
	return(Result);
}

// 
// NOTE(georgy): vec4
// 

inline vec4
operator+ (vec4 A, vec4 B)
{
	vec4 Result;

	Result.x = A.x + B.x;
	Result.y = A.y + B.y;
	Result.z = A.z + B.z;
	Result.w = A.w + B.w;

	return(Result);
}

inline vec4
operator- (vec4 A, vec4 B)
{
	vec4 Result;

	Result.x = A.x - B.x;
	Result.y = A.y - B.y;
	Result.z = A.z - B.z;
	Result.w = A.w - B.w;

	return(Result);
}

inline vec4
Hadamard(vec4 A, vec4 B)
{
	vec4 Result;

	Result.x = A.x * B.x;
	Result.y = A.y * B.y;
	Result.z = A.z * B.z;
	Result.w = A.w * B.w;

	return(Result);
}

inline vec4
operator* (vec4 A, float B)
{
	vec4 Result;

	Result.x = A.x * B;
	Result.y = A.y * B;
	Result.z = A.z * B;
	Result.w = A.w * B;

	return(Result);
}

inline vec4
operator* (float B, vec4 A)
{
	vec4 Result = A * B;

	return(Result);
}

inline vec4&
operator+= (vec4& A, vec4 B)
{
	A = A + B;

	return(A);
}

inline vec4&
operator-= (vec4& A, vec4 B)
{
	A = A - B;

	return(A);
}

inline vec4&
operator*= (vec4& A, float B)
{
	A = A * B;

	return(A);
}

inline vec4
operator- (vec4 A)
{
	vec4 Result;

	Result.x = -A.x;
	Result.y = -A.y;
	Result.z = -A.z;
	Result.w = -A.w;

	return(Result);
}

inline float
Dot(vec4 A, vec4 B)
{
	float Result = A.x * B.x + A.y * B.y + A.z * B.z + A.w * B.w;

	return(Result);
}

inline float
LengthSq(vec4 A)
{
	float Result = Dot(A, A);

	return(Result);
}

inline float
Length(vec4 A)
{
	float Result = sqrtf(LengthSq(A));

	return(Result);
}

inline vec4
Normalize(vec4 A)
{
	vec4 Result = A * (1.0f / Length(A));

	return(Result);
}

inline vec4
Lerp(vec4 A, vec4 B, float t)
{
	vec4 Result = A + (B - A) * t;

	return(Result);
}

// 
// NOTE(georgy): mat2
// 

static mat2
Rotation2x2(float Angle)
{
	mat2 Result;

	float Rad = Radians(Angle);
	vec2 XAxis = vec2(cosf(Rad), sinf(Rad));
	vec2 YAxis = Perp(XAxis);

	Result.a11 = XAxis.x;
	Result.a21 = XAxis.y;
	Result.a12 = YAxis.x;
	Result.a22 = YAxis.y;

	return(Result);
}

inline vec2
operator*(mat2 M, vec2 V)
{
	vec2 Result;

	Result.x = M.a11*V.x + M.a12*V.y;
	Result.y = M.a21*V.x + M.a22*V.y;

	return(Result);
}

// 
// NOTE(georgy): mat3
// 

inline mat3
Identity3x3(float Diagonal = 1.0f)
{
	mat3 Result = {};

	Result.a11 = Diagonal;
	Result.a22 = Diagonal;
	Result.a33 = Diagonal;

	return(Result);
}

inline mat3
Scaling3x3(float Scale)
{
	mat3 Result = {};

	Result.a11 = Scale;
	Result.a22 = Scale;
	Result.a33 = Scale;

	return(Result);
}

inline mat3
Scaling3x3(vec3 Scale)
{
	mat3 Result = {};

	Result.a11 = Scale.x;
	Result.a22 = Scale.y;
	Result.a33 = Scale.z;

	return(Result);
}

static mat3
Rotation3x3(float Angle, vec3 Axis)
{
	mat3 Result;

	float Rad = Radians(Angle);
	float Sine = sinf(Rad);
	float Cosine = cosf(Rad);

	Axis = NOZ(Axis);
	if ((Axis.x == 0.0f) && (Axis.y == 0.0f) && (Axis.z == 0.0f))
	{
		Result = Identity3x3();
	}
	else
	{
		Result.a11 = Axis.x * Axis.x * (1.0f - Cosine) + Cosine;
		Result.a21 = Axis.x * Axis.y * (1.0f - Cosine) + Axis.z * Sine;
		Result.a31 = Axis.x * Axis.z * (1.0f - Cosine) - Axis.y * Sine;

		Result.a12 = Axis.x * Axis.y * (1.0f - Cosine) - Axis.z * Sine;
		Result.a22 = Axis.y * Axis.y * (1.0f - Cosine) + Cosine;
		Result.a32 = Axis.y * Axis.z * (1.0f - Cosine) + Axis.x * Sine;

		Result.a13 = Axis.x * Axis.z * (1.0f - Cosine) + Axis.y * Sine;
		Result.a23 = Axis.y * Axis.z * (1.0f - Cosine) - Axis.x * Sine;
		Result.a33 = Axis.z * Axis.z * (1.0f - Cosine) + Cosine;
	}

	return(Result);
}

static mat3
Transpose3x3(const mat3& M)
{
	mat3 Result;

	Result.a11 = M.a11;
	Result.a21 = M.a12;
	Result.a31 = M.a13;

	Result.a12 = M.a21;
	Result.a22 = M.a22;
	Result.a32 = M.a23;

	Result.a13 = M.a31;
	Result.a23 = M.a32;
	Result.a33 = M.a33;

	return(Result);
}

static mat3
Inverse3x3(const mat3& M)
{
	mat3 InverseMatrix = {};

	float Determinant = M.a11 * M.a22 * M.a33 + M.a12 * M.a23 * M.a31 + M.a13 * M.a21 * M.a32 -
		(M.a31 * M.a22 * M.a13 + M.a32 * M.a23 * M.a11 + M.a33 * M.a21 * M.a12);
	if (Determinant > Epsilon)
	{
		float OneOverDeterminant = 1.0f / Determinant;

		InverseMatrix.a11 = (M.a22 * M.a33 - M.a32 * M.a23) * OneOverDeterminant;
		InverseMatrix.a12 = (-(M.a21 * M.a33 - M.a31 * M.a23)) * OneOverDeterminant;
		InverseMatrix.a13 = (M.a21 * M.a32 - M.a31 * M.a22) * OneOverDeterminant;
		InverseMatrix.a21 = (-(M.a12 * M.a33 - M.a32 * M.a13)) * OneOverDeterminant;
		InverseMatrix.a22 = (M.a11 * M.a33 - M.a31 * M.a13) * OneOverDeterminant;
		InverseMatrix.a23 = (-(M.a11 * M.a32 - M.a31 * M.a12)) * OneOverDeterminant;
		InverseMatrix.a31 = (M.a12 * M.a23 - M.a22 * M.a13) * OneOverDeterminant;
		InverseMatrix.a32 = (-(M.a11 * M.a23 - M.a21 * M.a13)) * OneOverDeterminant;
		InverseMatrix.a33 = (M.a11 * M.a22 - M.a21 * M.a12) * OneOverDeterminant;

		InverseMatrix = Transpose3x3(InverseMatrix);
	}

	return(InverseMatrix);
}

static mat3
operator*(mat3 A, mat3 B)
{
	mat3 Result;

	for (uint32_t Row = 0;
		Row < 3;
		Row++)
	{
		for (uint32_t Column = 0;
			Column < 3;
			Column++)
		{
			float Sum = 0.0f;
			for (uint32_t E = 0;
				E < 3;
				E++)
			{
				Sum += A.E[Row + E * 3] * B.E[Column * 3 + E];
			}
			Result.E[Row + Column * 3] = Sum;
		}
	}

	return(Result);
}

static vec3
operator*(mat3 A, vec3 B)
{
	vec3 Result;

	Result.x = A.a11 * B.x + A.a12 * B.y + A.a13 * B.z;
	Result.y = A.a21 * B.x + A.a22 * B.y + A.a23 * B.z;
	Result.z = A.a31 * B.x + A.a32 * B.y + A.a33 * B.z;

	return(Result);
}

// 
// NOTE(georgy): mat4
// 

inline mat4
Identity(float Diagonal = 1.0f)
{
	mat4 Result = {};

	Result.a11 = Diagonal;
	Result.a22 = Diagonal;
	Result.a33 = Diagonal;
	Result.a44 = Diagonal;

	return(Result);
}

inline mat4
Translation(vec3 Translate)
{
	mat4 Result = Identity(1.0f);

	Result.a14 = Translate.x;
	Result.a24 = Translate.y;
	Result.a34 = Translate.z;

	return(Result);
}

inline mat4
Scaling(float Scale)
{
	mat4 Result = {};

	Result.a11 = Scale;
	Result.a22 = Scale;
	Result.a33 = Scale;
	Result.a44 = 1.0f;

	return(Result);
}

inline mat4
Scaling(vec3 Scale)
{
	mat4 Result = {};

	Result.a11 = Scale.x;
	Result.a22 = Scale.y;
	Result.a33 = Scale.z;
	Result.a44 = 1.0f;

	return(Result);
}

static mat4
Rotation(float Angle, vec3 Axis)
{
	mat4 Result;

	float Rad = Radians(Angle);
	float Sine = sinf(Rad);
	float Cosine = cosf(Rad);

	Axis = Normalize(Axis);

	Result.a11 = Axis.x * Axis.x * (1.0f - Cosine) + Cosine;
	Result.a21 = Axis.x * Axis.y * (1.0f - Cosine) + Axis.z * Sine;
	Result.a31 = Axis.x * Axis.z * (1.0f - Cosine) - Axis.y * Sine;
	Result.a41 = 0.0f;

	Result.a12 = Axis.x * Axis.y * (1.0f - Cosine) - Axis.z * Sine;
	Result.a22 = Axis.y * Axis.y * (1.0f - Cosine) + Cosine;
	Result.a32 = Axis.y * Axis.z * (1.0f - Cosine) + Axis.x * Sine;
	Result.a42 = 0.0f;

	Result.a13 = Axis.x * Axis.z * (1.0f - Cosine) + Axis.y * Sine;
	Result.a23 = Axis.y * Axis.z * (1.0f - Cosine) - Axis.x * Sine;
	Result.a33 = Axis.z * Axis.z * (1.0f - Cosine) + Cosine;
	Result.a43 = 0.0f;

	Result.a14 = 0.0f;
	Result.a24 = 0.0f;
	Result.a34 = 0.0f;
	Result.a44 = 1.0f;

	return(Result);
}

static mat4
Mat4(const mat3& M)
{
	mat4 Result;

	Result.a11 = M.a11;
	Result.a21 = M.a21;
	Result.a31 = M.a31;
	Result.a12 = M.a12;
	Result.a22 = M.a22;
	Result.a32 = M.a32;
	Result.a13 = M.a13;
	Result.a23 = M.a23;
	Result.a33 = M.a33;
	Result.a41 = Result.a42 = Result.a43 = Result.a14 = Result.a24 = Result.a34 = 0.0f;
	Result.a44 = 1.0f;

	return(Result);
}

static mat4
LookAt(vec3 From, vec3 Target, vec3 UpAxis = vec3(0.0f, 1.0f, 0.0f))
{
	vec3 Forward = Normalize(From - Target);
	vec3 Right = Normalize(Cross(UpAxis, Forward));
	vec3 Up = Cross(Forward, Right);

	mat4 Result;

	Result.a11 = Right.x;
	Result.a21 = Up.x;
	Result.a31 = Forward.x;
	Result.a41 = 0.0f;

	Result.a12 = Right.y;
	Result.a22 = Up.y;
	Result.a32 = Forward.y;
	Result.a42 = 0.0f;

	Result.a13 = Right.z;
	Result.a23 = Up.z;
	Result.a33 = Forward.z;
	Result.a43 = 0.0f;

	Result.a14 = -Dot(Right, From);
	Result.a24 = -Dot(Up, From);
	Result.a34 = -Dot(Forward, From);
	Result.a44 = 1.0f;

	return(Result);
}

static mat4
ViewRotationMatrixFromDirection(vec3 Dir, vec3 UpAxis = vec3(0.0f, 1.0f, 0.0f))
{
	vec3 Forward = Normalize(-Dir);
	vec3 Right = Normalize(Cross(UpAxis, Forward));
	vec3 Up = Cross(Forward, Right);

	mat4 Result;

	Result.a11 = Right.x;
	Result.a21 = Up.x;
	Result.a31 = Forward.x;
	Result.a41 = 0.0f;

	Result.a12 = Right.y;
	Result.a22 = Up.y;
	Result.a32 = Forward.y;
	Result.a42 = 0.0f;

	Result.a13 = Right.z;
	Result.a23 = Up.z;
	Result.a33 = Forward.z;
	Result.a43 = 0.0f;

	Result.a14 = 0.0f;
	Result.a24 = 0.0f;
	Result.a34 = 0.0f;
	Result.a44 = 1.0f;

	return(Result);
}

static mat4
Ortho(float Bottom, float Top, float Left, float Right, float Near, float Far)
{
	mat4 Result = {};

	Result.a11 = 2.0f / (Right - Left);
	Result.a22 = 2.0f / (Top - Bottom);
	Result.a33 = -2.0f / (Far - Near);
	Result.a14 = -(Right + Left) / (Right - Left);
	Result.a24 = -(Top + Bottom) / (Top - Bottom);
	Result.a34 = -(Far + Near) / (Far - Near);
	Result.a44 = 1.0f;

	return(Result);
}

static mat4
Perspective(float FoV, float AspectRatio, float Near, float Far)
{
	float Scale = tanf(Radians(FoV) * 0.5f) * Near;
	float Top = Scale;
	float Bottom = -Top;
	float Right = Scale * AspectRatio;
	float Left = -Right;

	mat4 Result = {};

	Result.a11 = 2.0f * Near / (Right - Left);
	Result.a22 = 2.0f * Near / (Top - Bottom);
	Result.a13 = (Right + Left) / (Right - Left);
	Result.a23 = (Top + Bottom) / (Top - Bottom);
	Result.a33 = -(Far + Near) / (Far - Near);
	Result.a43 = -1.0f;
	Result.a34 = -(2.0f * Far * Near) / (Far - Near);

	return(Result);
}

static mat4
operator*(mat4 A, mat4 B)
{
	mat4 Result;

	for (uint32_t Row = 0;
		Row < 4;
		Row++)
	{
		for (uint32_t Column = 0;
			Column < 4;
			Column++)
		{
			float Sum = 0.0f;
			for (uint32_t E = 0;
				E < 4;
				E++)
			{
				Sum += A.E[Row + E * 4] * B.E[Column * 4 + E];
			}
			Result.E[Row + Column * 4] = Sum;
		}
	}

	return(Result);
}
