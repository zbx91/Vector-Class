
#include <iostream>

// Constants

#define M_PI = 3.14159265f
#define M_FLT_EPSILON FLT_MIN

// Angle Conversion Functions

inline float Rad2Deg(const float &f)
{
	return f * (180.f / M_PI);
}

inline float Deg2Rad(const float &f)
{
	return f * (M_PI / 180.f);
}

// Trigonometry

inline void SinCos(const float &f, float &s, float &c)
{
	s = (float)sin(f);
	c = (float)cos(f);
}

inline float ATan2(const float &y, const float &x)
{
	return (float)atan2(y, x);
}

// Rounding Functions

inline float Floor(const float &f)
{
	__m128 mm1 = _mm_set_ss(f);

	__m128 mm2 = _mm_round_ps(mm1, _MM_FROUND_FLOOR);

	return _mm_cvtss_f32(mm2);
}

// Sqrt Function

inline float Sqrt(const float &f)
{
	__m128 mm1 = _mm_set_ss(f);

	__m128 mm2 = _mm_sqrt_ss(mm1);

	return _mm_cvtss_f32(mm2);
}

// Angle Normalization

inline void _NormalizeAngle(float &a)
{
	if (!(a > 180.f || a < -180.f))
		return;

	float r = Floor((-a + 180.f) / 360.f);

	r *= 360;

	a += r;
}