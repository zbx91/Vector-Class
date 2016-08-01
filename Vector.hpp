#pragma once

#include "Math.h"

// CVector Two Dimensional

class CVector2
{
public:
	float	x, y;

	CVector2() {};
	CVector2(float a, float b) : x(a), y(b) {};

	float			&operator[](const int);
	float			operator[](const int) const;
	CVector2		&operator=(const CVector2 &);

	CVector2		&operator+=(const CVector2 &);
	CVector2		&operator-=(const CVector2 &);
	CVector2		&operator*=(const float);
	CVector2		&operator/=(const float);

	CVector2		operator+(const CVector2 &) const;
	CVector2		operator-(const CVector2 &) const;
	CVector2		operator*(const float) const;
	CVector2		operator/(const float) const;

	bool			operator==(const CVector2 &) const;
	bool			operator!=(const CVector2 &) const;

	void			Rotate(const float);

	float			Length() const;
	float			LengthSqr() const;

	float			DistTo(const CVector2 &) const;
	float			DistToSqr(const CVector2 &) const;
	bool			WithinAARect(const CVector2 &, const CVector2 &) const;

	float			Dot(const CVector2 &) const;

	float			Normalize();

	float			ToAngle() const;

	void			Random(const float, const float);
	void			Negate();
	void			Clear();

	void			CopyToArray(float *) const;
	float*			Base();
	float const*	Base() const;

	bool			IsValid() const;
	bool			IsZero(const float) const;
};

inline float &CVector2::operator[](const int i)
{
	return ((float *)this)[i];
}

inline float CVector2::operator[](const int i) const
{
	return ((float *)this)[i];
}

FORCEINLINE CVector2 &CVector2::operator=(const CVector2 &v)
{
	x = v.x;
	y = v.y;

	return *this;
}

FORCEINLINE CVector2 &CVector2::operator+=(const CVector2 &v)
{
	x += v.x;
	y += v.y;

	return *this;
}

FORCEINLINE CVector2 &CVector2::operator-=(const CVector2 &v)
{
	x -= v.x;
	y -= v.y;

	return *this;
}

FORCEINLINE CVector2 &CVector2::operator*=(const float f)
{
	x *= f;
	y *= f;

	return *this;
}

FORCEINLINE CVector2 &CVector2::operator/=(const float f)
{
	x /= f + M_EPSILON;
	y /= f + M_EPSILON;

	return *this;
}

inline CVector2 CVector2::operator+(const CVector2 &v) const
{
	return CVector2(x + v.x, y + v.y);
}

inline CVector2 CVector2::operator-(const CVector2 &v) const
{
	return CVector2(x - v.x, y - v.y);
}

inline CVector2 CVector2::operator*(const float f) const
{
	return CVector2(x * f, y * f);
}

inline CVector2 CVector2::operator/(const float f) const
{
	return CVector2(x / (f + M_EPSILON),
		y / (f + M_EPSILON));
}

inline bool CVector2::operator==(const CVector2 &v) const
{
	return (v.x == x && v.y == y);
}

inline bool CVector2::operator!=(const CVector2 &v) const
{
	return (v.x != x || v.y != y);
}

inline float CVector2::Length() const
{
	return Sqrt((x * x) + (y * y));
}

FORCEINLINE float CVector2::LengthSqr() const
{
	return (x * x) + (y * y);
}

inline float CVector2::DistTo(const CVector2 &v) const
{
	return (*this - v).Length();
}

FORCEINLINE float CVector2::DistToSqr(const CVector2 &v) const
{
	return (*this - v).LengthSqr();
}

inline bool CVector2::WithinAARect(const CVector2 &min, const CVector2 &max) const
{
	return ((x > min.x) && (x < max.x) &&
		(y > min.y) && (y < max.y));
}

inline void CVector2::Rotate(const float f)
{
	float _x, _y;

	float s, c;

	SinCos(DEG2RAD(f), s, c);

	_x = x;
	_y = y;

	x = (_x * c) - (_y * s);
	y = (_x * s) + (_y * c);
}

FORCEINLINE float CVector2::Dot(const CVector2 &v) const
{
	return (x * v.x) + (y * v.y);
}

inline float CVector2::Normalize()
{
	float l = Length();
	float m = 1.f / (l + M_EPSILON);

	*this *= m;

	return l;
}

inline float CVector2::ToAngle() const
{
	return ATan2(y, x);
}

inline void CVector2::Random(const float min, const float max)
{
	x = min + ((float)rand() / RAND_MAX) * (max - min);
	y = min + ((float)rand() / RAND_MAX) * (max - min);
}

inline void CVector2::Negate()
{
	x = -x;
	y = -y;
}

inline void CVector2::Clear()
{
	x = y = 0.f;
}

inline void CVector2::CopyToArray(float *f) const
{
	f[0] = x;
	f[1] = y;
}

inline float* CVector2::Base()
{
	return (float*)this;
}

inline float const* CVector2::Base() const
{
	return (float const*)this;
}

inline bool CVector2::IsValid() const
{
	return IsFinite(x) && IsFinite(y);
}

inline bool CVector2::IsZero(const float tolerance = 0.01f) const
{
	return (x > -tolerance && x < tolerance &&
		y > -tolerance && y < tolerance);
}

// CVector3 Three-Dimensional (Angle & Vector)

class CVector3
{
public:
	float	x, y, z;

	CVector3() {};
	CVector3(float a, float b, float c) : x(a), y(b), z(c) {};

	float			&operator[](const int);
	float			operator[](const int) const;
	CVector3		&operator=(const CVector3 &);

	CVector3		&operator+=(const CVector3 &);
	CVector3		&operator-=(const CVector3 &);
	CVector3		&operator*=(const float);
	CVector3		&operator/=(const float);

	CVector3		operator+(const CVector3 &) const;
	CVector3		operator-(const CVector3 &) const;
	CVector3		operator*(const float) const;
	CVector3		operator/(const float) const;

	bool			operator==(const CVector3 &) const;
	bool			operator!=(const CVector3 &) const;

	float			Length() const;
	float			LengthSqr() const;
	float			Length2D() const;
	float			Length2DSqr() const;

	float			DistTo(const CVector3 &) const;
	float			DistToSqr(const CVector3 &) const;
	bool			WithinAABox(const CVector3 &, const CVector3 &) const;

	void			Rotate(const CVector3 &);
	void			Rotate2D(const float);

	CVector3		Cross(const CVector3 &) const;
	float			Dot(const CVector3 &) const;

	float			NormalizeVector();
	void			NormalizeAngle();

	CVector3		ToAngle() const;
	CVector3		ToVector() const;
	CVector2		ToCVector2() const;

	CVector3		Forward() const;
	CVector3		Right() const;
	CVector3		Up() const;

	void			Random(const float, const float);
	void			Negate();
	void			Clear();

	void			CopyToArray(float *) const;
	float*			Base();
	float const*	Base() const;


	bool			IsValid() const;
	bool			IsZero(const float) const;
};

inline float &CVector3::operator[](const int i)
{
	return ((float *)this)[i];
}

inline float CVector3::operator[](const int i) const
{
	return ((float *)this)[i];
}

inline CVector3 &CVector3::operator=(const CVector3 &e)
{
	x = e.x;
	y = e.y;
	z = e.z;

	return *this;
}

FORCEINLINE CVector3 &CVector3::operator+=(const CVector3 &e)
{
	x += e.x;
	y += e.y;
	z += e.z;

	return *this;
}

FORCEINLINE CVector3 &CVector3::operator-=(const CVector3 &e)
{
	x -= e.x;
	y -= e.y;
	z -= e.z;

	return *this;
}

FORCEINLINE CVector3 &CVector3::operator*=(const float f)
{
	x *= f;
	y *= f;
	z *= f;

	return *this;
}

FORCEINLINE CVector3 &CVector3::operator/=(const float f)
{
	x /= f + M_EPSILON;
	y /= f + M_EPSILON;
	z /= f + M_EPSILON;

	return *this;
}

inline CVector3 CVector3::operator+(const CVector3 &e) const
{
	return CVector3(x + e.x, y + e.y, z + e.z);
}

inline CVector3 CVector3::operator-(const CVector3 &e) const
{
	return CVector3(x - e.x, y - e.y, z - e.z);
}

inline CVector3 CVector3::operator*(const float f) const
{
	return CVector3(x * f, y * f, z * f);
}

inline CVector3 CVector3::operator/(const float f) const
{
	return CVector3(x / (f + M_EPSILON),
		y / (f + M_EPSILON),
		z / (f + M_EPSILON));
}

inline bool CVector3::operator==(const CVector3 &e) const
{
	return (e.x == x && e.y == y && e.z == z);
}

inline bool CVector3::operator!=(const CVector3 &e) const
{
	return (e.x != x || e.y != y || e.z != z);
}

inline float CVector3::Length() const
{
	return Sqrt((x * x) + (y * y) + (z * z));
}

FORCEINLINE float CVector3::LengthSqr() const
{
	return (x * x) + (y * y) + (z * z);
}

inline float CVector3::Length2D() const
{
	return Sqrt((x * x) + (y * y));
}

FORCEINLINE float CVector3::Length2DSqr() const
{
	return (x * x) + (y * y);
}

inline float CVector3::DistTo(const CVector3 &v) const
{
	return (*this - v).Length();
}

FORCEINLINE float CVector3::DistToSqr(const CVector3 &v) const
{
	return (*this - v).LengthSqr();
}

inline bool CVector3::WithinAABox(const CVector3 &min, const CVector3 &max) const
{
	return ((x > min.x) && (x < max.x) &&
		(y > min.y) && (y < max.y) &&
		(z > min.z) && (z < max.z));
}

inline void CVector3::Rotate(const CVector3 &a)
{
	float _x, _y, _z;

	float s, c;

	SinCos(DEG2RAD(a.x), s, c);

	_y = y;
	_z = z;

	y = (_y * c) - (_z * s);
	z = (_y * s) + (_z * c);

	SinCos(DEG2RAD(a.y), s, c);

	_x = x;
	_z = z;

	x = (_x * c) + (_z * s);
	z = (-_x * s) + (_z * c);

	SinCos(DEG2RAD(a.z), s, c);

	_x = x;
	_y = y;

	x = (_x * c) - (_y * s);
	y = (_x * s) + (_y * c);
}

inline void CVector3::Rotate2D(const float f)
{
	float _x, _y;

	float s, c;

	SinCos(DEG2RAD(f), s, c);

	_x = x;
	_y = y;

	x = (_x * c) - (_y * s);
	y = (_x * s) + (_y * c);
}

inline CVector3 CVector3::Cross(const CVector3 &v) const
{
	return CVector3((y * v.z) - (z * v.y),
		(z * v.x) - (x * v.z),
		(x * v.y) - (y * v.x));
}

FORCEINLINE float CVector3::Dot(const CVector3 &v) const
{
	return (x * v.x) + (y * v.y) + (z * v.z);
}

inline float CVector3::NormalizeVector()
{
	float l = Length();
	float m = 1.f / (l + M_EPSILON);

	*this *= m;

	return l;
}

inline void CVector3::NormalizeAngle()
{
	x = Mod(x + 180.f, 360.f) + 180.f;
	y = Mod(y + 180.f, 360.f) + 180.f;
	z = Mod(z + 180.f, 360.f) + 180.f;
}

inline CVector3 CVector3::ToAngle() const
{
	if (x == 0.f && y == 0.f)
		return CVector3(z > 0.f ? -90.f : 90.f, 0.f, 0.f);

	return CVector3(-RAD2DEG(ASin(z / (Length() + M_EPSILON))), RAD2DEG(ATan2(y, x)), 0.f);
}

inline CVector3 CVector3::ToVector() const
{
	float sx, cx, sy, cy;

	SinCos(DEG2RAD(x), sx, cx);
	SinCos(DEG2RAD(y), sy, cy);

	return CVector3(sx * cy, sx*sy, cx);
}

inline CVector2 CVector3::ToCVector2() const
{
	return CVector2(x, y);
}

inline CVector3 CVector3::Forward() const
{
	float sx, cx, sy, cy;

	SinCos(DEG2RAD(x), sx, cx);
	SinCos(DEG2RAD(y), sy, cy);

	return CVector3(sx * cy, sx * sy, cx);
}

inline CVector3 CVector3::Right() const
{
	float sx, cx, sy, cy, sz, cz;

	SinCos(DEG2RAD(x), sx, cx);
	SinCos(DEG2RAD(y), sy, cy);
	SinCos(DEG2RAD(z), sz, cz);

	return CVector3((-sx * cy * sz) + (cz * sy), (-sx * sy * sz) - (cy * cz), -cx * sz);
}

inline CVector3 CVector3::Up() const
{
	float sx, cx, sy, cy, sz, cz;

	SinCos(DEG2RAD(x), sx, cx);
	SinCos(DEG2RAD(y), sy, cy);
	SinCos(DEG2RAD(z), sz, cz);

	return CVector3((sx * cy * cz) + (sy * sz), (sx * sy * cz) - (cy * sz), cx * cz);
}

inline void CVector3::Random(const float min, const float max)
{
	x = min + ((float)rand() / RAND_MAX) * (max - min);
	y = min + ((float)rand() / RAND_MAX) * (max - min);
	z = min + ((float)rand() / RAND_MAX) * (max - min);
}

inline void CVector3::Negate()
{
	x = -x;
	y = -y;
	z = -z;
}

inline void CVector3::Clear()
{
	x = y = z = 0.f;
}

inline void CVector3::CopyToArray(float *f) const
{
	f[0] = x;
	f[1] = y;
	f[2] = z;
}

inline float* CVector3::Base()
{
	return (float*)this;
}

inline float const* CVector3::Base() const
{
	return (float const*)this;
}

inline bool CVector3::IsValid() const
{
	return IsFinite(x) && IsFinite(y) && IsFinite(z);
}

inline bool CVector3::IsZero(const float tolerance = 0.01f) const
{
	return (x > -tolerance && x < tolerance &&
		y > -tolerance && y < tolerance &&
		z > -tolerance && z < tolerance);
}

// Three-Dimensional Conversion Functions

void VectorToAngle(const CVector3 &v, CVector3 &a)
{
	if (v.x == 0.f && v.y == 0.f)
	{
		if (v.z > 0.f)
			a.x = -90.f;
		else
			a.x = 90.f;

		a.y = 0.f;
	}
	else
	{
		a.x = -RAD2DEG(asin(v.z / (v.Length() + M_EPSILON)));
		a.y = RAD2DEG(atan2(v.y, v.x));
	}

	a.z = 0.f;
}

void AngleToVector(const CVector3 &a, CVector3 *forward, CVector3 *right, CVector3 *up)
{
	float sx, cx, sy, cy, sz, cz;
	SinCos(DEG2RAD(a.x), sx, cx);
	SinCos(DEG2RAD(a.y), sy, cy);
	SinCos(DEG2RAD(a.z), sz, cz);

	if (forward)
	{
		forward->x = sx * cy;
		forward->y = sx * sy;
		forward->z = cx;
	}

	if (right)
	{
		right->x = cx*sy;
		right->y = (sx * sy * sz) + (cy * cz);
		right->z = (sx * sy * cz) - (cy * sz);
	}

	if (up)
	{
		up->x = -sx;
		up->y = cx * sz;
		up->z = cx * cz;
	}
}