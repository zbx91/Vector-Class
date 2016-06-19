#pragma once

#define FORCEINLINE __forceinline

#include "Math.hpp"

// Vector Two Dimensional

class Vector2D
{
public:
	float	x, y;

	Vector2D();
	Vector2D(float, float);

	float	&operator[](char);
	Vector2D	&operator=(const Vector2D &);

	FORCEINLINE Vector2D	&operator+=(const Vector2D &);
	FORCEINLINE Vector2D	&operator-=(const Vector2D &);
	FORCEINLINE Vector2D	&operator*=(const float);
	FORCEINLINE Vector2D	&operator/=(const float);

	Vector2D	operator+(const Vector2D &) const;
	Vector2D	operator-(const Vector2D &) const;
	Vector2D	operator*(const float) const;
	Vector2D	operator/(const float) const;

	bool	operator==(const Vector2D &) const;
	bool	operator!=(const Vector2D &) const;

	void	Rotate(const float);

	float	Length() const;
	FORCEINLINE float	LengthSqr() const;

	float	DistTo(const Vector2D &) const;
	FORCEINLINE float	DistToSqr(const Vector2D &) const;
	bool	WithinAARect(const Vector2D &, const Vector2D &) const;

	FORCEINLINE float	Dot(const Vector2D &) const;

	float	Normalize();

	float	ToAngle() const;

	void	Random(const float, const float);
	void	Negate();
	void	Clear();

	void	CopyToArray(float *) const;
	float	*Base() const;

	bool	IsValid() const;
	bool	IsZero(const float) const;
};

inline Vector2D::Vector2D()
{
	x = y = 0.f;
}

inline Vector2D::Vector2D(float _x, float _y)
{
	x = _x;
	y = _y;
}

inline float &Vector2D::operator[](char i)
{
	return ((float *)this)[i];
}

FORCEINLINE Vector2D &Vector2D::operator=(const Vector2D &v)
{
	x = v.x;
	y = v.y;

	return *this;
}

FORCEINLINE Vector2D &Vector2D::operator+=(const Vector2D &v)
{
	x += v.x;
	y += v.y;

	return *this;
}

FORCEINLINE Vector2D &Vector2D::operator-=(const Vector2D &v)
{
	x -= v.x;
	y -= v.y;

	return *this;
}

FORCEINLINE Vector2D &Vector2D::operator*=(const float f)
{
	x *= f;
	y *= f;

	return *this;
}

FORCEINLINE Vector2D &Vector2D::operator/=(const float f)
{
	x /= f + M_FLT_EPSILON;
	y /= f + M_FLT_EPSILON;

	return *this;
}

inline Vector2D Vector2D::operator+(const Vector2D &v) const
{
	return Vector2D(x + v.x, y + v.y);
}

inline Vector2D Vector2D::operator-(const Vector2D &v) const
{
	return Vector2D(x - v.x, y - v.y);
}

inline Vector2D Vector2D::operator*(const float f) const
{
	return Vector2D(x * f, y * f);
}

inline Vector2D Vector2D::operator/(const float f) const
{
	return Vector2D(x / (f + M_FLT_EPSILON),
		y / (f + M_FLT_EPSILON));
}

inline bool Vector2D::operator==(const Vector2D &v) const
{
	return (v.x == x && v.y == y);
}

inline bool Vector2D::operator!=(const Vector2D &v) const
{
	return (v.x != x || v.y != y);
}

inline float Vector2D::Length() const
{
	return Sqrt((x * x) + (y * y));
}

FORCEINLINE float Vector2D::LengthSqr() const
{
	return (x * x) + (y * y);
}

inline float Vector2D::DistTo(const Vector2D &v) const
{
	return (*this - v).Length();
}

FORCEINLINE float Vector2D::DistToSqr(const Vector2D &v) const
{
	return (*this - v).LengthSqr();
}

inline bool Vector2D::WithinAARect(const Vector2D &min, const Vector2D &max) const
{
	return ((x > min.x) && (x < max.x) &&
		(y > min.y) && (y < max.y));
}

inline void Vector2D::Rotate(const float f)
{
	float _x, _y;

	float s, c;

	SinCos(DEG2RAD(f), s, c);

	_x = x;
	_y = y;

	x = (_x * c) - (_y * s);
	y = (_x * s) + (_y * c);
}

FORCEINLINE float Vector2D::Dot(const Vector2D &v) const
{
	return (x * v.x) + (y * v.y);
}

inline float Vector2D::Normalize()
{
	float l = Length();
	float m = 1.f / (l + M_FLT_EPSILON);

	*this *= m;

	return l;
}

inline float Vector2D::ToAngle() const
{
	return ATan2(y, x);
}

inline void Vector2D::Random(const float min, const float max)
{
	x = min + ((float)rand() / RAND_MAX) * (max - min);
	y = min + ((float)rand() / RAND_MAX) * (max - min);
}

inline void Vector2D::Negate()
{
	x = -x;
	y = -y;
}

inline void Vector2D::Clear()
{
	x = y = 0.f;
}

inline void Vector2D::CopyToArray(float *f) const
{
	f[0] = x;
	f[1] = y;
}

inline float *Vector2D::Base() const
{
	return (float *)this;
}

inline bool Vector2D::IsValid() const
{
	using namespace std;

	return isfinite(x) && isfinite(y);
}

inline bool Vector2D::IsZero(const float tolerance = 0.01f) const
{
	return (x > -tolerance && x < tolerance &&
		y > -tolerance && y < tolerance);
}

// Euler Three-Dimensional (Angle & Vector)

typedef class Euler Angle, Vector;

class Euler
{
public:
	float	x, y, z;

	Euler();
	Euler(float, float, float);

	float		&operator[](char);
	Euler		&operator=(const Euler &);

	FORCEINLINE Euler		&operator+=(const Euler &);
	FORCEINLINE Euler		&operator-=(const Euler &);
	FORCEINLINE Euler		&operator*=(const float);
	FORCEINLINE Euler		&operator/=(const float);

	Euler	operator+(const Euler &) const;
	Euler	operator-(const Euler &) const;
	Euler	operator*(const float) const;
	Euler	operator/(const float) const;

	bool	operator==(const Euler &) const;
	bool	operator!=(const Euler &) const;

	float	Length() const;
	FORCEINLINE float	LengthSqr() const;
	float	Length2D() const;
	FORCEINLINE float	Length2DSqr() const;

	float	DistTo(const Vector &) const;
	FORCEINLINE float	DistToSqr(const Vector &) const;
	bool	WithinAABox(const Vector &, const Vector &) const;

	void	Rotate(const Angle &);
	void	Rotate2D(const float);

	Vector	Cross(const Vector &) const;
	FORCEINLINE float	Dot(const Vector &) const;

	float	NormalizeVector();
	void	NormalizeAngle();

	Angle	ToAngle() const;
	Vector	ToVector() const;
	Vector2D	ToVector2D() const;

	Vector	Forward() const;
	Vector	Right() const;
	Vector	Up() const;

	void	Random(const float, const float);
	void	Negate();
	void	Clear();

	void	CopyToArray(float *) const;
	float	*Base() const;

	bool	IsValid() const;
	bool	IsZero(const float) const;
};

inline Euler::Euler()
{
	x = y = z = 0.f;
}

inline Euler::Euler(float _x, float _y, float _z)
{
	x = _x;
	y = _y;
	z = _z;
}

inline float &Euler::operator[](char i)
{
	return ((float *)this)[i];
}

inline Euler &Euler::operator=(const Euler &e)
{
	x = e.x;
	y = e.y;
	z = e.z;

	return *this;
}

FORCEINLINE Euler &Euler::operator+=(const Euler &e)
{
	x += e.x;
	y += e.y;
	z += e.z;

	return *this;
}

FORCEINLINE Euler &Euler::operator-=(const Euler &e)
{
	x -= e.x;
	y -= e.y;
	z -= e.z;

	return *this;
}

FORCEINLINE Euler &Euler::operator*=(const float f)
{
	x *= f;
	y *= f;
	z *= f;

	return *this;
}

FORCEINLINE Euler &Euler::operator/=(const float f)
{
	x /= f + M_FLT_EPSILON;
	y /= f + M_FLT_EPSILON;
	z /= f + M_FLT_EPSILON;

	return *this;
}

inline Euler Euler::operator+(const Euler &e) const
{
	return Euler(x + e.x, y + e.y, z + e.z);
}

inline Euler Euler::operator-(const Euler &e) const
{
	return Euler(x - e.x, y - e.y, z - e.z);
}

inline Euler Euler::operator*(const float f) const
{
	return Euler(x * f, y * f, z * f);
}

inline Euler Euler::operator/(const float f) const
{
	return Euler(x / (f + M_FLT_EPSILON),
		y / (f + M_FLT_EPSILON),
		z / (f + M_FLT_EPSILON));
}

inline bool Euler::operator==(const Euler &e) const
{
	return (e.x == x && e.y == y && e.z == z);
}

inline bool Euler::operator!=(const Euler &e) const
{
	return (e.x != x || e.y != y || e.z != z);
}

inline float Euler::Length() const
{
	return Sqrt((x * x) + (y * y) + (z * z));
}

FORCEINLINE float Euler::LengthSqr() const
{
	return (x * x) + (y * y) + (z * z);
}

inline float Euler::Length2D() const
{
	return Sqrt((x * x) + (y * y));
}

FORCEINLINE float Euler::Length2DSqr() const
{
	return (x * x) + (y * y);
}

inline float Euler::DistTo(const Vector &v) const
{
	return (*this - v).Length();
}

FORCEINLINE float Euler::DistToSqr(const Vector &v) const
{
	return (*this - v).LengthSqr();
}

inline bool Euler::WithinAABox(const Vector &min, const Vector &max) const
{
	return ((x > min.x) && (x < max.x) &&
			(y > min.y) && (y < max.y) &&
			(z > min.z) && (z < max.z));
}

inline void Euler::Rotate(const Angle &a)
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

inline void Euler::Rotate2D(const float f)
{
	float _x, _y;

	float s, c;

	SinCos(DEG2RAD(f), s, c);

	_x = x;
	_y = y;

	x = (_x * c) - (_y * s);
	y = (_x * s) + (_y * c);
}

inline Vector Euler::Cross(const Vector &v) const
{
	return Vector((y * v.z) - (z * v.y),
		(z * v.x) - (x * v.z),
		(x * v.y) - (y * v.x));
}

FORCEINLINE float Euler::Dot(const Vector &v) const
{
	return (x * v.x) + (y * v.y) + (z * v.z);
}

inline float Euler::NormalizeVector()
{
	float l = Length();
	float m = 1.f / (l + M_FLT_EPSILON);

	*this *= m;

	return l;
}

inline void Euler::NormalizeAngle()
{
	x = _NormalizeAngle(x);
	y = _NormalizeAngle(y);
	z = _NormalizeAngle(z);
}

inline Angle &Euler::ToAngle() const
{
	if (x == 0.f && y == 0.f)
		return Angle(z > 0.f ? -90.f : 90.f, 0.f, 0.f);

	return Angle(-RAD2DEG(ASin(z / (Length() + M_FLT_EPSILON))), RAD2DEG(ATan2(y, x)), 0.f);
}

inline Vector Euler::ToVector() const
{
	float sx, cx, sy, cy;

	SinCos(DEG2RAD(x), sx, cx);
	SinCos(DEG2RAD(y), sy, cy);

	return Vector(sx * cy, sx*sy, cx);
}

inline Vector2D Euler::ToVector2D() const
{
	return Vector2D(x, y);
}

inline Vector Euler::Forward() const
{
	float sx, cx, sy, cy;

	SinCos(DEG2RAD(x), sx, cx);
	SinCos(DEG2RAD(y), sy, cy);

	return Vector(sx * cy, sx * sy, cx);
}

inline Vector Euler::Right() const
{
	float sx, cx, sy, cy, sz, cz;

	SinCos(DEG2RAD(x), sx, cx);
	SinCos(DEG2RAD(y), sy, cy);
	SinCos(DEG2RAD(z), sz, cz);

	return Vector((-sx * cy * sz) + (cz * sy), (-sx * sy * sz) - (cy * cz), -cx * sz);
}

inline Vector Euler::Up() const
{
	float sx, cx, sy, cy, sz, cz;

	SinCos(DEG2RAD(x), sx, cx);
	SinCos(DEG2RAD(y), sy, cy);
	SinCos(DEG2RAD(z), sz, cz);

	return Vector((sx * cy * cz) + (sy * sz), (sx * sy * cz) - (cy * sz), cx * cz);
}

inline void Euler::Random(const float min, const float max)
{
	x = min + ((float)rand() / RAND_MAX) * (max - min);
	y = min + ((float)rand() / RAND_MAX) * (max - min);
	z = min + ((float)rand() / RAND_MAX) * (max - min);
}

inline void Euler::Negate()
{
	x = -x;
	y = -y;
	z = -z;
}

inline void Euler::Clear()
{
	x = y = z = 0.f;
}

inline void Euler::CopyToArray(float *f) const
{
	f[0] = x;
	f[1] = y;
	f[2] = z;
}

inline float *Euler::Base() const
{
	return (float *)this;
}

inline bool Euler::IsValid() const
{
	using namespace std;

	return isfinite(x) && isfinite(y) && isfinite(z);
}

inline bool Euler::IsZero(const float tolerance = 0.01f) const
{
	return (x > -tolerance && x < tolerance &&
		y > -tolerance && y < tolerance &&
		z > -tolerance && z < tolerance);
}

// Three-Dimensional Conversion Functions

inline void VectorToAngle(const Vector &v, Angle &a)
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
		a.x = -RAD2DEG(asin(v.z / (v.Length() + M_FLT_EPSILON)));
		a.y = RAD2DEG(atan2(v.y, v.x));
	}

	a.z = 0.f;
}

inline void AngleToVector(const Angle &a, Vector *forward,  Vector *right, Vector *up)
{
	float sx, cx, sy, cy, sz, cz;
	SinCos(Deg2Rad(a.x), sx, cx);
	SinCos(Deg2Rad(a.y), sy, cy);
	SinCos(Deg2Rad(a.z), sz, cz);

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