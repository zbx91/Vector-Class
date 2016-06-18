#pragma once

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

//Rounding Functions

inline float Floor(const float &f)
{
	__m128 mm1 = _mm_set_ss(f);

	__m128 mm2 = _mm_round_ps(mm1, _MM_FROUND_FLOOR);

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

// Vector Two Dimensional

class Vector2D
{
public:
	float	x, y;

	Vector2D();
	Vector2D(float, float);

	float		&operator[](char);
	Vector2D	&operator=(const Vector2D &);

	Vector2D	&operator+=(const Vector2D &);
	Vector2D	&operator-=(const Vector2D &);
	Vector2D	&operator*=(const float &);
	Vector2D	&operator/=(const float &);

	Vector2D	operator+(const Vector2D &) const;
	Vector2D	operator-(const Vector2D &) const;
	Vector2D	operator*(const float &) const;
	Vector2D	operator/(const float &) const;

	bool		operator==(const Vector2D &) const;
	bool		operator!=(const Vector2D &) const;

	void		Rotate(const float &);

	float		Length() const;
	float		LengthSqr() const;

	float		DistTo(const Vector2D &) const;
	float		DistToSqr(const Vector2D &) const;
	bool		WithinAARect(const Vector2D &, const Vector2D &) const;

	float		Dot(const Vector2D &) const;

	float		Normalize();

	float		ToAngle() const;

	void		Random(const float &, const float &);
	void		Negate();
	void		Clear();

	void		CopyToArray(float *) const;
	float		*Base() const;

	bool		IsValid() const;
	bool		IsZero(const float &) const;
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

inline Vector2D &Vector2D::operator=(const Vector2D &v)
{
	x = v.x;
	y = v.y;

	return *this;
}

inline Vector2D &Vector2D::operator+=(const Vector2D &v)
{
	x += v.x;
	y += v.y;

	return *this;
}

inline Vector2D &Vector2D::operator-=(const Vector2D &v)
{
	x -= v.x;
	y -= v.y;

	return *this;
}

inline Vector2D &Vector2D::operator*=(const float &f)
{
	x *= f;
	y *= f;

	return *this;
}

inline Vector2D &Vector2D::operator/=(const float &f)
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

inline Vector2D Vector2D::operator*(const float &f) const
{
	return Vector2D(x * f, y * f);
}

inline Vector2D Vector2D::operator/(const float &f) const
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
	return sqrt((x * x) + (y * y));
}

inline float Vector2D::LengthSqr() const
{
	return (x * x) + (y * y);
}

inline float Vector2D::DistTo(const Vector2D &v) const
{
	return (*this - v).Length();
}

inline float Vector2D::DistToSqr(const Vector2D &v) const
{
	return (*this - v).LengthSqr();
}

inline bool Vector2D::WithinAARect(const Vector2D &m, const Vector2D &M) const // Rect Minimum, Rect Maximum
{
	return ((x > m.x) && (x < M.x) &&
		(y > m.y) && (y < M.y));
}

inline void Vector2D::Rotate(const float &f)
{
	float _x, _y;

	float s, c;

	SinCos(Deg2Rad(f), s, c);

	_x = x;
	_y = y;

	x = (_x * c) - (_y * s);
	y = (_x * s) + (_y * c);
}

inline float Vector2D::Dot(const Vector2D &v) const
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

inline void Vector2D::Random(const float &m, const float &M) // Minimum, Maximum
{
	x = m + ((float)rand() / RAND_MAX) * (M - m);
	y = m + ((float)rand() / RAND_MAX) * (M - m);
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

inline bool Vector2D::IsZero(const float &t = 0.01f) const // Tolerance
{
	return (x > -t && x < t &&
		y > -t && y < t);
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

	Euler		&operator+=(const Euler &);
	Euler		&operator-=(const Euler &);
	Euler		&operator*=(const float &);
	Euler		&operator/=(const float &);

	Euler		operator+(const Euler &) const;
	Euler		operator-(const Euler &) const;
	Euler		operator*(const float &) const;
	Euler		operator/(const float &) const;

	bool		operator==(const Euler &) const;
	bool		operator!=(const Euler &) const;

	float		Length() const;
	float		LengthSqr() const;
	float		Length2D() const;
	float		Length2DSqr() const;

	float		DistTo(const Vector &) const;
	float		DistToSqr(const Vector &) const;
	bool		WithinAABox(const Vector &, const Vector &) const;

	void		Rotate(const Angle &);
	void		Rotate2D(const float &);

	Vector		Cross(const Vector &) const;
	float		Dot(const Vector &) const;

	float		NormalizeVector();
	void		NormalizeAngle();

	Angle		&ToAngle() const;
	Vector		&ToVector() const;
	Vector2D	&ToVector2D() const;

	Vector		&Forward() const;
	Vector		&Right() const;
	Vector		&Up() const;

	void		Random(const float &, const float &);
	void		Negate();
	void		Clear();

	void		CopyToArray(float *) const;
	float		*Base() const;

	bool		IsValid() const;
	bool		IsZero(const float &) const;
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

inline Euler &Euler::operator+=(const Euler &e)
{
	x += e.x;
	y += e.y;
	z += e.z;

	return *this;
}

inline Euler &Euler::operator-=(const Euler &e)
{
	x -= e.x;
	y -= e.y;
	z -= e.z;

	return *this;
}

inline Euler &Euler::operator*=(const float &f)
{
	x *= f;
	y *= f;
	z *= f;

	return *this;
}

inline Euler &Euler::operator/=(const float &f)
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

inline Euler Euler::operator*(const float &f) const
{
	return Euler(x * f, y * f, z * f);
}

inline Euler Euler::operator/(const float &f) const
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
	return sqrt((x * x) + (y * y) + (z * z));
}

inline float Euler::LengthSqr() const
{
	return (x * x) + (y * y) + (z * z);
}

inline float Euler::Length2D() const
{
	return sqrt((x * x) + (y * y));
}

inline float Euler::Length2DSqr() const
{
	return (x * x) + (y * y);
}

inline float Euler::DistTo(const Vector &v) const
{
	return (*this - v).Length();
}

inline float Euler::DistToSqr(const Vector &v) const
{
	return (*this - v).LengthSqr();
}

inline bool Euler::WithinAABox(const Vector &m, const Vector &M) const // Box Minimum, Box Maximum
{
	return ((x > m.x) && (x < M.x) &&
			(y > m.y) && (y < M.y) &&
			(z > m.z) && (z < M.z));
}

inline void Euler::Rotate(const Angle &a)
{
	float _x, _y, _z;

	float s, c;

	SinCos(Deg2Rad(a.x), s, c);

	_y = y;
	_z = z;

	y = (_y * c) - (_z * s);
	z = (_y * s) + (_z * c);

	SinCos(Deg2Rad(a.y), s, c);

	_x = x;
	_z = z;

	x = (_x * c) + (_z * s);
	z = (-_x * s) + (_z * c);

	SinCos(Deg2Rad(a.z), s, c);

	_x = x;
	_y = y;

	x = (_x * c) - (_y * s);
	y = (_x * s) + (_y * c);
}

inline void Euler::Rotate2D(const float &f)
{
	float _x, _y;

	float s, c;

	SinCos(Deg2Rad(f), s, c);

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

inline float Euler::Dot(const Vector &v) const
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
	_NormalizeAngle(x);
	_NormalizeAngle(y);
	_NormalizeAngle(z);
}

inline Angle &Euler::ToAngle() const
{
	if (x == 0.f && y == 0.f)
		return Angle(z > 0.f ? -90.f : 90.f, 0.f, 0.f);

	return Angle(-Rad2Deg(ASin(z / (Length() + M_FLT_EPSILON))), Rad2Deg(ATan2(y, x)), 0.f);
}

inline Vector &Euler::ToVector() const
{
	static float sx, cx, sy, cy;

	SinCos(Deg2Rad(x), sx, cx);
	SinCos(Deg2Rad(y), sy, cy);

	return Vector(sx * cy, sx*sy, cx);
}

inline Vector2D &Euler::ToVector2D() const
{
	return Vector2D(x, y);
}

inline Vector &Euler::Forward() const
{
	static float sx, cx, sy, cy;

	SinCos(Deg2Rad(x), sx, cx);
	SinCos(Deg2Rad(y), sy, cy);

	return Vector(sx * cy, sx * sy, cx);
}

inline Vector &Euler::Right() const
{
	static float sx, cx, sy, cy, sz, cz;

	SinCos(Deg2Rad(x), sx, cx);
	SinCos(Deg2Rad(y), sy, cy);
	SinCos(Deg2Rad(z), sz, cz);

	return Vector((-sx * cy * sz) + (cz * sy), (-sx * sy * sz) - (cy * cz), -cx * sz);
}

inline Vector &Euler::Up() const
{
	static float sx, cx, sy, cy, sz, cz;

	SinCos(Deg2Rad(x), sx, cx);
	SinCos(Deg2Rad(y), sy, cy);
	SinCos(Deg2Rad(z), sz, cz);

	return Vector((sx * cy * cz) + (sy * sz), (sx * sy * cz) - (cy * sz), cx * cz);
}

inline void Euler::Random(const float &m, const float &M) // Minimum, Maximum
{
	x = m + ((float)rand() / RAND_MAX) * (M - m);
	y = m + ((float)rand() / RAND_MAX) * (M - m);
	z = m + ((float)rand() / RAND_MAX) * (M - m);
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

inline bool Euler::IsZero(const float &t = 0.01f) const // Tolerance
{
	return (x > -t && x < t &&
		y > -t && y < t &&
		z > -t && z < t);
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
		a.x = -Rad2Deg(asin(v.z / (v.Length() + M_FLT_EPSILON)));
		a.y = Rad2Deg(atan2(v.y, v.x));
	}

	a.z = 0.f;
}

inline void AngleToVector(const Angle &a, Vector *f,  Vector *r, Vector *u) // Angle, Forward, Right, Up
{
	float sx, cx, sy, cy, sz, cz;
	SinCos(Deg2Rad(a.x), sx, cx);
	SinCos(Deg2Rad(a.y), sy, cy);
	SinCos(Deg2Rad(a.z), sz, cz);

	if (f)
	{
		f->x = sx * cy;
		f->y = sx * sy;
		f->z = cx;
	}

	if (r)
	{
		r->x = cx*sy;
		r->y = (sx * sy * sz) + (cy * cz);
		r->z = (sx * sy * cz) - (cy * sz);
	}

	if (u)
	{
		u->x = -sx;
		u->y = cx * sz;
		u->z = cx * cz;
	}
}