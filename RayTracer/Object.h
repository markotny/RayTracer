#pragma once
#include "pch.h"
#include "geometry.h"

static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<> dis(0, 1);

enum material_type
{
	diffuse,
	glass,
	metal
};

class object
{
public:
	material_type material;
	Vec3f surface_color;
	float transparency, reflection;
	Vec3f albedo = 0.18; // reflected/received ratio
	float Kd = 0.8; // phong model diffuse weight
	float Ks = 0.2; // phong model specular weight
	float n = 10;   // phong specular exponent
	object(
		const Vec3f &surf_col,
		const float &refl,
		const float &transp,
		const material_type mat)
	: material(mat), surface_color(surf_col), transparency(transp), reflection(refl)
	{
		if (mat == glass)
		{
			Kd = 0.4;
			Ks = 0.1;
			n = 20;
		}
		else if (mat == metal)
		{
			Kd = 0.8;
			Ks = 0.15;
			n = 100;
		}
	}

	~object() = default;

	/**
	 * \brief compute the intersection of the object with a ray
	 * \param orig the ray origin
	 * \param dir the ray direction
	 * \param t0 distance from ray origin to first intersection point
	 * \param t0 distance from ray origin to second intersection point
	 * \return true if an intersection was found, false otherwise
	 */
	virtual bool intersect(const Vec3f &orig, const Vec3f &dir, float& t0, float& t1) const = 0;

	/**
	 * \brief 
	 * \param point_hit point ont the surface we want to get data on
	 * \param normal_hit normal at point_hit
	 */
	virtual void get_surface_data(const Vec3f &point_hit, Vec3f &normal_hit) const = 0;
};

class sphere : public object
{
public:
	float radius, radius2;
	Vec3f center;

	sphere(
		const Vec3f &c,
		const float &r,
		const Vec3f &surf_col,
		const material_type mat = glass,
		const float &refl = 0,
		const float &transp = 0)
		: object(surf_col, refl, transp, mat),
		  radius(r), radius2(r * r), center(c)
	{}

	bool intersect(const Vec3f& orig, const Vec3f& dir, float& t0, float& t1) const override;

	void get_surface_data(const Vec3f& point_hit, Vec3f& normal_hit) const override;
};

class plane : public object
{
public:
	Vec3f p0, normal;	//plane defined by a point and a normal

	plane(
		const Vec3f &p,
		const Vec3f &n,
		const Vec3f &surf_col,
		const material_type mat = diffuse,
		const float &refl = 0,
		const float &transp = 0)
		: object(surf_col, refl, transp, mat),
		p0(p), normal(n)
	{}

	bool intersect(const Vec3f& orig, const Vec3f& dir, float& t0, float& t1) const override;

	void get_surface_data(const Vec3f& point_hit, Vec3f& normal_hit) const override;
};

class cylinder : public object
{
public:
	Vec3f p0, pk, d0; //starting and ending points, d0 direction vector
	float radius, radius2;

	cylinder(
		const Vec3f &p,
		const Vec3f &k,
		const float &r,
		const Vec3f &surf_col,
		const material_type mat = metal,
		const float &refl = 0,
		const float &transp = 0)
		: object(surf_col, refl, transp, mat),
		p0(p), pk(k), radius(r), radius2(r * r)
	{
		d0 = pk - p0;
		d0.normalize();
	}

	bool intersect(const Vec3f& p1, const Vec3f& d1, float& t0, float& t1) const override;

	void get_surface_data(const Vec3f& point_hit, Vec3f& normal_hit) const override;
};
