#pragma once
#include "geometry.h"


class light
{
public:
	explicit light(const Vec3f &c = 1, const float &i = 1) :  color(c), intensity(i) {}
	virtual ~light() = default;
	virtual void illuminate(const Vec3f &P, Vec3f &, Vec3f &, float &) const = 0;
	Vec3f color;
	float intensity;
};

class distant_light : public light
{
public:
	Vec3f dir;

	explicit distant_light(const Vec3f &d, const Vec3f &c = 1, const float &i = 1) : light(c, i), dir(d)
	{
		dir.normalize();
	}
	void illuminate(const Vec3f &P, Vec3f &light_dir, Vec3f &light_intensity, float &distance) const override
	{
		light_dir = dir;
		light_intensity = color * intensity;
		distance = kInfinity;
	}
};

class point_light : public light
{
public:
	Vec3f pos;

	explicit point_light(const Vec3f p, const Vec3f &c = 1, const float &i = 1) : light(c, i), pos(p)
	{}

	// P: is the shaded point
	void illuminate(const Vec3f &P, Vec3f &light_dir, Vec3f &light_intensity, float &distance) const override
	{
		light_dir = (P - pos);
		auto r2 = light_dir.norm();
		distance = sqrt(r2);
		light_dir.x /= distance, light_dir.y /= distance, light_dir.z /= distance;
		// avoid division by 0
		light_intensity = color * (intensity / (4 * M_PI * r2));
	}
};