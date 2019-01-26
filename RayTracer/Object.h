#pragma once
#include "pch.h"
#include "geometry.h"

const float kInfinity = std::numeric_limits<float>::max();
static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<> dis(0, 1);

class object
{
public:
	Vec3f color;
	object() : color(dis(gen), dis(gen), dis(gen)){}
	~object(){}

	/**
	 * \brief compute the intersection of the object with a ray
	 * \param orig the ray origin
	 * \param dir the ray direction
	 * \param dist distance from ray origin to intersection point
	 * \return true if an intersection was found, false otherwise
	 */
	virtual bool intersect(const Vec3f &orig, const Vec3f &dir, float &dist) const = 0;

	/**
	 * \brief 
	 * \param point_hit point ont the surface we want to get data on
	 * \param normal_hit normal at point_hit
	 * \param tex texture coordinates at point_hit
	 */
	virtual void get_surface_data(const Vec3f &point_hit, Vec3f &normal_hit, Vec2f &tex) const = 0;
};

class sphere : public object
{
	float radius_, radius2_;
	Vec3f center_;
public:
	sphere(const Vec3f &c, const float &r) : radius_(r), radius2_(r * r), center_(c) {}

	bool intersect(const Vec3f& orig, const Vec3f& dir, float& dist) const override;

	void get_surface_data(const Vec3f& point_hit, Vec3f& normal_hit, Vec2f& tex) const override;
};