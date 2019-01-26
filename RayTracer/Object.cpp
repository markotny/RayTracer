#include "pch.h"
#include "Object.h"

bool solve_quadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
	const auto delta = b * b - 4 * a * c;
	if (delta < 0) return false;
	
	if (delta == 0) {
		x0 = x1 = -0.5 * b / a;
	}
	else {
		float q = (b > 0) ?
			-0.5 * (b + sqrt(delta)) :
			-0.5 * (b - sqrt(delta));
		x0 = q / a;
		x1 = c / q;
	}
	return true;
}

bool sphere::intersect(const Vec3f & orig, const Vec3f & dir, float & dist) const
{
	float x0, x1; // solutions for dist if the ray intersects
	auto L = orig - center_;
	auto a = dir.dotProduct(dir);
	auto b = 2 * dir.dotProduct(L);
	auto c = L.dotProduct(L) - radius2_;

	if (!solve_quadratic(a, b, c, x0, x1)) return false;
	if (x0 > x1) std::swap(x0, x1);

	if (x0 < 0) {
		x0 = x1; // if t0 is negative, let's use t1 instead
		if (x0 < 0) return false; // both t0 and t1 are negative
	}

	dist = x0;

	return true;
}

void sphere::get_surface_data(const Vec3f& point_hit, Vec3f& normal_hit, Vec2f& tex) const
{
	normal_hit = point_hit - center_;
	normal_hit.normalize();
	// In this particular case, the normal is similar to a point on a unit sphere
	// centered around the origin. We can thus use the normal coordinates to compute
	// the spherical coordinates of Phit.
	// atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
	// acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
	tex.x = (1 + atan2(normal_hit.z, normal_hit.x) / M_PI) * 0.5;
	tex.y = acosf(normal_hit.y) / M_PI;
}


