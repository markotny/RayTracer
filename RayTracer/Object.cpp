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

/**
 * \brief Compute a ray-sphere intersection using the geometric solution
 */
bool sphere::intersect(const Vec3f & orig, const Vec3f & dir, float& t0, float& t1) const
{
	auto l = center - orig;
	auto tca = l.dotProduct(dir);
	if (tca < 0) return false;

	auto d2 = l.dotProduct(l) - tca * tca;
	if (d2 > radius2) return false;

	auto thc = sqrt(radius2 - d2);
	t0 = tca - thc;
	t1 = tca + thc;

	return true;
}

void sphere::get_surface_data(const Vec3f& point_hit, Vec3f& normal_hit) const
{
	normal_hit = point_hit - center;
	normal_hit.normalize();
}

bool plane::intersect(const Vec3f& orig, const Vec3f& dir, float& t0, float& t1) const
{
	// assuming vectors are all normalized
	const auto denom = dir.dotProduct(normal);
	if (denom > 1e-6 || denom < -1e-6) {
		auto p0_orig = p0 - orig;
		t0 = p0_orig.dotProduct(normal) / denom;
		t1 = t0;
		return t0 >= 0;
	}

	return false;
}

void plane::get_surface_data(const Vec3f& point_hit, Vec3f& normal_hit) const
{
	normal_hit = normal;
}

bool cylinder::intersect(const Vec3f& p1, const Vec3f& d1, float& t0, float& t1) const
{
	auto AB = pk - p0;
	auto AO = p1 - p0;
	auto AOxAB = AO.crossProduct(AB);
	auto VxAB = d1.crossProduct(AB);
	auto ab2 = AB.dotProduct(AB);
	auto a = VxAB.dotProduct(VxAB);
	auto b = 2 * VxAB.dotProduct(AOxAB);
	auto c = AOxAB.dotProduct(AOxAB) - radius2 * ab2;
	if (!solve_quadratic(a, b, c, t0, t1))
		return false;

	if (t0 > t1)
		std::swap(t0, t1);

	const auto dt0 = (p1 + t0 * d1 - p0).dotProduct(d0) * d0;
	const auto dt1 = (p1 + t1 * d1 - p0).dotProduct(d0) * d0;

	if (dt0.length() > (pk - p0).length())
	{
		if (dt1.length() > (pk - p0).length())
			return false;

		auto denom = d1.dotProduct(d0);
		auto p_p1 = pk - p1;
		t0 = p_p1.dotProduct(d0) / denom;

		if (dt1.dotProduct(d0) < 0)
		{
			denom = d1.dotProduct(-d0);
			p_p1 = p0 - p1;
			t1 = p_p1.dotProduct(-d0) / denom;
		}
	}

	if (dt0.dotProduct(d0) < 0)
	{
		if (dt1.dotProduct(d0) < 0)
			return false;

		auto denom = d1.dotProduct(-d0);
		auto p_p1 = p0 - p1;
		t0 = p_p1.dotProduct(-d0) / denom;

		if (dt1.length() > (pk-p0).length())
		{
			denom = d1.dotProduct(d0);
			p_p1 = pk - p1;
			t1 = p_p1.dotProduct(d0) / denom;
		}
	}
	return true;
}

void cylinder::get_surface_data(const Vec3f& point_hit, Vec3f& normal_hit) const
{
	const auto c0 = p0 + (point_hit - p0).dotProduct(d0) * d0;
	if ((c0 - p0).length() < 1e-4)	// point on bottom or top of cylinder
		normal_hit = -d0;
	else if ((c0 - pk).length() < 1e-4)
		normal_hit = d0;
	else
	{
		normal_hit = point_hit - c0;
		normal_hit.normalize();
	}
}


