// RayTracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <cstdlib>
#include <cmath>
#include <GL/glut.h>
#include <GL/gl.h>
#include <vector>
#include "Object.h"
#include "light.h"

#define MAX_RAY_DEPTH 5
std::vector<std::unique_ptr<object>> objects;
std::vector<std::unique_ptr<light>> lights;
Vec3f background_color(0.8);
float ambient_intensity(0.4);
float bias = 1e-4; // add some bias to the point from which we will be tracing

int im_size = 800;

int image_width = 800;
int image_height = 800;
float fov = 50;

float viewport_size = 3.0;
GLubyte pixel[1][1][3];          // składowe koloru rysowanego piksela

bool trace(const Vec3f& orig, const Vec3f& dir, float& t_nearest, const object* &hit_object)
{
	t_nearest = kInfinity;

	for (auto &object : objects)
	{
		auto t0 = kInfinity, t1 = kInfinity;
		if (object->intersect(orig, dir, t0, t1)) {
			if (t0 < 0) t0 = t1;
			if (t0 < t_nearest) {
				hit_object = object.get();
				t_nearest = t0;
			}
		}
	}
	return (hit_object != nullptr);
}

inline
Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mix_val)
{
	return a * (1 - mix_val) + b * mix_val;
}

inline
float mix(const float &a, const float &b, const float &mix_val)
{
	return a * (1 - mix_val) + b * mix_val;
}

inline
float clamp(const float &lo, const float &hi, const float &v)
{
	return std::max(lo, std::min(hi, v));
}

Vec3f reflect(const Vec3f& dir, const Vec3f& n)
{
	return dir - 2 * dir.dotProduct(n) * n;
}

Vec3f phong(const object* &hit_object, const Vec3f &point_hit, const Vec3f normal_hit, const Vec3f &dir)
{
	Vec3f diffuse = hit_object->surface_color * ambient_intensity, specular = 0;
	for (auto& light : lights)
	{
		Vec3f light_dir, light_intensity;
		float dist;
		light->illuminate(point_hit, light_dir, light_intensity, dist);

		const object * covering_obj = nullptr;
		float dist_to_covering_obj;
		const auto visible = !trace(point_hit + normal_hit * bias, -light_dir, dist_to_covering_obj, covering_obj);

		// compute the diffuse component
		if (visible)
		{
			diffuse += hit_object->surface_color * hit_object->albedo * light_intensity * std::max(0.f, normal_hit.dotProduct(-light_dir));

			// compute the specular component
			// what would be the ideal reflection direction for this light ray
			auto reflected = reflect(light_dir, normal_hit);
			specular += light_intensity * std::pow(std::max(0.f, reflected.dotProduct(-dir)), hit_object->n);
		}
	}
	return diffuse * hit_object->Kd + specular * hit_object->Ks;
}

Vec3f cast_ray(const Vec3f& orig, const Vec3f& dir, const int& depth)
{
	auto hit_color = background_color;
	const object *hit_object = nullptr; // this is a pointer to the hit object
	float dist; // this is the intersection distance from the ray origin to the hit point
	
	if (trace(orig, dir, dist, hit_object)) {
		const auto point_hit = orig + dir * dist;
		Vec3f normal_hit;
		hit_object->get_surface_data(point_hit, normal_hit);
		
		auto inside = false;
		if (dir.dotProduct(normal_hit) > 0)
			normal_hit = -normal_hit, inside = true;

		hit_color += phong(hit_object, point_hit, normal_hit, dir);

		if (hit_object->material != diffuse)
		{
			if ((hit_object->transparency > 0 || hit_object->reflection > 0) && depth < MAX_RAY_DEPTH)
			{
				const auto facing_ratio = -dir.dotProduct(normal_hit);

				// change the mix value to tweak the effect
				const auto fresnel_effect = mix(pow(1 - facing_ratio, 3), 1, 0.1);

				// compute reflection direction (not need to normalize because all vectors
				// are already normalized)
				auto refl_dir = reflect(dir, normal_hit);
				refl_dir.normalize();
				const auto reflection = cast_ray(point_hit + normal_hit * bias, refl_dir, depth + 1);

				Vec3f refraction = 0;
				// if the sphere is also transparent compute refraction ray (transmission)
				if (hit_object->transparency) {
					float ior = 1.4, eta = inside ? ior : 1 / ior; // are we inside or outside the surface?
					auto cosi = -normal_hit.dotProduct(dir);
					auto k = 1 - eta * eta * (1 - cosi * cosi);
					auto refr_dir = dir * eta + normal_hit * (eta *  cosi - sqrt(k));
					refr_dir.normalize();
					refraction = cast_ray(point_hit - normal_hit * bias, refr_dir, depth + 1);
				}

				// the result is a mix of reflection and refraction (if the sphere is transparent)
				hit_color = (
					reflection * fresnel_effect +
					refraction * (1 - fresnel_effect) * hit_object->transparency
					) * hit_object->surface_color;
			}
		}

		return hit_color;
	}

	return hit_color;
}


void display()
{
	int  x, y;               // pozycja rysowanego piksela "całkowitoliczbowa"
	float x_fl, y_fl;      // pozycja rysowanego piksela "zmiennoprzecinkowa"
	int im_size_2;       // połowa rozmiaru obrazu w pikselach


	im_size_2 = im_size / 2;    // obliczenie połowy rozmiaru obrazu w pikselach

	glClear(GL_COLOR_BUFFER_BIT);

	glFlush();

	auto aspect_ratio = image_width / static_cast<float>(image_height); // assuming width > height
	
	auto framebuffer = new Vec3f[image_width * image_height];
	auto pixs = framebuffer;

	// rysowanie pikseli od lewego górnego narożnika do prawego dolnego narożnika

	uint32_t j = 0, i = 0;
	for (y = im_size_2; y > -im_size_2; y--, j++)
	{
		for (x = -im_size_2; x < im_size_2; x++, i++)
		{

			x_fl = (float)x / (im_size / viewport_size);
			y_fl = (float)y / (im_size / viewport_size);
			
			float Px = x_fl * tan(fov / 2 * M_PI / 180) * aspect_ratio;
			float Py = y_fl * tan(fov / 2 * M_PI / 180);
			Vec3f ray_origin(0);
			auto ray_direction = Vec3f(Px, Py, -1) - ray_origin; // note that this just equal to Vec3f(Px, Py, -1);
			ray_direction.normalize();	// it's a direction so don't forget to normalize 

			auto pix = cast_ray(Vec3f(0), ray_direction, 0);

			*(pixs++) = pix;
			
			pixel[0][0][0] = 255 * clamp(0, 1, pix[0]);
			pixel[0][0][1] = 255 * clamp(0, 1, pix[1]);
			pixel[0][0][2] = 255 * clamp(0, 1, pix[2]);

			glRasterPos3f(x_fl, y_fl, 0);
			glDrawPixels(1, 1, GL_RGB, GL_UNSIGNED_BYTE, pixel);
		}
	}
	glFlush();

	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::ofstream ofs("./out.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << im_size << " " << im_size << "\n255\n";
	for (uint32_t i = 0; i < im_size * im_size; ++i) {
		char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
		char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
		char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
		ofs << r << g << b;
	}

	ofs.close();
	std::cout << "Wrote to file.\n";
	delete[] framebuffer;
}

void setup_scene1()
{
	// generate a scene made of random spheres
	const uint32_t num_spheres = 20;
	gen.seed(0);
	//lights
	lights.push_back(std::unique_ptr<light>(std::make_unique<distant_light>(
		Vec3f(-0.5, -0.5, 0), 1, 5)));

	lights.push_back(std::unique_ptr<light>(std::make_unique<point_light>(
		Vec3f(0, 4, -1), 1, 1000)));

	for (uint32_t i = 0; i < num_spheres; ++i) {
		Vec3f rand_pos((0.5 - dis(gen)) * 10, (0.5 - dis(gen)) * 10, -1 - dis(gen) * 10);
		float rand_radius = (0.5 + dis(gen) * 0.5);
		float refl = dis(gen);
		objects.push_back(std::unique_ptr<object>(
			std::make_unique<sphere>(rand_pos, rand_radius, Vec3f(dis(gen), dis(gen), dis(gen)), glass, refl, 1-refl)));
	}
}

void setup_scene2()
{
	gen.seed();

	//lights
	lights.push_back(std::unique_ptr<light>(std::make_unique<distant_light>(
		Vec3f(0, -1, 0), 1, 0.7)));

	lights.push_back(std::unique_ptr<light>(std::make_unique<point_light>(
		Vec3f(0, 4, -1), 1, 1000)));

	lights.push_back(std::unique_ptr<light>(std::make_unique<point_light>(
		Vec3f(-2, 10, -13), 1, 1000)));

	//main sphere
	objects.push_back(std::unique_ptr<object>(std::make_unique<sphere>(
		Vec3f(0, 1, -7), 3.0,
		Vec3f(1, 0.9, 0.1), glass, 0.1, 0.9)));
	
	//inside
	objects.push_back(std::unique_ptr<object>(std::make_unique<sphere>(
		Vec3f(0, 1, -7), 1,
		Vec3f(1,0.84,0),glass)));
	
	//reflective spheres
	objects.push_back(std::unique_ptr<object>(std::make_unique<sphere>(
		Vec3f(-3, -5, -3), 4,
		Vec3f(0.5, 0.7, 1),glass,1)));

	objects.push_back(std::unique_ptr<object>(std::make_unique<sphere>(
		Vec3f(18, -12, -25), 25,
		Vec3f(0.4, 0.5, 0.8), glass, 1)));

	//triangle
	objects.push_back(std::unique_ptr<object>(std::make_unique<cylinder>(
		Vec3f(-2, -3, -5), Vec3f(5, 5, -9), 0.5,
		Vec3f(0.8, 0.1, 0.1),glass,1,1)));

	objects.push_back(std::unique_ptr<object>(std::make_unique<cylinder>(
		Vec3f(-2, -3, -5), Vec3f(-5, 4, -9), 0.5,
		Vec3f(0.8, 0.1, 0.1), glass, 1, 1)));

	objects.push_back(std::unique_ptr<object>(std::make_unique<cylinder>(
		Vec3f(-5, 4, -9), Vec3f(5, 5, -9), 0.5,
		Vec3f(0.8, 0.1, 0.1), glass, 1, 1)));

	//planes
	objects.push_back(std::unique_ptr<object>(std::make_unique<plane>(
		Vec3f(0, -4, 0), Vec3f(0.2, 1, 0.3).normalize(),
		Vec3f(0.8, 1, 1), metal, 1)));

	objects.push_back(std::unique_ptr<object>(std::make_unique<plane>(
		Vec3f(0, 0, -20), Vec3f(-0.2, 1, 1).normalize(),
		Vec3f(0, 1, 1), metal, 1)));
}

void my_init(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-viewport_size / 2, viewport_size / 2, -viewport_size / 2, viewport_size / 2, -viewport_size / 2, viewport_size / 2);
	glMatrixMode(GL_MODELVIEW);
	
	setup_scene2();
}

int main()
{
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	glutInitWindowSize(im_size, im_size);
	glutCreateWindow("Ray Casting");
	my_init();

	glutDisplayFunc(display);
	glutMainLoop();
}
