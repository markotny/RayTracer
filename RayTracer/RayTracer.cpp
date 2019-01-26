// RayTracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <cstdlib>
#include <cmath>
#include <GL/glut.h>
#include <GL/gl.h>
#include <vector>
#include "Object.h"

std::vector<std::unique_ptr<object>> objects;


typedef float point[3];
//
///**
// * \brief Funkcja obliczająca punkt przecięcia promienia i powierzchni sfery
// * \param p 
// * \param v 
// * \return 
// */
//int trace(float *p, float *v);

/**
 * \brief Funkcja obliczająca oświetlenie punktu na powierzchni sfery według modelu Phonga
 */
void phong(void);

/**
 * \brief Funkcja obliczająca iloczyn skalarny dwóch wektorów
 * \param p1 
 * \param p2 
 * \return 
 */
float dot_product(point p1, point p2);

/**
 * \brief Funkcja normalizująca wektor
 * \param p 
 */
void normalization(point p);

/**
 * \brief Rozmiar obrazu w pikselach (obraz jest kwadratem)
 */
int im_size = 400;

int image_width = 400;
int image_height = 400;
float fov = 50;

/**
 * \brief Rozmiar okna obserwatora
 */
float viewport_size = 3.0;


// Położenie i parametry źródła światła

float    light_position[] = { 3.0, 2.5, 5.0 };
float    light_specular[] = { 1.0, 1.0, 0.0 };
float    light_diffuse[] = { 0.0, 1.0, 1.0 };
float    light_ambient[] = { 0.0, 0.0, 0.0 };


// Promień i parametry rysowanej sfery

float    sphere_radius = 1.0;
float    sphere_specular[] = { 0.8, 0.8, 0.8 };
float    sphere_diffuse[] = { 0.6, 0.7, 0.8 };
float    sphere_ambient[] = { 1.0, 1.0, 1.0 };
float    sphere_specular_shininess = 30.0;


/**
 * \brief Parametry światła rozproszonego
 */
float global_a[] = { 0.25, 0.15, 0.1 };


// Parametry "śledzonego" promienia

/**
 * \brief punkt, z którego wychodzi promień
 */
float starting_point[3];

/**
 * \brief wektor opisujący kierunek promienia
 */
float starting_directions[] = { 0.0, 0.0, -1.0 };


// Zmienne pomocnicze

/**
 * \brief współrzędne (x,y,z) punktu przecięcia promienia i sfery
 */
float inter[3];

/**
 * \brief zmienna określająca, czy sfera została przecięta przez
 */
int inters;

/**
 * \brief składowe koloru dla oświetlonego punktu na powierzchni sfery
 */
float inters_c[3];

/**
 * \brief składowe koloru rysowanego piksela
 */
GLubyte pixel[1][1][3];

///**
// * \brief Funkcja oblicza punkt przecięcia promienia i powierzchni sfery
// * \param p punkt początkowy promienia
// * \param v wektor opisujący kierunek biegu promienia
// * \return 1 jeśli promień przecina sferę, 0 gdy nie przecina
// */
//int trace(float *p, float *v)
//{
//	float a, b, c, d, r;
//
//	a = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
//	b = 2 * (p[0] * v[0] + p[1] * v[1] + p[2] * v[2]);
//	c = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 1.0;
//
//	d = b * b - 4 * a*c;
//
//	if (d >= 0)                              // jest co najmniej jeden punkt przecięcia
//	{
//		r = (-b - sqrt(d)) / (2 * a);     // parametr dla bliższego punktu przecięcia
//
//		inter[0] = p[0] + r * v[0];    // współrzędne punktu przecięcia
//		inter[1] = p[1] + r * v[1];
//		inter[2] = p[2] + r * v[2];
//
//		return 1;                         // jest punkt przecięcia
//	}
//	else                                    // promień nie przecina sfery
//		return 0;
//}


/**
 * \brief Funkcja oblicza oświetlenie punktu na powierzchni sfery używając modelu Phonga
 */
void phong()
{
	float normal_vector[3];                      // wektor normalny do powierzchni
	float light_vec[3];                             // wektor wskazujący źródeł
	float reflection_vector[3];                  // wektor kierunku odbicia światła
	float viewer_v[3] = { 0.0, 0.0, 1.0 };   // wektor kierunku obserwacji
	float n_dot_l, v_dot_r;                      // zmienne pomocnicze

	normal_vector[0] = inter[0];         // wektor normalny do powierzchni sfery       
	normal_vector[1] = inter[1];
	normal_vector[2] = inter[2];

	light_vec[0] = light_position[0] - inter[0]; // wektor wskazujący kierunek źródła
	light_vec[1] = light_position[1] - inter[1];
	light_vec[2] = light_position[2] - inter[2];

	normalization(light_vec);                        // normalizacja wektora kierunku świecenia źródła           

	n_dot_l = dot_product(light_vec, normal_vector);

	reflection_vector[0] = 2 * (n_dot_l)*normal_vector[0] - light_vec[0];
	reflection_vector[1] = 2 * (n_dot_l)*normal_vector[1] - light_vec[1];
	reflection_vector[2] = 2 * (n_dot_l)*normal_vector[2] - light_vec[2];

	// obliczenie wektora opisującego kierunek światła odbitego od punktu na powierzchni sfery
	normalization(reflection_vector); // normalizacja wektora kierunku światła odbitego

	v_dot_r = dot_product(reflection_vector, viewer_v);

	if (v_dot_r < 0)                        // obserwator nie widzi oświetlanego punktu

		v_dot_r = 0;

	// sprawdzenie czy punkt na powierzchni sfery jest oświetlany przez źródło
	if (n_dot_l > 0)     // punkt jest oświetlany,oświetlenie wyliczane jest ze wzorów dla modelu Phonga
	{
		inters_c[0] = (sphere_diffuse[0] * light_diffuse[0] * n_dot_l + sphere_specular[0] * light_specular[0] * pow(double(v_dot_r), 20.0)) + sphere_ambient[0] * light_ambient[0] + sphere_ambient[0] * global_a[0];
		inters_c[1] = (sphere_diffuse[1] * light_diffuse[1] * n_dot_l + sphere_specular[1] * light_specular[1] * pow(double(v_dot_r), 20.0)) + sphere_ambient[1] * light_ambient[1] + sphere_ambient[1] * global_a[1];
		inters_c[2] = (sphere_diffuse[2] * light_diffuse[2] * n_dot_l + sphere_specular[2] * light_specular[2] * pow(double(v_dot_r), 20.0)) + sphere_ambient[2] * light_ambient[2] + sphere_ambient[2] * global_a[2];
	}
	else                // punkt nie jest oświetlany, uwzględniane jest tylko światło rozproszone  
	{
		inters_c[0] = sphere_ambient[0] * global_a[0];
		inters_c[1] = sphere_ambient[1] * global_a[1];
		inters_c[2] = sphere_ambient[2] * global_a[2];
	}
}

/**
 * \brief Funkcja przeprowadza normalizację wektora
 * \param p wektor do znormalizowania
 */
void normalization(point p)
{
	float d = 0.0;
	int i;

	for (i = 0; i < 3; i++)
		d += p[i] * p[i];

	d = sqrt(d);

	if (d > 0.0)
		for (i = 0; i < 3; i++)
			p[i] /= d;
}

/**
 * \brief Funkcja oblicza iloczyn skalarny wektorów
 * \param p1 wektor1
 * \param p2 wektor2
 * \return iloczyn skalarny
 */
float dot_product(point p1, point p2)
{
	return (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]);
}

bool trace(const Vec3f& orig, const Vec3f& dir, float& dist_nearest, const object* &hit_object)
{
	dist_nearest = kInfinity;
	std::vector<std::unique_ptr<object>>::const_iterator iter = objects.begin();
	for (; iter != objects.end(); ++iter) {
		auto t = kInfinity;
		if ((*iter)->intersect(orig, dir, t) && t < dist_nearest) {
			hit_object = iter->get();
			dist_nearest = t;
		}
	}
	return (hit_object != nullptr);
}

inline
Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mixValue)
{
	return a * (1 - mixValue) + b * mixValue;
}

Vec3f cast_ray(const Vec3f& orig, const Vec3f& dir)
{
	Vec3f hit_color = 0;
	const object *hit_object = nullptr; // this is a pointer to the hit object
	float dist; // this is the intersection distance from the ray origin to the hit point
	
	if (trace(orig, dir, dist, hit_object)) {
		auto point_hit = orig + dir * dist;
		Vec3f normal_hit;
		Vec2f tex;
		hit_object->get_surface_data(point_hit, normal_hit, tex);
		// Use the normal and texture coordinates to shade the hit point.
		// The normal is used to compute a simple facing ratio and the texture coordinate
		// to compute a basic checker board pattern
		float scale = 4;
		float pattern = (fmodf(tex.x * scale, 1) > 0.5) ^ (fmodf(tex.y * scale, 1) > 0.5);
		hit_color = std::max(0.f, normal_hit.dotProduct(-dir)) * mix(hit_object->color, hit_object->color * 0.8, pattern);
	}

	return hit_color;
}

inline
float clamp(const float &lo, const float &hi, const float &v)
{
	return std::max(lo, std::min(hi, v));
}

/**
 * \brief Funkcja rysująca obraz oświetlonej sceny
 */
void display()
{
	int  x, y;               // pozycja rysowanego piksela "całkowitoliczbowa"
	float x_fl, y_fl;      // pozycja rysowanego piksela "zmiennoprzecinkowa"
	int im_size_2;       // połowa rozmiaru obrazu w pikselach


	im_size_2 = im_size / 2;    // obliczenie połowy rozmiaru obrazu w pikselach

	glClear(GL_COLOR_BUFFER_BIT);

	glFlush();

	auto aspect_ratio = image_width / static_cast<float>(image_height); // assuming width > height

	/*unsigned width = im_size, height = im_size;
	auto inv_width = 1 / float(width), inv_height = 1 / float(height);
	float fov = 30, aspect_ratio = width / float(height);
	float angle = tan(M_PI * 0.5 * fov / 180.);*/

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

			/*float xx = (2 * ((x + 0.5) * inv_width) - 1) * angle * aspect_ratio;
			float yy = (1 - 2 * ((y + 0.5) * inv_height)) * angle;
			Vec3f ray_dir(xx, yy, -1);
			ray_dir.normalize();*/


			/*float Px = (2 * (x + 0.5) / (float)image_width - 1) * tan(fov / 2 * M_PI / 180) * aspect_ratio;
			float Py = (1 - 2 * (y + 0.5) / (float)image_height) * tan(fov / 2 * M_PI / 180);*/

			float Px = x_fl * tan(fov / 2 * M_PI / 180) * aspect_ratio;
			float Py = y_fl * tan(fov / 2 * M_PI / 180);
			Vec3f ray_origin(0);
			auto ray_direction = Vec3f(Px, Py, -1) - ray_origin; // note that this just equal to Vec3f(Px, Py, -1);
			ray_direction.normalize();	// it's a direction so don't forget to normalize 

			// przeliczenie pozycji(x,y) w pikselach na pozycję "zmiennoprzecinkową" w oknie obserwatora

			/*starting_point[0] = x_fl;
			starting_point[1] = y_fl;
			starting_point[2] = viewport_size;*/

			//Vec3f orig{ x_fl, y_fl, viewport_size };
			//Vec3f dir{ 0.0, 0.0, -1.0 };
			//auto pix = cast_ray(orig, dir);
			auto pix = cast_ray(Vec3f(0), ray_direction);

			*(pixs++) = pix;
			

			pixel[0][0][0] = 255 * clamp(0, 1, pix[0]);
			pixel[0][0][1] = 255 * clamp(0, 1, pix[1]);
			pixel[0][0][2] = 255 * clamp(0, 1, pix[2]);

			glRasterPos3f(x_fl, y_fl, 0);

			// inkrementacja pozycji rastrowej dla rysowania piksela

			glDrawPixels(1, 1, GL_RGB, GL_UNSIGNED_BYTE, pixel);

			// narysowanie kolejnego piksela na ekranie
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

/**
 * \brief Funkcja inicjalizująca definiująca sposób rzutowania
 */
void my_init(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-viewport_size / 2, viewport_size / 2, -viewport_size / 2, viewport_size / 2, -viewport_size / 2, viewport_size / 2);
	glMatrixMode(GL_MODELVIEW);
	
	// generate a scene made of random spheres
	const uint32_t num_spheres = 32;
	gen.seed(0);

	for (uint32_t i = 0; i < num_spheres; ++i) {
		Vec3f rand_pos( (0.5 - dis(gen)) * 10, (0.5 - dis(gen)) * 10, -1 - dis(gen) * 10);
		float rand_radius = (0.5 + dis(gen) * 0.5);
		objects.push_back(std::unique_ptr<object>(new sphere(rand_pos, rand_radius)));
	}
	//objects.push_back(std::unique_ptr<object>(new sphere(Vec3f(0,0,-2), 1.0)));
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
