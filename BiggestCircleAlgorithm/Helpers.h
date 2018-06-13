#pragma once
#include <iostream>
#include <ctime>
#include <string>
#include <random>
#include <set>
#include <cmath>
#include <limits>
#include <tuple>

#define MIN 0
#define MAX 800

constexpr double epsilon = (double)std::numeric_limits<float>::epsilon();
inline double A(std::tuple<double, double, double> params) { return std::get<0>(params); }
inline double B(std::tuple<double, double, double> params) { return std::get<1>(params); }
inline double C(std::tuple<double, double, double> params) { return std::get<2>(params); }

inline bool fuzzyEquals(double d1, double d2)
{
	return abs(d1 - d2) < epsilon;
}

typedef struct _point
{
	double x;
	double y;
	_point() : x(0.0), y(0.0) {}
	_point(double xp, double yp) : x(xp), y(yp) {}
} point;


double randf();
void random_points(point ps[], int count, double min, double max, bool integer = false);
point random_point_on_segment(const point &s, const point &e);
bool to_left(const point &s, const point &e, const point &p);
void swap(point ps[], int i, int j);
void space_partion_rec(point ps[], int count, int l, int r);
void space_partition(point ps[], int count);

/*int main(int argc, char *argv[])
{
	srand(time(NULL));
	if (argc != 2)
	{
		std::cerr << "Not enough arguments" << std::endl;
		return -1;
	}

	int count = atoi(argv[1]);

	point *ps = new point[count];
	random_points(ps, count, MIN, MAX, true);
	space_partition(ps, count);


	std::cout << count << std::endl;
	for (int i = 0; i < count; i++)
		std::cout << (int)ps[i].x << " " << (int)ps[i].y << std::endl;

	delete[] ps;
	return 0;
}*/

