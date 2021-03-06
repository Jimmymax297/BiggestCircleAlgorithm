#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <chrono>
#include "BruteGeometry.h"
#include "Segment.h"
#include "Point.h"
#ifdef _WIN32
#include "getopt.h"
#endif

using namespace std::chrono;

void generate(std::vector<Point>& v, int n)
{
	v.resize(n);

	point *ps = new point[n];
	random_points(ps, n, MIN, MAX, true);
	space_partition(ps, n);

	std::transform(ps, ps + n, v.begin(), [](point p)->Point {return Point(p.x, p.y);});

	delete[] ps;
}

int main(int argc, char **argv)
{
	BruteGeometry bg;

	int mode = 0, n = 0;
	int k = 0, s = 0, r = 0;

	for (;;)
	{
		int opt = getopt(argc, argv, "m:n:k:s:r:");
		if (opt == -1)
			break;

		switch (opt)
		{
		case 'm':
			try
			{
				mode = std::stoi(std::string(optarg));
			}
			catch (std::exception) 
			{
				std::cerr << "Invalid argument for option 'm'\n";
				return -1;
			}
			break;
		case 'n':
			try
			{
				n = std::stoi(std::string(optarg));
			}
			catch (std::exception)
			{
				std::cerr << "Invalid argument for option 'n'\n";
				return -1;
			}
			break;
		case 'k':
			try
			{
				k = std::stoi(std::string(optarg));
			}
			catch (std::exception)
			{
				std::cerr << "Invalid argument for option 'k'\n";
				return -1;
			}
			break;
		case 's':
			try
			{
				s = std::stoi(std::string(optarg));
			}
			catch (std::exception)
			{
				std::cerr << "Invalid argument for option 's'\n";
				return -1;
			}
			break;
		case 'r':
			try
			{
				r = std::stoi(std::string(optarg));
			}
			catch (std::exception)
			{
				std::cerr << "Invalid argument for option 'r'\n";
				return -1;
			}
			break;
		case '?':
			std::cerr << "Unknown option\n";
			return -1;
		}
	}
	
	std::vector<milliseconds> times;
	switch (mode)
	{
	case 1:
		std::cin >> n;

		if (n < 3)
		{
			std::cerr << "Invalid number of vertices\n";
			return -1;
		}

		bg.input.resize(n);

		for (int i = 0; i < n; ++i)
		{
			int x, y;
			std::cin >> x >> y;
			bg.input[i].setX(x);
			bg.input[i].setY(y);
		}
		break;
	case 2:
		if (n < 3)
		{
			std::cerr << "Invalid number of vertices\n";
			return -1;
		}
		srand(time(NULL));
		generate(bg.input, n);
		break;
	case 3:
	{
		if (n < 3)
		{
			std::cerr << "Invalid number of vertices\n";
			return -1;
		}
		if (k <= 0)
		{
			std::cerr << "Invalid number of steps\n";
			return -1;
		}
		if (s <= 0)
		{
			std::cerr << "Invalid step value\n";
			return -1;
		}
		if (r <= 0)
		{
			std::cerr << "Invalid number of problems per step\n";
			return -1;
		}
		srand(time(NULL));

		times.resize(k);
		for (int i = 0; i < k; ++i)
		{
			std::vector<std::vector<Point>> problems(r);
			for (int j = 0; j < r; ++j)
			{
				generate(problems[j], n);
			}
			milliseconds start = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
			for (int j = 0; j < r; ++j)
			{
				bg.input = problems[j];
				bg.process();
				bg.max();
			}
			milliseconds stop = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
			times[i] = (stop - start) / r;
			n += s;
		}
		n -= s * k;
		int med = k / 2;
		int n_med = n + s * med;
		double q = (double)(n_med * n_med * n_med * n_med) / times[med].count();
		std::cout << "n  t(n)  q(n)\n";
		for (int i = 0; i < k; ++i)
			std::cout << n + s * i << ' ' << times[i].count()  << ' ' << times[i].count() * q / ((n + s * i)  * (n + s * i) * (n + s * i) * (n + s * i)) << std::endl;
		return 0;
	}
	default:
		std::cerr << "Unsupported mode";
		break;
	}
	bg.process();
	//std::for_each(bg.output.begin(), bg.output.end(), [](Circle& c)->void{cout << c.getCenter().getX() << ' ' << c.getCenter().getY() << ' ' << c.getRadius() << endl;});
	Circle c = bg.max();
	std::cout << c.getCenter().getX() << ' ' << c.getCenter().getY() << ' ' << c.getRadius() << std::endl;
	std::cout << bg.input.size() << std::endl;
	for (auto& p: bg.input)
		std::cout << p.getX() << ' ' << p.getY() << std::endl;
	return 0;
}

