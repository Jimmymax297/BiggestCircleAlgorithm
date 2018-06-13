#pragma once
#include "Point.h"
#include "Segment.h"
#include "Circle.h"
#include "Helpers.h"
#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <vector>
#include <list>

class BruteGeometry
{
public:
	std::vector<Point> in;
	std::vector<Point>& input;
	std::list<Circle> output;

	BruteGeometry(): in(), input(in), output() {}

	void process()
	{
		output.clear();
		#define ind(x) ((x + input.size()) % input.size())
		for (int i = 0; i < input.size(); ++i)
			for (int j = i + 1; j < input.size(); ++j)
				for (int k = j + 1; k < input.size(); ++k)
				{
					circleCenter(input[ind(i)], input[ind(j)], input[ind(k)]);
					circleCenter(input[ind(i)], input[ind(j)], Segment(input[ind(k)], input[ind(k + 1)]));
					circleCenter(input[ind(i)], input[ind(k)], Segment(input[ind(j)], input[ind(j + 1)]));
					circleCenter(input[ind(j)], input[ind(k)], Segment(input[ind(i)], input[ind(i + 1)]));
					circleCenter(Segment(input[ind(j)], input[ind(j + 1)]), Segment(input[ind(k)], input[ind(k + 1)]), input[ind(i)]);
					circleCenter(Segment(input[ind(i)], input[ind(i + 1)]), Segment(input[ind(k)], input[ind(k + 1)]), input[ind(j)]);
					circleCenter(Segment(input[ind(i)], input[ind(i + 1)]), Segment(input[ind(j)], input[ind(j + 1)]), input[ind(k)]);
					circleCenter(Segment(input[ind(i)], input[ind(i + 1)]), Segment(input[ind(j)], input[ind(j + 1)]), Segment(input[ind(k)], input[ind(k + 1)]));
				}
		//std::for_each(output.begin(), output.end(), [](Circle& c)->void {std::cout << c.getCenter().getX() << ' ' << c.getCenter().getY() << ' ' << c.getRadius() << std::endl;});
		//std::cout << std::endl;

		auto it = std::remove_if(output.begin(), output.end(), [this](Circle & c) {return !containsPoint(c.getCenter());});
		output.erase(it, output.end());

		//std::for_each(output.begin(), output.end(), [](Circle& c)->void {std::cout << c.getCenter().getX() << ' ' << c.getCenter().getY() << ' ' << c.getRadius() << std::endl;});
		//std::cout << std::endl;
		it = std::remove_if(output.begin(), output.end(), [this](Circle & c) {
			for (int i = 0; i < input.size(); ++i)
			{
				Segment s(input[ind(i)], input[ind(i + 1)]);
				if (!fuzzyEquals(distanceFromSegment(c.getCenter(), s), c.getRadius()) && distanceFromSegment(c.getCenter(), s) < c.getRadius())
					return true;
			}
			return false;
		});
		output.erase(it, output.end());
#undef ind(x)
		//std::for_each(output.begin(), output.end(), [](Circle& c)->void {std::cout << c.getCenter().getX() << ' ' << c.getCenter().getY() << ' ' << c.getRadius() << std::endl;});
		//std::cout << std::endl;
	}

	Circle max()
	{
		Circle maxCircle(Point(0, 0), 0);
		for (auto &c : output)
			if (c.getRadius() > maxCircle.getRadius())
				maxCircle = c;
		//std::cout << maxCircle.getCenter().getX() << ' ' << maxCircle.getCenter().getY() << ' ' << maxCircle.getRadius() << std::endl;
		return maxCircle;
	}

	static void intersect_line(const Point& p1, const Point& p2, const Point& pos, int *winding)
	{
		double x1 = p1.getX();
		double y1 = p1.getY();
		double x2 = p2.getX();
		double y2 = p2.getY();
		double y = pos.getY();

		int dir = 1;

		if (fuzzyEquals(y1, y2))
			// ignore horizontal lines according to scan conversion rule
			return;
		else if (y2 < y1)
		{
			double x_tmp = x2; x2 = x1; x1 = x_tmp;
			double y_tmp = y2; y2 = y1; y1 = y_tmp;
			dir = -1;
		}
		if (y >= y1 && y < y2)
		{
			double x = x1 + ((x2 - x1) / (y2 - y1)) * (y - y1);
			// count up the winding number if we're
			if (x <= pos.getX())
				(*winding) += dir;
		}
	}

	bool containsPoint(const Point& pt) const
	{
		if (input.empty())
			return false;

		int winding_number = 0;

		Point last_pt = input[0];
		Point last_start = input[0];
		for (int i = 1; i < input.size(); ++i)
		{
			const Point &e = input[i];
			intersect_line(last_pt, e, pt, &winding_number);
			last_pt = e;
		}

		// implicitly close last subpath
		/*if (last_pt != last_start)*/
		intersect_line(last_pt, last_start, pt, &winding_number);

		return (winding_number % 2) != 0;
	}

	static inline double distanceFromLine(const Point& p, double a, double b, double c)
	{
		return (abs(a * p.getX() + b * p.getY() + c) / sqrt(a * a + b * b));
	}

	static inline double distanceFromLine(const Point& p, std::tuple<double, double, double> line)
	{
		return distanceFromLine(p, A(line), B(line), C(line));
	}

	static double distanceFromSegment(const Point& p, const Segment& segment)
	{
		/*double lineDistance = distanceFromLine(p, Segment::lineSegment(segment));
		double distanceP1 = Point::distance(p, segment.getPointA());
		double distanceP2 = Point::distance(p, segment.getPointB());
		return std::min(std::max(lineDistance, distanceP1), std::max(lineDistance, distanceP2));*/
		auto line = Segment::lineSegment(segment);
		Point q = intersectionBLines(line, perpendicularLine(p, segment));
		if ((segment.getPointA().getX() <= q.getX() && q.getX() <= segment.getPointB().getX() && segment.getPointA().getY() <= q.getY() && q.getY() <= segment.getPointB().getY())
			|| (segment.getPointA().getX() >= q.getX() && q.getX() >= segment.getPointB().getX() && segment.getPointA().getY() <= q.getY() && q.getY() <= segment.getPointB().getY())
			|| (segment.getPointA().getX() <= q.getX() && q.getX() <= segment.getPointB().getX() && segment.getPointA().getY() >= q.getY() && q.getY() >= segment.getPointB().getY())
			|| (segment.getPointA().getX() >= q.getX() && q.getX() >= segment.getPointB().getX() && segment.getPointA().getY() >= q.getY() && q.getY() >= segment.getPointB().getY()))
			return distanceFromLine(p, line);
		return std::min(Point::distance(p, segment.getPointA()), Point::distance(p, segment.getPointB()));
	}

	static std::tuple<double, double, double> perpendicularLine(const Point& p, const Segment& s)
	{
		auto line = Segment::lineSegment(s);
		double c = B(line) * p.getX() - A(line)*p.getY();
		return std::make_tuple(-B(line), A(line), c);
	}

	static double distanceFromParallelLines(const Segment& segmentA, const Segment& segmentB)
	{
		auto lineA = Segment::lineSegment(segmentA);
		auto lineB = Segment::lineSegment(segmentB);

		double a1 = A(lineA);
		double a2 = A(lineB);
		double b1 = B(lineA);
		double b2 = B(lineB);
		double c1 = C(lineA);
		double c2 = C(lineB);

		if (b1 == 0)
		{
			return abs(c1 / a1 - c2 / a2);
		}

		Point test = Point(1, -a1 / b1 - c1 / b1);

		return distanceFromLine(test, lineB);
	}

	static inline bool collinearPoints(const Point& p1, const Point& p2, const Point& p3)
	{
		return ((p1.getY() - p2.getY()) / (p1.getX() - p2.getX())) == ((p1.getY() - p3.getY()) / (p1.getX() - p3.getX()));
	}

	static inline bool parallelSegments(const Segment& segmentA, const Segment& segmentB)
	{
		auto lineA = Segment::lineSegment(segmentA);
		auto lineB = Segment::lineSegment(segmentB);
		if ((A(lineA) == 0 && B(lineB) == 0) || (B(lineA) == 0 && A(lineB) == 0))
			return false;
		return (A(lineA) * B(lineB) == A(lineB) * B(lineA));
	}

	static inline Point intersectionBLines(const std::tuple<double, double, double>& line1, const std::tuple<double, double, double>& line2)
	{
		double a1 = A(line1);
		double a2 = A(line2);
		double b1 = B(line1);
		double b2 = B(line2);
		double c1 = C(line1);
		double c2 = C(line2);

		double w = a1 * b2 - a2 * b1;
		double wx = -c1 * b2 + c2 * b1;
		double wy = -a1 * c2 + a2 * c1;

		return Point(wx / w, wy / w);
	}

	static std::tuple<double, double, double> bisector(Segment segmentA, Segment segmentB)
	{
		if (parallelSegments(segmentA, segmentB))
		{
			throw std::exception();
		}
		double a1, a2, b1, b2, c1, c2;

		auto tupleA = Segment::lineSegment(segmentA);
		auto tupleB = Segment::lineSegment(segmentB);

		double sq1 = sqrt(A(tupleA) * A(tupleA) + B(tupleA) * B(tupleA));
		double sq2 = sqrt(A(tupleB) * A(tupleB) + B(tupleB) * B(tupleB));

		//first bisector
		a1 = A(tupleA) * sq2 - A(tupleB) * sq1;
		b1 = B(tupleA) * sq2 - B(tupleB) * sq1;
		c1 = C(tupleA) * sq2 - C(tupleB) * sq1;
		
		//second bisector
		a2 = A(tupleA) * sq2 + A(tupleB) * sq1;
		b2 = B(tupleA) * sq2 + B(tupleB) * sq1;
		c2 = C(tupleA) * sq2 + C(tupleB) * sq1;

		Point inter = Segment::intersect(segmentA, segmentB);

		double angle1 = Point::angle(inter, Point::distance(inter, segmentA.getPointA()) > Point::distance(inter, segmentA.getPointB()) ? segmentA.getPointA() : segmentA.getPointB());
		double angle2 = Point::angle(inter, Point::distance(inter, segmentB.getPointA()) > Point::distance(inter, segmentB.getPointB()) ? segmentB.getPointA() : segmentB.getPointB());
		double angleBe = (angle1 + angle2) / 2;

		//checking which bisector is the correct one
		if (fuzzyEquals(a1, 0))
		{
			Point inc = Point(inter.getX() + 1, inter.getY());
			double angleInc = Point::angle(inter, inc);
			//if (fuzzyEquals(distanceFromLine(inc, A(tupleA), B(tupleA), C(tupleA)), distanceFromLine(inc, A(tupleB), B(tupleB), C(tupleB))))
			if (fuzzyEquals(angleBe, angleInc) || fuzzyEquals(angleBe, angleInc - 180) || fuzzyEquals(angleBe, angleInc + 180))
				return std::tuple<double, double, double>(a1, b1, c1);
			else
				return std::tuple<double, double, double>(a2, b2, c2);
		}
		else if (fuzzyEquals(a2, 0))
		{
			Point inc = Point(inter.getX() + 1, inter.getY());
			double angleInc = Point::angle(inter, inc);
			//if (fuzzyEquals(distanceFromLine(inc, A(tupleA), B(tupleA), C(tupleA)), distanceFromLine(inc, A(tupleB), B(tupleB), C(tupleB))))
			if (fuzzyEquals(angleBe, angleInc) || fuzzyEquals(angleBe, angleInc - 180) || fuzzyEquals(angleBe, angleInc + 180))
				return std::tuple<double, double, double>(a2, b2, c2);
			else
				return std::tuple<double, double, double>(a1, b1, c1);
		}
		else
		{
			Point inc = Point(inter.getX() + 1 / a1, inter.getY() - 1 / b1);
			double dis1 = distanceFromLine(inc, tupleA);
			double dis2 = distanceFromLine(inc, tupleB);
			double angleInc = Point::angle(inter, inc);
			//if (fuzzyEquals(distanceFromLine(inc, tupleA), distanceFromLine(inc, tupleB)))
			if (fuzzyEquals(angleBe, angleInc) || fuzzyEquals(angleBe, angleInc - 180) || fuzzyEquals(angleBe, angleInc + 180))
				return std::tuple<double, double, double>(a1, b1, c1);
			else
				return std::tuple<double, double, double>(a2, b2, c2);
		}
	}

	//triple point center
	void circleCenter(const Point& a, const Point& b, const Point& c)
	{
		//all points are on one line
		if (collinearPoints(a, b, c))
			return;

		auto tuple1 = Segment::lineSegment(Segment(a, b));
		auto tuple2 = Segment::lineSegment(Segment(a, c));
		if (fuzzyEquals(A(tuple1), A(tuple2)) && fuzzyEquals(B(tuple1), B(tuple2)) && fuzzyEquals(C(tuple1), C(tuple2)))
			return;
		//powers
		double xa2 = a.getX() * a.getX();
		double ya2 = a.getY() * a.getY();
		double xb2 = b.getX() * b.getX();
		double yb2 = b.getY() * b.getY();
		double xc2 = c.getX() * c.getX();
		double yc2 = c.getY() * c.getY();
		//x and y of circle center
		double xs, ys;

		//X argument of center
		xs = ((xc2 - xa2) * (b.getY() - a.getY()) - ((xb2 - xa2) * (c.getY() - a.getY())) + (c.getY() - b.getY()) * (b.getY() - a.getY()) * (c.getY() - a.getY())) /
			(2 * ((c.getX() - a.getX()) * (b.getY() - a.getY()) - (b.getX() - a.getX()) * (c.getY() - a.getY())));

		//Y argument of center
		if (!(b.getY() == a.getY()))
			ys = (a.getX() - b.getX()) / (b.getY() - a.getY()) * (xs - (a.getX() + b.getX()) / 2) + (a.getY() + b.getY()) / 2;
		else
			ys = (b.getX() - c.getX()) / (c.getY() - b.getY()) * (xs - (b.getX() + c.getX()) / 2) + (b.getY() + c.getY()) / 2;

		Point center = Point(xs, ys);
		double radius = Point::distance(center, a);
		
		Circle circle = Circle(center, radius);

		output.push_back(circle);

		return;
	}

	//two segments and one point circle
	void circleCenter(const Segment& segmentA, const Segment& segmentB, const Point& p)
	{
		//atributes of Point p
		double q = p.getX();
		double r = p.getY();
		double m, n, o, delta;
		//lines are parallel

		if (Segment::containsPoint(p, segmentA) || Segment::containsPoint(p, segmentB))
			return;

		if (parallelSegments(segmentA, segmentB))
		{
			double d1 = distanceFromLine(p, Segment::lineSegment(segmentA));
			double d2 = distanceFromLine(p, Segment::lineSegment(segmentB));
			double pDistance = distanceFromParallelLines(segmentA, segmentB);
			if((d1 + d2) != distanceFromParallelLines(segmentA, segmentB))
				return;

			Point middle = Point((segmentA.getPointA().getX() + segmentB.getPointA().getX()) / 2, (segmentA.getPointA().getY() + segmentB.getPointA().getY()) / 2);
			auto lineA = Segment::lineSegment(segmentA);
			double k = -A(lineA) * middle.getX() - B(lineA) * middle.getY();
			double h = A(lineA);
			double i = B(lineA);
			std::tuple<double, double, double> middleLine = { h, i, k };
			double distance = distanceFromLine(segmentA.getPointA(), middleLine);
			double a = A(middleLine);
			double b = B(middleLine);
			double c = C(middleLine);
			//middle line is paralle to OX
			if (b == 0)
			{
				double w = -c / a;
				m = 1;
				n = -2 * r;
				o = r * r + (c / a + q) *(c / a + q) - distance * distance;

				delta = sqrt(n * n - 4 * m * o);

				double v1 = (-n - delta) / (2 * m);
				double v2 = (-n + delta) / (2 * m);

				Point center1 = Point(w, v1);
				double radius1 = distanceFromLine(center1, Segment::lineSegment(segmentA));

				Circle circle1 = Circle(center1, radius1);

				output.push_back(circle1);

				Point center2 = Point(w, v2);
				double radius2 = distanceFromLine(center2, Segment::lineSegment(segmentA));

				Circle circle2 = Circle(center2, radius2);

				output.push_back(circle2);

				return;
			}

			m = ((a * a) / (b * b)) + 1;
			n = 2 * ((c / b + r) * a / b - q);
			o = q * q + (c / b + r)* (c / b + r) - distance * distance;

			delta = sqrt(n * n - 4 * m * o);

			double w1 = (-n - delta) / (2 * m);
			double w2 = (-n + delta) / (2 * m);
			double v1 = -a / b * w1 - c / b;
			double v2 = -a / b * w2 - c / b;

			Point center1 = Point(w1, v1);
			double radius1 = distanceFromLine(center1, Segment::lineSegment(segmentA));

			Circle circle1 = Circle(center1, radius1);

			output.push_back(circle1);

			Point center2 = Point(w2, v2);
			double radius2 = distanceFromLine(center2, Segment::lineSegment(segmentA));

			Circle circle2 = Circle(center2, radius2);

			output.push_back(circle2);

			return;
		}

		//bisector of two straights that contain segment lineA and segment lineB
		auto tupleBis = bisector(segmentA, segmentB);
		double a = A(tupleBis);
		double b = B(tupleBis);
		double c = C(tupleBis);
		//point of intersection of straights that contain segment lineA and segment lineB
		Point inter = Segment::intersect(segmentA, segmentB);

		auto lineA = Segment::lineSegment(segmentA);
		double d = A(lineA);
		double e = B(lineA);
		double f = C(lineA);

		double sq = d * d + e * e;

		if (fuzzyEquals(b, 0))
		{
			double w = -c / a;
			
			m = e * e - sq;
			n = 2 * ((f - d * c / a)*e + r * sq);
			o = (f - d * c / a) * (f - d * c / a) - sq * (r*r + (q + c / a) * (q + c / a));

			delta = sqrt(n * n - 4 * m * o);
			double v1 = (-n - delta) / (2 * m);
			double v2 = (-n + delta) / (2 * m);

			Point center1 = Point(w, v1);
			double radius1 = distanceFromLine(center1, Segment::lineSegment(segmentA));

			Circle circle1 = Circle(center1, radius1);

			output.push_back(circle1);

			Point center2 = Point(w, v2);
			double radius2 = distanceFromLine(center2, Segment::lineSegment(segmentA));

			Circle circle2 = Circle(center2, radius2);

			output.push_back(circle2);

			return;
		}

		//simple mats
		
		m = (d - e * a / b) * (d - e * a / b) - sq * (1 + (a*a) / (b*b));
		n = 2 * ((f - e * c / b) * (d - e * a / b) - sq * ((a / b) * (c / b + r) - q));
		o = (f - e * c / b) * (f - e * c / b) - sq * (q*q + (c / b + r) * (c / b + r));

		delta = sqrt(n * n - 4 * m * o);
		double w1 = (-n - delta) / (2 * m);
		double w2 = (-n + delta) / (2 * m);
		double v1 = -a / b * w1 - c / b;
		double v2 = -a / b * w2 - c / b;

		Point center1 = Point(w1, v1);
		double radius1 = distanceFromLine(center1, Segment::lineSegment(segmentA));

		Circle circle1 = Circle(center1, radius1);

		output.push_back(circle1);

		Point center2 = Point(w2, v2);
		double radius2 = distanceFromLine(center2, Segment::lineSegment(segmentA));

		Circle circle2 = Circle(center2, radius2);

		output.push_back(circle2);
	}

	//triple Segment circle
	void circleCenter(const Segment& segmentA, const Segment& segmentB, const Segment& segmentC)
	{
		double a1, a2, b1, b2, c1, c2;
		//all three are parallel
		if (parallelSegments(segmentA, segmentB) && parallelSegments(segmentA, segmentC))
		{
			return;
		}

		Point inter;
		//one pair is parallel
		if (parallelSegments(segmentA, segmentB))
		{
			auto bise1 = bisector(segmentA, segmentC);
			auto bise2 = bisector(segmentC, segmentB);
		
			Point center = intersectionBLines(bise1, bise2);
			double radius = distanceFromLine(center, Segment::lineSegment(segmentA));

			Circle circle = Circle(center, radius);

			output.push_back(circle);

			return;
		}
		if (parallelSegments(segmentA, segmentC))
		{
			auto bise1 = bisector(segmentA, segmentB);
			auto bise2 = bisector(segmentB, segmentC);

			Point center = intersectionBLines(bise1, bise2);
			double radius = distanceFromLine(center, Segment::lineSegment(segmentA));

			Circle circle = Circle(center, radius);

			output.push_back(circle);

			return;
		}
		if (parallelSegments(segmentB, segmentC))
		{
			auto bise1 = bisector(segmentA, segmentB);
			auto bise2 = bisector(segmentC, segmentA);

			Point center = intersectionBLines(bise1, bise2);
			double radius = distanceFromLine(center, Segment::lineSegment(segmentA));

			Circle circle = Circle(center, radius);

			output.push_back(circle);

			return;
		}
		//tabs for bisectors of segments
		auto bise1 = bisector(segmentA, segmentB);
		auto bise2 = bisector(segmentA, segmentC);

		//point of intersection of two bisectors

		//return center of circle
		Point center = intersectionBLines(bise1, bise2);
		double radius = distanceFromLine(center, Segment::lineSegment(segmentA));

		Circle circle = Circle(center, radius);

		output.push_back(circle);
	}

	//double Point and one Segment circle
	void circleCenter(const Point& p, const Point& q, const Segment& line)
	{
		if (Segment::containsPoint(p, line) || Segment::containsPoint(q, line))
			return;

		//all elements are collinear
		if (collinearPoints(p, q, line.getPointA()) && collinearPoints(p, q, line.getPointB()))
			return;

		auto straight = Segment::lineSegment(line);
		double a = A(straight);
		double b = B(straight);
		double c = C(straight);

		//points are not on the plane cut by the segment
		if ((a*p.getX() + b * p.getY() + c) * (a*q.getX() + b * q.getY() + c) <= 0)
			return;

		auto mid = Segment::lineSegment(Segment(p, q));
		double e = A(mid);
		Point middle = Point(((p.getX() + q.getX()) / 2), ((p.getY() + q.getY()) / 2));
		double r = middle.getX();
		double s = middle.getY();

		double f = B(mid) * r - A(mid) * s;

		auto sym = std::make_tuple(-B(mid), A(mid), f);
		double d = A(sym);

		double sq = a * a + b * b;
		double x, y, z, delta;

		double m = p.getX();
		double n = p.getY();

		if (parallelSegments(Segment(p, q), line))
		{
			Point inter = intersectionBLines(straight, sym);
			
			circleCenter(p, q, inter);
			return;
		}

		if (e == 0)
		{
			double w = -f / d;
			x = b * b - sq;
			y = 2 * ((c - a * f / d)*b + sq * n);
			z = (c - a * f / d)*(c - a * f / d) - sq * (n*n + (f / d + m)*(f / d + m));

			delta = sqrt(y * y - 4 * x * z);
			double v1 = (-y - delta) / (2 * x);
			double v2 = (-y + delta) / (2 * x);

			Point center1 = Point(w, v1);
			double radius1 = Point::distance(center1, p);

			Circle circle1 = Circle(center1, radius1);

			output.push_back(circle1);

			Point center2 = Point(w, v2);
			double radius2 = Point::distance(center2, p);

			Circle circle2 = Circle(center2, radius2);
			
			output.push_back(circle2);

			return;
		}

		x = (a - d * b / e) * (a - d * b / e) - sq * (1 + (d*d) / (e*e));
		y = 2 * ((a - d * b / e) * (c - f * b / e) - sq * (d / e * (f / e + n) - m));
		z = (c - f * b / e)*(c - f * b / e) - sq * (m*m + (f / e + n)*(f / e + n));
		delta = sqrt(y * y - 4 * x * z);
		double w1 = (-y - delta) / (2 * x);
		double w2 = (-y + delta) / (2 * x);
		double v1 = -d / e * w1 - f / e;
		double v2 = -d / e * w2 - f / e;

		Point center1 = Point(w1, v1);
		double radius1 = Point::distance(center1, p);

		Circle circle1 = Circle(center1, radius1);

		output.push_back(circle1);

		Point center2 = Point(w2, v2);
		double radius2 = Point::distance(center2, p);

		Circle circle2 = Circle(center2, radius2);

		output.push_back(circle2);
	}
};