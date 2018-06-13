#pragma once
#include "Point.h"
#include <utility>
#include <tuple>

class Segment 
{
private:
	Point a, b;

public:
	Segment() {}

	Segment(double x1, double x2, double y1, double y2);

	Segment(Point a, Point b);

	static std::tuple<double, double, double> lineSegment(Segment a);
	const Point& getPointA() const;
	const Point& getPointB() const;

	bool parallel(Segment a, Segment b);

	static Point intersect(Segment a, Segment b);
	static bool containsPoint(const Point& p, const Segment& s);
};