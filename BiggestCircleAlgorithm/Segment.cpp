#include "Segment.h"
#include "BruteGeometry.h"
#include <tuple>

Segment::Segment(double x1, double x2, double y1, double y2)
{
	this->a.setX(x1);
	this->a.setY(y1);
	this->b.setX(x2);
	this->b.setY(y2);
}

Segment::Segment(Point a, Point b)
{
	this->a = Point(a.getX(), a.getY());
	this->b = Point(b.getX(), b.getY());
}

std::tuple<double, double, double> Segment::lineSegment(Segment seg)
{
	double a = seg.getPointA().getY() - seg.getPointB().getY();
	double b = seg.getPointB().getX() - seg.getPointA().getX();
	double c = seg.getPointA().getX() * (seg.getPointB().getY() - seg.getPointA().getY()) + seg.getPointA().getY() * (seg.getPointA().getX() - seg.getPointB().getX());
	return std::tuple<double, double, double>(a, b, c);
}

const Point& Segment::getPointA() const
{
	return this->a;
}

const Point& Segment::getPointB() const
{
	return this->b;
}

bool Segment::parallel(Segment a, Segment b)
{
	auto t1 = lineSegment(a);
	auto t2 = lineSegment(b);
	return abs(A(t1) * B(t2) - A(t2) * B(t1)) < epsilon;
}

Point Segment::intersect(Segment a, Segment b)
{
	auto tupleA = lineSegment(a);
	auto tupleB = lineSegment(b);

	if (BruteGeometry::parallelSegments(a, b))
	{
		throw std::exception();
	}

	double a1, a2, b1, b2, c1, c2;
	a1 = A(tupleA);
	a2 = A(tupleB);
	b1 = B(tupleA);
	b2 = B(tupleB);
	c1 = C(tupleA);
	c2 = C(tupleB);

	double w = a1 * b2 - a2 * b1;
	double wx = -c1 * b2 + c2 * b1;
	double wy = -a1 * c2 + a2 * c1;

	return Point(wx/w, wy/w);
}

bool Segment::containsPoint(const Point & p, const Segment & s)
{
	return ((s.getPointA().getX() == p.getX() && s.getPointA().getY() == p.getY()) || (s.getPointB().getX() == p.getX() && s.getPointB().getY() == p.getY() ));
}
