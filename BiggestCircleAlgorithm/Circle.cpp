#include "Circle.h"


Circle::Circle(Point center = Point(0, 0), double radius = 0) : center(center)
{
	this->radius = radius < 0 ? 0 : radius;
}

const Point & Circle::getCenter() const
{
	return this->center;
}

const double & Circle::getRadius() const
{
	return this->radius;
}

void Circle::setCenter(const Point & p)
{
	this->center = p;
}

void Circle::setRadius(double radius)
{
	if (radius < 0) return;
	this->radius = radius;
}
