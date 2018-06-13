#include <cmath>
#include "Point.h"
#include "Helpers.h"


Point::Point(double x, double y)
{
	this->x = x;
	this->y = y;
}

double Point::getX() const
{
	return this->x;
}

double Point::getY() const
{
	return this->y;
}

void Point::setX(double x)
{
	this->x = x;
}

void Point::setY(double y)
{
	this->y = y;
}

double Point::distance(const Point& a, const Point& b)
{
	return sqrt((b.getX() - a.getX()) * (b.getX() - a.getX()) + (b.getY() - a.getY()) * (b.getY() - a.getY()));
}

double Point::angle(const Point & p1, const Point & p2)
{
	double dx = p2.getX() - p1.getX();
	double dy = p2.getY() - p1.getY();

	double theta = std::atan2(dy, dx) * 180 / (std::atan(1) * 4);

	double theta_normalized = theta < 0 ? theta + 360 : theta;
	theta_normalized = 360 - theta_normalized;

	if (fuzzyEquals(theta_normalized, 360))
		return 0;
	else
		return theta_normalized;
}
