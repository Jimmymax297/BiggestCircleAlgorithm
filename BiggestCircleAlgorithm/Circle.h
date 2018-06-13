#pragma once
#include "Point.h"
#include "Segment.h"

class Circle
{
private:
	Point center;
	double radius;
public:
	Circle(Point center, double radius);
	const Point& getCenter() const;
	const double& getRadius() const;
	void setCenter(const Point& p);
	void setRadius(double radius);
	//boolean isValid(std::list<Segment> Segments);
};