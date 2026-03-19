#pragma once
#include "stafx.h"
class AlgEnergyCalculation
{
public:
	AlgEnergyCalculation() {}
	double calEnergy(const Handle(Geom_BSplineCurve)& curve, int type = 2);
	double calEnergy(const Handle(Geom_BSplineSurface)& surface, int type = 2);


private:
	double BsplineSecondDerivative(double u);
	double BsplineSecondDerivative(double u, double v);


	Handle(Geom_BSplineCurve) m_curve;
	Handle(Geom_BSplineSurface) m_surface;

	int m_type{ 2 };

};

