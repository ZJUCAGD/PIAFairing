#pragma once
#include "stafx.h"

class fairFunctional
{
public:
	fairFunctional(const std::vector<double>& knotVec, int degree, int r);
	fairFunctional(const TColStd_Array1OfReal& knots, const TColStd_Array1OfInteger& mults, int degree, int r);
	double getValue(int i1, int i2);
	void getAllValue(std::vector<std::vector<double>>& value);

private:
	double func(int r, int i1, int i2, double u);
	double derOneBasisFun(int r, int i, double u);
	double integration(double a, double b, int r, int i1, int i2);
	double calValue(int r, int i1, int i2);
private:
	int m_degree;
	std::vector<double> m_knot;
	std::vector<double> m_noRepeatKnotVec;
	const std::vector<double> w = { 0.2369269, 0.4786287,0.5688889, 0.4786287, 0.2369269 };
	const std::vector<double> x = { -0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985 };
	std::vector<std::vector<double>> m_value;
};
