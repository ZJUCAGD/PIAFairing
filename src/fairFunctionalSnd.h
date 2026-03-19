#pragma once
#include "stafx.h"

//用于计算曲面的能量项
class fairFunctionalSnd
{
public:
	fairFunctionalSnd(const TColStd_Array1OfReal& UKnots, 
					const TColStd_Array1OfReal& VKnots, 
					const TColStd_Array1OfInteger& UMults, 
					const TColStd_Array1OfInteger& VMults, 
					int uDegree, 
					int vDegree);

	double getValue(int i1, int i2);
	void getAllValue(std::vector<std::vector<double>>& value);

private:
	double func(double u, double v, int i1, int i2);
	double calValue(int i1, int i2);
	double integration(double a, double b, double c, double d, int i1, int i2);
	std::vector<double> derOneBasisFun(int r, int i, double u, bool isU);
	void calParamValue(int row, double u , bool isU);
	double packageKnot(int indexU, int indexV, int i1, int i2);
	double packageParam(int indexU, int indexV, int i1, int i2);

private:
	int m_uDegree;
	int m_vDegree;

	int m_rowNum;
	int m_colNum;

	std::vector<double> m_uKnot;
	std::vector<double> m_vKnot;

	TColStd_Array1OfReal m_noRepeatUKnotVec;
	TColStd_Array1OfReal m_noRepeatVKnotVec;

	TColStd_Array1OfReal m_flatUKnotVec;
	TColStd_Array1OfReal m_flatVKnotVec;

	const std::vector<double> w = { 0.2369269, 0.4786287,0.5688889, 0.4786287, 0.2369269 };
	const std::vector<double> x = { -0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985 };
	
	std::vector<std::vector<double>> m_value;

	std::vector<std::vector<double>> m_derUValues;
	std::vector<std::vector<double>> m_derVValues;

	std::vector<std::vector<gp_Vec>> m_derUBasis;
	std::vector<std::vector<gp_Vec>> m_derVBasis;

};

