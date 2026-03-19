#pragma once
#include "stafx.h"
//B样条曲线拟合数据点
//通过最小二乘方法
class AlgPointFitCurve
{
public:
	AlgPointFitCurve() {}
	bool Init(const std::vector<gp_Pnt>& data, int poleNum, int degree);
	Handle(Geom_BSplineCurve) getResult();

private:
	bool averageVec(int num, int firstNum, int endNum, std::vector<int>& avegVec);
	bool initPoles(int poleNum);
	void setParameter();
	void setKnots();
	void creatBSplineCurve();
	void executeEquation();

	int FindSpan(double u);
	void BasisFuns(int i, double u, std::vector<double>& N);
	double oneBasisFun(int i, double u);

private:
	Handle(Geom_BSplineCurve) myCurve;
	TColgp_Array1OfPnt myPoleArray;
	TColStd_Array1OfReal myParameterArray;
	TColStd_Array1OfReal myKnotArray;
	TColStd_Array1OfInteger myMultArray;
	Standard_Integer myDegree{ 3 };

	std::vector<gp_Pnt> myOriginFitData;
	std::vector<double> myKnotVec;
	std::vector<int> myInitialPoleIndex;
};

