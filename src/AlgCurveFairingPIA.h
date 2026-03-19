#pragma once
#include "fairFunctional.h"


class AlgCurveFairingPIA
{
public:
	AlgCurveFairingPIA() {}

	//无初始曲线
	bool Init(const std::vector<gp_Pnt>& data, int poleNum, int degree, int der = 2);

	//有初始曲线
	void Init(const std::vector<gp_Pnt>& data, const Handle(Geom_BSplineCurve)& aCurve, int der = 2);

	//全局Fairing-PIA
	void executeGobalPoles(double maxTimes, double error);

	//局部Fairing-PIA，仅调整部分控制顶点
	void executeLocalPoles(std::vector<int> indexList, double maxTimes, double error);

	//局部Fairing-PIA，仅调整部分数据点对应的曲线
	void executeLocalPoles(int dataIndexBegin, int dataIndexEnd, double maxTimes, double error);

	//用解方程的方式去做
	void executeEquation();


	//设置控制顶点的权重
	void setWeights(std::vector<int> indexList, std::vector<double> weights);

	Handle(Geom_BSplineCurve) getResult();

private:
	bool averageVec(int num, int firstNum, int endNum, std::vector<int>& avegVec);
	bool initPoles(int poleNum);
	void setParameter();
	void setKnots();
	void creatBSplineCurve();
	int FindSpan(double u);
	void BasisFuns(int i, double u, std::vector<double>& N);
	void updateGobalPoles();
	void updateLocalPoles();
	double getFittingError();
	double oneBasisFun(int i, double u);
	void generateParamMu();

private:

	std::vector<gp_Pnt> myOriginFitData;
	Handle(Geom_BSplineCurve) myCurve;
	TColgp_Array1OfPnt myPoleArray;
	TColStd_Array1OfReal myParameterArray;
	TColStd_Array1OfReal myKnotArray;
	TColStd_Array1OfReal myFlatKnots;
	TColStd_Array1OfInteger myMultArray;
	Standard_Integer myDegree{ 3 };
	std::vector<double> myKnotVec;
	std::vector<int> myInitialPoleIndex;
	int myDerivative{ 2 };
	std::vector<double> myWeights;
	std::vector<int> myAdjustPoleIndex;
	std::vector<double> muParamMu;
};

