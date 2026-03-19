#pragma once
#include "fairFunctional.h"

class AlgCurFPIAByAdjustingControlPoint
{
public:
	AlgCurFPIAByAdjustingControlPoint() {}

	~AlgCurFPIAByAdjustingControlPoint();

	//初始化
	bool Init(const Handle(Geom_BSplineCurve)& curve, int energyType);

	//全局
	void execute(int maxTimes, double error);

	//局部
	void execute(const std::vector<size_t>& indexList, const std::vector<double>& weights, int maxTimes);

	//设置控制顶点的权重
	void setWeights(const std::vector<size_t>& indexList, const std::vector<double>& weights);

	Handle(Geom_BSplineCurve) getResult();

	void equationCheck();

private:
	void calDiffVec(std::vector<gp_Vec>& diffVec);
	void calFairVec(std::vector<gp_Vec>& fairVec);
	void checkVec(const gp_Vec& aVec);
	void checkPnt(const gp_Pnt& aPnt);
	void iterative();
	void localIterative(const std::vector<size_t>& indexList);

	void checkWeightInterval(Eigen::MatrixXd aMatrix);
	void calMuParam();
	double iterativeDiff(const TColgp_Array1OfPnt& lastPoles, const TColgp_Array1OfPnt& curPoles);
	void poleDistance(const TColgp_Array1OfPnt& poles1, const TColgp_Array1OfPnt& poles2);
private:
	TColgp_Array1OfPnt m_originPoles;
	Handle(Geom_BSplineCurve) m_curve;
	int m_energyType{ 2 };
	std::vector<double> m_weights;
	std::vector<double> m_mus;
	fairFunctional* m_fairFunctional;
	
};

