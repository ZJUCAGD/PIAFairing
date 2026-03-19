#pragma once
#include "fairFunctionalSnd.h"

class AlgSurfFPIAByAdjustingControlPoint
{
public:
	AlgSurfFPIAByAdjustingControlPoint() {}

	~AlgSurfFPIAByAdjustingControlPoint();

	//初始化
	bool Init(const Handle(Geom_BSplineSurface)& surface);

	//全局
	void execute(double maxTimes, double error);

	//局部
	void execute(const std::vector<size_t>& indexList, const std::vector<double>& weights, int maxTimes);

	//设置控制顶点的权重
	void setWeights(const std::vector<size_t>& indexList, const std::vector<double>& weights);

	Handle(Geom_BSplineSurface) getResult();

	//void equationCheck();

private:
	void calDiffVec(std::vector<std::vector<gp_Vec>>& diffVec);
	void calFairVec(std::vector<std::vector<gp_Vec>>& fairVec);
	void checkVec(const gp_Vec& aVec);
	void checkPnt(const gp_Pnt& aPnt);
	void iterative();
	void localIterative(const std::vector<size_t>& indexList);
	//void checkWeightInterval(Eigen::MatrixXd aMatrix);
	void calMuParam();
	double iterativeDiff(const TColgp_Array2OfPnt& lastPoles, const TColgp_Array2OfPnt& curPoles);

private:
	TColgp_Array2OfPnt m_originPoles;
	Handle(Geom_BSplineSurface) m_surface;
	std::vector<std::vector<double>> m_weights;
	std::vector<std::vector<double>> m_mus;
	//fairFunctional* m_uFairFunctional;
	//fairFunctional* m_vFairFunctional;

	std::vector<std::vector<double>> m_FairFuncValue; 

};

