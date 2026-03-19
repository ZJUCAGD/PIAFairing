#pragma once
#include "stafx.h"
//B样条曲面拟合数据点
//通过最小二乘方法

class AlgPointFitSurface
{
public:
	
	AlgPointFitSurface() {}
	bool Init(const std::vector<std::vector<gp_Pnt>>& data, int poleNum1, int poleNum2, int uDegree, int vDegree);
	Handle(Geom_BSplineSurface) getResult();






private:
	void setParameter();
	void setKnots();
	void executeEquation();
	void creatBSplineCurve();


private:
	std::vector<std::vector<gp_Pnt>> myOriginFitData;
	TColStd_Array2OfReal myUParameterArray;
	TColStd_Array2OfReal myVParameterArray;
	
	Handle(Geom_BSplineSurface) mySurface;

	TColgp_Array2OfPnt myControlPoints;

	Standard_Integer myUDegree{ 3 };
	Standard_Integer myVDegree{ 3 };

	TColStd_Array1OfReal myUKnots;
	TColStd_Array1OfReal myVKnots;
	TColStd_Array1OfReal myUFlatKnots;
	TColStd_Array1OfReal myVFlatKnots;

	TColStd_Array1OfInteger myUMults;
	TColStd_Array1OfInteger myVMults;
};

