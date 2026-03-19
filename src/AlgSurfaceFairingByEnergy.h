#pragma once
#include "stafx.h"
class AlgSurfaceFairingByEnergy
{
public:
	AlgSurfaceFairingByEnergy() {};
	void Init(const std::vector<std::vector<gp_Pnt>>& data, int uNum, int vNum, int uDegree, int vDegree);
	void Init(const Handle(Geom_BSplineSurface)& surface/*, int uNum, int vNum, int uDegree, int vDegree*/);
	void execute(double weight);
	void getResult(Handle(Geom_BSplineSurface)& surface);

private:
	void setParameter();
	void setKnots();
	void executeEquation();
	void creatBSplineSurface();

private:
	Handle(Geom_BSplineSurface) m_surface;
	std::vector<std::vector<gp_Pnt>> m_data;
	int m_uPoleNum{ 2 };
	int m_vPoleNum{ 2 };
	int m_uDegree{ 3 };
	int m_vDegree{ 3 };
	TColgp_Array2OfPnt m_controlPoints;
	TColStd_Array2OfReal m_uParameterArray;
	TColStd_Array2OfReal m_vParameterArray;

	TColStd_Array1OfReal m_uKnots;
	TColStd_Array1OfReal m_vKnots;
	TColStd_Array1OfReal m_uFlatKnots;
	TColStd_Array1OfReal m_vFlatKnots;

	TColStd_Array1OfInteger m_uMults;
	TColStd_Array1OfInteger m_vMults;

	double m_weight{ 1e-7 };
	std::vector<std::vector<double>> m_FairFuncValue;

	

};

