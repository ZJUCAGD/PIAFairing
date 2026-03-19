#pragma once
#include "fairFunctional.h"

class AlgCurveIterativeFairing
{
public:
	AlgCurveIterativeFairing(){}
	bool Init(const Handle(Geom_BSplineCurve)& theCurve, const std::vector<double>& distances);
	void execute();
	void getResult(Handle(Geom_BSplineCurve)& aCurve);

private:
	void getDiffVec(std::vector<gp_Vec>& aVec);
	void getFairVec(std::vector<gp_Vec>& aVec);
	void Iterative();
private:
	Handle(Geom_BSplineCurve) myCurve;
	std::vector<double> myDistances;
	TColgp_Array1OfPnt myOriginPoles;
	std::vector<double> myFlatKnots;
	std::vector<double> myLambdas;
	std::vector<double> myVis;


};

