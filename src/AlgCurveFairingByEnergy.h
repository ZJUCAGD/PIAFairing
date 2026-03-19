#pragma once
#include "fairFunctional.h"

class AlgCurveFairingByEnergy
{
public:
	AlgCurveFairingByEnergy() {}
	bool Init(const Handle(Geom_BSplineCurve)& curve, int energyType);
	void execute(double weight);
	Handle(Geom_BSplineCurve) getResult();


private:
	Handle(Geom_BSplineCurve) m_curve;
	fairFunctional* m_fairFunctional;
	void poleDistance(const TColgp_Array1OfPnt& poles1, const TColgp_Array1OfPnt& poles2);
};

