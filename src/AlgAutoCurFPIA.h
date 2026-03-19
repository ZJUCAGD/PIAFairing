#pragma once
#include "stafx.h"

class AlgAutoCurFPIA
{
public:
	AlgAutoCurFPIA(const Handle(Geom_BSplineCurve)& curve, int par) : m_curve(curve), m_r(par){};
	bool Fairing(int fairingNum, const std::vector<double> weights);
	Handle(Geom_BSplineCurve) getResult() { return m_curve; }

private:
	void sortByEnergy();


private:
	Handle(Geom_BSplineCurve) m_curve;
	std::vector<size_t> m_orderedIndex;
	std::vector<std::vector<double>> m_parBasisFun;
	int m_r{ 2 };
};

