#pragma once
#include "stafx.h"

class AlgAutoSurfFPIA
{
public:
	AlgAutoSurfFPIA(const Handle(Geom_BSplineSurface)& surface, int par) : m_surface(surface), m_r(par) {};
	bool Fairing(int fairingNum,const std::vector<double> weights);
	Handle(Geom_BSplineSurface) getResult() { return m_surface; }

private:
	void sortByEnergy();


private:
	Handle(Geom_BSplineSurface) m_surface;
	std::vector<size_t> m_orderedIndex;
	std::vector<std::vector<double>> m_parBasisFun;
	int m_r{ 2 };
};

