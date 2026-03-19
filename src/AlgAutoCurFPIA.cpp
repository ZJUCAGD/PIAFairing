#include "AlgAutoCurFPIA.h"
#include "fairFunctional.h"
#include "AlgEnergyCalculation.h"
#include "AlgCurFPIAByAdjustingControlPoint.h"

bool AlgAutoCurFPIA::Fairing(int fairingNum, const std::vector<double> weights)
{
	if (m_curve.IsNull()) return false;
	int num = m_curve->NbPoles();
	if (fairingNum > num) return false;
	sortByEnergy();

	//TODO:只动前fairingNum个控制顶点
	std::cout << "一共有" << num << "个控制顶点，其中以下" << fairingNum << "为首选：\n";
	for (int i = 0; i < fairingNum; i++)
	{
		std::cout << m_orderedIndex[i] << std::endl;
	}

	AlgCurFPIAByAdjustingControlPoint algPIA;
	algPIA.Init(m_curve, m_r);
	std::vector<size_t> tmpo(m_orderedIndex.begin(), m_orderedIndex.begin() + fairingNum);
	//std::vector<double> weights(fairingNum, 5e-5);
	
	algPIA.execute(tmpo, weights, 100);
	m_curve = algPIA.getResult();
	if(m_curve.IsNull())
		std::cout << "curve is null";
	return true;
}

void AlgAutoCurFPIA::sortByEnergy()
{
	if (m_parBasisFun.size() == 0)
	{
		//TODO:计算基函数导数积分
		fairFunctional alg(m_curve->Knots(), m_curve->Multiplicities(),m_curve->Degree(),m_r);
		alg.getAllValue(m_parBasisFun);
	}

	//TODO:计算新的控制顶点
	const auto poles = m_curve->Poles();
	auto newPoles = poles;
	int num = poles.Size();
	for (int i = 1; i <= num; i++)
	{
		gp_XYZ newPole(0,0,0);
		for (int j = 1; j <= num; j++)
		{
			if (j == i) continue;
			newPole += m_parBasisFun[i-1][j-1] * poles.Value(j).Coord();
		}
		newPole /= (-m_parBasisFun[i - 1][i - 1]);
		newPoles.SetValue(i, gp_Pnt(newPole));
	}

	//TODO:计算新控制顶点带给曲线的光顺提升
	AlgEnergyCalculation alg;
	double originEnergy = alg.calEnergy(m_curve, m_r);
	std::vector<double> adjustEnergy(num);
	for (int i = 1; i <= num; i++)
	{
		Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(m_curve->Poles(), m_curve->Knots(), m_curve->Multiplicities(), m_curve->Degree());
		curve->SetPole(i, newPoles.Value(i));
		AlgEnergyCalculation tmpAlg;
		adjustEnergy[i - 1] = originEnergy - (tmpAlg.calEnergy(curve, m_r));
	}

	//TODO:对控制顶点进行排序
	auto sortDescendingWithIndices = [](const std::vector<double>& arr) -> std::vector<size_t> {
		std::vector<size_t> indices(arr.size());
		std::iota(indices.begin(), indices.end(), 0);

		std::sort(indices.begin(), indices.end(),
			[&arr](size_t i, size_t j) { return arr[i] > arr[j]; });

		return indices;
		};

	m_orderedIndex = sortDescendingWithIndices(adjustEnergy);
}
