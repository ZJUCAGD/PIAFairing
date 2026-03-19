#include "AlgAutoSurfFPIA.h"
#include "fairFunctionalSnd.h"
#include "AlgSurfFPIAByAdjustingControlPoint.h"
#include "AlgEnergyCalculation.h"

bool AlgAutoSurfFPIA::Fairing(int fairingNum, const std::vector<double> weights)
{
	if (m_surface.IsNull()) return false;
	int uNum = m_surface->NbUPoles(); 
	int vNum = m_surface->NbVPoles();
	if (fairingNum > uNum * vNum) return false;

	sortByEnergy();

	AlgSurfFPIAByAdjustingControlPoint algPIA;
	algPIA.Init(m_surface);

	std::vector<size_t> tmpo = m_orderedIndex;
	//std::vector<double> weights(tmpo.size(), 0.1);
	algPIA.execute(tmpo, weights, 1000);

	m_surface = algPIA.getResult();
	if (m_surface.IsNull())
		std::cout << "surface is null\n";
	return true;
}

void AlgAutoSurfFPIA::sortByEnergy()
{
	if (m_parBasisFun.size() == 0)
	{
		//计算基函数导数积分
		fairFunctionalSnd alg(m_surface->UKnots(), m_surface->VKnots(), m_surface-> UMultiplicities(), m_surface->VMultiplicities(), m_surface->UDegree(), m_surface->VDegree());
		alg.getAllValue(m_parBasisFun);
	}

	//计算新的控制顶点
	const auto poles = m_surface->Poles();
	auto newPoles = poles;
	int num = poles.Size();
	int uNum = m_surface->NbUPoles();
	int vNum = m_surface->NbVPoles();
	for (int i = 1; i <= num; i++)
	{
		gp_XYZ newPole(0, 0, 0);
		for (int k = 1; k <= uNum; k++)
		{
			for (int q = 1; q <= vNum; q++)
			{
				int index = (k - 1) * vNum + q;
				if (index == i) continue;
				newPole += m_parBasisFun[i - 1][index - 1] * poles.Value(k,q).Coord();
			}
		}
		
		newPole /= (-m_parBasisFun[i - 1][i - 1]);
		newPoles.SetValue((i-1)/ vNum + 1, (i-1)%vNum+1, gp_Pnt(newPole));
	}

	//计算新控制顶点带给曲线的光顺提升
	AlgEnergyCalculation alg;
	double originEnergy = alg.calEnergy(m_surface, m_r);
	std::vector<double> adjustEnergy(num);
	for (int i = 1; i <= num; i++)
	{
		Handle(Geom_BSplineSurface) surface = new Geom_BSplineSurface(m_surface->Poles(), m_surface->UKnots(), m_surface->VKnots(),m_surface->UMultiplicities(), m_surface->VMultiplicities(), m_surface->UDegree(), m_surface->VDegree());
		surface->SetPole((i - 1) / vNum + 1, (i - 1) % vNum + 1, newPoles.Value((i - 1) / vNum + 1, (i - 1) % vNum + 1));
		AlgEnergyCalculation tmpAlg;
		adjustEnergy[i - 1] = originEnergy - (tmpAlg.calEnergy(surface, m_r));
	}

	//对控制顶点进行排序
	auto sortDescendingWithIndices = [](const std::vector<double>& arr) -> std::vector<size_t> {
		std::vector<size_t> indices(arr.size());
		std::iota(indices.begin(), indices.end(), 0);

		std::sort(indices.begin(), indices.end(),
			[&arr](size_t i, size_t j) { return arr[i] > arr[j]; });

		return indices;
		};

	m_orderedIndex = sortDescendingWithIndices(adjustEnergy);
}
