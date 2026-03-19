#include "AlgSurfaceFairingByEnergy.h"
#include "fairFunctionalSnd.h"
void AlgSurfaceFairingByEnergy::Init(const std::vector<std::vector<gp_Pnt>>& data, int uNum, int vNum, int uDegree, int vDegree)
{
	m_data = data;
	m_uPoleNum = uNum;
	m_vPoleNum = vNum;
	m_uDegree = uDegree;
	m_vDegree = vDegree;
	m_controlPoints = TColgp_Array2OfPnt(1, m_uPoleNum, 1, m_vPoleNum);

	setParameter();
	setKnots();
	
	auto uFairFunctional1 = new fairFunctionalSnd(m_uKnots, m_vKnots, m_uMults, m_vMults, m_uDegree, m_vDegree);
	uFairFunctional1->getAllValue(m_FairFuncValue);

}

void AlgSurfaceFairingByEnergy::Init(const Handle(Geom_BSplineSurface)& surface)
{
	m_surface = surface;
	m_uPoleNum = m_surface->NbUPoles();
	m_vPoleNum = m_surface->NbVPoles();

	m_uDegree = m_surface->UDegree();
	m_vDegree = m_surface->VDegree();

	TColgp_Array2OfPnt m_controlPoints = TColgp_Array2OfPnt(m_surface->Poles());
	
	TColStd_Array1OfReal m_uKnots = TColStd_Array1OfReal(m_surface->UKnots());
	TColStd_Array1OfReal m_vKnots = TColStd_Array1OfReal(m_surface->VKnots());

	TColStd_Array1OfInteger m_uMults = TColStd_Array1OfInteger(m_surface->UMultiplicities());
	TColStd_Array1OfInteger m_vMults = TColStd_Array1OfInteger(m_surface->VMultiplicities());

	auto uFairFunctional1 = new fairFunctionalSnd(m_uKnots, m_vKnots, m_uMults, m_vMults, m_uDegree, m_vDegree);
	uFairFunctional1->getAllValue(m_FairFuncValue);

}

void AlgSurfaceFairingByEnergy::execute(double weight)
{
	m_weight = weight;
	executeEquation();
	creatBSplineSurface();
}

void AlgSurfaceFairingByEnergy::getResult(Handle(Geom_BSplineSurface)& surface)
{
	surface = m_surface;
}

void AlgSurfaceFairingByEnergy::setParameter()
{
	int dataNum1 = m_data.size();
	int dataNum2 = m_data[0].size();
	m_uParameterArray = TColStd_Array2OfReal(1, dataNum1, 1, dataNum2);
	m_vParameterArray = TColStd_Array2OfReal(1, dataNum1, 1, dataNum2);


	//每行
	for (int i = 0; i < dataNum1; i++)
	{
		double TotalD = 0;
		double firstParam = 0;
		double lastParam = 1;
		double paramLength = lastParam - firstParam;
		//每列
		for (int j = 0; j < dataNum2 - 1; j++)
		{
			TotalD += m_data[i][j].Distance(m_data[i][j + 1]);
		}

		for (int j = 1; j <= dataNum2; j++)
		{
			if (j == 1) m_uParameterArray.SetValue(i + 1, j, firstParam);
			else if (j == dataNum2) m_uParameterArray.SetValue(i + 1, j, lastParam);
			else
			{
				double temp = m_data[i][j - 1].Distance(m_data[i][j - 2]) / TotalD * paramLength;
				m_uParameterArray.SetValue(i + 1, j, m_uParameterArray.Value(i + 1, j - 1) + temp);
			}
		}
	}

	//每列
	for (int i = 0; i < dataNum2; i++)
	{
		double TotalD = 0;
		double firstParam = 0;
		double lastParam = 1;
		double paramLength = lastParam - firstParam;
		//每行
		for (int j = 0; j < dataNum1 - 1; j++)
		{
			TotalD += m_data[j][i].Distance(m_data[j + 1][i]);
		}

		for (int j = 1; j <= dataNum1; j++)
		{
			if (j == 1) m_vParameterArray.SetValue(j, i + 1, firstParam);
			else if (j == dataNum1) m_vParameterArray.SetValue(j, i + 1, lastParam);
			else
			{
				double temp = m_data[j - 1][i].Distance(m_data[j - 2][i]) / TotalD * paramLength;
				m_vParameterArray.SetValue(j, i + 1, m_vParameterArray.Value(j - 1, i + 1) + temp);
			}
		}
	}
}

void AlgSurfaceFairingByEnergy::setKnots()
{
	int uKnotNum = m_controlPoints.UpperRow() - m_uDegree + 1;
	int vKnotNum = m_controlPoints.UpperCol() - m_vDegree + 1;

	m_uKnots = TColStd_Array1OfReal(1, uKnotNum);
	m_vKnots = TColStd_Array1OfReal(1, vKnotNum);

	m_uMults = TColStd_Array1OfInteger(1, uKnotNum);
	m_vMults = TColStd_Array1OfInteger(1, vKnotNum);

	double uSpace = 1.0 / (double)(uKnotNum - 1);
	double vSpace = 1.0 / (double)(vKnotNum - 1);
	for (int i = 0; i < uKnotNum; i++)
	{
		double temp = i * uSpace;
		m_uKnots.SetValue(i + 1, temp);
	}
	for (int j = 0; j < vKnotNum; j++)
	{
		double temp = j * vSpace;
		m_vKnots.SetValue(j + 1, temp);
	}

	m_uMults.Init(1);
	m_uMults.SetValue(1, m_uDegree + 1);
	m_uMults.SetValue(uKnotNum, m_uDegree + 1);

	m_vMults.Init(1);
	m_vMults.SetValue(1, m_vDegree + 1);
	m_vMults.SetValue(vKnotNum, m_vDegree + 1);


	std::vector<double> tempFlatUKnot, tempFlatVKnot;
	for (int i = 1; i <= m_uKnots.Size(); i++)
	{
		for (int j = 1; j <= m_uMults.Value(i); j++)
		{
			tempFlatUKnot.push_back(m_uKnots.Value(i));
		}
	}
	for (int i = 1; i <= m_vKnots.Size(); i++)
	{
		for (int j = 1; j <= m_vMults.Value(i); j++)
		{
			tempFlatVKnot.push_back(m_vKnots.Value(i));
		}
	}

	m_uFlatKnots.Resize(1, tempFlatUKnot.size(), false);
	m_vFlatKnots.Resize(1, tempFlatVKnot.size(), false);
	for (int i = 0; i < tempFlatUKnot.size(); i++)
	{
		m_uFlatKnots.SetValue(i + 1, tempFlatUKnot[i]);
	}
	for (int i = 0; i < tempFlatVKnot.size(); i++)
	{
		m_vFlatKnots.SetValue(i + 1, tempFlatVKnot[i]);
	}
}

void AlgSurfaceFairingByEnergy::executeEquation()
{
	int dataNum = m_data.size() * m_data[0].size();
	int poleNum = m_controlPoints.NbColumns() * m_controlPoints.NbRows();
	
	Eigen::MatrixXd aMatrixDr = Eigen::MatrixXd::Zero(poleNum, poleNum);
	for (int i = 0; i < poleNum; i++)
	{
		for (int j = 0; j < poleNum; j++)
		{
			aMatrixDr(i, j) = m_FairFuncValue[i][j];
		}
	}

	Eigen::MatrixXd aMatrixI = Eigen::MatrixXd::Identity(poleNum, poleNum);
	Eigen::MatrixXd aMatrixA = (1 - m_weight) * aMatrixI + m_weight * aMatrixDr;

	Eigen::VectorXd Qx(poleNum);
	Eigen::VectorXd Qy(poleNum);
	Eigen::VectorXd Qz(poleNum);

	for (int i = 0; i < m_controlPoints.NbRows(); i++)
	{
		for (int j = 0; j < m_controlPoints.NbColumns(); j++)
		{
			gp_Pnt point = m_controlPoints.Value(i + 1, j + 1);
			int row = i * m_controlPoints.NbColumns() + j;
			Qx(row) = point.X();
			Qy(row) = point.Y();
			Qz(row) = point.Z();
		}
	}

	Eigen::MatrixXd Bx = (1 - m_weight) * Qx;
	Eigen::MatrixXd By = (1 - m_weight) * Qy;
	Eigen::MatrixXd Bz = (1 - m_weight) * Qz;

	Eigen::VectorXd Px = aMatrixA.colPivHouseholderQr().solve(Bx);
	Eigen::VectorXd Py = aMatrixA.colPivHouseholderQr().solve(By);
	Eigen::VectorXd Pz = aMatrixA.colPivHouseholderQr().solve(Bz);

	for (int i = 1; i <= m_controlPoints.Size(); i++)
	{
		int row = (i - 1) / m_controlPoints.UpperCol();
		int col = (i - 1) % m_controlPoints.UpperCol();
		m_controlPoints.SetValue(row + 1, col + 1, gp_Pnt(Px(i - 1), Py(i - 1), Pz(i - 1)));
	}
}

void AlgSurfaceFairingByEnergy::creatBSplineSurface()
{
	try
	{
		m_surface = new Geom_BSplineSurface(m_controlPoints, m_uKnots, m_vKnots, m_uMults, m_vMults, m_uDegree, m_vDegree);
	}
	catch (...)
	{
		std::cout << "Generate surface failed.";
	}
}
