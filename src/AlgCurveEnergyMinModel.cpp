#include "AlgCurveEnergyMinModel.h"

bool AlgCurveEnergyMinModel::Init(const std::vector<gp_Pnt>& data, int poleNum, int degree, int der)
{
	myOriginFitData = data;
	myDegree = degree;
	myDerivative = der;

	if (!initPoles(poleNum))
	{
		return false;
	}

	setParameter();
	setKnots();
	creatBSplineCurve();

	return true;
}

bool AlgCurveEnergyMinModel::Execute()
{
	generateMatrixN();
	generateMatrixDr();
	generateMatrixA();

	Eigen::VectorXd Qx(myOriginFitData.size());
	Eigen::VectorXd Qy(myOriginFitData.size());
	Eigen::VectorXd Qz(myOriginFitData.size());

	for (int i = 0; i < myOriginFitData.size(); i++)
	{
		gp_Pnt point = myOriginFitData[i];
		Qx(i) = point.X();
		Qy(i) = point.Y();
		Qz(i) = point.Z();
	}
	Eigen::MatrixXd B = (1 - myOmega) * m_matrixN.transpose();

	Eigen::MatrixXd Bx = B * Qx;
	Eigen::MatrixXd By = B * Qy;
	Eigen::MatrixXd Bz = B * Qz;

	

	Eigen::VectorXd Px = m_matrixA.colPivHouseholderQr().solve(Bx);
	Eigen::VectorXd Py = m_matrixA.colPivHouseholderQr().solve(By);
	Eigen::VectorXd Pz = m_matrixA.colPivHouseholderQr().solve(Bz);

	//std::cout << "bx\n";
	//std::cout << Pz << std::endl;


	for (int i = 1; i <= myPoleArray.Size(); i++)
	{
		myPoleArray.SetValue(i, gp_Pnt(Px(i - 1), Py(i - 1), Pz(i - 1)));
	}
	creatBSplineCurve();
}

Handle(Geom_BSplineCurve) AlgCurveEnergyMinModel::getResult()
{
	return myCurve;
}

bool AlgCurveEnergyMinModel::averageVec(int num, int firstNum, int endNum, std::vector<int>& avegVec)
{
	if (endNum - firstNum < num)
	{
		return false;
	}
	double stp = (double)(endNum - firstNum) / (double)(num - 1);
	avegVec.push_back(firstNum);
	for (int i = 1; i <= num - 2; i++)
	{
		avegVec.push_back((int)(i * stp));
	}
	avegVec.push_back(endNum);
	return true;
}

bool AlgCurveEnergyMinModel::initPoles(int poleNum)
{
	int dataNum = myOriginFitData.size();

	if (poleNum > dataNum)
	{
		std::cout << "The number of poles is bigger than that of data points.\nPlease input a number belong to (1, " << dataNum << ")." << std::endl;
		return false;
	}

	if (poleNum < 1)
	{
		std::cout << "The number of poles is too small to generate the PIA curve.\nPlease input a number belong to (1, " << dataNum << ")." << std::endl;
		return false;
	}

	std::vector<int> vec;
	if (!averageVec(poleNum, 1, dataNum, vec))
	{
		return false;
	}

	myPoleArray = TColgp_Array1OfPnt(1, poleNum);
	for (int i = 0; i < vec.size(); i++)
	{
		//std::cout << vec[i] << std::endl;
		myPoleArray.SetValue(i + 1, myOriginFitData[vec[i] - 1]);
	}

	myInitialPoleIndex = vec;

	return true;
}

void AlgCurveEnergyMinModel::setParameter()
{
	int dataNum = myOriginFitData.size();
	myParameterArray = TColStd_Array1OfReal(1, dataNum);

	double TotalD = 0;
	for (int i = 0; i < dataNum - 1; i++)
	{
		TotalD += myOriginFitData[i].Distance(myOriginFitData[i + 1]);
	}

	for (int i = 1; i <= dataNum; i++)
	{
		if (i == 1) myParameterArray.SetValue(i, 0);
		else if (i == dataNum) myParameterArray.SetValue(i, 1);
		else
		{
			double temp = myOriginFitData[i - 1].Distance(myOriginFitData[i - 2]) / TotalD;
			myParameterArray.SetValue(i, myParameterArray.Value(i - 1) + temp);
		}
		//std::cout << m_parameterArray.Value(i) << std::endl;
	}

}

void AlgCurveEnergyMinModel::setKnots()
{
	int knotNum = myPoleArray.Size() - myDegree + 1;//m_Poles.Size() - m_degree - 1+2
	myKnotArray = TColStd_Array1OfReal(1, knotNum);
	myMultArray = TColStd_Array1OfInteger(1, knotNum);

	myKnotVec.resize(myDegree, 0);
	for (int i = 1; i <= knotNum; i++)
	{
		if (i == 1) myKnotArray.SetValue(i, 0);
		else if (i == knotNum) myKnotArray.SetValue(i, 1);
		else
		{
			double temp = 0;
			for (int j = i; j < i + myDegree; j++)
			{
				temp += myParameterArray.Value(myInitialPoleIndex[j - 1]);
			}
			myKnotArray.SetValue(i, temp / (double)myDegree);
		}
		myKnotVec.push_back(myKnotArray.Value(i));
	}

	for (int i = 0; i < myDegree; i++)
	{
		myKnotVec.push_back(1);
	}

	myMultArray.Init(1);
	myMultArray.SetValue(1, myDegree + 1);
	myMultArray.SetValue(knotNum, myDegree + 1);
}

void AlgCurveEnergyMinModel::creatBSplineCurve()
{
	try
	{
		myCurve = new Geom_BSplineCurve(myPoleArray, myKnotArray, myMultArray, myDegree);
	}
	catch (...)
	{
		std::cout << "Generate curve failed.";
	}
}

int AlgCurveEnergyMinModel::FindSpan(double u)
{
	std::vector<double> U = myKnotVec;
	int nu = U.size() - myDegree - 2;
	int p = myDegree;
	if (u == U[nu + 1]) return (nu);
	if (u > U[nu + 1]) return -1;
	if (u < U[0]) return -1;
	int low = p, high = nu + 1;
	int mid = (low + high) / 2;
	while (u < U[mid] || u >= U[mid + 1])
	{
		if (u < U[mid])
		{
			high = mid;
		}
		else
		{
			low = mid;
		}
		mid = (low + high) / 2;
	}
	return (mid);
}

void AlgCurveEnergyMinModel::BasisFuns(int i, double u, std::vector<double>& N)
{
	int p = myDegree;
	std::vector<double> U = myKnotVec;
	N.push_back(1.0);
	std::vector<double> left(p + 1, 0);
	std::vector<double> right(p + 1, 0);

	for (int j = 1; j <= p; j++)
	{
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		double saved = 0.0;
		for (int r = 0; r < j; r++)
		{
			double temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N.push_back(saved);
	}
}

void AlgCurveEnergyMinModel::generateMatrixN()
{
	m_matrixN = Eigen::MatrixXd::Zero(myParameterArray.Size(), myPoleArray.Size());
	for (int i = 0; i < myParameterArray.Size(); i++)
	{
		double temp = myParameterArray.Value(i + 1);
		int index = FindSpan(temp);
		std::vector<double> nonZeroBasis;
		BasisFuns(index, temp, nonZeroBasis);
		for (int j = 0; j < nonZeroBasis.size(); j++)
		{
			int basisIndex = index + j - myDegree;
			m_matrixN(i, basisIndex) = nonZeroBasis[j];
		}
	}
}

void AlgCurveEnergyMinModel::generateMatrixDr()
{
	m_matrixDr = Eigen::MatrixXd::Zero(myPoleArray.Size(), myPoleArray.Size());
	fairFunctional IntegEngine(myKnotVec, myDegree);

	for (int i = 0; i < myPoleArray.Size(); i++)
	{
		for (int j = 0; j < myPoleArray.Size(); j++)
		{
			m_matrixDr(i, j) = IntegEngine.getValue(myDerivative, i, j);
		}
	}
}

//void AlgCurveEnergyMinModel::generateMatrixOmega()
//{
//	m_matrixOmega = Eigen::MatrixXd::Zero(myPoleArray.Size(), myPoleArray.Size());
//	for (int i = 0; i < myPoleArray.Size(); i++)
//	{
//		m_matrixOmega(i, i) = 1e-6;
//	}
//}

void AlgCurveEnergyMinModel::generateMatrixA()
{
	//Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(myPoleArray.Size(), myPoleArray.Size());
	m_matrixA = (1-myOmega) * m_matrixN.transpose() * m_matrixN + myOmega * m_matrixDr;
}

//void AlgCurveEnergyMinModel::generateMatrixMu()
//{
//	Eigen::MatrixXd Mu = Eigen::MatrixXd::Zero(myPoleArray.Size(), myPoleArray.Size());
//	for (int i = 0; i < myPoleArray.Size(); i++)
//	{
//		double sum = 0;
//		for (int j = 0; j < myPoleArray.Size(); j++)
//		{
//			sum += abs(m_matrixA(i, j));
//		}
//		Mu(i, i) = 1.0 / sum;
//	}
//
//	Eigen::MatrixXd final = Eigen::MatrixXd::Identity(myPoleArray.Size(), myPoleArray.Size()) - Mu * m_matrixA;
//	Eigen::EigenSolver<Eigen::MatrixXd> es(final);
//}
