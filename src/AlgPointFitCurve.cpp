#include "AlgPointFitCurve.h"

bool AlgPointFitCurve::Init(const std::vector<gp_Pnt>& data, int poleNum, int degree)
{
	myOriginFitData = data;
	myDegree = degree;

	if (!initPoles(poleNum))
	{
		return false;
	}

	setParameter();
	setKnots();
	executeEquation();
	creatBSplineCurve();
}

Handle(Geom_BSplineCurve) AlgPointFitCurve::getResult()
{
	return myCurve;
}

void AlgPointFitCurve::creatBSplineCurve()
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

void AlgPointFitCurve::executeEquation()
{
	//N
	Eigen::MatrixXd aMatrixN = Eigen::MatrixXd::Zero(myParameterArray.Size(), myPoleArray.Size());
	for (int i = 0; i < myParameterArray.Size(); i++)
	{
		double temp = myParameterArray.Value(i + 1);
		int index = FindSpan(temp);
		std::vector<double> nonZeroBasis;
		BasisFuns(index, temp, nonZeroBasis);
		for (int j = 0; j < nonZeroBasis.size(); j++)
		{
			int basisIndex = index + j - myDegree;
			aMatrixN(i, basisIndex) = nonZeroBasis[j];
		}
	}

	//A
	Eigen::MatrixXd aMatrixA = aMatrixN.transpose() * aMatrixN;

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
	Eigen::MatrixXd B = aMatrixN.transpose();

	Eigen::MatrixXd Bx = B * Qx;
	Eigen::MatrixXd By = B * Qy;
	Eigen::MatrixXd Bz = B * Qz;

	Eigen::VectorXd Px = aMatrixA.colPivHouseholderQr().solve(Bx);
	Eigen::VectorXd Py = aMatrixA.colPivHouseholderQr().solve(By);
	Eigen::VectorXd Pz = aMatrixA.colPivHouseholderQr().solve(Bz);

	for (int i = 1; i <= myPoleArray.Size(); i++)
	{
		myPoleArray.SetValue(i, gp_Pnt(Px(i - 1), Py(i - 1), Pz(i - 1)));
	}
	
}

int AlgPointFitCurve::FindSpan(double u)
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

void AlgPointFitCurve::BasisFuns(int i, double u, std::vector<double>& N)
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

double AlgPointFitCurve::oneBasisFun(int i, double u)
{
	std::vector<double> U = myKnotVec;
	int p = myDegree;
	double nip;
	int m = U.size() - 1;
	if (i == 0 && u == U[0] || i == m - p - 1 && u == U[m])
	{
		nip = 1.0;
		return nip;
	}
	if (u < U[i] || u >= U[i + p + 1])
	{
		nip = 0.0;
		return nip;
	}
	std::vector<double> N(p + 1, 0);
	for (int j = 0; j <= p; j++)
	{
		if (u >= U[i + j] && u < U[i + j + 1])
		{
			N[j] = 1.0;
		}
		else
		{
			N[j] = 0.0;
		}
	}
	for (int k = 1; k <= p; k++)
	{
		double saved;
		if (N[0] == 0.0)
		{
			saved = 0.0;
		}
		else
		{
			saved = ((u - U[i]) * N[0]) / (U[i + k] - U[i]);
		}
		for (int j = 0; j < p - k + 1; j++)
		{
			double Uleft = U[i + j + 1];
			double Uright = U[i + j + k + 1];
			if (N[j + 1] == 0.0)
			{
				N[j] = saved;
				saved = 0.0;
			}
			else
			{
				double temp = N[j + 1] / (Uright - Uleft);
				N[j] = saved + (Uright - u) * temp;
				saved = (u - Uleft) * temp;
			}
		}
	}
	nip = N[0];
	return nip;
}

bool AlgPointFitCurve::averageVec(int num, int firstNum, int endNum, std::vector<int>& avegVec)
{
	if (endNum - firstNum < num)
	{
		return false;
	}
	double stp = (double)(endNum - firstNum) / (double)(num - 1);
	avegVec.push_back(firstNum);
	for (int i = 1; i <= num - 2; i++)
	{
		avegVec.push_back(firstNum + (int)(i * stp));
	}
	avegVec.push_back(endNum);
	return true;
}

bool AlgPointFitCurve::initPoles(int poleNum)
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
	/*for (int i = 0; i < vec.size(); i++)
	{
		myPoleArray.SetValue(i + 1, myOriginFitData[vec[i] - 1]);
	}*/

	myInitialPoleIndex = vec;
	
	return true;
}

void AlgPointFitCurve::setParameter()
{
	int dataNum = myOriginFitData.size();
	myParameterArray = TColStd_Array1OfReal(1, dataNum);

	double TotalD = 0;
	for (int i = 0; i < dataNum - 1; i++)
	{
		TotalD += myOriginFitData[i].Distance(myOriginFitData[i + 1]);
	}

	double firstParam = 0;
	double lastParam = 1;
	if (!myCurve.IsNull())
	{
		firstParam = myCurve->FirstParameter();
		lastParam = myCurve->LastParameter();
	}

	double paramLength = lastParam - firstParam;
	for (int i = 1; i <= dataNum; i++)
	{
		if (i == 1) myParameterArray.SetValue(i, firstParam);
		else if (i == dataNum) myParameterArray.SetValue(i, lastParam);
		else
		{
			double temp = myOriginFitData[i - 1].Distance(myOriginFitData[i - 2]) / TotalD * paramLength;
			myParameterArray.SetValue(i, myParameterArray.Value(i - 1) + temp);
		}
		//std::cout << myParameterArray.Value(i) << std::endl;
	}
}

void AlgPointFitCurve::setKnots()
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

	/*myFlatKnots.Resize(1, myKnotVec.size(), false);
	for (int i = 0; i < myKnotVec.size(); i++)
	{
		myFlatKnots.SetValue(i + 1, myKnotVec[i]);
	}*/

	myMultArray.Init(1);
	myMultArray.SetValue(1, myDegree + 1);
	myMultArray.SetValue(knotNum, myDegree + 1);
}
