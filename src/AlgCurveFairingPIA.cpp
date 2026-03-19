#include "AlgCurveFairingPIA.h"

bool AlgCurveFairingPIA::Init(const std::vector<gp_Pnt>& data, int poleNum, int degree, int der)
{
	myOriginFitData = data;
	myDegree = degree;
	myDerivative = der;

	if (!initPoles(poleNum))
	{
		return false;
	}
	myWeights.resize(poleNum, 1e-4);

	setParameter();
	setKnots();
	creatBSplineCurve();
	generateParamMu();


	return true;
}

void AlgCurveFairingPIA::Init(const std::vector<gp_Pnt>& data, const Handle(Geom_BSplineCurve)& aCurve, int der)
{
	myCurve = aCurve;
	myOriginFitData = data;
	myDerivative = der;

	myDegree = aCurve->Degree();

	myPoleArray = TColgp_Array1OfPnt(1, aCurve->Poles().Size());
	myPoleArray.Assign(aCurve->Poles());

	myKnotArray = TColStd_Array1OfReal(1, aCurve->Knots().Size());
	myKnotArray.Assign(aCurve->Knots());

	myMultArray = TColStd_Array1OfInteger(1, aCurve->Multiplicities().Size());
	myMultArray.Assign(aCurve->Multiplicities());

	myWeights.resize(myPoleArray.Size(), 1e-4);

	for (int i=1;i<= myKnotArray.Size();i++)
	{
		int num = myMultArray.Value(i);
		double knot = myKnotArray.Value(i);

		for (int j = 1; j <= num; j++)
		{
			myKnotVec.push_back(knot);
		}
		
	}

	myFlatKnots.Resize(1, myKnotVec.size(), false);
	for (int i = 0; i < myKnotVec.size(); i++)
	{
		myFlatKnots.SetValue(i + 1, myKnotVec[i]);
		//std::cout << myKnotVec[i] << std::endl;

	}

	setParameter();
	generateParamMu();

}

void AlgCurveFairingPIA::executeGobalPoles(double maxTimes, double error)
{
	int times = 0;
	double error1 = getFittingError();
	double error2 = error1;
	while (times < maxTimes)
	{
		updateGobalPoles();
		creatBSplineCurve();
		//std::cout << error1 << " " << error2 << std::endl;
		
		error1 = getFittingError();
		if (error1 < error || abs(error1- error2)< 1e-8)
		{
			return;
		}
		//std::cout << error1 << " " << error2 << std::endl;
		error2 = error1;
		times++;
	}
}

void AlgCurveFairingPIA::executeLocalPoles(std::vector<int> indexList, double maxTimes, double error)
{
	myAdjustPoleIndex = indexList;

	int times = 0;
	double error1 = getFittingError();
	double error2 = error1;
	while (times < maxTimes)
	{
		updateLocalPoles();
		creatBSplineCurve();

		error1 = getFittingError();
		//std::cout << error1 << " " << error2 << std::endl;

		if (error1 < error || abs(error1 - error2) < 1e-8)
		{
			return;
		}
		
		error2 = error1;
		times++;
	}
}

void AlgCurveFairingPIA::executeLocalPoles(int dataIndexBegin, int dataIndexEnd, double maxTimes, double error)
{
	TColStd_Array1OfReal flatKnots(1,myKnotVec.size());
	for (int i = 0; i < myKnotVec.size(); i++)
	{
		flatKnots.SetValue(i+1, myKnotVec[i]);
	}
	std::set<int> adjustIndex;

	for (int i = dataIndexBegin; i <= dataIndexEnd; i++)
	{
		double uParam = myParameterArray.Value(i);
		Standard_Integer theFirstIndex;//以1开头
		math_Matrix BsplineBasisU(1, 1, 1, myDegree + 1);//以1开头
		BSplCLib::EvalBsplineBasis(0, myDegree + 1, flatKnots, uParam, theFirstIndex, BsplineBasisU);
		for (int j=0;j<= myDegree;j++)
		{
			adjustIndex.insert(theFirstIndex+j);
		}
	}

	myAdjustPoleIndex.clear();
	for (auto& index : adjustIndex)
	{
		myAdjustPoleIndex.push_back(index);
	}

	//可以计算的更少
	int times = 0;
	double error1 = getFittingError();
	double error2 = error1;
	while (times < maxTimes)
	{
		updateLocalPoles();
		creatBSplineCurve();

		error1 = getFittingError();
		if (error1 < error || abs(error1 - error2) < 1e-8)
		{
			return;
		}
		error2 = error1;
		times++;
	}
}

void AlgCurveFairingPIA::executeEquation()
{
	//Omega
	Eigen::MatrixXd aMatrixOmega = Eigen::MatrixXd::Zero(myPoleArray.Size(), myPoleArray.Size());
	for (int i = 0; i < myPoleArray.Size(); i++)
	{
		aMatrixOmega(i, i) = myWeights[i];
	}

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

	//Dr
	Eigen::MatrixXd aMatrixDr = Eigen::MatrixXd::Zero(myPoleArray.Size(), myPoleArray.Size());
	fairFunctional IntegEngine(myKnotVec, myDegree);

	for (int i = 0; i < myPoleArray.Size(); i++)
	{
		for (int j = 0; j < myPoleArray.Size(); j++)
		{
			aMatrixDr(i, j) = IntegEngine.getValue(myDerivative, i, j);
		}
	}

	//A
	Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(myPoleArray.Size(), myPoleArray.Size());
	Eigen::MatrixXd aMatrixA = (identity - aMatrixOmega) * aMatrixN.transpose() * aMatrixN + aMatrixOmega * aMatrixDr;

	//Mu
	Eigen::MatrixXd Mu = Eigen::MatrixXd::Zero(myPoleArray.Size(), myPoleArray.Size());
	for (int i = 0; i < myPoleArray.Size(); i++)
	{
		double sum = 0;
		for (int j = 0; j < myPoleArray.Size(); j++)
		{
			sum += abs(aMatrixA(i, j));
		}
		muParamMu.push_back(1.0 / sum);
	}

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
	Eigen::MatrixXd B = (identity - aMatrixOmega) * aMatrixN.transpose();

	Eigen::MatrixXd Bx = B * Qx;
	Eigen::MatrixXd By = B * Qy;
	Eigen::MatrixXd Bz = B * Qz;

	


	Eigen::VectorXd Px = aMatrixA.colPivHouseholderQr().solve(Bx);
	Eigen::VectorXd Py = aMatrixA.colPivHouseholderQr().solve(By);
	Eigen::VectorXd Pz = aMatrixA.colPivHouseholderQr().solve(Bz);

	//std::cout << "aMatrixN\n";

	//std::cout << aMatrixOmega << std::endl;


	//std::cout << aMatrixN.row(1) << std::endl;
	//std::cout << Pz << std::endl;


	for (int i = 1; i <= myPoleArray.Size(); i++)
	{
		myPoleArray.SetValue(i, gp_Pnt(Px(i - 1), Py(i - 1), Pz(i - 1)));
	}
	creatBSplineCurve();
}

void AlgCurveFairingPIA::setWeights(std::vector<int> indexList, std::vector<double> weights)
{
	for (int i = 0; i < indexList.size(); i++)
	{
		myWeights[indexList[i]] = weights[i];
	}
}

Handle(Geom_BSplineCurve) AlgCurveFairingPIA::getResult()
{
	return myCurve;
}

bool AlgCurveFairingPIA::averageVec(int num, int firstNum, int endNum, std::vector<int>& avegVec)
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

bool AlgCurveFairingPIA::initPoles(int poleNum)
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

void AlgCurveFairingPIA::setParameter()
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

void AlgCurveFairingPIA::setKnots()
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

	myFlatKnots.Resize(1, myKnotVec.size(), false);
	for (int i = 0; i < myKnotVec.size(); i++)
	{
		myFlatKnots.SetValue(i + 1, myKnotVec[i]);
	}

	myMultArray.Init(1);
	myMultArray.SetValue(1, myDegree + 1);
	myMultArray.SetValue(knotNum, myDegree + 1);
}

void AlgCurveFairingPIA::creatBSplineCurve()
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

int AlgCurveFairingPIA::FindSpan(double u)
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

void AlgCurveFairingPIA::BasisFuns(int i, double u, std::vector<double>& N)
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

void AlgCurveFairingPIA::updateGobalPoles()
{
	std::vector<gp_Vec> diffVec, fairVec;
	int pointNum = myOriginFitData.size();
	int poleNum = myPoleArray.Size();

	//计算每个数据点差向量
	for (int i = 1; i <= pointNum; i++)
	{
		gp_Pnt temp;
		myCurve->D0(myParameterArray.Value(i), temp);
		auto aData = myOriginFitData[i - 1];
		gp_Vec dif(aData.XYZ() - temp.XYZ());
		diffVec.push_back(dif);
	}

	//计算每个控制点差向量
	std::vector<gp_Vec> fitVec(poleNum, gp_Vec(0, 0, 0));
	for (int i = 1; i <= pointNum; i++)
	{
		double aParam = myParameterArray.Value(i);
		Standard_Integer theFirstIndex;//以1开头
		math_Matrix BsplineBasisU(1, 1, 1, myDegree + 1);//以1开头
		BSplCLib::EvalBsplineBasis(0, myDegree + 1, myFlatKnots, aParam, theFirstIndex, BsplineBasisU);

		for (int j = 0; j <= myDegree; j++)
		{
			gp_Vec tempVec = BsplineBasisU(1, 1 + j) * diffVec[i - 1];
			fitVec[theFirstIndex + j - 1] += tempVec;
		}
	}

	fairFunctional IntegEngine(myKnotVec, myDegree);
	for (int i = 1; i <= poleNum; i++)
	{
		gp_Vec tempPoint(0, 0, 0);
		for (int j = 1; j <= poleNum; j++)
		{
			double para = IntegEngine.getValue(myDerivative, i - 1, j - 1);
			gp_Pnt aPole = myPoleArray.Value(j);
			tempPoint += para * aPole.XYZ();
		}
		fairVec.push_back(tempPoint);
	}

	for (int i = 1; i <= poleNum; i++)
	{
		double omega = myWeights[i - 1];
		gp_Vec iterVec = (1 - omega) * fitVec[i - 1] - omega * fairVec[i - 1];
		gp_Pnt newPole = gp_Pnt(myPoleArray.Value(i).XYZ() + muParamMu[i - 1] * iterVec.XYZ());
		myPoleArray.SetValue(i, newPole);
	}
}

void AlgCurveFairingPIA::updateLocalPoles()
{
	/*int pointNum = myOriginFitData.size();
	int poleNum = myPoleArray.Size();

	std::vector<gp_Vec> fitVecList(myAdjustPoleIndex.size(), gp_Vec(0, 0, 0));
	for (int i = 1; i <= pointNum; i++)
	{
		double aParam = myParameterArray.Value(i);
		auto aPoint = myOriginFitData[i];

		for (int j=0;j< myAdjustPoleIndex.size();j++)
		{
			int indexBasisFun = myAdjustPoleIndex[j] - 1;
			double basisFun = oneBasisFun(indexBasisFun, aParam);
			if (abs(basisFun) > 1e-8)
			{
				gp_Pnt temp;
				myCurve->D0(aParam, temp);
				gp_Vec dif(aPoint.XYZ() - temp.XYZ());
				fitVecList[j] += dif * basisFun;
			}
		}
	}

	fairFunctional IntegEngine(myKnotVec, myDegree);
	std::vector<gp_Vec> fairVecList(myAdjustPoleIndex.size(), gp_Vec(0, 0, 0));

	for (int i = 0; i < myAdjustPoleIndex.size(); i++)
	{
		int index = myAdjustPoleIndex[i];
		for (int j = 1; j <= poleNum; j++)
		{
			double para = IntegEngine.getValue(myDerivative, index - 1, j - 1);
			gp_Vec vec(myPoleArray.Value(j).XYZ());
			fairVecList[i] += para * vec;
		}
	}
	
	for (int i = 0; i < myAdjustPoleIndex.size(); i++)
	{
		int index = myAdjustPoleIndex[i];
		double omega = myWeights[index - 1];

		gp_Vec difVec = (1 - omega) * fitVecList[i] - omega * fairVecList[i];
		gp_Pnt newPole = gp_Pnt(myPoleArray.Value(index).XYZ() + muParamMu[index - 1] * difVec.XYZ());

		myPoleArray.SetValue(index, newPole);
	}
	
	creatBSplineCurve();*/

	std::vector<gp_Vec> diffVec, fairVec;
	int pointNum = myOriginFitData.size();
	int poleNum = myPoleArray.Size();

	//计算每个数据点差向量
	for (int i = 1; i <= pointNum; i++)
	{
		gp_Pnt temp;
		myCurve->D0(myParameterArray.Value(i), temp);
		auto aData = myOriginFitData[i - 1];
		gp_Vec dif(aData.XYZ() - temp.XYZ());
		diffVec.push_back(dif);
	}

	//计算每个控制点差向量
	std::vector<gp_Vec> fitVec(poleNum, gp_Vec(0, 0, 0));
	for (int i = 1; i <= pointNum; i++)
	{
		double aParam = myParameterArray.Value(i);
		Standard_Integer theFirstIndex;//以1开头
		math_Matrix BsplineBasisU(1, 1, 1, myDegree + 1);//以1开头
		BSplCLib::EvalBsplineBasis(0, myDegree + 1, myFlatKnots, aParam, theFirstIndex, BsplineBasisU);

		for (int j = 0; j <= myDegree; j++)
		{
			gp_Vec tempVec = BsplineBasisU(1, 1 + j) * diffVec[i - 1];
			fitVec[theFirstIndex + j - 1] += tempVec;
		}
	}

	fairFunctional IntegEngine(myKnotVec, myDegree);
	for (int i = 1; i <= poleNum; i++)
	{
		gp_Vec tempPoint(0, 0, 0);
		for (int j = 1; j <= poleNum; j++)
		{
			double para = IntegEngine.getValue(myDerivative, i - 1, j - 1);
			gp_Pnt aPole = myPoleArray.Value(j);
			tempPoint += para * aPole.XYZ();
		}
		fairVec.push_back(tempPoint);
	}

	for (int i = 1; i <= poleNum; i++)
	{
		double omega = myWeights[i - 1];
		gp_Vec iterVec = (1 - omega) * fitVec[i - 1] - omega * fairVec[i - 1];
		gp_Pnt newPole = gp_Pnt(myPoleArray.Value(i).XYZ() + muParamMu[i - 1] * iterVec.XYZ());
		myPoleArray.SetValue(i, newPole);
	}

	for (int i = 0; i < myAdjustPoleIndex.size(); i++)
	{
		int index = myAdjustPoleIndex[i];
		double omega = myWeights[index - 1];

		gp_Vec iterVec = (1 - omega) * fitVec[index - 1] - omega * fairVec[index - 1];
		gp_Pnt newPole = gp_Pnt(myPoleArray.Value(index).XYZ() + muParamMu[index - 1] * iterVec.XYZ());
		myPoleArray.SetValue(index, newPole);
	}
}

double AlgCurveFairingPIA::getFittingError()
{
	double res{ 0.0 };
	int num = myParameterArray.Size();
	for (int i = 1; i < num; i++)
	{
		gp_Pnt temp;
		gp_Pnt data = myOriginFitData[i-1];
		myCurve->D0(myParameterArray.Value(i), temp);
		gp_Vec dif(data.X() - temp.X(), data.Y() - temp.Y(), data.Z() - temp.Z());
		res += dif.Dot(dif);
	}
	res /= (double)num;
	return sqrt(res);
}

double AlgCurveFairingPIA::oneBasisFun(int i, double u)
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

void AlgCurveFairingPIA::generateParamMu()
{
	//Omega
	Eigen::MatrixXd aMatrixOmega = Eigen::MatrixXd::Zero(myPoleArray.Size(), myPoleArray.Size());
	for (int i = 0; i < myPoleArray.Size(); i++)
	{
		aMatrixOmega(i, i) = myWeights[i];
	}

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

	//Dr
	Eigen::MatrixXd aMatrixDr = Eigen::MatrixXd::Zero(myPoleArray.Size(), myPoleArray.Size());
	fairFunctional IntegEngine(myKnotVec, myDegree);

	for (int i = 0; i < myPoleArray.Size(); i++)
	{
		for (int j = 0; j < myPoleArray.Size(); j++)
		{
			aMatrixDr(i, j) = IntegEngine.getValue(myDerivative, i, j);
		}
	}

	//A
	Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(myPoleArray.Size(), myPoleArray.Size());
	Eigen::MatrixXd aMatrixA = (identity - aMatrixOmega) * aMatrixN.transpose() * aMatrixN + aMatrixOmega * aMatrixDr;

	//Mu
	Eigen::MatrixXd Mu = Eigen::MatrixXd::Zero(myPoleArray.Size(), myPoleArray.Size());
	for (int i = 0; i < myPoleArray.Size(); i++)
	{
		double sum = 0;
		for (int j = 0; j < myPoleArray.Size(); j++)
		{
			sum += abs(aMatrixA(i, j));
		}
		muParamMu.push_back(1.0 / sum);
		//std::cout << 1.0 / sum << std::endl;
	}
}

