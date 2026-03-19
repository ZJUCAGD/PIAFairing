#include "AlgPointFitSurface.h"

bool AlgPointFitSurface::Init(const std::vector<std::vector<gp_Pnt>>& data, int poleNum1, int poleNum2, int uDegree, int vDegree)
{
    myOriginFitData = data;
	myUDegree = uDegree;
	myVDegree = vDegree;

	myControlPoints = TColgp_Array2OfPnt(1, poleNum1,1, poleNum2);
	setParameter();
	setKnots();
	executeEquation();
	creatBSplineCurve();

    return true;
}

Handle(Geom_BSplineSurface) AlgPointFitSurface::getResult()
{
    return mySurface;
}

void AlgPointFitSurface::setParameter()
{
	int dataNum1 = myOriginFitData.size();
	int dataNum2 = myOriginFitData[0].size();
	myUParameterArray = TColStd_Array2OfReal(1, dataNum1,1, dataNum2);
	myVParameterArray = TColStd_Array2OfReal(1, dataNum1, 1, dataNum2);

	
	//每行
	for (int i = 0; i < dataNum1; i++)
	{
		double TotalD = 0;
		double firstParam = 0;
		double lastParam = 1;
		double paramLength = lastParam - firstParam;
		//每列
		for (int j = 0; j < dataNum2-1; j++)
		{
			TotalD += myOriginFitData[i][j].Distance(myOriginFitData[i][j+1]);
		}

		for (int j = 1; j <= dataNum2; j++)
		{
			if (j == 1) myUParameterArray.SetValue(i+1,j, firstParam);
			else if (j == dataNum2) myUParameterArray.SetValue(i + 1, j, lastParam);
			else
			{
				double temp = myOriginFitData[i][j - 1].Distance(myOriginFitData[i][j - 2]) / TotalD * paramLength;
				myUParameterArray.SetValue(i+1, j, myUParameterArray.Value(i+1, j - 1) + temp);
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
			TotalD += myOriginFitData[j][i].Distance(myOriginFitData[j+1][i]);
		}

		for (int j = 1; j <= dataNum1; j++)
		{
			if (j == 1) myVParameterArray.SetValue(j,i+1, firstParam);
			else if (j == dataNum1) myVParameterArray.SetValue(j,i+1, lastParam);
			else
			{
				double temp = myOriginFitData[j-1][i].Distance(myOriginFitData[j-2][i]) / TotalD * paramLength;
				myVParameterArray.SetValue(j,i+1, myVParameterArray.Value(j - 1,i + 1) + temp);
			}
		}
	}
}

void AlgPointFitSurface::setKnots()
{
	int uKnotNum = myControlPoints.UpperRow() - myUDegree + 1;
	int vKnotNum = myControlPoints.UpperCol() - myVDegree + 1;

	myUKnots = TColStd_Array1OfReal(1, uKnotNum);
	myVKnots = TColStd_Array1OfReal(1, vKnotNum);

	myUMults = TColStd_Array1OfInteger(1, uKnotNum);
	myVMults = TColStd_Array1OfInteger(1, vKnotNum);

	double uSpace = 1.0 / (double)(uKnotNum - 1);
	double vSpace = 1.0 / (double)(vKnotNum - 1);
	for (int i = 0; i < uKnotNum; i++)
	{
		double temp = i * uSpace;
		myUKnots.SetValue(i+1, temp);
	}
	for (int j = 0; j< vKnotNum; j++)
	{
		double temp = j * vSpace;
		myVKnots.SetValue(j + 1, temp);
	}

	myUMults.Init(1);
	myUMults.SetValue(1, myUDegree + 1);
	myUMults.SetValue(uKnotNum, myUDegree + 1);

	myVMults.Init(1);
	myVMults.SetValue(1, myVDegree + 1);
	myVMults.SetValue(vKnotNum, myVDegree + 1);


	std::vector<double> tempFlatUKnot, tempFlatVKnot;
	for (int i = 1; i <= myUKnots.Size(); i++)
	{
		for (int j = 1; j <= myUMults.Value(i); j++)
		{
			tempFlatUKnot.push_back(myUKnots.Value(i));
		}
	}
	for (int i = 1; i <= myVKnots.Size(); i++)
	{
		for (int j = 1; j <= myVMults.Value(i); j++)
		{
			tempFlatVKnot.push_back(myVKnots.Value(i));
		}
	}

	myUFlatKnots.Resize(1, tempFlatUKnot.size(), false);
	myVFlatKnots.Resize(1, tempFlatVKnot.size(), false);
	for (int i = 0; i < tempFlatUKnot.size(); i++)
	{
		myUFlatKnots.SetValue(i + 1, tempFlatUKnot[i]);
	}
	for (int i = 0; i < tempFlatVKnot.size(); i++)
	{
		myVFlatKnots.SetValue(i + 1, tempFlatVKnot[i]);
	}
}

void AlgPointFitSurface::executeEquation()
{
	int dataNum = myOriginFitData.size() * myOriginFitData[0].size();
	int poleNum = myControlPoints.NbColumns() * myControlPoints.NbRows();
	Eigen::MatrixXd aMatrixN = Eigen::MatrixXd::Zero(dataNum, poleNum);
	//std::cout << dataNum << std::endl;
	//std::cout << poleNum << std::endl;

	for (int i = 0; i < myOriginFitData.size(); i++)
	{
		for (int j = 0; j < myOriginFitData[i].size(); j++)
		{
			double uParam = myUParameterArray.Value(i + 1, j + 1);
			double vParam = myVParameterArray.Value(i + 1, j + 1);

			Standard_Integer theFirstUIndex;//以1开头
			math_Matrix BsplineBasisU(1, 1, 1, myUDegree + 1);//以1开头
			BSplCLib::EvalBsplineBasis(0, myUDegree + 1, myUFlatKnots, uParam, theFirstUIndex, BsplineBasisU);
			
			Standard_Integer theFirstVIndex;//以1开头
			math_Matrix BsplineBasisV(1, 1, 1, myVDegree + 1);//以1开头
			BSplCLib::EvalBsplineBasis(0, myVDegree + 1, myVFlatKnots, vParam, theFirstVIndex, BsplineBasisV);

			int row = i * myOriginFitData[i].size() + j;
			for (int k = 0; k <= myUDegree; k++)
			{
				for (int p = 0; p <= myVDegree; p++)
				{
					int col = (theFirstUIndex + k - 1) * (myControlPoints.UpperCol()) + theFirstVIndex + p - 1;
					aMatrixN(row, col) = BsplineBasisU(1, 1 + k) * BsplineBasisV(1, 1 + p);
				}
			}
		}
		
	}

	//A
	Eigen::MatrixXd aMatrixA = aMatrixN.transpose() * aMatrixN;

	Eigen::VectorXd Qx(dataNum);
	Eigen::VectorXd Qy(dataNum);
	Eigen::VectorXd Qz(dataNum);

	for (int i = 0; i < myOriginFitData.size(); i++)
	{
		for (int j = 0; j < myOriginFitData[i].size(); j++)
		{
			gp_Pnt point = myOriginFitData[i][j];
			int row = i * myOriginFitData[i].size() + j;
			Qx(row) = point.X();
			Qy(row) = point.Y();
			Qz(row) = point.Z();
		}
	}

	Eigen::MatrixXd B = aMatrixN.transpose();
	Eigen::MatrixXd Bx = B * Qx;
	Eigen::MatrixXd By = B * Qy;
	Eigen::MatrixXd Bz = B * Qz;

	Eigen::VectorXd Px = aMatrixA.colPivHouseholderQr().solve(Bx);
	Eigen::VectorXd Py = aMatrixA.colPivHouseholderQr().solve(By);
	Eigen::VectorXd Pz = aMatrixA.colPivHouseholderQr().solve(Bz);

	for (int i = 1; i <= myControlPoints.Size(); i++)
	{
		int row = (i-1) / myControlPoints.UpperCol();
		int col = (i-1) % myControlPoints.UpperCol();
		myControlPoints.SetValue(row + 1, col + 1, gp_Pnt(Px(i - 1), Py(i - 1), Pz(i - 1)));
	}
}

void AlgPointFitSurface::creatBSplineCurve()
{
	try
	{
		mySurface = new Geom_BSplineSurface(myControlPoints, myUKnots, myVKnots, myUMults, myVMults, myUDegree, myVDegree);
	}
	catch (...)
	{
		std::cout << "Generate surface failed.";
	}
}
