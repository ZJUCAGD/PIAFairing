#pragma once
#include "fairFunctional.h"


class AlgCurveEnergyMinModel
{
public:

	AlgCurveEnergyMinModel() {}

	//Ä¬ÈÏÓÃ¶þ½×
	bool Init(const std::vector<gp_Pnt>& data, int poleNum, int degree, int der = 2);
	bool Execute();
	Handle(Geom_BSplineCurve) getResult();

private:
	bool averageVec(int num, int firstNum, int endNum, std::vector<int>& avegVec);
	bool initPoles(int poleNum);
	void setParameter();
	void setKnots();
	void creatBSplineCurve();
	int FindSpan(double u);
	

	void BasisFuns(int i, double u, std::vector<double>& N);
	

private:
	void generateMatrixN();
	void generateMatrixDr();
	//void generateMatrixOmega();
	void generateMatrixA();
	

private:
	std::vector<gp_Pnt> myOriginFitData;

	Handle(Geom_BSplineCurve) myCurve;

	TColgp_Array1OfPnt myPoleArray;
	TColStd_Array1OfReal myParameterArray;
	TColStd_Array1OfReal myKnotArray;
	TColStd_Array1OfInteger myMultArray;
	Standard_Integer myDegree{ 3 };
	std::vector<double> myKnotVec;
	std::vector<int> myInitialPoleIndex;
	
	int myDerivative{ 2 };

private:
	Eigen::MatrixXd m_matrixN;
	Eigen::MatrixXd m_matrixDr;
	Eigen::MatrixXd m_matrixA;
	//Eigen::MatrixXd m_matrixOmega;
	double myOmega{ 5e-4 };
};

