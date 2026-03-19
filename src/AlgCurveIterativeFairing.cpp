#include "AlgCurveIterativeFairing.h"

bool AlgCurveIterativeFairing::Init(const Handle(Geom_BSplineCurve)& theCurve, const std::vector<double>& distances)
{
	if (theCurve.IsNull())
	{
		std::cout << "曲线为空\n";
		return false;
	}
	myCurve = theCurve;
	myDistances = distances;
	auto poles = myCurve->Poles();
	int poleNum = poles.Size();
	myOriginPoles = TColgp_Array1OfPnt(1, poleNum);
	for (int i = 1; i <= poleNum; i++)
	{
		myOriginPoles.SetValue(i, poles.Value(i));
	}

	auto knots = myCurve->Knots();
	auto mults = myCurve->Multiplicities();
	for (int i = 1; i <= knots.Size(); i++)
	{
		auto temp = knots.Value(i);
		for (int j = 1; j <= mults.Value(i); j++)
		{
			myFlatKnots.push_back(temp);
		}
	}
	myLambdas.resize(myDistances.size(), 1e-3);
	myVis = myDistances;
	for (int i = 0; i < myVis.size(); i++)
	{
		double vis = myVis[i];
		myVis[i] = vis / 100.0;
	}
	return true;
}

void AlgCurveIterativeFairing::execute()
{
	//迭代停止的条件是什么
	double dist1 = 0, dist2 = 1;
	for (int i = 0; i < 600; i++)
	{
		Iterative();
		auto poles = myCurve->Poles();
		double maxDist = 0;
		for (int j = 1; j <= poles.Size(); j++)
		{
			double dis = sqrt(poles.Value(j).Distance(myOriginPoles.Value(j)));
			if (dis > maxDist)
			{
				maxDist = dis;
			}
		}
		dist1 = maxDist;
		if (abs(dist1 - dist2) < 1e-6)
		{
			std::cout << i << std::endl;
			break;
		}
		dist2 = dist1;
		std::cout << maxDist << std::endl;
	}
}

void AlgCurveIterativeFairing::getResult(Handle(Geom_BSplineCurve)& aCurve)
{
	aCurve = myCurve;
}

void AlgCurveIterativeFairing::getDiffVec(std::vector<gp_Vec>& aVec)
{
	aVec.resize(myOriginPoles.Size(), gp_Vec{0,0,0});
	const TColgp_Array1OfPnt currentPoles = myCurve->Poles();
	for (int i = 1; i <= currentPoles.Size(); i++)
	{
		auto p1 = currentPoles.Value(i);
		auto p2 = myOriginPoles.Value(i);
		gp_Vec diff = p1.XYZ() - p2.XYZ();
		aVec[i - 1] = diff;
	}
}

void AlgCurveIterativeFairing::getFairVec(std::vector<gp_Vec>& aVec)
{
	aVec.resize(myOriginPoles.Size(), gp_Vec{ 0,0,0 });
	const TColgp_Array1OfPnt currentPoles = myCurve->Poles();
	
	fairFunctional alg(myFlatKnots,myCurve->Degree());
	for (int i = 1; i <= currentPoles.Size(); i++)
	{
		gp_Vec temp{ 0,0,0 };
		for (int j = 1; j <= currentPoles.Size(); j++)
		{
			gp_Vec pole = gp_Vec(currentPoles.Value(j).XYZ());
			temp += alg.getValue(2, i-1, j-1) * pole;
		}
		aVec[i - 1] = temp;
	}
}

void AlgCurveIterativeFairing::Iterative()
{
	//迭代终止的条件是什么？
	const TColgp_Array1OfPnt currentPoles = myCurve->Poles();
	const std::vector<double> currentLambdas = myLambdas;

	std::vector<gp_Vec> aDiffVec, aFairVec;
	getDiffVec(aDiffVec);
	getFairVec(aFairVec);

	std::vector<gp_Pnt> newPnts;

	for (int i = 1; i <= currentPoles.Size(); i++)
	{
		//TODO:如何确定alpha
		double alpha = 0.0001;
		auto curPole = currentPoles.Value(i);
		gp_Vec aVec = curPole.XYZ() + alpha * (2 * aFairVec[i - 1].XYZ() + 2 * currentLambdas[i - 1] * aDiffVec[i - 1].XYZ());
		newPnts.push_back(gp_Pnt(aVec.XYZ()));
	}

	
	for (int i = 0; i < myLambdas.size(); i++)
	{
		//TODO:如何确定alpha
		double alpha = 3e-3;
		double dist = currentPoles.Value(i + 1).Distance(myOriginPoles.Value(i + 1));
		double lambda = currentLambdas[i] - alpha * (dist * dist - myDistances[i] * myDistances[i]);
		myLambdas[i] = lambda;
	}

	for (int i = 1; i <= currentPoles.Size(); i++)
	{
		
		myCurve->SetPole(i, newPnts[i-1]);
	}
}
