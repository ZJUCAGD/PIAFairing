#pragma once
#include "AlgCurFPIAByAdjustingControlPoint.h"
#include "AlgSurfFPIAByAdjustingControlPoint.h"

#include "AlgPointFitCurve.h"
#include "AlgPointFitSurface.h"

#include "AlgCurveFairingByEnergy.h"
#include "AlgSurfaceFairingByEnergy.h"

class CurveAndSurfaceTestCases
{
public:
	CurveAndSurfaceTestCases() {}
	void surfaceCase3_Car();


	void curveCase1_G();
	void curveCase2_curve();
	void curveCase3();

	void curveAutoCase();
	void surfaceAutoCase();

private:
	void SurfaceToIgs(const Handle(Geom_BSplineSurface)& surface, Standard_CString theFileName);
    void CurveToIgs(const Handle(Geom_BSplineCurve)& fairingCurve, Standard_CString theFileName);
    
	void PointsToPLY(const std::vector<gp_Pnt>& points, const std::string& filename);
    bool readFileToData(const Standard_CString& theFileName, std::vector<std::vector<gp_Pnt>>& theData);
    void readFileToCurve(std::string fileName, Handle(Geom_BSplineCurve)& theCurve);
    //¶ÁÈ¡IgsÎÄ¼þ
	bool readIGESFileToCurve(const Standard_CString& theFileName, Handle(Geom_BSplineCurve)& theCurve);
	bool readIGESFileToSurface(const Standard_CString& theFileName, Handle(Geom_BSplineSurface)& theSurface);
    void createSpiralCurve(Handle(Geom_BSplineCurve)& theCurve);

	std::vector<gp_Pnt> addGaussianNoiseToPoint(const std::vector<gp_Pnt>& data, int begin, int end, double mean = 0.0, double stddev = 1.0);
	std::vector<std::vector<gp_Pnt>> addGaussianNoiseToSurfacePoint(const std::vector<std::vector<gp_Pnt>>& data, int begin1, int end1, int begin2, int end2, double mean = 0.0, double stddev = 1.0);
	void checkCurve(const Handle(Geom_BSplineCurve)& theCurve);
	void checkSurface(const Handle(Geom_BSplineSurface)& theSurface);
};

