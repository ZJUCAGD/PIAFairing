#include "CurveAndSurfaceTestCases.h"

#pragma comment(lib, "TKernel.lib")
#pragma comment(lib, "TKMath.lib")
#pragma comment(lib, "TKBRep.lib")
#pragma comment(lib, "TKTopAlgo.lib")
#pragma comment(lib, "TKDEIGES.lib")

int main()
{
    //curve PIA by adjust control Points
    CurveAndSurfaceTestCases aCase;
    aCase.curveCase1_G();
    //aCase.curveCase2_curve();

    //auto PIA by adjust control Points
    //aCase.curveAutoCase();
    return 0;
}