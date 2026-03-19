#include "AlgEnergyCalculation.h"

// 计算 B 样条曲线的二阶导数
double AlgEnergyCalculation::BsplineSecondDerivative(double u) 
{
    double value = 0;
    if (m_type == 1)
    {
        gp_Pnt p1;
        gp_Vec v1;
        m_curve->D1(u, p1, v1);
        value = v1.SquareMagnitude();
    }
    else if (m_type == 2)
    {
        gp_Pnt p1;
        gp_Vec v1, v2;
        m_curve->D2(u, p1, v1, v2);
        value = v2.SquareMagnitude();
    }
    else
    {
        gp_Pnt p1;
        gp_Vec v1, v2, v3;
        m_curve->D3(u, p1, v1, v2, v3);
        value = v3.SquareMagnitude();
    }

    return value;
}

double AlgEnergyCalculation::BsplineSecondDerivative(double u, double v)
{
    gp_Pnt P;
    gp_Vec u1, v1, u2,v2,uv;
    m_surface->D2(u,v,P, u1, v1, u2, v2, uv);

    return u2.SquareMagnitude() + v2.SquareMagnitude() + 2*uv.SquareMagnitude();
}

double AlgEnergyCalculation::calEnergy(const Handle(Geom_BSplineCurve)& curve, int type)
{
    m_curve = curve;
    m_type = type;
    int sampleNum = 2000;
    //auto knots = curve->Knots();

    double energy = 0.0;
    double uMin = curve->FirstParameter();
    double uMax = curve->LastParameter();
    double deltaU = (uMax - uMin) / (sampleNum - 1);

    for (int k = 0; k < sampleNum; ++k) {
        double u = uMin + k * deltaU;
        double secondDeriv = BsplineSecondDerivative(u);
        energy += secondDeriv;
    }

    return energy * deltaU;
}

double AlgEnergyCalculation::calEnergy(const Handle(Geom_BSplineSurface)& surface, int type)
{
    m_surface = surface;
    m_type = type;
    int sampleNum = 2000;
    //auto knots = curve->Knots();

    double energy = 0.0;
    if (surface.IsNull())
    {
        std::cout << "isnull\n";
    }
    auto knots = surface->UKnots();
    std::cout << knots.Size() << std::endl;

    double uMin = surface->UKnot(1);
    double uMax = surface->UKnot(surface->UKnots().Size());

    double vMin = surface->VKnot(1);
    double vMax = surface->VKnot(surface->VKnots().Size());

    double deltaU = (uMax - uMin) / (sampleNum - 1);
    double deltaV = (vMax - vMin) / (sampleNum - 1);

    for (int i = 0; i < sampleNum; i++)
    {
        for (int k = 0; k < sampleNum; ++k) {
            double u = uMin + i * deltaU;
            double v = vMin + k * deltaV;

            double secondDeriv = BsplineSecondDerivative(u, v);
            energy += secondDeriv;
        }
    }
    

    return energy * deltaU * deltaV;
}
