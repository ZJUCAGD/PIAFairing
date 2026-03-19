#include "AlgCurveFairingByEnergy.h"

bool AlgCurveFairingByEnergy::Init(const Handle(Geom_BSplineCurve)& curve, int energyType)
{
    if (curve.IsNull())
    {
        return false;
    }
    m_curve = curve;
    m_fairFunctional = new fairFunctional(curve->Knots(), curve->Multiplicities(), curve->Degree(), energyType);
    return true;
}

void AlgCurveFairingByEnergy::execute(double weight)
{
    int poleNum = m_curve->NbPoles();
    const auto poles = m_curve->Poles();

    Eigen::MatrixXd aMatrixDr = Eigen::MatrixXd::Zero(poleNum, poleNum);
    for (int i = 0; i < poleNum; i++)
    {
        for (int j = 0; j < poleNum; j++)
        {
            aMatrixDr(i, j) = m_fairFunctional->getValue(i, j);
        }
    }

    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(poleNum, poleNum);

    Eigen::MatrixXd aMatrixA = (1- weight) * identity + weight * aMatrixDr;

    Eigen::VectorXd Qx(poleNum);
    Eigen::VectorXd Qy(poleNum);
    Eigen::VectorXd Qz(poleNum);

    for (int i = 0; i < poleNum; i++)
    {
        gp_Pnt pole = poles.Value(i + 1);
        Qx(i) = pole.X();
        Qy(i) = pole.Y();
        Qz(i) = pole.Z();
    }

    Eigen::MatrixXd Bx = (1 - weight) * Qx;
    Eigen::MatrixXd By = (1 - weight) * Qy;
    Eigen::MatrixXd Bz = (1 - weight) * Qz;

    Eigen::VectorXd Px = aMatrixA.colPivHouseholderQr().solve(Bx);
    Eigen::VectorXd Py = aMatrixA.colPivHouseholderQr().solve(By);
    Eigen::VectorXd Pz = aMatrixA.colPivHouseholderQr().solve(Bz);

    for (int i = 1; i <= poleNum; i++)
    {
        auto newPole = gp_Pnt(Px(i - 1), Py(i - 1), Pz(i - 1));
        m_curve->SetPole(i, newPole);
    }

    poleDistance(poles, m_curve->Poles());

}

Handle(Geom_BSplineCurve) AlgCurveFairingByEnergy::getResult()
{
    return m_curve;
}

void AlgCurveFairingByEnergy::poleDistance(const TColgp_Array1OfPnt& poles1, const TColgp_Array1OfPnt& poles2)
{
    if (poles1.Size() != poles2.Size())
    {
        std::cout << "¿ØÖÆµãÓÐ´íÎó£¡\n";
    }
    
    int num = poles1.Size();
    double maxDist = -1;
    double avgDist = 0;
    for (int i = 1; i <= num; i++)
    {
        auto p1 = poles1.Value(i);
        auto p2 = poles2.Value(i);
        double dist = p1.Distance(p2);
        if (dist > maxDist)
        {
            maxDist = dist;
        }
        avgDist += p1.SquareDistance(p2);
    }

    avgDist /= (double)num;
    avgDist = sqrt(avgDist);

}
