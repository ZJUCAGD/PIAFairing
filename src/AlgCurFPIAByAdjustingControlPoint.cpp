#include "AlgCurFPIAByAdjustingControlPoint.h"

AlgCurFPIAByAdjustingControlPoint::~AlgCurFPIAByAdjustingControlPoint()
{
    delete m_fairFunctional;
}

bool AlgCurFPIAByAdjustingControlPoint::Init(const Handle(Geom_BSplineCurve)& curve, int energyType)
{
    if (curve.IsNull())
    {
        return false;
    }
    m_curve = curve;
    m_energyType = energyType;
    m_weights.resize(curve->NbPoles(), 0);

    m_originPoles = TColgp_Array1OfPnt(curve->Poles());
    m_fairFunctional = new fairFunctional(curve->Knots(), curve->Multiplicities(), curve->Degree(), m_energyType);

    //calMuParam();
    return true;
}

void AlgCurFPIAByAdjustingControlPoint::execute(int maxTimes, double error)
{
    calMuParam();
    int time = 0;
    auto lastPoles = m_curve->Poles();
    //AlgEnergyCalculation algEnergy;
    //std::cout << "第一个" << algEnergy.calEnergy(m_curve, 1) << std::endl;

    while (time < maxTimes)
    {
        iterative();

        auto curPoles = m_curve->Poles();
        
        double diff = iterativeDiff(lastPoles, curPoles);
        if (diff < 1e-6)
        {
            std::cout << time << std::endl;
            return;
        }
        lastPoles = TColgp_Array1OfPnt(curPoles);
        time++;
    }

    //poleDistance(lastPoles, m_curve->Poles());
    /*auto poles = m_curve->Poles();
    std::cout << "==============result===============\n";
    for (auto& p : poles)
    {
        checkPnt(p);
    }*/
}

void AlgCurFPIAByAdjustingControlPoint::execute(const std::vector<size_t>& indexList, const std::vector<double>& weights, int maxTimes)
{
    if (weights.size() < indexList.size())
    {
        std::cout << "请输入正确的权重.\n";
        return;
    }

    for (int i = 0; i < indexList.size(); i++)
    {
        m_weights[indexList[i]] = weights[i];
    }

    calMuParam();

    int time = 0;
    auto lastPoles = m_curve->Poles();
    
    while (time < maxTimes)
    {
        localIterative(indexList);

        auto curPoles = m_curve->Poles();

        double diff = iterativeDiff(lastPoles, curPoles);
        if (diff < 1e-6)
        {
            std::cout << time << std::endl;
            return;
        }
        lastPoles = TColgp_Array1OfPnt(curPoles);
        time++;
    }
}

void AlgCurFPIAByAdjustingControlPoint::setWeights(const std::vector<size_t>& indexList, const std::vector<double>& weights)
{
    if (weights.size() < indexList.size())
    {
        std::cout << "请输入正确的权重.\n";
        return;
    }

    for (int i = 0; i < indexList.size(); i++)
    {
        m_weights[indexList[i]] = weights[i];
    }
    calMuParam();
}

Handle(Geom_BSplineCurve) AlgCurFPIAByAdjustingControlPoint::getResult()
{
    return m_curve;
}

void AlgCurFPIAByAdjustingControlPoint::equationCheck()
{
    //检测矩阵形式
    int poleNum = m_originPoles.Size();

    //Omega
    Eigen::MatrixXd aMatrixOmega = Eigen::MatrixXd::Zero(poleNum, poleNum);
    for (int i = 0; i < poleNum; i++)
    {
        aMatrixOmega(i, i) = m_weights[i];
    }

    //Dr
    Eigen::MatrixXd aMatrixDr = Eigen::MatrixXd::Zero(poleNum, poleNum);
    for (int i = 0; i < poleNum; i++)
    {
        for (int j = 0; j < poleNum; j++)
        {
            aMatrixDr(i, j) = m_fairFunctional->getValue(i, j);
        }
    }
    checkWeightInterval(aMatrixDr);

    //I
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(poleNum, poleNum);
    
    //A
    Eigen::MatrixXd aMatrixA = identity - aMatrixOmega + aMatrixOmega * aMatrixDr;

    //mu
    m_mus.resize(poleNum, 0);
    for (int i = 0; i < poleNum; i++)
    {
        double temp = 0;
        for (int j = 0; j < poleNum; j++)
        {
            temp += abs(aMatrixA(i, j));
        }
        m_mus[i] = (1.0/temp);
        std::cout << m_mus[i] << std::endl;
    }

    Eigen::MatrixXd aMatrixB = identity - aMatrixOmega;

    Eigen::VectorXd Qx(poleNum);
    Eigen::VectorXd Qy(poleNum);
    Eigen::VectorXd Qz(poleNum);

    for (int i = 0; i < poleNum; i++)
    {
        gp_Pnt pole = m_originPoles.Value(i + 1);
        Qx(i) = pole.X();
        Qy(i) = pole.Y();
        Qz(i) = pole.Z();
    }

    Eigen::MatrixXd Bx = aMatrixB * Qx;
    Eigen::MatrixXd By = aMatrixB * Qy;
    Eigen::MatrixXd Bz = aMatrixB * Qz;

    Eigen::VectorXd Px = aMatrixA.colPivHouseholderQr().solve(Bx);
    Eigen::VectorXd Py = aMatrixA.colPivHouseholderQr().solve(By);
    Eigen::VectorXd Pz = aMatrixA.colPivHouseholderQr().solve(Bz);

   /* for (int i = 1; i <= poleNum; i++)
    {
       auto newPole = gp_Pnt(Px(i - 1), Py(i - 1), Pz(i - 1));
       checkPnt(newPole);
    }*/
   // creatBSplineCurve();
}

//计算控制顶点偏离值
void AlgCurFPIAByAdjustingControlPoint::calDiffVec(std::vector<gp_Vec>& diffVec)
{
    int poleNum = m_curve->NbPoles();
    diffVec.resize(poleNum, gp_Vec(0, 0, 0));

    for (int i = 1; i <= poleNum; i++)
    {
        auto oriPole = m_originPoles.Value(i);
        auto curPole = m_curve->Pole(i);
        gp_Vec tempVec(oriPole.XYZ()-curPole.XYZ());
        diffVec[i - 1] = tempVec;
    }
}

void AlgCurFPIAByAdjustingControlPoint::calFairVec(std::vector<gp_Vec>& fairVec)
{
    int poleNum = m_curve->NbPoles();
    fairVec.resize(poleNum, gp_Vec(0, 0, 0));
    for (int i = 0; i < poleNum; i++)
    {
        gp_Vec tempVec(0, 0, 0);
        for (int j = 0; j < poleNum; j++)
        {
            auto curPole = m_curve->Pole(j + 1);
            double fairParam = m_fairFunctional->getValue(i, j);
            tempVec += gp_Vec(curPole.XYZ() * fairParam);
        }
        fairVec[i] = tempVec;
    }
}

void AlgCurFPIAByAdjustingControlPoint::checkVec(const gp_Vec& aVec)
{
    std::cout << "[" << aVec.X() << ", " << aVec.Y() << ", " << aVec.Z() << "]\n";
}

void AlgCurFPIAByAdjustingControlPoint::checkPnt(const gp_Pnt& aPnt)
{
    std::cout << "(" << aPnt.X() << ", " << aPnt.Y() << ", " << aPnt.Z() << ")\n";
}

void AlgCurFPIAByAdjustingControlPoint::iterative()
{
    const auto curPoles = m_curve->Poles();
    int poleNum = m_curve->NbPoles();
    std::vector<gp_Vec> diffVec, fairVce;
    calDiffVec(diffVec);
    calFairVec(fairVce);
    std::vector<gp_Pnt> newPoles(poleNum, gp_Pnt(0, 0, 0));

    //首尾两个不动
    for (int i = 1; i < poleNum-1; i++)
    {
        double mu = m_mus[i];
        double weight = m_weights[i];
        auto curPole = curPoles.Value(i + 1);
        gp_Vec addDiffVec = (1 - weight) * diffVec[i] - weight * fairVce[i];
        gp_Pnt newPole = gp_Pnt(curPole.XYZ() + mu * addDiffVec.XYZ());
        newPoles[i] = newPole;
    }

    for (int i = 1; i < poleNum-1; i++)
    {
        m_curve->SetPole(i + 1, newPoles[i]);
    }
}

void AlgCurFPIAByAdjustingControlPoint::localIterative(const std::vector<size_t>& indexList)
{
    const auto curPoles = m_curve->Poles();
    int poleNum = m_curve->NbPoles();
    std::vector<gp_Vec> diffVec, fairVce;
    calDiffVec(diffVec);
    calFairVec(fairVce);
    std::vector<gp_Pnt> newPoles(poleNum, gp_Pnt(0, 0, 0));

    //首尾两个不动
    for (int i = 0; i < indexList.size(); i++)
    {
        double mu = m_mus[indexList[i]];
        double weight = m_weights[indexList[i]];
        auto curPole = curPoles.Value(indexList[i] + 1);
        gp_Vec addDiffVec = (1 - weight) * diffVec[indexList[i]] - weight * fairVce[indexList[i]];
        gp_Pnt newPole = gp_Pnt(curPole.XYZ() + mu * addDiffVec.XYZ());
        newPoles[indexList[i]] = newPole;
    }

    for (int i = 0; i < indexList.size(); i++)
    {
        m_curve->SetPole(indexList[i] + 1, newPoles[indexList[i]]);
    }
}

void AlgCurFPIAByAdjustingControlPoint::checkWeightInterval(Eigen::MatrixXd aMatrix)
{
    int row = aMatrix.rows();
    int col = aMatrix.cols();
    std::vector<double> maxWeight(row, 0);

    for (int i = 0; i < row; i++)
    {
        double maxEntity = 0;
        for (int j = 0; j < col; j++)
        {
            if (abs(aMatrix(i,j)> maxEntity))
            {
                maxEntity = abs(aMatrix(i, j));
            }
        }
        maxWeight[i] = 1.0 / 2.0 / row / maxEntity;
        //std::cout << maxWeight[i] << std::endl;
    }
}

void AlgCurFPIAByAdjustingControlPoint::calMuParam()
{
    int poleNum = m_originPoles.Size();

    //Omega
    Eigen::MatrixXd aMatrixOmega = Eigen::MatrixXd::Zero(poleNum, poleNum);
    for (int i = 0; i < poleNum; i++)
    {
        aMatrixOmega(i, i) = m_weights[i];
    }
    //std::cout << aMatrixOmega << std::endl;

    //Dr
    Eigen::MatrixXd aMatrixDr = Eigen::MatrixXd::Zero(poleNum, poleNum);
    for (int i = 0; i < poleNum; i++)
    {
        for (int j = 0; j < poleNum; j++)
        {
            aMatrixDr(i, j) = m_fairFunctional->getValue(i, j);
        }
    }

    //I
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(poleNum, poleNum);

    //A
    Eigen::MatrixXd aMatrixA = identity - aMatrixOmega + aMatrixOmega * aMatrixDr;
    //std::cout << "============mu================" << std::endl;

    //mu
    m_mus.resize(poleNum, 0);
    for (int i = 0; i < poleNum; i++)
    {
        double temp = 0;
        for (int j = 0; j < poleNum; j++)
        {
            temp += abs(aMatrixA(i, j));
        }
        m_mus[i] = (1.0 / temp);
        //std::cout << m_mus[i] << std::endl;
    }
}

double AlgCurFPIAByAdjustingControlPoint::iterativeDiff(const TColgp_Array1OfPnt& lastPoles, const TColgp_Array1OfPnt& curPoles)
{
    double diff = 0;
    for (int i = 1; i <= lastPoles.Size(); i++)
    {
        double dist = lastPoles.Value(i).Distance(curPoles.Value(i));
        if (diff < dist)
        {
            diff = dist;
        }
    }
    //std::cout << diff;
    return diff;
}

void AlgCurFPIAByAdjustingControlPoint::poleDistance(const TColgp_Array1OfPnt& poles1, const TColgp_Array1OfPnt& poles2)
{
    if (poles1.Size() != poles2.Size())
    {
        std::cout << "控制点有错误！\n";
    }

    int num = poles1.Size();
    double maxDist = -1;
    double avgDist = 0;
    for (int i = 1; i <= num; i++)
    {
        auto p1 = poles1.Value(i);
        auto p2 = poles2.Value(i);
        double dist = p1.Distance(p2);
        //std::cout << dist << std::endl;
        if (dist > maxDist)
        {
            maxDist = dist;
        }
        avgDist += p1.SquareDistance(p2);
    }

    avgDist /= (double)num;
    avgDist = sqrt(avgDist);

}
