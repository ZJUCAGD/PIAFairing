#include "AlgSurfFPIAByAdjustingControlPoint.h"

AlgSurfFPIAByAdjustingControlPoint::~AlgSurfFPIAByAdjustingControlPoint()
{
}

bool AlgSurfFPIAByAdjustingControlPoint::Init(const Handle(Geom_BSplineSurface)& surface)
{
    if (surface.IsNull())
    {
        return false;
    }
    m_surface = surface;
    m_weights.resize(surface->NbUPoles(),std::vector<double>(surface->NbVPoles(), 0.0));

    m_originPoles = TColgp_Array2OfPnt(surface->Poles());

    auto uKnots = surface->UKnots();
    auto uMults = surface->UMultiplicities();

    auto vKnots = surface->VKnots();
    auto vMults = surface->VMultiplicities();

    auto uFairFunctional1 = new fairFunctionalSnd(uKnots, vKnots, uMults, vMults, surface->UDegree(), surface->VDegree());
    uFairFunctional1->getAllValue(m_FairFuncValue);

    return true;
}

void AlgSurfFPIAByAdjustingControlPoint::execute(double maxTimes, double error)
{
    int time = 0;
    auto lastPoles = m_surface->Poles();
    //auto curPoles = lastPoles;
    double diff = 0;
    while (time < maxTimes)
    {
        iterative();
        //curPoles = m_surface->Poles();
        //diff = iterativeDiff(lastPoles, curPoles);
        //std::cout << diff << std::endl;
        //if (diff < 1e-12)
        //{
        //    /*std::cout << diff << std::endl;
        //    std::cout << time << std::endl;

        //    double sum = 0;
        //    for (int i = 1; i <= lastPoles.UpperRow(); i++)
        //    {
        //        for (int j = 1; j <= lastPoles.UpperCol(); j++)
        //        {
        //            sum += lastPoles.Value(i, j).SquareDistance(curPoles.Value(i, j));
        //        }
        //    }

        //    sum /= (double)(lastPoles.UpperRow() * lastPoles.UpperCol());
        //    std::cout <<"平均" << sqrt(sum) << std::endl;*/

        //    //auto poles = m_surface->Poles();
        //    //std::cout << "==============result===============\n";
        //    /*for (auto& p : poles)
        //    {
        //        checkPnt(p);
        //    }*/
        //    return;
        //}
        //lastPoles = TColgp_Array2OfPnt(curPoles);
        //std::cout << time << std::endl;
        time++;
    }

    auto curPoles = m_surface->Poles();
    double sum = 0;
    double maxtemp = -1;
    for (int i = 1; i <= lastPoles.UpperRow(); i++)
    {
        for (int j = 1; j <= lastPoles.UpperCol(); j++)
        {
            double temp = lastPoles.Value(i, j).SquareDistance(curPoles.Value(i, j));
            sum += temp;
            if (maxtemp < temp)
            {
                maxtemp = temp;
            }
        }
    }

    sum /= (double)(lastPoles.UpperRow() * lastPoles.UpperCol());
    //std::cout << "平均" << sqrt(sum) << std::endl;
    //std::cout << "最大" << sqrt(maxtemp) << std::endl;
}

void AlgSurfFPIAByAdjustingControlPoint::execute(const std::vector<size_t>& indexList, const std::vector<double>& weights, int maxTimes)
{
    if (weights.size() != indexList.size())
    {
        std::cout << "请输入正确的权重.\n";
        return;
    }

    int uNum = m_surface->NbUPoles();
    int vNum = m_surface->NbVPoles();

    for (int i = 0; i < indexList.size(); i++)
    {
        int row = indexList[i] / vNum;
        int col = indexList[i] % vNum;
        m_weights[row][col] = weights[i];
    }

    calMuParam();

    int time = 0;
    auto lastPoles = m_surface->Poles();

    while (time < maxTimes)
    {
        localIterative(indexList);

        auto curPoles = m_surface->Poles();

        double diff = iterativeDiff(lastPoles, curPoles);
        if (diff < 1e-6)
        {
            //std::cout << time << std::endl;
            return;
        }
        lastPoles = TColgp_Array2OfPnt(curPoles);
        time++;
    }
}

void AlgSurfFPIAByAdjustingControlPoint::setWeights(const std::vector<size_t>& indexList, const std::vector<double>& weights)
{
    if (weights.size() != indexList.size())
    {
        std::cout << "请输入正确的权重.\n";
        return;
    }

    int uNum = m_surface->NbUPoles();
    int vNum = m_surface->NbVPoles();

    for (int i = 0; i < indexList.size(); i++)
    {
        int row = indexList[i] / vNum;
        int col = indexList[i] % vNum;
        m_weights[row][col] = weights[i];
    }

    calMuParam();
}

Handle(Geom_BSplineSurface) AlgSurfFPIAByAdjustingControlPoint::getResult()
{
	return m_surface;
}

void AlgSurfFPIAByAdjustingControlPoint::calDiffVec(std::vector<std::vector<gp_Vec>>& diffVec)
{
    int rowNum = m_originPoles.UpperRow();
    int colNum = m_originPoles.UpperCol();
    diffVec.resize(rowNum, std::vector<gp_Vec>(colNum, gp_Vec(0, 0, 0)));

    for (int i = 1; i <= rowNum; i++)
    {
        for (int j = 1; j <= colNum; j++)
        {
            auto oriPole = m_originPoles.Value(i,j);
            auto curPole = m_surface->Pole(i, j);
            gp_Vec tempVec(oriPole.XYZ() - curPole.XYZ());
            diffVec[i - 1][j - 1] = tempVec;
        }
    }
}

void AlgSurfFPIAByAdjustingControlPoint::calFairVec(std::vector<std::vector<gp_Vec>>& fairVec)
{
    int rowNum = m_originPoles.UpperRow();
    int colNum = m_originPoles.UpperCol();
    int allNum = rowNum * colNum;

    fairVec.resize(rowNum, std::vector<gp_Vec>(colNum, gp_Vec(0, 0, 0)));
    for (int i = 0; i < allNum; i++)
    {
        //对每个控制点
        gp_Vec tempVec(0, 0, 0);
        int row = i / colNum;
        int col = i % colNum;

        for (int j = 0; j < allNum; j++)
        {
            int tempRow = j / colNum;
            int tempCol = j % colNum;

            auto curPole = m_surface->Pole(tempRow + 1, tempCol + 1);
            double fairParam = m_FairFuncValue[i][j];
            tempVec += gp_Vec(curPole.XYZ() * fairParam);
        }
        fairVec[row][col] = tempVec;
    }
}

void AlgSurfFPIAByAdjustingControlPoint::checkVec(const gp_Vec& aVec)
{
    std::cout << "(" << aVec.X() << ", " << aVec.Y() << ", " << aVec.Z() << ")\n";
}

void AlgSurfFPIAByAdjustingControlPoint::checkPnt(const gp_Pnt& aPnt)
{
    std::cout << "(" << aPnt.X() << ", " << aPnt.Y() << ", " << aPnt.Z() << ")\n";
}

void AlgSurfFPIAByAdjustingControlPoint::iterative()
{
    const auto curPoles = m_surface->Poles();
    int rowNum = m_originPoles.UpperRow();
    int colNum = m_originPoles.UpperCol();

    std::vector<std::vector<gp_Vec>> diffVec, fairVce;
    calDiffVec(diffVec);
    calFairVec(fairVce);
    std::vector<std::vector<gp_Pnt>> newPoles(rowNum, std::vector<gp_Pnt>(colNum, gp_Pnt(0, 0, 0)));

    //首尾两个不动
    for (int i = 1; i < rowNum - 1; i++)
    {    
        for (int j = 1; j < colNum - 1; j++)
        {
            double mu = m_mus[i][j];
            double weight = m_weights[i][j];
            auto curPole = curPoles.Value(i + 1, j+1);
            gp_Vec addDiffVec = (1 - weight) * diffVec[i][j] - weight * fairVce[i][j];
            //std::cout << addDiffVec.X() << std::endl;
            gp_Pnt newPole = gp_Pnt(curPole.XYZ() + mu * addDiffVec.XYZ());
            newPoles[i][j] = newPole;
        }
    }

    for (int i = 1; i < rowNum - 1; i++)
    {
        for (int j = 1; j < colNum - 1; j++)
        {
            m_surface->SetPole(i + 1,j + 1, newPoles[i][j]);
        }
    }
}

void AlgSurfFPIAByAdjustingControlPoint::localIterative(const std::vector<size_t>& indexList)
{
    const auto curPoles = m_surface->Poles();
    int rowNum = m_originPoles.UpperRow();
    int colNum = m_originPoles.UpperCol();

    std::vector<std::vector<gp_Vec>> diffVec, fairVce;
    calDiffVec(diffVec);
    calFairVec(fairVce);
    std::vector<std::vector<gp_Pnt>> newPoles(rowNum, std::vector<gp_Pnt>(colNum, gp_Pnt(0, 0, 0)));

    for (int i = 0; i < indexList.size(); i++)
    {
        int row = indexList[i] / colNum;
        int col = indexList[i] % colNum;
        double mu = m_mus[row][col];
        double weight = m_weights[row][col];
        auto curPole = curPoles.Value(row + 1, col + 1);
        gp_Vec addDiffVec = (1 - weight) * diffVec[row][col] - weight * fairVce[row][col];
        gp_Pnt newPole = gp_Pnt(curPole.XYZ() + mu * addDiffVec.XYZ());
        newPoles[row][col] = newPole;
    }

    for (int i = 0; i < indexList.size(); i++)
    {
        int row = indexList[i] / colNum;
        int col = indexList[i] % colNum;
        m_surface->SetPole(row + 1, col + 1, newPoles[row][col]);
    }
}

void AlgSurfFPIAByAdjustingControlPoint::calMuParam()
{
    int poleNum = m_originPoles.UpperRow() * m_originPoles.UpperCol();

    //Omega
    Eigen::MatrixXd aMatrixOmega = Eigen::MatrixXd::Zero(poleNum, poleNum);
    int index = 0;
    for (int i = 0; i < m_weights.size(); i++)
    {
        for (int j = 0; j < m_weights[i].size(); j++)
        {
            aMatrixOmega(index, index) = m_weights[i][j];
            index++;
        }
    }

    //Dr
    Eigen::MatrixXd aMatrixDr = Eigen::MatrixXd::Zero(poleNum, poleNum);
    for (int i = 0; i < poleNum; i++)
    {
        for (int j = 0; j < poleNum; j++)
        {
            aMatrixDr(i, j) = m_FairFuncValue[i][j];
            //std::cout << aMatrixDr(i, j) << " ";
        }
        //std::cout << "\n";
    }

    //I
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(poleNum, poleNum);

    //A
    Eigen::MatrixXd aMatrixA = identity - aMatrixOmega + aMatrixOmega * aMatrixDr;

    //mu
    std::vector<double> mus(poleNum,0);
    for (int i = 0; i < poleNum; i++)
    {
        double temp = 0;
        for (int j = 0; j < poleNum; j++)
        {
            temp += abs(aMatrixA(i, j));
        }
        mus[i] = (1.0 / temp);
    }

    m_mus.resize(m_originPoles.UpperRow(), std::vector<double>(m_originPoles.UpperCol(), 1e-3));
    index = 0;
    for (int i = 0; i < m_mus.size(); i++)
    {
        for (int j = 0; j < m_mus[i].size(); j++)
        {
            m_mus[i][j] = mus[index];
            index++;
        }
    }
}

double AlgSurfFPIAByAdjustingControlPoint::iterativeDiff(const TColgp_Array2OfPnt& lastPoles, const TColgp_Array2OfPnt& curPoles)
{
    double maxDiff = 0;
    for (int i = 1; i <= lastPoles.UpperRow(); i++)
    {
        for (int j = 1; j <= lastPoles.UpperCol(); j++)
        {
            double dist = lastPoles.Value(i, j).Distance(curPoles.Value(i, j));
            if (maxDiff < dist)
            {
                maxDiff = dist;
            }
        }
    }
    return maxDiff;
}
