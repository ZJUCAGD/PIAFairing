#include "fairFunctionalSnd.h"
#include<ctime>

fairFunctionalSnd::fairFunctionalSnd(const TColStd_Array1OfReal& UKnots, const TColStd_Array1OfReal& VKnots, const TColStd_Array1OfInteger& UMults, const TColStd_Array1OfInteger& VMults, int uDegree, int vDegree)
    : m_noRepeatUKnotVec(UKnots), m_noRepeatVKnotVec(VKnots), m_uDegree(uDegree), m_vDegree(vDegree)
{
    int flatUNum = 0;
    for (auto& u : UMults)
    {
        flatUNum += u;
    }

    int flatVNum = 0;
    for (auto& v : VMults)
    {
        flatVNum += v;
    }

    m_flatUKnotVec = TColStd_Array1OfReal(1, flatUNum);
    m_flatVKnotVec = TColStd_Array1OfReal(1, flatVNum);

    int index = 0;
    for (int i = 1; i <= UKnots.Size(); i++)
    {
        int num = UMults.Value(i);
        double value = UKnots.Value(i);
        for (int j = 1; j <= num; j++)
        {
            m_flatUKnotVec.SetValue(++index, value);
        }
    }

    index = 0;
    for (int i = 1; i <= VKnots.Size(); i++)
    {
        int num = VMults.Value(i);
        double value = VKnots.Value(i);
        for (int j = 1; j <= num; j++)
        {
            m_flatVKnotVec.SetValue(++index, value);
        }
    }

    m_rowNum = flatUNum - m_uDegree - 1;
    m_colNum = flatVNum - m_vDegree - 1;

    int allPoleNum = m_rowNum * m_colNum;
    m_value.resize(allPoleNum, std::vector<double>(allPoleNum, 0));



    m_derUValues.resize((m_noRepeatUKnotVec.Size() - 1) * 5);
    m_derVValues.resize((m_noRepeatVKnotVec.Size() - 1) * 5);

    //m_derUBasis.resize(m_rowNum, std::vector<gp_Vec>((m_noRepeatUKnotVec.Size() - 1) * 5, gp_Vec(0, 0, 0)));
    //m_derVBasis.resize(m_colNum, std::vector<gp_Vec>((m_noRepeatVKnotVec.Size() - 1) * 5, gp_Vec(0, 0, 0)));

    //for u
    int indexForU = 0;
    for (int i = 1; i < m_noRepeatUKnotVec.Size(); i++)
    {
        double a = m_noRepeatUKnotVec.Value(i);
        double b = m_noRepeatUKnotVec.Value(i + 1);

        for (int j = 0; j < 5; j++) {
            double param = (b + a) / 2.0 + x[j] * (b - a) / 2.0;
            //std::cout << indexForU << std::endl;
            calParamValue(indexForU++, param, true);
        }
    }

    //for v
    int indexForV = 0;
    for (int i = 1; i < m_noRepeatVKnotVec.Size(); i++)
    {
        double a = m_noRepeatVKnotVec.Value(i);
        double b = m_noRepeatVKnotVec.Value(i + 1);

        for (int j = 0; j < 5; j++) {
            double param = (b + a) / 2.0 + x[j] * (b - a) / 2.0;
            //std::cout << indexForV << std::endl;
            calParamValue(indexForV++, param, false);
        }
    }

    /*for (int i = 0; i < m_derUBasis.size(); i++)
    {
        for (int j = 0; j < m_derUBasis[i].size(); j++)
        {
            std::cout << m_derUBasis[i][j].X() << " ";
        }
        std::cout << "\n";
    }*/

    //auto astart = clock();
    for (int i = 0; i < allPoleNum; i++)
    {
        for (int j = i; j < allPoleNum; j++)
        {
            m_value[i][j] = calValue(i, j);
            m_value[j][i] = m_value[i][j];

        }
        //std::cout << i << std::endl;

        /*if (i == 10)
        {
            auto anend = clock();
            std::cout << (anend - astart)/1000.0 << std::endl;
        }*/
    }
}

double fairFunctionalSnd::getValue(int i1, int i2)
{
	return m_value[i1][i2];
}

void fairFunctionalSnd::getAllValue(std::vector<std::vector<double>>& value)
{
    value = m_value;
}

double fairFunctionalSnd::func(double u, double v, int i1, int i2)
{
    //Standard_Integer theFirstUIndex;
    //math_Matrix BsplineBasisU(1, 3, 1, m_uDegree + 1);
    //BSplCLib::EvalBsplineBasis(2, m_uDegree + 1, m_flatUKnotVec, u, theFirstUIndex, BsplineBasisU);

    //Standard_Integer theFirstVIndex;
    //math_Matrix BsplineBasisV(1, 3, 1, m_vDegree + 1);
    //BSplCLib::EvalBsplineBasis(2, m_vDegree + 1, m_flatVKnotVec, v, theFirstVIndex, BsplineBasisV);

    //for i1
    int rowIndex1 = i1 / m_colNum;
    int colIndex1 = i1 % m_colNum;

    //int col1 = rowIndex1 - theFirstUIndex + 1;
    //int col2 = colIndex1 - theFirstVIndex + 1;

    //if (col1 < 0 || col1 > m_uDegree || col2 < 0 || col2 > m_vDegree)
    //{
    //    return 0;
    //}

    //for i2
    int rowIndex2 = i2 / m_colNum;
    int colIndex2 = i2 % m_colNum;

    //int col3 = rowIndex2 - theFirstUIndex + 1;
    //int col4 = colIndex2 - theFirstVIndex + 1;

    //if (col3 < 0 || col3 > m_uDegree || col4 < 0 || col4 > m_vDegree)
    //{
    //    return 0;
    //}

    //col1++;
    //col2++;
    //col3++;
    //col4++;

    /*double N1UU = BsplineBasisU(3, col1) * BsplineBasisV(1, col2);
    double N1UV = BsplineBasisU(2, col1) * BsplineBasisV(2, col2);
    double N1VV = BsplineBasisU(1, col1) * BsplineBasisV(3, col2);

    double N2UU = BsplineBasisU(3, col3) * BsplineBasisV(1, col4);
    double N2UV = BsplineBasisU(2, col3) * BsplineBasisV(2, col4);
    double N2VV = BsplineBasisU(1, col3) * BsplineBasisV(3, col4);*/

    auto derU1 = derOneBasisFun(2, rowIndex1, u, true);
    auto derV1 = derOneBasisFun(2, colIndex1, v, false);

    auto derU2 = derOneBasisFun(2, rowIndex2, u, true);
    auto derV2 = derOneBasisFun(2, colIndex2, v, false);

    double N1UU = derU1[2] * derV1[0];
    double N1UV = derU1[1] * derV1[1];
    double N1VV = derU1[0] * derV1[2];

    double N2UU = derU2[2] * derV2[0];
    double N2UV = derU2[1] * derV2[1];
    double N2VV = derU2[0] * derV2[2];

    //N1_uu(u,v)N2_uu(u,v) + 2 * N1_uv(u,v)N2_uv(u,v) + N1_vv(u,v)N2_vv(u,v)
    return N1UU * N2UU + 2 * N1UV * N2UV + N1VV * N2VV;
}

double fairFunctionalSnd::calValue(int i1, int i2)//下标从0开始
{
    double sum = 0;

    for (int i = 1; i < m_noRepeatUKnotVec.Size(); i++)
    {
        double ab = m_noRepeatUKnotVec.Value(i + 1) - m_noRepeatUKnotVec.Value(i);
        for (int j = 1; j < m_noRepeatVKnotVec.Size(); j++)
        {
            //sum += integration(m_noRepeatUKnotVec.Value(i), m_noRepeatUKnotVec.Value(i + 1), m_noRepeatVKnotVec.Value(j), m_noRepeatVKnotVec.Value(j + 1), i1, i2);
            double cd = m_noRepeatVKnotVec.Value(j + 1) - m_noRepeatVKnotVec.Value(j);
            sum += (packageKnot(i-1, j-1, i1, i2) * ab * cd);
        }
    }

    return sum;
}

double fairFunctionalSnd::integration(double a, double b, double c, double d, int i1, int i2)
{
    std::vector<double> t1 = { 0,0,0,0,0 };
    std::vector<double> t2 = { 0,0,0,0,0 };

    for (int i = 0; i < 5; i++) {
        t1[i] = (b + a) / 2. + x[i] * (b - a) / 2.;
        t2[i] = (d + c) / 2. + x[i] * (d - c) / 2.;
    }

    double s = 0.;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            double fun = func(t1[i], t2[j], i1, i2);
            //std::cout << t1[i] <<" " << t2[j] << std::endl;
            s = s + w[i] * w[j] * fun;
        }
    }
    //std::cout << "============\n\n";

    return s * (b - a) * (d - c) / 4.;
}

double fairFunctionalSnd::packageKnot(int indexU, int indexV, int i1, int i2)
{
    double s = 0.;
    for (int i = 0; i < 5; i++) {
        int index1 = indexU * 5 + i;
        for (int j = 0; j < 5; j++) {
            int index2 = indexV * 5 + j;
            double fun = packageParam(index1, index2, i1, i2);
            s = s + w[i] * w[j] * fun;
        }
    }
    return s/ 4.0;
}

double fairFunctionalSnd::packageParam(int indexU, int indexV, int i1, int i2)
{
    ////for i1
    //int rowIndex1 = i1 / m_colNum;
    //int colIndex1 = i1 % m_colNum;

    //auto tempDerUBasis1 = m_derUBasis[rowIndex1];
    //auto tempDerVBasis1 = m_derVBasis[colIndex1];

    ////for i2
    //int rowIndex2 = i2 / m_colNum;
    //int colIndex2 = i2 % m_colNum;

    //auto tempDerUBasis2 = m_derUBasis[rowIndex2];
    //auto tempDerVBasis2 = m_derVBasis[colIndex2];

    //double N1UU = tempDerUBasis1[indexU].X() * tempDerVBasis1[indexV].X();
    //double N1UV = tempDerUBasis1[indexU].Y() * tempDerVBasis1[indexV].Y();
    //double N1VV = tempDerUBasis1[indexU].Z() * tempDerVBasis1[indexV].Z();

    //double N2UU = tempDerUBasis2[indexU].X() * tempDerVBasis2[indexV].X();
    //double N2UV = tempDerUBasis2[indexU].Y() * tempDerVBasis2[indexV].Y();
    //double N2VV = tempDerUBasis2[indexU].Z() * tempDerVBasis2[indexV].Z();

    //return N1UU * N2UU + 2 * N1UV * N2UV + N1VV * N2VV;

    auto paramUList = m_derUValues[indexU];
    auto paramVList = m_derVValues[indexV];

    //for i1
    int rowIndex1 = i1 / m_colNum;
    int colIndex1 = i1 % m_colNum;

    int col1 = rowIndex1 + 1 - paramUList.back();
    int col2 = colIndex1 + 1 - paramVList.back();

    if (col1 < 0 || col1 > m_uDegree || col2 < 0 || col2 > m_vDegree)
    {
        return 0;
    }

    //for i2
    int rowIndex2 = i2 / m_colNum;
    int colIndex2 = i2 % m_colNum;

    int col3 = rowIndex2 + 1 - paramUList.back();
    int col4 = colIndex2 + 1 - paramVList.back();

    if (col3 < 0 || col3 > m_uDegree || col4 < 0 || col4 > m_vDegree)
    {
        return 0;
    }

    int uDe = m_uDegree + 1;
    int vDe = m_vDegree + 1;

    double N1UU = paramUList[col1 + 2 * uDe]  *  paramVList[col2];
    double N1UV = paramUList[col1 + uDe] * paramVList[col2 + vDe];
    double N1VV = paramUList[col1]  *  paramVList[col2 + 2 * vDe];

    double N2UU = paramUList[col3 + 2 * uDe] * paramVList[0 * (m_vDegree + 1) + col4];
    double N2UV = paramUList[1 * (m_uDegree + 1) + col3] * paramVList[1 * (m_vDegree + 1) + col4];
    double N2VV = paramUList[0 * (m_uDegree + 1) + col3] * paramVList[2 * (m_vDegree + 1) + col4];

    //N1_uu(u,v)N2_uu(u,v) + 2 * N1_uv(u,v)N2_uv(u,v) + N1_vv(u,v)N2_vv(u,v)
    return N1UU * N2UU + 2 * N1UV * N2UV + N1VV * N2VV;
}

std::vector<double> fairFunctionalSnd::derOneBasisFun(int r, int i, double u, bool isU)
{
    std::vector<double> ders;
    ders.resize(r + 1, 0);

    int p = 0;
    std::vector<double> knots;

    if (isU)
    {
        p = m_uDegree;
        knots.resize(m_flatUKnotVec.Size(), 0);
        for (int k = 0; k < m_flatUKnotVec.Size(); k++)
        {
            knots[k] = m_flatUKnotVec.Value(i + 1);
        }
    }
    else
    {
        p = m_vDegree;
        knots.resize(m_flatVKnotVec.Size(),0);
        for (int k = 0; k < m_flatVKnotVec.Size(); k++)
        {
            knots[k] = m_flatVKnotVec.Value(i + 1);
        }
    }
    
    if (u < knots[i] || u >= knots[i + p + 1])   /* Local property */
    {
        for (int k = 0; k <= r; k++)
        {
            ders[k] = 0.0;
        }
        return ders;
    }

    std::vector<std::vector<double>> N(p + 1, std::vector<double>(p + 1, 0));

    for (int j = 0; j <= p; j++)   /* Initialize zeroth-degree functs */
    {
        if (u >= knots[i + j] && u < knots[i + j + 1])
        {
            N[j][0] = 1.0;
        }
        else
        {
            N[j][0] = 0.0;
        }
    }

    double Uleft = 0.0, Uright = 0.0, saved = 0.0;
    for (int k = 1; k <= p; k++)   /* Compute full triangular table */
    {
        if (N[0][k - 1] == 0.0)
        {
            saved = 0.0;
        }
        else
        {
            saved = ((u - knots[i]) * N[0][k - 1]) / (knots[i + k] - knots[i]);
        }
        for (int j = 0; j < p - k + 1; j++)
        {
            Uleft = knots[i + j + 1];
            Uright = knots[i + j + k + 1];
            if (N[j + 1][k - 1] == 0.0)
            {
                N[j][k] = saved;
                saved = 0.0;
            }
            else
            {
                double temp = N[j + 1][k - 1] / (Uright - Uleft);
                N[j][k] = saved + (Uright - u) * temp;
                saved = (u - Uleft) * temp;
            }
        }
    }
    ders[0] = N[0][p];    /* The function value */
    std::vector<double> ND(r + 1, 0);
    for (int k = 1; k <= r; k++) /* Compute the derivatives */
    {
        for (int j = 0; j <= k; j++)   /* Load appropriate column */
        {
            ND[j] = N[j][p - k];
        }
        for (int jj = 1; jj <= k; jj++)   /* Compute table of width k */
        {
            if (ND[0] == 0.0) saved = 0.0;
            else saved = ND[0] / (knots[i + p - k + jj] - knots[i]);
            for (int j = 0; j < k - jj + 1; j++)
            {
                Uleft = knots[i + j + 1];
                Uright = knots[i + j + p + jj - k + 1];
                if (ND[j + 1] == 0.0)
                {
                    ND[j] = (p - k + jj) * saved;
                    saved = 0.0;
                }
                else
                {
                    double temp = ND[j + 1] / (Uright - Uleft);
                    ND[j] = (p - k + jj) * (saved - temp);
                    saved = temp;
                }
            }
        }
        ders[k] = ND[0]; /* kth derivative */
    }
    return ders;
}

void fairFunctionalSnd::calParamValue(int row, double param, bool isU)
{
    if (isU)
    {
        Standard_Integer theFirstUIndex;
        math_Matrix BsplineBasisU(1, 3, 1, m_uDegree + 1);
        BSplCLib::EvalBsplineBasis(2, m_uDegree + 1, m_flatUKnotVec, param, theFirstUIndex, BsplineBasisU);
        std::vector <double> tempValue(3 * (m_uDegree + 1) + 1);
        
        int index = 0;
        for (int i = 1; i <= 3; i++)
        {
            for (int j = 1; j <= m_uDegree + 1; j++)
            {
                tempValue[index++] = BsplineBasisU(i, j);
            }
        }
        tempValue[index] = theFirstUIndex;
        m_derUValues[row] = tempValue;

        //for (int i = 0; i < m_uDegree + 1; i++)
        //{
        //    int basisIndex = theFirstUIndex + i - 1;//第basisIndex个基函数
        //    m_derUBasis[basisIndex][row] = gp_Vec(BsplineBasisU(3, i+1), BsplineBasisU(2,i+1), BsplineBasisU(1,i+1));//U方向为反序
        //}


    }
    else
    {
        Standard_Integer theFirstVIndex;
        math_Matrix BsplineBasisV(1, 3, 1, m_vDegree + 1);
        BSplCLib::EvalBsplineBasis(2, m_vDegree + 1, m_flatVKnotVec, param, theFirstVIndex, BsplineBasisV);
        std::vector <double> tempValue(3 * (m_vDegree + 1) + 1);

        int index = 0;
        for (int i = 1; i <= 3; i++)
        {
            for (int j = 1; j <= m_vDegree + 1; j++)
            {
                tempValue[index++] = BsplineBasisV(i, j);
            }
        }
        tempValue[index] = theFirstVIndex;
        m_derVValues[row] = tempValue;

        //for (int i = 0; i < m_vDegree + 1; i++)
        //{
        //    int basisIndex = theFirstVIndex + i - 1;//第basisIndex个基函数
        //    m_derVBasis[basisIndex][row] = gp_Vec(BsplineBasisV(1, i + 1), BsplineBasisV(2, i + 1), BsplineBasisV(3, i + 1));//U方向为反序
        //}
    }

}




