#include "fairFunctional.h"

fairFunctional::fairFunctional(const std::vector<double>& knotVec, int degree, int r)
{
	for (int i = 0; i < knotVec.size() - 1; i++)
	{
		if (knotVec[i] != knotVec[i + 1])
		{
			m_noRepeatKnotVec.push_back(knotVec[i]);
		}
	}
	m_noRepeatKnotVec.push_back(knotVec.back());
	
	m_knot = knotVec;
	m_degree = degree;
	int poleNum = knotVec.size() - degree - 1;
	m_value.resize(poleNum, std::vector<double>(poleNum,0));

	for (int i = 0; i < poleNum; i++)
	{
		for (int j = 0; j < poleNum; j++)
		{
			m_value[i][j] = calValue(r, i, j);
		}
	}
}

fairFunctional::fairFunctional(const TColStd_Array1OfReal& knots, const TColStd_Array1OfInteger& mults, int degree, int r)
{
	for (int i = 1; i <= knots.Size(); i++)
	{
		double knot = knots.Value(i);
		int length = mults.Value(i);
		m_noRepeatKnotVec.push_back(knot);
		for (int j = 1; j <= length; j++)
		{
			m_knot.push_back(knot);
		}
	}
	/*std::cout << "=====konts=====\n";
	for (auto& n : m_noRepeatKnotVec)
	{
		std::cout << n << std::endl;
	}*/

	m_degree = degree;
	int poleNum = m_knot.size() - degree - 1;
	m_value.resize(poleNum, std::vector<double>(poleNum, 0));

	for (int i = 0; i < poleNum; i++)
	{
		for (int j = 0; j < poleNum; j++)
		{
			m_value[i][j] = calValue(r, i, j);
			//std::cout << m_value[i][j] << std::endl;
		}
	}
}

double fairFunctional::func(int r, int i1, int i2, double u)
{
	return derOneBasisFun(r, i1, u) * derOneBasisFun(r, i2, u);
}

double fairFunctional::derOneBasisFun(int r, int i, double u)
{
	std::vector<double> ders;
	ders.resize(r + 1, 0);
	int p = m_degree;
	if (u < m_knot[i] || u >= m_knot[i + p + 1])   /* Local property */
	{
		for (int k = 0; k <= r; k++)
		{
			ders[k] = 0.0;
		}
		return ders[r];
	}

	std::vector<std::vector<double>> N(p + 1, std::vector<double>(p + 1, 0));

	for (int j = 0; j <= p; j++)   /* Initialize zeroth-degree functs */
	{
		if (u >= m_knot[i + j] && u < m_knot[i + j + 1])
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
			saved = ((u - m_knot[i]) * N[0][k - 1]) / (m_knot[i + k] - m_knot[i]);
		}
		for (int j = 0; j < p - k + 1; j++)
		{
			Uleft = m_knot[i + j + 1];
			Uright = m_knot[i + j + k + 1];
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
			else saved = ND[0] / (m_knot[i + p - k + jj] - m_knot[i]);
			for (int j = 0; j < k - jj + 1; j++)
			{
				Uleft = m_knot[i + j + 1];
				Uright = m_knot[i + j + p + jj - k + 1];
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
	return ders[r];
}

double fairFunctional::getValue(int i1, int i2)
{
	return m_value[i1][i2];
}

void fairFunctional::getAllValue(std::vector<std::vector<double>>& value)
{
	value = m_value;
}

double fairFunctional::integration(double a, double b, int r, int i1, int i2)
{
	double sum = 0.0;
	for (int i = 0; i < 5; i++)
	{
		sum += w[i] * func(r, i1, i2, (x[i] * (b - a) + b + a) / 2);
	}
	double integral = ((b - a) / 2) * sum;
	return integral;
}

double fairFunctional::calValue(int r, int i1, int i2)
{
	double sum = 0;
	for (int i = 0; i < m_noRepeatKnotVec.size() - 1; i++)
	{
		sum += integration(m_noRepeatKnotVec[i], m_noRepeatKnotVec[i + 1], r, i1, i2);
	}
	return sum;
}
