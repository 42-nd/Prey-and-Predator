#include "Cauchy_Problem.h"
#include <iomanip>
#include <iostream>

namespace Com_Methods
{
	/// <summary>
	/// реализация явных схем Рунге-Кутты
	/// FD_Type - тип разностной схемы
	/// Time_Begin - начальный момент
	/// Time_End - конечный момент
	/// Initial_Conditions - начальные условия
	/// F - правая часть системы ДУ
	/// </summary>
	std::vector<double> Cauchy_Problem<double>::ERKs_Scheme_Start(Cauchy_Problem<double>::Difference_Scheme_Type FD_Type,
														          double Time_Begin,
														          double Time_End,
														          const std::vector<double>& Initial_Conditions,
														          const std::function<std::vector<double>(const std::vector<double>& U, const Point& P)>& F)
	{
		//число шагов явной разностной схемы Рунге-Кутты
        int Stages_Number = 4;

		//число неизвестных
		int Variable_Number = Initial_Conditions.size();

		//шаг по времени
		double h = Time_End - Time_Begin;

		//вектор результата инициализируется начальными условиями
        std::vector<double> Solution(Variable_Number);
        for (int i = 0; i < Variable_Number; i++)
            Solution[i] = Initial_Conditions[i];
		
        // выделение памяти под параметры таблицы Бутчера для методов Рунге - Кутты
        std::vector<double> b, c; 
        std::vector<std::vector<double>> a;

        switch (FD_Type)
        {
            //явная схема Эйлера(одношаговый явный метод)
            case Difference_Scheme_Type::ERK1:
            {
                Stages_Number = 1;
                a = std::vector<std::vector<double>>{ {0} };
                b = std::vector<double > {1};
                c = std::vector<double > {0.5};
                break;
            }
            //двухшаговый явный метод (ERK2)
            case Difference_Scheme_Type::ERK2:
            {
                Stages_Number = 2;
                a = std::vector<std::vector<double>>{ {0.0,       0.0}, 
                                                      {2.0 / 3.0, 0.0} };
                b = std::vector<double >{ 1.0 / 4.0, 3.0 / 4.0 };
                c = std::vector<double >{ 0.0,       2.0 / 3.0 };
                break;
            }
            //четырёхшаговый явный метод (ERK4)
            default:
            {
                Stages_Number = 4;
                a = std::vector<std::vector<double>>{ {0,   0,   0,   0}, 
                                                      {0.5, 0,   0,   0}, 
                                                      {0,   0.5, 0,   0}, 
                                                      {0,   0,   1.0, 0}, };
                b = std::vector<double >{ 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 };
                c = std::vector<double >{ 0.0, 0.5, 0.5, 1.0};
                break;
            }
        }  

        //реализация схем ERKs по таблице Бутчера
        std::vector<std::vector<double>> w(Stages_Number, std::vector<double>(Variable_Number));

        for (int k = 0; k < Stages_Number; k++)
        {
            std::vector<double> adjustment1(Variable_Number);
            for (int l = 0; l < k; l++)
            {
                for (int i = 0; i < Variable_Number; i++)
                    adjustment1[i] += a[k][l] * w[l][i];
            }

            std::vector<double> U(Variable_Number);
            for (int i = 0; i < Variable_Number; i++)
                U[i] = Solution[i] + h * adjustment1[i];
            Point P(Time_Begin + h * c[k], 0, 0);
            w[k] = F(U, P);
        }
        
        std::vector<double> adjustment2(Variable_Number);
        for (int k = 0; k < Stages_Number; k++)
        {
            for (int i = 0; i < Variable_Number; i++)
            {
                adjustment2[i] = adjustment2[i] + b[k] * w[k][i];
            }
        }

        for (int i = 0; i < Variable_Number; i++)
            Solution[i] += h * adjustment2[i];

		return Solution;
	}
}

//------------------------------------------------------------------------------------------------------------------------------

//ниже представлено численное решение системы уравнений
// u' = v
// v' = exp(x)
// u(0) = 1, v(0) = 1

int main()
{
	Com_Methods::Cauchy_Problem<double> Problem;

    //число неизвестных
    int Variable_Number = 2;
    
    //правая часть
    std::function<std::vector<double>(const std::vector<double>& U, const Com_Methods::Point& P)> Source_Term = [&](const std::vector<double>& U, const Com_Methods::Point& P)
    {
        std::vector<double> F(Variable_Number);
        F[0] = U[1];
        F[1] = exp(P.x());
        return F;
    };

    //начальные условия
    std::function<std::vector<double>(const Com_Methods::Point& P)> Initial_Solution = [&](const Com_Methods::Point& P)
    {
        std::vector<double> U(Variable_Number);
        U[0] = exp(P.x());
        U[1] = exp(P.x());
        return U;
    };

    //аналитическое решение
    std::function<std::vector<double>(const Com_Methods::Point& P)> Analitical_Solution = [&](const Com_Methods::Point& P)
    {
        std::vector<double> U(Variable_Number);
        U[0] = exp(P.x());
        U[1] = exp(P.x());
        return U;
    };

    //отрезок, на котором отыскивается численное решение
    double X0 = 0, Xn = 1;

    //мелкость разбиения отрезка и шаг
    int n = 10;
    double h = (Xn - X0) / n;

    for (int i = 1; i <= n; i++)
    {
        Com_Methods::Point P_Begin((i - 1) * h, 0, 0), P_End(i * h, 0, 0);

        //начальные условия
        auto U0 = Initial_Solution(P_Begin);

        //аналитическое решение задачи
        auto U = Analitical_Solution(P_End);

        //решение задачи Коши по схеме Рунге-Кутты 4
        auto Result = Problem.ERKs_Scheme_Start(Com_Methods::Cauchy_Problem<double>::Difference_Scheme_Type::ERK4, P_Begin.x(), P_End.x(), U0, Source_Term);

        std::cout << std::scientific << std::left << "\n\nCurrent time: " << P_End.x();

        std::cout << "\nAnalitical solution" << " | Numerical solution | " << "Error";

        for (int j = 0; j < Variable_Number; j++)
            std::cout << "\n" << std::setw(19) << U[j] << " | " << std::setw(18) << Result[j] << " | " << fabs(U[j] - Result[j]);
    }
}
