//#include "Cauchy_Problem.h"
#include <iostream>
#include "Cauchy_Problem.h"
using namespace std;


// реализация явных схем Рунге-Кутты
// FD_Type - тип разностной схемы
// Time_Begin - начальный момент
// Time_End - конечный момент
// Initial_Conditions - начальные условия
// F - правая часть системы ДУ
	
	vector<double> Cauchy_Problem<double>::ERKs_Scheme_Start(Cauchy_Problem<double>::Difference_Scheme_Type FD_Type,
														          double Time_Begin,
														          double Time_End,
														          const vector<double>& Initial_Conditions,
														          const function<vector<double>(const vector<double>& U, const Point& P)>& F)
	{
		//число шагов явной разностной схемы Рунге-Кутты
        int Stages_Number = 4;

		//число неизвестных
		int Variable_Number = Initial_Conditions.size();

		//шаг по времени
		double h = Time_End - Time_Begin;

		//вектор результата инициализируется начальными условиями
        vector<double> Solution(Variable_Number);
        for (int i = 0; i < Variable_Number; i++) {
            Solution[i] = Initial_Conditions[i];
        }
		
        // выделение памяти под параметры таблицы Бутчера для методов Рунге - Кутты
        vector<double> b, c; 
        vector<vector<double>> a;

        switch (FD_Type){
            //явная схема Эйлера(одношаговый явный метод)
            case Difference_Scheme_Type::ERK1:{
                Stages_Number = 1;
                a = vector<vector<double>>{ {0} };
                b = vector<double > {1};
                c = vector<double > {0.5};
                break;
            }
            //двухшаговый явный метод (ERK2)
            case Difference_Scheme_Type::ERK2:{
                Stages_Number = 2;
                a = vector<vector<double>>{ {0.0, 0.0}, {2.0 / 3.0, 0.0} };
                b = vector<double >{ 1.0 / 4.0, 3.0 / 4.0 };
                c = vector<double >{ 0.0, 2.0 / 3.0 };
                break;
            }
            //четырёхшаговый явный метод (ERK4)
            default:{
                Stages_Number = 4;
                a = vector<vector<double>>{ {0,   0,   0,   0}, 
                                                      {0.5, 0,   0,   0}, 
                                                      {0,   0.5, 0,   0}, 
                                                      {0,   0,   1.0, 0}, };
                b = vector<double >{ 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 };
                c = vector<double >{ 0.0, 0.5, 0.5, 1.0};
                break;
            }
        }  

        //реализация схем ERKs по таблице Бутчера
        vector<vector<double>> w(Stages_Number, vector<double>(Variable_Number));

        for (int k = 0; k < Stages_Number; k++){
            vector<double> adjustment1(Variable_Number);
            for (int l = 0; l < k; l++){
                for (int i = 0; i < Variable_Number; i++) {
                    adjustment1[i] += a[k][l] * w[l][i];
                }
            }

            vector<double> U(Variable_Number);
            for (int i = 0; i < Variable_Number; i++) {
                U[i] = Solution[i] + h * adjustment1[i];
            }
            Point P(Time_Begin + h * c[k], 0, 0);
            w[k] = F(U, P);
        }
        
        vector<double> adjustment2(Variable_Number);
        for (int k = 0; k < Stages_Number; k++){
            for (int i = 0; i < Variable_Number; i++){
                adjustment2[i] = adjustment2[i] + b[k] * w[k][i];
            }
        }

        for (int i = 0; i < Variable_Number; i++) {
            Solution[i] += h * adjustment2[i];
        }

		return Solution;
	}

    vector<vector<double>> Cauchy_Problem<double>::RungeKutta_Start(double Time_Begin, double h, int steps,
        const vector<double>& Initial_Conditions,
        const function<vector<double>(const vector<double>& U, const Point& P)>& F) {

        // Инициализация значений для первых шагов
        vector<vector<double>> Precomputed_Steps(steps, Initial_Conditions);
        vector<double> Solution = Initial_Conditions;

        for (int step = 0; step < steps; ++step) {
            Solution = ERKs_Scheme_Start(Difference_Scheme_Type::ERK4, Time_Begin + step * h, Time_Begin + (step + 1) * h, Solution, F);
            Precomputed_Steps[step] = Solution;
        }

        return Precomputed_Steps;
    }

    vector<double> Cauchy_Problem<double>::Adams_Bashforth_Scheme_Start(Difference_Scheme_Type FD_Type,
        double Time_Begin,
        double Time_End,
        const vector<vector<double>>& Initial_Conditions,

        const function<vector<double>(const vector<double>& U, const Point& P)>& F) {

        int Variable_Number = 2;
        double h = (Time_End - Time_Begin);
        vector<double> Next_Solution(Variable_Number);
        Point P(Time_Begin, 0, 0);

        if (FD_Type == Difference_Scheme_Type::AdamsBashforth2) {
            for (int i = 0; i < Variable_Number; ++i) {
                
                Next_Solution[i] = Initial_Conditions[3][i] + h * (3.0 / 2.0 * F(Initial_Conditions[3], P)[i] - 1.0 / 2.0 * F(Initial_Conditions[2], P)[i]);
            }
        }
        else if (FD_Type == Difference_Scheme_Type::AdamsBashforth3) {
            for (int i = 0; i < Variable_Number; ++i) {
                Next_Solution[i] = Initial_Conditions[3][i] + h * (23.0 / 12.0 * F(Initial_Conditions[3], P)[i] - 16.0 / 12.0 * F(Initial_Conditions[2], P)[i] + 5.0 / 12.0 * F(Initial_Conditions[1], P)[i]);
            }
        }
        else if (FD_Type == Difference_Scheme_Type::AdamsBashforth4) {
            for (int i = 0; i < Variable_Number; ++i) {
                Next_Solution[i] = Initial_Conditions[3][i] + h * (55.0 / 24.0 * F(Initial_Conditions[3], P)[i] - 59.0 / 24.0 * F(Initial_Conditions[2], P)[i] + 37.0 / 24.0 * F(Initial_Conditions[1], P)[i] - 9.0 / 24.0 * F(Initial_Conditions[0], P)[i]);
            }
        }
        else if (FD_Type == Difference_Scheme_Type::AdamsMoulton3) {
            // Итеративное решение для неявной схемы Адамса-Мултона 3
            const int max_iterations = 10;
            const double tolerance = 1e-6;

            vector<double> Previous_Solution = Initial_Conditions[2];
            vector<double> Current_Solution = Initial_Conditions[3];

            for (int iter = 0; iter < max_iterations; ++iter) {
                Point P_next(Time_Begin + h, 0, 0);

                for (int i = 0; i < Variable_Number; ++i) {
                    Current_Solution[i] = Initial_Conditions[3][i] +h / 12.0 * (5.0 * F(Current_Solution, P_next)[i] +8.0 * F(Initial_Conditions[3], P)[i] -F(Initial_Conditions[2], P)[i]);
                }

                double error = 0.0;
                for (int i = 0; i < Variable_Number; ++i) {
                    error += abs(Current_Solution[i] - Previous_Solution[i]);
                }
                if (error < tolerance) {
                    break;
                }

                Previous_Solution = Current_Solution;
            }
            Next_Solution = Current_Solution;
        }

        return Next_Solution;
    }

