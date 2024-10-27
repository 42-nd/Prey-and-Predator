#include <iomanip>
#include <iostream>
#include "Cauchy_Problem.h"
using namespace std;

//------------------------------------------------------------------------------------------------------------------------------

// Нейтральное равновесие хищников и жертв
// u' = α * u - β * u * v
// v' = δ * u * v - γ * v
// Начальные условия: u(0) = 4212, v(0) = 552

int main() {
    Cauchy_Problem<double> Problem;

    // Количество переменных
    int Variable_Number = 2;

    // Параметры модели
    double alpha = 4.01; // Коэффициент рождаемости жертв
    double beta = 0.01; // Коэффициент успешной охоты
    double gamma = 19; // Коэффициент естественной убыли хищников
    double delta = 0.01; // Коэффициент воспроизводства хищников

    // Правая часть
    function<vector<double>(const vector<double>& U, const Point& P)> Source_Term = [&](const vector<double>& U, const Point& P) {
        vector<double> F(Variable_Number);
        F[0] = alpha * U[0] - beta * U[0] * U[1]; // Изменение количества жертв
        F[1] = delta * U[0] * U[1] - gamma * U[1]; // Изменение количества хищников
        return F;
    };

    // Начальные условия
    vector<double> U0_erk = { 4212, 552 }; // Начальные условия: x0 = 4212, y0 = 552
    vector<double> U0_adams = { 4212, 552 };
    // Отрезок, на котором отыскивается численное решение
    double X0 = 0, Xn = 10; // Можете увеличить длину отрезка для более долгосрочного моделирования

    // Мелкость разбиения отрезка и шаг
    int n = 100; // Увеличьте количество шагов для лучшей точности
    double h = (Xn - X0) / n;
    for (int i = 1; i <= n; i++) {
        Point P_Begin((i - 1) * h, 0, 0), P_End(i * h, 0, 0);

        // Решение задачи Коши по схеме Рунге-Кутты 4
        auto Result_ERK = Problem.ERKs_Scheme_Start(
            Cauchy_Problem<double>::Difference_Scheme_Type::ERK4,
            P_Begin.x(), P_End.x(), U0_erk, Source_Term
        );
        U0_erk = Result_ERK;
        // Решение задачи Коши по схеме Адамса-Башфорта 4-го порядка
        auto Result_Adams = Problem.Adams_Bashforth_Scheme_Start(
            Cauchy_Problem<double>::Difference_Scheme_Type::AdamsBashforth4,
            P_Begin.x(), P_End.x(), U0_adams, Source_Term
        );
        U0_adams = Result_Adams;
        // Вывод значений решений
        cout << left << "\n\nCurrent time: " << P_End.x();
        cout << "\nERK solution " << setw(18) << "| Adams solution";

        cout << "\nPrey:     " << setw(18) << Result_ERK[0]
            << " | " << setw(18) << Result_Adams[0]
            <<  "\nPredator: " << setw(18) << Result_ERK[1] 
            << " | " << setw(18) << Result_Adams[1];
    }

    return 0;
}
