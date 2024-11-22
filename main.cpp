#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib> 
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
    ofstream dataFile_ERK("data_erk.csv");
    ofstream dataFile_Adams("data_adams.csv");

    // Заголовки для файлов
    dataFile_ERK << "Time,Prey,Predator\n";
    dataFile_Adams << "Time,Prey,Predator\n";
    // Параметры модели
    //double alpha = 4.01; // Коэффициент рождаемости жертв
    //double beta = 0.01; // Коэффициент успешной охоты
    //double gamma = 10; // Коэффициент естественной убыли хищников
    //double delta = 0.01; // Коэффициент воспроизводства хищников
    double alpha = 0.3; // Коэффициент рождаемости жертв
    double beta = 0.06; // Коэффициент успешной охоты
    double gamma = 0.7; // Коэффициент естественной убыли хищников
    double delta = 0.035; // Коэффициент воспроизводства хищников

    // Правая часть
    function<vector<double>(const vector<double>& U, const Point& P)> Source_Term = [&](const vector<double>& U, const Point& P) {
        vector<double> F(Variable_Number);
        F[0] = alpha * U[0] - beta * U[0] * U[1]; 
        F[1] = delta * U[0] * U[1] - gamma * U[1];
        //vector<double> F(Variable_Number);
        //F[0] = U[1];
        //F[1] = exp(P.x());
        return F;
    };

    // Начальные условия
    //vector<double> U0_erk = { 4212, 552 }; // Начальные условия: x0 = 4212, y0 = 552
    //vector<double> U0_adams = { 4212, 552 };
    vector<double> U0_erk = { 5, 2 }; 
    vector<double> U0_adams = { 5, 2 };
    double X0 = 0, Xn = 365;

    // Мелкость разбиения отрезка и шаг
    int n = 365;
    double h = (Xn - X0) / n;
    vector<vector<double>> adams_starter = Problem.RungeKutta_Start(0, h/10, 3, U0_adams, Source_Term);
    vector<vector<double>> combined_vector;
    combined_vector.push_back(U0_adams);
    combined_vector.insert(combined_vector.end(), adams_starter.begin(), adams_starter.end());
    // Отрезок, на котором отыскивается численное решение

    for (int i = 1; i <= n; i++) {
        Point P_Begin((i - 1) * h, 0, 0), P_End(i * h, 0, 0);
        // Решение задачи Коши по схеме Рунге-Кутты 4
        auto Result_ERK = Problem.ERKs_Scheme_Start(
            Cauchy_Problem<double>::Difference_Scheme_Type::ERK4,
            P_Begin.x(), P_End.x(), U0_erk, Source_Term
        );
        dataFile_ERK << P_End.x() << "," << Result_ERK[0] << "," << Result_ERK[1] << "\n";
        U0_erk = Result_ERK;
        cout << endl;
        for (auto it : combined_vector) {
            for (auto i : it) {
                cout << i << endl;

            }
            cout << endl;
        }

        auto Result_Adams = Problem.Adams_Bashforth_Scheme_Start(
            Cauchy_Problem<double>::Difference_Scheme_Type::AdamsBashforth2,
            P_Begin.x(), P_End.x(), combined_vector, Source_Term
        );
        dataFile_Adams << P_End.x() << "," << Result_Adams[0] << "," << Result_Adams[1] << "\n";
        combined_vector.erase(combined_vector.begin());
        combined_vector.push_back(Result_Adams);
        //U0_adams[0] = Result_Adams[0];
        //U0_adams[1] = Result_Adams[0];
        //U0_adams = Result_ERK;
        // Вывод значений решений
        //cout << left << "\n\nCurrent time: " << P_End.x();
        //cout << "\n ERK Analitical solution" << " | Numerical solution | ";
        //cout << "\n" << setw(19) << U0_erk[0] << " | " << setw(18) << Result_ERK[0];
        //cout << "\n" << setw(19) << U0_erk[1] << " | " << setw(18) << Result_ERK[1];
        //cout << "\n ADAMS Analitical solution" << " | Numerical solution | ";
        //cout << "\n" << setw(19) << combined_vector[0][0] << " | " << setw(18) << Result_Adams[0];
        //cout << "\n" << setw(19) << combined_vector[0][1] << " | " << setw(18) << Result_Adams[1];
        cout << left << "\n\nCurrent time: " << P_End.x();
        cout << "\nERK solution " << setw(18) << "| Adams solution";

        cout << "\nPrey:     " << setw(18) << Result_ERK[0]
            << " | " << setw(18) << Result_Adams[0]
            <<  "\nPredator: " << setw(18) << Result_ERK[1] 
            << " | " << setw(18) << Result_Adams[1];
    }
    dataFile_ERK.close();
    dataFile_Adams.close(); 

    int result_ERK = system("python3 plot.py data_erk.csv");
    int result_Adams = system("python3 plot.py data_adams.csv");
    return 0;
}
