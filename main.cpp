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
void update_parameters(int month, double& alpha, double& beta, double& gamma, double& delta) {
    if (month >= 3 && month <= 4) {
        // Весна: средний рост зайцев, умеренная активность лис
        alpha = 0.3;
        beta = 0.001;
        gamma = 0.15;
        delta = 0.00008;
    } 
    else if (month >= 5 && month <= 8) {
        // Лето: высокий рост зайцев, высокая активность лис
        alpha = 0.5;
        beta = 0.0012;
        gamma = 0.17;
        delta = 0.0001;
    }
    else if (month >= 9 && month <= 10) {
        // Осень: умеренный спад зайцев, активность лис снижается
        alpha = 0.25;
        beta = 0.0009;
        gamma = 0.16;
        delta = 0.00007;
    }
    else {
        // Зима: сильное снижение численности зайцев, низкая активность лис
        alpha = 0.15;
        beta = 0.0005;
        gamma = 0.2;
        delta = 0.00005;
    }
}
int main() {
    Cauchy_Problem<double> Problem;

    // Количество переменных
    int Variable_Number = 2;
    ofstream dataFile_ERK("data_erk.csv");
    ofstream dataFile_Adams("data_adams.csv");
    // Заголовки для файлов
    dataFile_ERK << "Time,Prey,Predator\n";
    dataFile_Adams << "Time,Prey,Predator\n";
     //Параметры модели
    double alpha = 0.401; // Коэффициент рождаемости жертв
    double beta = 0.001; // Коэффициент успешной охоты
    double gamma = 0.19; // Коэффициент естественной убыли хищников
    double delta = 0.0001; // Коэффициент воспроизводства хищников

    // Правая часть
    function<vector<double>(const vector<double>& U, const Point& P)> Source_Term = [&](const vector<double>& U, const Point& P) {
        vector<double> F(Variable_Number);
        F[0] = (alpha - beta  * U[1]) * U[0];
        F[1] = (-gamma + delta * U[0]) * U[1];

        return F;
    };


    double X0 = 0, Xn = 365*3;
    int n = 365*3;
    double h = (Xn - X0) / n;

    // Начальные условия
    vector<double> U0_erk = { 4212, 552 }; // Начальные условия: x0 = 4212, y0 = 552
    vector<double> U0_adams = { 4212, 552 };

    vector<vector<double>> adams_starter = Problem.RungeKutta_Start(0, h, 3, U0_adams, Source_Term);
    vector<vector<double>> combined_vector_adams;
    combined_vector_adams.push_back(U0_adams);
    combined_vector_adams.insert(combined_vector_adams.end(), adams_starter.begin(), adams_starter.end());


    for (int i = 270; i <= n; i++) {
        int month = ((i * (int)h) / 30) % 12 + 1; 
        update_parameters(month, alpha, beta, gamma, delta);
        cout << left << "\n\nCurrent month: " << month << endl;
        cout << left << "Params:\n" << "Alpha: " << alpha << "\nBeta" << beta << endl;
        cout << left << "Gamma: " << gamma << "\nDelta" << delta;
        Point P_Begin((i - 1) * h, 0, 0), P_End(i * h, 0, 0);
        // Решение задачи Коши по схеме Рунге-Кутты 4
        auto Result_ERK = Problem.ERKs_Scheme_Start(
            Cauchy_Problem<double>::Difference_Scheme_Type::ERK4,
            P_Begin.x(), P_End.x(), U0_erk, Source_Term
        );

        dataFile_ERK << P_End.x() << "," << Result_ERK[0] << "," << Result_ERK[1] << "\n";
        U0_erk = Result_ERK;

        auto Result_Adams = Problem.Adams_Predictor_Corrector_Scheme(
            Cauchy_Problem<double>::Difference_Scheme_Type::Adams,
            P_Begin.x(), P_End.x(), combined_vector_adams, Source_Term
        );
        dataFile_Adams << P_End.x() << "," << Result_Adams[0] << "," << Result_Adams[1] << "\n";
        combined_vector_adams.erase(combined_vector_adams.begin());
        combined_vector_adams.push_back(Result_Adams);


        cout << left << "\n\nCurrent time: " << P_End.x();
        cout << "\nERK solution " << setw(18) << "| Adams solution";

        cout << "\nPrey:     " << setw(18) << Result_ERK[0]
            << " | " << setw(18) << Result_Adams[0]
            << "\nPredator: " << setw(18) << Result_ERK[1]
            << " | " << setw(18) << Result_Adams[1];
    }
    dataFile_ERK.close();
    dataFile_Adams.close(); 
    int result_ERK = system("python3 plot.py data_erk.csv");
    int result_Adams = system("python3 plot.py data_adams.csv");
    return 0;
}
