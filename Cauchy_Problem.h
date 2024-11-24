#pragma once
#ifndef Cauchy_Problem_h
#define Cauchy_Problem_h
#include <vector>
#include <functional>
#include "Point.h"
using namespace std;
// задача Коши
template<class Type>
class Cauchy_Problem {
public:


	// доступные конечно-разностные схемы

	enum class Difference_Scheme_Type {
		ERK1 = 1,
		ERK2,
		ERK4,
		Adams,
	};


	// реализация явных схем Рунге-Кутты
	// FD_Type - тип разностной схемы
	// Time_Begin - начальный момент
	// Time_End - конечный момент
	// Initial_Conditions - начальные условия
	// F - правая часть системы ДУ

	vector<Type> Adams_Predictor_Corrector_Scheme(Difference_Scheme_Type FD_Type,
		double Time_Begin,
		double Time_End,
		const vector<vector<double>>& Initial_Conditions,
		const function<vector<Type>(const vector<Type>& U, const Point& P)>& F);

	vector<Type> ERKs_Scheme_Start(Difference_Scheme_Type FD_Type,
		double Time_Begin,
		double Time_End,
		const vector<Type>& Initial_Conditions,
		const function<vector<Type>(const vector<Type>& U, const Point& P)>& F);

	vector<vector<Type>> RungeKutta_Start(double Time_Begin, double h, int steps,
		const vector<Type>& Initial_Conditions,
		const function<vector<Type>(const vector<double>& U, const Point& P)>& F);


};

#endif 

