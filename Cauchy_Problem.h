#pragma once
#ifndef Cauchy_Problem_h
#define Cauchy_Problem_h
#include <vector>
#include <functional>
#include "Point.h"
using namespace std;
// ������ ����
template<class Type>
class Cauchy_Problem {
public:


	// ��������� �������-���������� �����

	enum class Difference_Scheme_Type {
		ERK1 = 1,
		ERK2,
		ERK4,
		AdamsBashforth2,
		AdamsBashforth3,
		AdamsBashforth4
	};


	// ���������� ����� ���� �����-�����
	// FD_Type - ��� ���������� �����
	// Time_Begin - ��������� ������
	// Time_End - �������� ������
	// Initial_Conditions - ��������� �������
	// F - ������ ����� ������� ��

	vector<Type> Adams_Bashforth_Scheme_Start(Difference_Scheme_Type FD_Type,
		double Time_Begin,
		double Time_End,
		const vector<Type>& Initial_Conditions,
		const function<vector<Type>(const vector<Type>& U, const Point& P)>& F);

	vector<Type> ERKs_Scheme_Start(Difference_Scheme_Type FD_Type,
		double Time_Begin,
		double Time_End,
		const vector<Type>& Initial_Conditions,
		const function<vector<Type>(const vector<Type>& U, const Point& P)>& F);



};

#endif 
