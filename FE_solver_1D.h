/*************************************************
Copyright:杨昌鹄
Author:杨昌鹄
Date:2021-01-15
Description:一维有限元求解器 派生类
**************************************************/

#pragma once
#include<iostream>
#include<string>
#include <Eigen/Dense>
using namespace Eigen;     // 改成这样亦可 using Eigen::MatrixXd; 
#include<iostream>
using namespace std;
//using namespace cv;
#include"FE_solver.h"

//一维求解器类
class FE_solver_1D:public FE_solver
{
public:

	//FE_solver_1D构造函数
	FE_solver_1D(int a_,int b_,int n_,int gauss_type_,double ga_,double gb_,int basis_type_trial_,int basis_type_test_,int boundary_,double qbub_=0);

	//计算P、Pb_trial、Pb_test:
	virtual void Generate_PT();            //子类重写父类的虚函数或者纯虚函数，virtual关键字可删除也可不删除

	virtual void Generate_PT(int mesh_type);  //空实现，在一维中无意义

	//设定边界条件
	virtual void Generate_BoundaryNodes();

	//组装A矩阵
	virtual void Assemble_matrix_A(bool T=false);

	//组装b向量
	virtual void Assemble_b(bool T=false);

	//处理边界条件
	virtual void Treat_Boundary();

	//求出Uh
	virtual void Solution();

	//计算最大误差
	virtual void Compute_Error();

	//trial基函数
	 double FE_basis_local_fun_trial(double x, int basis_index, int basis_der_x);

	//virtual double FE_basis_local_fun_trial(double x, double y, int basis_index, int basis_der_x, int basis_der_y) ;

	//test基函数
	 double FE_basis_local_fun_test(double x, int basis_index, int basis_der_x);

	//virtual double FE_basis_local_fun_test(double x, double y, int basis_index, int basis_der_x, int basis_der_y) ;

	//计算高斯积分的权重和节点
	virtual void Compute_Gauss(int n);

	//显示基本信息
	virtual void Print_message_normal();



	//c(x)
	double  Cx(double x);

	//f(x)
	double fx(double x);

	//计算trial_test高斯积分
	double Gauss_qual_trial_test( int alpha, int belta) ;

	//计算fx_test高斯积分
	double Gauss_qual_fx_test( int belta);

	//真实u（x）值
	double Real_Ux(double x);

	//网格数组p_;
	RowVectorXd p_;

	//单元数组pb_test
	RowVectorXd pb_test_;

	//trial单元数组 pb_trial
	RowVectorXd pb_trial_;



};
