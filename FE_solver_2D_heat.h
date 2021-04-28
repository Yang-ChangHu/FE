#pragma once
#include<iostream>
#include<string>
#include <Eigen/Dense>
using namespace Eigen;
#include"FE_solver_2D.h"
using namespace std;
class  FE_solver_2D_heat :public FE_solver_2D
{
public:
	//非稳态c(x,y)，f（x,y,t），u（x,y,t)
	FE_solver_2D_heat(double start, double end, double dt, int N1_, int N2_, int mesh_type_, double a_x, double a_y, double b_x, double b_y, int basis_type_trial_, int basis_type_test_);

	virtual void autoRun();

	/*组装A矩阵
	调用:
			Caculate_vertices(n);网格单元节点坐标
			Compute_Gauss(); //计算高斯点权重(local)
			this->Gauss_qual_trial_test_2D;//用高斯积分法计算积分
					Cx(x,y)
					this->FE_basis_local_fun_trial(x,y, alpha, r, s) 计算偏导*/
	virtual void Assemble_matrix_A();

	/*计算第n个单元的trial_test高斯积分(对于M矩阵
	输入：
			r,s,p,q:trial基函数和test基函数分别对x、y的偏导数，这里对各项偏导都等于0
			n：第n个单元
	调用:
		   Compute_Gauss（）
			FE_basis_local_fun_trial：∂ψ(r+s)/∂x(r)∂y(s)
		   FE_basis_local_fun_test：∂ψ(p+q)/∂x(p)∂y(q)  目前test基函数==trial基函数
	输出:
			integral（1*ψtrial*ψtest）   */
	double Gauss_qual_trial_test_M_2D(int alpha, int belta, int n, int r, int s, int p, int q);

	/*计算第n个网格单元的fx_test高斯积分(非稳态)
	输入：
			p,q:test基函数分别对x、y的偏导数 在这里p=q=0
	调用:
		   Compute_Gauss（）
			fx(xi,yi,t)
			FE_basis_local_fun_test：ψtest(xi,yi)
	输出:
			integral（fx（xi,yi,t）*ψ_nβ(xi,yi)）*/
	double Gauss_qual_fxyt_test_2D(int belta, int n, int p, int q, double t);
	void IterateInTime(double theta);

	//组装b（tm）.b(tm+1)
	void Assemble_b(double tm, double tmp1);


	/*处理Dirichlet边界条件
调用：
		g_boundary:计算边界有限元的真实值*/
	virtual void Treat_Boundary_Dirichlet();

	/*处理neumann边界条件
	调用：
			Gauss_qual_neumann_test_2D:计算在第k条边界边上的c(x,y)*p(x,y)*local_basis_function_test(x,y)积分值*/
	virtual void Treat_Boundary_Neumann();

	//处理Robin边界条件
	virtual void Treat_Boundary_Robin();

	void Solution_heat();

	//计算最大误差
	virtual void Compute_Error();

	//c(x,y,t)
	double  Cx(double x, double y);

	//f(x,y,t)
	double fx(double x, double y, double t);

	//Neumann边界的cp（x，y）
	double  cp(double x, double y);

	//非稳态的u（x，y，t）
	double Real_Ux(double x, double y, double t);

	//边界有限元节点的值
	double g_boundary(double x, double y);

	//非稳态参数
	//
	double theta;
	//开始时间
	double start_;

	//结束时间
	double end_;

	//时间步长
	double dt_;

	//质量矩阵
	MatrixXd m_matrix_;

	//初始时间的X0
	MatrixXd X_init;

	//
	MatrixXd bm_vector_;

	MatrixXd bmp1_vector_;

	MatrixXd a_tilde_;
	MatrixXd a_fixed_;
	MatrixXd b_tilde_;

	//计算X0
	void Assemble_X_init();
};
