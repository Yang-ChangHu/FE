#include "FE_solver_2D_heat.h"
#include<math.h>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
#include<iostream>

//非稳态
FE_solver_2D_heat::FE_solver_2D_heat(double start, double end, double dt, int N1_, int N2_, int mesh_type_, double a_x, double a_y, double b_x, double b_y, int basis_type_trial_, int basis_type_test_)
{
	this->start_ = start;
	this->end_ = end;
	this->dt_ = dt;
	this->N1_ = N1_;
	this->N2_ = N2_;

	this->n_m_ = (N1_ + 1) * (N2_ + 1);
	this->mesh_type = mesh_type_;

	this->a_ = ArrayXd::Zero(2, 1);
	this->a_ << a_x,
		a_y;

	this->b_ = ArrayXd::Zero(2, 1);
	this->b_ << b_x,
		b_y;

	this->basis_type_trial_ = basis_type_trial_;         //201 二维线性trial基函数   202 二维二次trial基函数
	this->basis_type_test_ = basis_type_test_;           //201 二维线性test基函数   202 二维二次test基函数

	this->number_of_gauss_points = 9;    //高斯插值点个数为9
}

void FE_solver_2D_heat::Assemble_matrix_A()
{
	//组装刚度矩阵a_matrix_，循环遍历每一个网格单元，计算aij

	this->a_matrix_ = MatrixXd::Zero(this->nb_test_, this->nb_trial_);
	MatrixXd a_matrix_1 = MatrixXd::Zero(this->nb_test_, this->nb_trial_);
	MatrixXd a_matrix_2 = MatrixXd::Zero(this->nb_test_, this->nb_trial_);

	this->m_matrix_ = MatrixXd::Zero(this->nb_test_, this->nb_trial_);

	for (int n = 0; n < this->n_; n++)
	{
		for (int alpha = 0; alpha < this->number_of_local_basis_trial_; alpha++)
		{
			for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
			{
				a_matrix_1(this->tb_test_(belta, n) - 1, this->tb_trial_(alpha, n) - 1) += this->Gauss_qual_trial_test_2D(alpha, belta, n, 1, 0, 1, 0); // 对x的偏微分
				a_matrix_2(this->tb_test_(belta, n) - 1, this->tb_trial_(alpha, n) - 1) += this->Gauss_qual_trial_test_2D(alpha, belta, n, 0, 1, 0, 1); //对y的偏微分

				this->m_matrix_(this->tb_test_(belta, n) - 1, this->tb_trial_(alpha, n) - 1) += this->Gauss_qual_trial_test_M_2D(alpha, belta, n, 0, 0, 0, 0);	//计算M矩阵
			}
			a_matrix_ = a_matrix_1 + a_matrix_2;
		}
	}
	cout << "\tAssemble_matrix_A()" << endl;
	cout << this->a_matrix_ << endl;
	cout << "\tAssemble_matrix_M()" << endl;
	cout << this->m_matrix_ << endl;
}

void FE_solver_2D_heat::Assemble_X_init()
{
	cout << "*********************************************" << endl;
	cout << "\tAssemble_X_init" << endl;

	this->X_init = MatrixXd::Zero(this->nb_test_, 1);
	double x;
	double y;
	for (int i = 0; i < this->nb_test_; i++)
	{
		x = this->pb_test_(0, i);
		y = this->pb_test_(1, i);

		this->X_init(i, 0) = Real_Ux(x, y, 0);
	}

	cout << this->X_init << endl;
}

void  FE_solver_2D_heat::IterateInTime(double theta)
{
	int numOfTimeStep = (this->end_ - this->start_) / (this->dt_);

	double tm = 0;
	double tmp1 = 0;
	a_tilde_ = (1 / dt_) * m_matrix_ + theta * a_matrix_;
	a_fixed_ = (1 / dt_) * m_matrix_ - (1 - theta) * a_matrix_;

	MatrixXd X_m = this->X_init;

	for (int m = 0; m < numOfTimeStep; m++)
	{
		tm = m * (this->dt_);
		tmp1 = (m + 1) * (this->dt_);
		Assemble_b(tm, tmp1);
		b_tilde_ = theta * bmp1_vector_ + (1 - theta) * bm_vector_ + a_fixed_ * X_m;
		Treat_Boundary_Neumann();
		Treat_Boundary_Robin();
		Treat_Boundary_Dirichlet();
		Solution_heat();
		X_m = this->solution_;
	}
}

void FE_solver_2D_heat::Assemble_b(double tm, double tmp1)
{
	cout << "*********************************************" << endl;
	cout << "\tAssemble_bm_bmp1()" << endl;
	this->bm_vector_ = MatrixXd::Zero(this->nb_test_, 1);
	this->bmp1_vector_ = MatrixXd::Zero(this->nb_test_, 1);

	for (int n = 0; n < this->n_; n++)
	{
		Caculate_vertices(n);
		Compute_Gauss(this->number_of_gauss_points); //计算高斯点权重

		for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
		{
			//integral（c（x,y）*▽ψ_nα*▽ψ_nβ）
			double gauss_quad_2D_trial_test_tm = 0;
			double gauss_quad_2D_trial_test_tmp1 = 0;

			gauss_quad_2D_trial_test_tm = Gauss_qual_fxyt_test_2D(belta, n, 0, 0, tm);
			gauss_quad_2D_trial_test_tmp1 = Gauss_qual_fxyt_test_2D(belta, n, 0, 0, tmp1);

			this->bm_vector_(this->tb_test_(belta, n) - 1, 0) += gauss_quad_2D_trial_test_tm;
			this->bmp1_vector_(this->tb_test_(belta, n) - 1, 0) += gauss_quad_2D_trial_test_tmp1;
		}
	}

	cout << bm_vector_ << endl;
	cout << "*************" << endl;
	cout << bmp1_vector_ << endl;
}

void FE_solver_2D_heat::Solution_heat()
{
	this->solution_ = MatrixXd::Zero(this->nb_test_, 1);

	this->solution_ = this->a_tilde_.inverse() * (this->b_tilde_);
	cout << "solution_:" << endl;
	cout << this->solution_ << endl;
}

void FE_solver_2D_heat::Compute_Error()
{
	MatrixXd real_val = MatrixXd::Zero(this->nb_test_, 1);
	for (int belta = 0; belta < this->nb_test_; belta++)
	{
		real_val(belta, 0) = Real_Ux(this->pb_test_(0, belta), this->pb_test_(1, belta), 1);
	}
	MatrixXd error;
	error = this->solution_ - real_val;
	error = error.cwiseAbs();
	MatrixXd::Index maxRow, maxCol;
	double max = error.maxCoeff(&maxRow, &maxCol);

	cout << "Max abs error:" << endl;
	cout << max << endl;
}

double FE_solver_2D_heat::Cx(double x, double y)
{
	return 2;
}

double  FE_solver_2D_heat::fx(double x, double y, double t)
{
	return -3 * exp(x + y + t);
}

double  FE_solver_2D_heat::cp(double x, double y)
{
	return 0;//？还没写
}

//计算t=0时刻dirichlet边界条件
double  FE_solver_2D_heat::g_boundary(double x, double y)
{
	if (x == 0)
	{
		return  exp(y);
	}
	else if (x == 2)
	{
		return exp(y + 2);
	}
	if (y == 0)
	{
		return  exp(x);
	}
	else if (y == 1)
	{
		return exp(x + 1);
	}
}

double FE_solver_2D_heat::Real_Ux(double x, double y, double t)
{
	return  exp(x + y + t);
}

double FE_solver_2D_heat::Gauss_qual_trial_test_M_2D(int alpha, int belta, int n, int r, int s, int p, int q)
{
	//

	double int_value = 0;

	//1.计算第n个网格单元点的高斯插值点以及权重，插值点个数为number_of_gauss_points
	MatrixXd local_Gauss_weights_point = Compute_Gauss(this->number_of_gauss_points, n);

	//2.用高斯插值法计算积分
	for (int k = 0; k < this->number_of_gauss_points; k++)   //
	{
		double fai_trial_x = 0;
		double fai_test_x = 0;

		fai_trial_x = this->FE_basis_local_fun_trial(local_Gauss_weights_point(1, k), local_Gauss_weights_point(2, k), n, alpha, r, s); //计算∂ψnα/∂x(xi,yi)
		fai_test_x = this->FE_basis_local_fun_test(local_Gauss_weights_point(1, k), local_Gauss_weights_point(2, k), n, belta, p, q);//计算∂ψβ/∂x(xi,yi)

		int_value += local_Gauss_weights_point(0, k) * fai_trial_x * fai_test_x;
	}

	return int_value;
}

double FE_solver_2D_heat::Gauss_qual_fxyt_test_2D(int belta, int n, int p, int q, double t)
{
	double int_value = 0;

	//计算高斯权重
	MatrixXd local_Gauss_weights_point = Compute_Gauss(this->number_of_gauss_points, n);

	for (int k = 0; k < this->number_of_gauss_points; k++)
	{
		double fx_value = 0;

		fx_value = this->fx(local_Gauss_weights_point(1, k), local_Gauss_weights_point(2, k), t); //计算fx（xi,yi）

		double fai_test_x = 0;
		fai_test_x = this->FE_basis_local_fun_test(local_Gauss_weights_point(1, k), local_Gauss_weights_point(2, k), n, belta, p, q);//计算ψ_testβ(xi,yi)

		int_value += local_Gauss_weights_point(0, k) * fx_value * fai_test_x;
	}

	return int_value;
}

void FE_solver_2D_heat::autoRun()
{
	Generate_PT();
	Generate_BoundaryNodes();
	Assemble_matrix_A();
	Assemble_X_init();
	IterateInTime(0.5);
	Compute_Error();
}