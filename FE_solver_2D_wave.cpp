//波方程，其实跟热方程大部分都一样，只是我懒得写虚函数了



#include "FE_solver_2D_wave.h"
#include<math.h>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
#include<iostream>



FE_solver_2D_wave::FE_solver_2D_wave(double start, double end, double dt, double u00, int N1_, int N2_, int mesh_type_, double a_x, double a_y, double b_x, double b_y, int basis_type_trial_, int basis_type_test_)
{
	this->start_ = start;
	this->end_ = end;
	this->dt_ = dt;
	this->u00_ = u00;
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

void FE_solver_2D_wave::Assemble_matrix_A()
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

void FE_solver_2D_wave::Assemble_X_init_X_second()
{
	cout << "*********************************************" << endl;
	cout << "\tAssemble_X_init_U00" << endl;

	this->X_init = MatrixXd::Zero(this->nb_test_, 1);
	double x;
	double y;
	for (int i = 0; i < this->nb_test_; i++)
	{
		x = this->pb_test_(0, i);
		y = this->pb_test_(1, i);

		this->X_init(i, 0) = Real_Ux(x, y, 0);
	}

	this->X_second = X_init.array() + dt_ * u00_;

	cout << "X0:" << endl;
	cout << this->X_init << endl;
	cout << "X1:" << endl;
	cout << this->X_second << endl;
}

void  FE_solver_2D_wave::IterateInTime(double theta)
{
	int numOfTimeStep = (this->end_ - this->start_) / (this->dt_);

	double tm = 0;
	a_tilde_ = (1 / pow(dt_, 2)) * m_matrix_ + 0.25 * a_matrix_;		//p79/83
	a_fixed_ = (2 / pow(dt_, 2)) * m_matrix_ - 0.5 * a_matrix_;


	MatrixXd X_m_1 = this->X_init;
	MatrixXd X_m = this->X_second;

	for (int m = 1; m < numOfTimeStep; m++)
	{
		tm = m * (this->dt_);
		Assemble_b(tm);
		b_tilde_ = this->bm_vector_+a_fixed_*X_m+a_tilde_*X_m_1;
		Treat_Boundary_Neumann();
		Treat_Boundary_Robin();
		Treat_Boundary_Dirichlet();
		Solution_heat();
		X_m_1 = X_m;
		X_m = this->solution_;
	}
}

void FE_solver_2D_wave::Assemble_b(double tm)
{
	cout << "*********************************************" << endl;
	cout << "\tAssemble_bm()" << endl;
	this->bm_vector_ = MatrixXd::Zero(this->nb_test_, 1);

	for (int n = 0; n < this->n_; n++)
	{
		Caculate_vertices(n);
		Compute_Gauss(this->number_of_gauss_points); //计算高斯点权重

		for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
		{
			//integral（c（x,y）*▽ψ_nα*▽ψ_nβ）
			double gauss_quad_2D_trial_test_tm = 0;
			gauss_quad_2D_trial_test_tm = Gauss_qual_fxyt_test_2D(belta, n, 0, 0, tm);
			this->bm_vector_(this->tb_test_(belta, n) - 1, 0) += gauss_quad_2D_trial_test_tm;
		}
	}

	cout << bm_vector_ << endl;
	cout << "*************" << endl;
}

void FE_solver_2D_wave::Solution_heat()
{
	this->solution_ = MatrixXd::Zero(this->nb_test_, 1);

	this->solution_ = this->a_tilde_.inverse() * (this->b_tilde_);
	cout << "solution_:" << endl;
	cout << this->solution_ << endl;
}

void FE_solver_2D_wave::Compute_Error()
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

double FE_solver_2D_wave::Cx(double x, double y)
{
	return 2;
}

double  FE_solver_2D_wave::fx(double x, double y, double t)
{
	return -3 * exp(x + y + t);
}

double  FE_solver_2D_wave::cp(double x, double y)
{
	return 0;//？还没写
}

//计算t=0时刻dirichlet边界条件
double  FE_solver_2D_wave::g_boundary(double x, double y)
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

double FE_solver_2D_wave::Real_Ux(double x, double y, double t)
{
	return  exp(x + y + t);
}

double FE_solver_2D_wave::Gauss_qual_trial_test_M_2D(int alpha, int belta, int n, int r, int s, int p, int q)
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

double FE_solver_2D_wave::Gauss_qual_fxyt_test_2D(int belta, int n, int p, int q, double t)
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

//在边界有限元节点上有u=g,所以dirichlet边界用的是有限元节点
void FE_solver_2D_wave::Treat_Boundary_Dirichlet()
{
	for (int k = 0; k < this->nbn_; k++)
	{
		if (this->boundary_nodes_(0, k) == 1)
		{
			int i = this->boundary_nodes_(1, k);
			this->a_tilde_.row(i - 1).fill(0);             //把某一行赋值成0
			this->a_tilde_(i - 1, i - 1) = 1;
			this->b_tilde_(i - 1, 0) = g_boundary(this->pb_test_(0, i - 1), this->pb_test_(1, i - 1));   //这里应该是pb_test_
		}
	}
}

//neumann边界处理
void FE_solver_2D_wave::Treat_Boundary_Neumann()
{
	MatrixXd v = MatrixXd::Zero(this->nb_test_, 1);  //neumann

	for (int k = 0; k < this->nbn_; k++)  //第k条边（从0开始）
	{
		if (this->boundary_edges_(0, k) == 2)   //如果第k条边界边是neumann边
		{
			int nk = this->boundary_edges_(1, k) - 1;  //nk 第k条边界边对应的单元  第nk个单元（从0开始）
			for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
			{
				double r = this->Gauss_qual_neumann_test_2D(belta, nk, k, 0, 0);     //对于P96/138 中的r
				v(this->tb_test_(belta, nk) - 1, 0) = v(this->tb_test_(belta, nk) - 1, 0) + r;
			}
		}
	}
	this->b_tilde_ = this->b_tilde_.array() + v.array();
}

void FE_solver_2D_wave::Treat_Boundary_Robin()
{
	MatrixXd W = MatrixXd::Zero(this->nb_test_, 1);  //neumann
	MatrixXd R = MatrixXd::Zero(this->nb_test_, this->nb_trial_);

	for (int k = 0; k < this->nbn_; k++)
	{
		if (this->boundary_edges_(0, k) == 3)
		{
			int nk = this->boundary_edges_(1, k) - 1;  //nk  第nk个单元（从0开始）
			for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
			{
				double r1 = this->Gauss_qual_Robin_cq_2D(belta, nk, k, 0, 0);
				W(this->tb_test_(belta, nk) - 1, 0) = W(this->tb_test_(belta, nk) - 1, 0) + r1;
			}

			for (int alpha = 0; alpha < this->number_of_local_basis_trial_; alpha++)
			{
				for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
				{
					double r2 = this->Gauss_qual_Robin_cr_2D(alpha, belta, nk, k, 0, 0, 0, 0);
					R(this->tb_test_(belta, nk) - 1, this->tb_test_(alpha, nk) - 1) = R(this->tb_test_(belta, nk) - 1, this->tb_test_(alpha, nk) - 1) + r2;
				}
			}
		}
	}
	this->a_tilde_ = this->a_tilde_ + R;
	this->b_tilde_ = this->b_tilde_ + W;

	//cout << "*******************" << endl;
	//cout << "边界处理后：a_tilde_ 矩阵" << endl;
	//cout << this->a_tilde_ << endl;
	//cout << "边界处理后：b_tilde_ 矩阵" << endl;
	//cout << this->b_tilde_ << endl;
}

void FE_solver_2D_wave::autoRun()
{
	Generate_PT();
	Generate_BoundaryNodes();
	Assemble_matrix_A();
	Assemble_X_init_X_second();
	IterateInTime(0.5);
	Compute_Error();
}