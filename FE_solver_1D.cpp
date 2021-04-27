#include "FE_solver_1D.h"
#include<math.h>
#include <Eigen/Dense>
using namespace Eigen;     // 改成这样亦可 using Eigen::MatrixXd;
using namespace std;
#include<iostream>

//FE_solver_1D构造函数
FE_solver_1D::FE_solver_1D(int a_, int b_, int n_, int gauss_type_, double ga_, double gb_, int basis_type_trial_, int basis_type_test_, int boundary_, double qbub_)
{
	cout << "FE_solver_1D构造函数运行" << endl;
	this->a_ = ArrayXd::Zero(1);
	this->a_(0, 0) = a_;

	this->b_ = ArrayXd::Zero(1, 1);
	this->b_(0, 0) = b_;

	this->vertices_ = MatrixXd::Zero(2, 1);

	this->n_ = n_;
	this->n_m_ = n_ + 1;
	this->gauss_type_ = gauss_type_;
	this->ga_ = ga_;
	this->gb_ = gb_;
	this->basis_type_trial_ = basis_type_trial_;
	this->basis_type_test_ = basis_type_test_;
	this->boundary_ = boundary_;
	this->nbn_ = 2;//nbn：the number of boundary finite element nodes
	this->qbub_ = qbub_;
}

//计算P、Pb_trial、Pb_test:
void FE_solver_1D::Generate_PT()            //virtual关键字在这就不需要了
{
	/*************************************************
	Function:       // Generate_P
	Description:    // 生成网格矩阵P,训练函数单元矩阵Pb_trial,test函数单元矩阵Pb_test
	Input:          //
	Output:         // this->P_，this->pb_test_,this->pb_trial
	Return:         //
	Others:         // 其它说明
	*************************************************/
	this->p_.setLinSpaced(this->n_m_, this->a_(0, 0), this->b_(0, 0));

	RowVectorXi t1 = VectorXi::LinSpaced(this->n_, 0, this->n_ - 1);
	RowVectorXi t2 = VectorXi::LinSpaced(this->n_, 1, this->n_);
	this->t_ = MatrixXi::Zero(2, this->n_);
	this->t_ << t1, t2;

	//RowVectorXd t2 = t1+1;

	if (this->basis_type_test_ == 101)
	{
		this->number_of_local_basis_test_ = 2;
		this->nb_test_ = this->n_ + 1;
		this->pb_test_ = this->p_;
		this->tb_test_ = this->t_;
	}
	else if (this->basis_type_test_ == 102)
	{
		number_of_local_basis_test_ = 3;
		this->nb_test_ = 2 * this->n_ + 1;
		this->pb_test_.setLinSpaced(2 * this->n_ + 1, this->a_(0, 0), this->b_(0, 0));

		RowVectorXi tb1_test = VectorXi::LinSpaced(this->n_, 0, this->nb_test_ - 3);
		RowVectorXi tb2_test = VectorXi::LinSpaced(this->n_, 2, this->nb_test_ - 1);
		RowVectorXi tb3_test = VectorXi::LinSpaced(this->n_, 1, this->nb_test_ - 2);

		this->tb_test_ = MatrixXi::Zero(3, this->n_);
		this->tb_test_ << tb1_test, tb2_test, tb3_test;
	}

	if (this->basis_type_trial_ == 101)
	{
		this->number_of_local_basis_trial_ = 2;
		this->nb_trial_ = this->n_ + 1;
		this->pb_trial_ = this->p_;
		this->tb_trial_ = this->t_;
	}
	else if (this->basis_type_trial_ == 102)           //一维二次单元
	{
		this->number_of_local_basis_trial_ = 3;
		this->nb_trial_ = 2 * this->n_ + 1;
		this->pb_trial_.setLinSpaced(2 * this->n_ + 1, this->a_(0, 0), b_(0, 0));

		RowVectorXi tb1_trial = VectorXi::LinSpaced(this->n_, 0, this->nb_test_ - 3);
		RowVectorXi tb2_trial = VectorXi::LinSpaced(this->n_, 2, this->nb_test_ - 1);
		RowVectorXi tb3_trial = VectorXi::LinSpaced(this->n_, 1, this->nb_test_ - 2);

		this->tb_trial_ = MatrixXi::Zero(3, this->n_);
		this->tb_trial_ << tb1_trial, tb2_trial, tb3_trial;
	}
}

//设定边界条件
void FE_solver_1D::Generate_BoundaryNodes()
{
	//Dirichlet:1  Neumann :2   Robin:3
	this->boundary_nodes_ = MatrixXi::Zero(3, 2);
	this->boundary_nodes_ << this->boundary_ / 10, this->boundary_ % 10,
		0, this->nb_test_ - 1,
		-1, 1;
	cout << "boundary_nodes_" << endl;
	cout << this->boundary_nodes_ << endl;
}
//组装A矩阵
void FE_solver_1D::Assemble_matrix_A()
{
	/*************************************************
	Function:       // Assemble_matrix_A
	Description:    // 组装A矩阵
	Input:          //
	Output:         // A
	Return:         //
	Others:         // 其它说明
	*************************************************/
	cout << "\tAssemble_matrix_A()" << endl;

	this->a_matrix_ = MatrixXd::Zero(this->nb_test_, this->nb_trial_);

	for (int n = 0; n < this->n_; n++)
	{
		this->vertices_(0, 0) = this->p_(this->t_(0, n));
		this->vertices_(1, 0) = this->p_(this->t_(1, n));
		Compute_Gauss(this->number_of_gauss_points); //计算高斯点权重

		for (int alpha = 0; alpha < this->number_of_local_basis_trial_; alpha++)
		{
			for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
			{
				double gauss_quad_1D_trial_test = 0;
				gauss_quad_1D_trial_test = this->Gauss_qual_trial_test(alpha, belta);
				this->a_matrix_(this->tb_test_(belta, n), this->tb_trial_(alpha, n)) += gauss_quad_1D_trial_test;
			}
		}
	}

	cout << this->a_matrix_ << endl;
}
//组装b向量
void FE_solver_1D::Assemble_b()
{
	/*************************************************
	Function:       // Assemble_matrix_b
	Description:    // 组装b向量
	Input:          //
	Output:         //  b_vertor_
	Return:         //
	Others:         // 其它说明
	*************************************************/
	cout << "\tAssemble_matrix_b()" << endl;

	this->b_vector_ = MatrixXd::Zero(this->nb_test_, 1);

	for (int n = 0; n < this->n_; n++)
	{
		this->vertices_(0, 0) = this->p_(this->t_(0, n));
		this->vertices_(1, 0) = this->p_(this->t_(1, n));
		Compute_Gauss(this->number_of_gauss_points); //计算高斯点权重
		for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
		{
			double gauss_quad_1D_fx_test = 0;
			gauss_quad_1D_fx_test = this->Gauss_qual_fx_test(belta);
			this->b_vector_(this->tb_test_(belta, n), 0) += gauss_quad_1D_fx_test;
		}
	}
	cout << this->b_vector_ << endl;
}

//处理边界条件
void FE_solver_1D::Treat_Boundary_Dirichlet()
{
	//double g[2] = { this->ga_,this->gb_ };
	Vector2d g_boundary(this->ga_, this->gb_);

	double ab[2] = { this->a_(0,0),this->b_(0,0) };
	for (int i = 0; i < 2; i++)
	{
		if (this->boundary_nodes_(0, i) == 1)      //Dirichlet
		{
			for (int belta = 0; belta < this->nb_test_; belta++)
			{
				this->a_matrix_(this->boundary_nodes_(1, i), belta) = 0;
			}
			this->a_matrix_(this->boundary_nodes_(1, i), this->boundary_nodes_(1, i)) = 1;
			//this->b_vector_(this->boundary_nodes_(1, i), 0) = g[i];
			this->b_vector_(this->boundary_nodes_(1, i), 0) = g_boundary(i);
		}
		else if (this->boundary_nodes_(0, i) == 2)           //NeuMann
		{
			this->b_vector_(nb_test_ - 1, 0) += this->boundary_nodes_(2, i) * g_boundary(i) * Cx(ab[i]);
		}
		else if (this->boundary_nodes_(0, i) == 3)             //Robin
		{
			MatrixXd tmp = MatrixXd::Ones(this->nb_test_, this->nb_trial_);
			this->a_matrix_ += (this->qbub_ * tmp);
			this->b_vector_(nb_test_ - 1, 0) += this->boundary_nodes_(2, i) * g_boundary(i) * Cx(ab[i]);
		}
	}

	cout << "边界处理后" << endl;
	cout << "A:" << endl;
	cout << this->a_matrix_ << endl;
	cout << "b:" << endl;
	cout << this->b_vector_ << endl;
}

//处理neumann边界条件
void FE_solver_1D::Treat_Boundary_Neumann()
{
	Vector2d g_boundary(this->ga_, this->gb_);

	double ab[2] = { this->a_(0,0),this->b_(0,0) };
	for (int i = 0; i < 2; i++)
	{
		if (this->boundary_nodes_(0, i) == 2)           //NeuMann
		{
			this->b_vector_(nb_test_ - 1, 0) += this->boundary_nodes_(2, i) * g_boundary(i) * Cx(ab[i]);
		}
	}

	cout << "Neumann边界处理后" << endl;
	cout << "A:" << endl;
	cout << this->a_matrix_ << endl;
	cout << "b:" << endl;
	cout << this->b_vector_ << endl;
}

//处理Robin边界条件
void FE_solver_1D::Treat_Boundary_Robin()
{
	Vector2d g_boundary(this->ga_, this->gb_);

	double ab[2] = { this->a_(0,0),this->b_(0,0) };
	for (int i = 0; i < 2; i++)
	{
		if (this->boundary_nodes_(0, i) == 3)             //Robin
		{
			MatrixXd tmp = MatrixXd::Ones(this->nb_test_, this->nb_trial_);
			this->a_matrix_ += (this->qbub_ * tmp);
			this->b_vector_(nb_test_ - 1, 0) += this->boundary_nodes_(2, i) * g_boundary(i) * Cx(ab[i]);
		}
	}

	cout << "Robin边界处理后" << endl;
	cout << "A:" << endl;
	cout << this->a_matrix_ << endl;
	cout << "b:" << endl;
	cout << this->b_vector_ << endl;
}

//求出Uh
void FE_solver_1D::Solution()
{
	this->solution_ = MatrixXd::Zero(this->nb_test_, 1);

	this->solution_ = this->a_matrix_.inverse() * (this->b_vector_);
	cout << "solution_:" << endl;
	cout << this->solution_ << endl;
}
//计算最大误差
void FE_solver_1D::Compute_Error()
{
	MatrixXd real_val = MatrixXd::Zero(this->nb_test_, 1);
	for (int belta = 0; belta < this->nb_test_; belta++)
	{
		real_val(belta, 0) = Real_Ux(this->pb_test_(0, belta));
	}
	MatrixXd error;
	error = this->solution_ - real_val;
	error = error.cwiseAbs();
	MatrixXd::Index maxRow, maxCol;
	double max = error.maxCoeff(&maxRow, &maxCol);

	cout << "Max abs error:" << endl;
	cout << max << endl;
}
//计算高斯积分的权重和节点
void FE_solver_1D::Compute_Gauss(int n)
{
	this->number_of_gauss_points = n;
	this->gauss_weight_nodes = MatrixXd::Zero(2, this->number_of_gauss_points);
	MatrixXd bias = MatrixXd::Ones(1, this->number_of_gauss_points);
	switch (this->number_of_gauss_points)
	{
	case 4:
	{
		bias = bias * (this->vertices_.mean());
		this->gauss_weight_nodes << 0.3478548451, 0.3478548451, 0.6521451549, 0.6521451549,
			0.8611363116, -0.8611363116, 0.3399810436, -0.3399810436;//高斯权重
		break;
	}
	default:
		break;
	}

	//高斯节点赋值
	cout << this->gauss_weight_nodes << endl;
	gauss_weight_nodes.row(1) = (gauss_weight_nodes.row(1) * (this->vertices_(1, 0) - this->vertices_(0, 0)) / 2.0) + bias;
}
//计算trial_test高斯积分
double FE_solver_1D::Gauss_qual_trial_test(int alpha, int belta)
{
	double int_value = 0;
	for (int k = 0; k < this->number_of_gauss_points; k++)
	{
		double cx = 0;
		cx = this->Cx(gauss_weight_nodes(1, k));

		double fai_trial_x = 0;
		fai_trial_x = this->FE_basis_local_fun_trial(gauss_weight_nodes(1, k), alpha, 1);

		double fai_test_x = 0;//目前暂时用trial基函数代替
		fai_test_x = this->FE_basis_local_fun_trial(gauss_weight_nodes(1, k), belta, 1);

		int_value += this->gauss_weight_nodes(0, k) * cx * fai_trial_x * fai_test_x;
	}
	int_value = int_value * (this->vertices_(1, 0) - this->vertices_(0, 0)) / 2.0;
	return int_value;
}

//c(x)
double  FE_solver_1D::Cx(double x)
{
	return exp(x);
}
double FE_solver_1D::Real_Ux(double x)
{
	return x * cos(x);
}
//f(x)
double FE_solver_1D::fx(double x)
{
	return -(cos(x) - 2 * sin(x) - x * cos(x) - x * sin(x)) * exp(x);
}
//计算fx_test高斯积分
double FE_solver_1D::Gauss_qual_fx_test(int belta)
{
	double int_value = 0;
	double fx = 0;
	double fai_test_x0 = 0;
	for (int k = 0; k < this->number_of_gauss_points; k++)
	{
		fx = this->fx(gauss_weight_nodes(1, k));
		fai_test_x0 = this->FE_basis_local_fun_trial(gauss_weight_nodes(1, k), belta, 0);
		int_value += this->gauss_weight_nodes(0, k) * fx * fai_test_x0;
	}
	int_value = int_value * (this->vertices_(1, 0) - this->vertices_(0, 0)) / 2.0;
	return int_value;
}

//trial基函数
double FE_solver_1D::FE_basis_local_fun_trial(double x, int basis_index, int basis_der_x)
{
	double h = (this->vertices_(1, 0) - this->vertices_(0, 0));
	if (this->basis_type_trial_ == 101)     //101：一维线性
	{
		if (basis_index == 0)
		{
			if (basis_der_x == 0)
			{
				return (this->vertices_(1, 0) - x) / h;
			}
			else if (basis_der_x == 1)
			{
				return -1 / h;
			}
			else if (basis_der_x == 2)
			{
				return 0;
			}
			else
			{
				cout << "你搞错了" << endl;
			}
		}
		else if (basis_index == 1)
		{
			if (basis_der_x == 0)
			{
				return (x - this->vertices_(0, 0)) / h;
			}
			else if (basis_der_x == 1)
			{
				return 1 / h;
			}
			else if (basis_der_x == 2)
			{
				return 0;
			}
			else
			{
				cout << "你搞错了" << endl;
			}
		}
	}
	else if (this->basis_type_trial_ == 102)     //102：一维二次
	{
		if (basis_index == 0)
		{
			if (basis_der_x == 0)
			{
				return 2 * pow(((x - this->vertices_(0, 0)) / h), 2) - 3 * ((x - this->vertices_(0, 0)) / h) + 1;
			}
			else if (basis_der_x == 1)
			{
				return 4 * (x - this->vertices_(0, 0)) / (pow(h, 2)) - 3 / h;
			}
			else if (basis_der_x == 2)
			{
				return 4 / (pow(h, 2));
			}
			else if (basis_der_x > 2)
			{
				return 0;
			}
			else
			{
				cout << "导数阶搞错" << endl;
			}
		}
		else if (basis_index == 1)
		{
			if (basis_der_x == 0)
			{
				return 2 * pow(((x - this->vertices_(0, 0)) / h), 2) - ((x - this->vertices_(0, 0)) / h);
			}
			else if (basis_der_x == 1)
			{
				return 4 * (x - this->vertices_(0, 0)) / (pow(h, 2)) - 1 / h;
			}
			else if (basis_der_x == 2)
			{
				return 4 / (pow(h, 2));
			}
			else if (basis_der_x > 2)
			{
				return 0;
			}
			else
			{
				cout << "导数阶搞错" << endl;
			}
		}
		else if (basis_index == 2)
		{
			if (basis_der_x == 0)
			{
				return -4 * pow(((x - this->vertices_(0, 0)) / h), 2) + 4 * ((x - this->vertices_(0, 0)) / h);
			}
			else if (basis_der_x == 1)
			{
				return -8 * (x - this->vertices_(0, 0)) / (pow(h, 2)) + 4 / h;
			}
			else if (basis_der_x == 2)
			{
				return -8 / (pow(h, 2));
			}
			else if (basis_der_x > 2)
			{
				return 0;
			}
			else
			{
				cout << "导数阶搞错" << endl;
			}
		}
	}
}
//test基函数
double  FE_solver_1D::FE_basis_local_fun_test(double x, int basis_index, int basis_der_x)
{
	//cout << "现在啥也没干" << endl;
	return 0;
}

//打印边界信息
void FE_solver_1D::Print_message_normal()
{
	cout << "\tPrint_message_normal()" << endl;
	cout << "边界为：" << "(" << this->a_ << "," << this->b_ << ")" << endl;
	cout << "边界值为：" << "(" << this->ga_ << "," << this->gb_ << ")" << "\t分成：" << this->n_ << "份" << endl;
	cout << "trial基函数类型为：" << this->basis_type_trial_ << "," << "\ttest基函数类型为：" << this->basis_type_test_ << endl;
	cout << "高斯类型:" << this->gauss_type_ << "\t边界类型为:" << this->boundary_ << endl;
	cout << "P向量" << endl;
	cout << this->p_ << endl;
	cout << "Pb_test向量" << endl;
	cout << this->pb_test_ << endl;
	cout << "Pb_trial向量" << endl;
	cout << this->pb_trial_ << endl;
	cout << "T向量" << endl;
	cout << this->t_ << endl;
	cout << "Tb_test向量" << endl;
	cout << this->tb_test_ << endl;
	cout << "Tb_trial向量" << endl;
	cout << this->tb_trial_ << endl;

	cout << "********************************" << endl;
}

void FE_solver_1D::autoRun()
{
}