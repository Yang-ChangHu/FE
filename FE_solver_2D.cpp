#include "FE_solver_2D.h"
#include<math.h>
#include <Eigen/Dense>
using namespace Eigen;    
using namespace std;
#include<iostream>




FE_solver_2D::FE_solver_2D(int N1_, int N2_, int gauss_type_, double a_x, double a_y, double b_x, double b_y, int basis_type_trial_, int basis_type_test_)
{
	this->N1_ = N1_;
	this->N2_ = N2_;

	this->gauss_type_ = gauss_type_;
	
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

void FE_solver_2D::Generate_PT(int mesh_type)  //mesh_type:网格类型 3 三角形网格  4 四边形
{
	switch (mesh_type)
	{
	case 3:
	{
		this->n_ = 2 * this->N1_ * this->N2_;
		this->n_m_ = (this->N1_ + 1) * (this->N2_ + 1);
		this->vertices_ = MatrixXd::Zero(2, 3);
		break;
	}
	case 4:
	{
		this->n_ =  N1_ * N2_;
		this->n_m_ = (N1_ + 1) * (N2_ + 1);
		this->vertices_ = MatrixXd::Zero(2, 4);
		break;
	}
	}

	//标准区间P矩阵
	this->p_ = MatrixXd::Zero(2, n_m_);
	RowVectorXd p2_tmp = VectorXd::LinSpaced(this->N2_ + 1, 0, 1);
	for (int j = 0; j < (this->N1_ + 1); j++)
	{
		this->p_.block(0,j*((this->N2_)+1),1,this->N2_+1)= (this->p_.block(0, j * ((this->N2_) + 1), 1, this->N2_ + 1)).array()+j/((float)(this->N1_));
		this->p_.block(1, j * ((this->N2_) + 1), 1, this->N2_ + 1)<< p2_tmp;

	}

	//变换到实际区间
	this->p_.row(0) = this->p_.row(0).array() * (this->b_(0,0) -this->a_(0,0))+this->a_(0,0);
	this->p_.row(1) = this->p_.row(1).array() * (this->b_(1, 0) - this->a_(1, 0)) + this->a_(1, 0);



	//构建T矩阵
	switch (mesh_type)//mesh_type:网格类型 3 三角形网格  4 四边形
	{
	case 3:          //三角形网格
	{
		this->t_ = MatrixXi::Zero(3, this->n_);

		//初始化前两个单元
		this->t_(0, 0) = 1;
		this->t_(1, 0) = this->N2_+2;
		this->t_(2, 0) = 2;

		this->t_(0, 1) = 2;
		this->t_(1, 1) = this->N2_ + 2;
		this->t_(2, 1) = this->N2_ + 3;

		//int i = 2;
		for (int i=2;i < this->n_;i++)
		{
			if(i/(2*this->N2_)==0)
			for (int j = 0; j < 3; j++)
			{
				this->t_(j, i) = this->t_(j, i - 2) + 1;
			}
			else
			{
				for (int j = 0; j < 3; j++)
				{
					this->t_(j, i) = this->t_(j, i - 2*(this->N2_)) + 1+this->N2_;
				}
			}
		}
	}
	case 4:          //四边形风格
	{

	}

	}
	
	
	if (this->basis_type_trial_ == 201)   //二维线性奇函数
	{
		this->pb_trial_ = this->p_;
		this->tb_trial_ = this->t_;

		this->number_of_local_basis_trial_ = 3;
		this->nb_trial_ = this->n_m_;
	}
	else if (this->basis_type_trial_ == 202)   //二维二次奇函数
	{
		//pb_trial还没写
		//tb_trail 也还没写
		this->number_of_local_basis_trial_ = 6;
		this->nb_trial_ = 2 * (this->n_m_);
	}



	if (this->basis_type_test_ == 201)
	{
		this->pb_test_ = this->p_;
		this->tb_test_ = this->t_;

		this->number_of_local_basis_test_ = 3;
		this->nb_test_ = this->n_m_;
	}
	else if (this->basis_type_test_ == 202)   //二维二次奇函数
	{

		//pb_test_还没写
		//tb_test_ 也还没写
		this->number_of_local_basis_test_ = 6;
		this->nb_test_ = 2 * (this->n_m_);
	}


	//输出信息
	cout << "\tP矩阵为" << endl;
	cout << this->p_ << endl;
	cout << "\tpb_trial_矩阵为" << endl;
	cout << this->pb_trial_ << endl;
	cout << "\tpb_test_矩阵为" << endl;
	cout << this->pb_test_ << endl;

	cout << "\t*****************" << endl;
	cout << "\tT矩阵为" << endl;
	cout << this->t_<< endl;
	cout << "\ttb_test_矩阵为" << endl;
	cout << this->tb_test_ << endl;
	cout << "\ttb_trial_矩阵为" << endl;
	cout << this->tb_trial_<<endl;
}
void  FE_solver_2D::Generate_PT()
{
	//空实现，二维不用这个
}

void FE_solver_2D::Generate_BoundaryNodes()
{
	//生产边界边矩阵
	Generate_boundary_edge();

	//生成边界点矩阵
	Generate_boundary_nodes();
}


//边界边矩阵 （有限元概念）
void FE_solver_2D::Generate_boundary_edge()
{
	//nbe:边界边的条数，nbe=2*N1*N2,网格概念
	//我们处理狄利克雷等边界条件，是在有限元节点上已知
	this->nbe_ = 2 * (this->N1_ + this->N2_);
	this->boundary_edges_ = MatrixXi::Zero(4, this->nbe_);

	//把第一行都设为1，即所有边界边都是Dirichelet边界
	this->boundary_edges_.row(0) = MatrixXi::Constant(1, this->nbe_, 1);
	this->boundary_edges_(1, 0) = 1;  //第一条边界边永远是1开始

	this->boundary_edges_(2, 0) = 1;
	this->boundary_edges_(3, 0) = 1 + this->N2_ + 1;

	for (int i = 1; i < this->nbe_; i++)
	{
		this->boundary_edges_(2, i) = this->boundary_edges_(3, i - 1);//第三行         边界边的第一个节点编号

		if (i < this->N1_)      //下面的边界边
		{
			this->boundary_edges_(1, i) = this->boundary_edges_(1, i - 1) + 2 * this->N2_;    //第二行 边界边的编号 
			this->boundary_edges_(3, i) = this->boundary_edges_(2, i) + this->N2_ + 1;   //第四行   边界边的第二个节点编号

		}
		else if (i == this->N1_)   //右边边界点第一个
		{
			this->boundary_edges_(1, i) = this->boundary_edges_(1, i - 1) + 1;
			this->boundary_edges_(3, i) = this->boundary_edges_(2, i) + 1;   //第四行   边界边的第二个节点编号

		}
		else if (i > this->N1_ && i < (this->N1_ + this->N2_))  //右边边界边
		{
			boundary_edges_(1, i) = boundary_edges_(1, i - 1) + 2;
			boundary_edges_(3, i) = boundary_edges_(2, i) + 1;
		}
		else if (i == (this->N1_ + this->N2_)) //上边边界点第一条
		{
			boundary_edges_(1, i) = boundary_edges_(1, i - 1);
			boundary_edges_(3, i) = boundary_edges_(2, i) - (this->N2_ + 1);
		}
		else if (i > (this->N1_ + this->N2_) && i < (2 * this->N1_ + this->N2_))  //上边边界
		{
			boundary_edges_(1, i) = boundary_edges_(1, i - 1) - 2 * this->N2_;
			boundary_edges_(3, i) = boundary_edges_(2, i) - (this->N2_ + 1);
		}
		else if (i == (2 * this->N1_ + this->N2_))                //右边边界最上的那条边
		{
			boundary_edges_(1, i) = boundary_edges_(1, i - 1) - 1;
			boundary_edges_(3, i) = boundary_edges_(2, i) - 1;
		}
		else
		{
			boundary_edges_(1, i) = boundary_edges_(1, i - 1) - 2;
			boundary_edges_(3, i) = boundary_edges_(2, i) - 1;
		}
	}


	cout << "boundary_edges_ matrix" << endl;
	cout << this->boundary_edges_ << endl;


}

//边界点矩阵（有限元概念）
void FE_solver_2D::Generate_boundary_nodes()
{
	//nbn_:边界节点数 nbn=2*(n1_+n2_)  有限元概念
	this->nbn_ = 2 * (this->N1_ + this->N2_);
	this->boundary_nodes_ = MatrixXi::Zero(2, this->nbn_);
	this->boundary_nodes_.row(0) = MatrixXi::Constant(1, this->nbn_, 1);  //Dirichlet边界条件
	this->boundary_nodes_.row(1) = boundary_edges_.row(2);

	cout << "\t boundary nodes matrix:" << endl;
	cout << this->boundary_nodes_ << endl;
}

void FE_solver_2D::Assemble_matrix_A(bool T=false)
{
	cout << "\tAssemble_matrix_A()" << endl;


	this->a_matrix_ = MatrixXd::Zero(this->nb_test_, this->nb_trial_);


	for (int n = 0; n < this->n_; n++)
	{
		Caculate_vertices(n);
		Compute_Gauss(this->number_of_gauss_points); //计算高斯点权重

		for (int alpha = 0; alpha < this->number_of_local_basis_trial_; alpha++)
		{
			for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
			{
				//integral（c（x,y）*▽ψ_nα*▽ψ_nβ）
				double gauss_quad_2D_trial_test = 0;

				//gauss_quad_2D_trial_test的第一项
				double gauss_quad_2D_trial_test_1 = 0;

				//gauss_quad_2D_trial_test的第二项
				double gauss_quad_2D_trial_test_2 = 0;
				
				gauss_quad_2D_trial_test_1 = this->Gauss_qual_trial_test_2D(alpha, belta,1,0,1,0,T);  //计算ppt36/138第一项

				gauss_quad_2D_trial_test_2 = this->Gauss_qual_trial_test_2D(alpha, belta, 0,1,0,1,T);  //计算ppt36/138第二项

				gauss_quad_2D_trial_test = gauss_quad_2D_trial_test_1 + gauss_quad_2D_trial_test_2;

				this->a_matrix_(this->tb_test_(belta, n)-1, this->tb_trial_(alpha, n)-1) += gauss_quad_2D_trial_test;
				
			}
		}
	}



	cout << this->a_matrix_ << endl;
}




void FE_solver_2D::Assemble_b(bool T=false)
{

	cout << "*********************************************" << endl;
	cout << "\tAssemble_b_vector()" << endl;

	this->b_vector_ = MatrixXd::Zero(this->nb_test_, 1);

	for (int n = 0; n < this->n_; n++)
	{
		Caculate_vertices(n);
		Compute_Gauss(this->number_of_gauss_points); //计算高斯点权重

		for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
		{
			//integral（c（x,y）*▽ψ_nα*▽ψ_nβ）
			double gauss_quad_2D_trial_test = 0;

			gauss_quad_2D_trial_test = Gauss_qual_fx_test_2D( belta, 0, 0, 0, 0,T);  //计算ppt36/138第一项
			this->b_vector_(this->tb_test_(belta, n) - 1,0) += gauss_quad_2D_trial_test;

		}

	}



	cout << this->b_vector_ << endl;
}

void FE_solver_2D::Treat_Boundary()
{
	
	for (int k = 0; k < this->nbn_; k++)
	{
		if (this->boundary_nodes_(0, k) == 1)  //dirichelet
		{
			int i = this->boundary_nodes_(1, k);
			this->a_matrix_.row(i - 1).fill(0);             //把某一行赋值成0
			this->a_matrix_(i - 1, i - 1) = 1;
			this->b_vector_(i - 1, 0) = g_boundary(this->pb_trial_(0, i-1), this->pb_trial_(1, i-1));
		}
		else if (this->boundary_nodes_(0, k) == 2)   //Neumann
		{
			MatrixXd v = MatrixXd::Zero(this->nb_test_, 1);
			for (int k = 0; k < this->nbe_; k++)
			{
				if (this->boundary_edges_(0, k) == 2)
				{
					int nk = this->boundary_edges_(1, k);
					for (int belta = 0; belta < this->number_of_local_basis_test_; belta++)
					{
						double int_value = this->Gauss_qual_neumann_test_2D(belta, 0, 0, 0, 0);
						v(this->tb_test_(belta, nk), 0) = v(this->tb_test_(belta,nk), 0) + int_value;
					}
				}
			}
			this->b_vector_ = this->b_vector_.array() + v.array();
		}
		else if (this->boundary_nodes_(0, k) == 3)
		{

		}
	}
	cout << "边界处理后：A矩阵" << endl;
	cout << this->a_matrix_ << endl;

	cout << "*******************" << endl;
	cout << "边界处理后：b矩阵" << endl;
	cout << this->b_vector_ << endl;

}


void FE_solver_2D::Compute_Error()
{
	MatrixXd real_val = MatrixXd::Zero(this->nb_test_, 1);
	for (int belta = 0; belta < this->nb_test_; belta++)
	{
		real_val(belta, 0) = Real_Ux(this->pb_test_(0, belta), this->pb_test_(1, belta));
	}
	MatrixXd error;
	error = this->solution_ - real_val;
	error = error.cwiseAbs();
	MatrixXd::Index maxRow, maxCol;
	double max = error.maxCoeff(&maxRow, &maxCol);


	cout << "Max abs error:" << endl;
	cout << max << endl;
}

void FE_solver_2D::Compute_Gauss(int n)
{
	MatrixXd reference_gauss_weight_nodes(3,3);	//参考高斯权重与节点
	switch (n)
	{
	case 3:
	{
		//标准高斯节点与权重
		this->gauss_weight_nodes = MatrixXd::Zero(3, 3);

		//Matrix3d reference_gauss_weight_nodes(3, 3);   //参考高斯权重与节点
		reference_gauss_weight_nodes << 1 / 6.0, 1 / 6.0, 1 / 6.0,
			1 / 2.0, 1 / 2.0, 0,
			0, 1 / 2.0, 1 / 2.0;
		break;
	}
	case 9:
	{
		reference_gauss_weight_nodes.resize(3, 9);
		this->gauss_weight_nodes = MatrixXd::Zero(3, 9);
		reference_gauss_weight_nodes << 8/ 81.0, 12.5 / 324.0 * (1 - sqrt(3 / 5.0)) , 12.5 / 324.0 * (1 - sqrt(3 / 5.0)) , 12.5 / 324.0 * (1 + sqrt(3 / 5.0)), 12.5 / 324.0 * (1 + sqrt(3 / 5.0)), \
			5 / 81.0, 5 / 81.0, 5/ 81.0 * (1 - sqrt(3 / 5.0)) , 5/ 81.0 * (1 + sqrt(3 / 5.0)),
			0.5, (1 + sqrt(3 / 5.0)) / 2.0, (1 + sqrt(3 / 5.0)) / 2.0, (1 - sqrt(3 / 5.0)) / 2.0, (1 - sqrt(3 / 5.0)) / 2.0, 0.5, 0.5, (1 + sqrt(3 / 5.0)) / 2.0, (1 - sqrt(3 / 5.0)) / 2.0,
			0.25, 0.1, (1 - sqrt(3 / 5.0))* (1 - sqrt(3 / 5.0)) / 4.0, (1 + sqrt(3 / 5.0))* (1 + sqrt(3 / 5.0)) / 4.0, 0.1,  (1 + sqrt(3 / 5.0)) / 4.0,  (1 - sqrt(3 / 5.0)) / 4.0, (1 - sqrt(3 / 5.0)) / 4.0, (1 + sqrt(3 / 5.0)) / 4.0;
	
	}
	}

		//cout << "reference gauss nodes" << endl;
		//cout << reference_gauss_weight_nodes << endl;



		//变换后的高斯节点与权重
		double x1 = this->vertices_(0, 0);
		double y1 = this->vertices_(1, 0);
		double x2 = this->vertices_(0, 1);
		double y2 = this->vertices_(1, 1);
		double x3 = this->vertices_(0, 2);
		double y3 = this->vertices_(1, 2);

		double Jacobi = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
		//cout << Jacobi << endl;
		this->gauss_weight_nodes.row(0) = reference_gauss_weight_nodes.row(0).array() * Jacobi;
		this->gauss_weight_nodes.row(1) = x1 + (x2 - x1) * (reference_gauss_weight_nodes.row(1).array()) + (x3 - x1) * (reference_gauss_weight_nodes.row(2).array());
		this->gauss_weight_nodes.row(2) = y1 + (y2 - y1) * (reference_gauss_weight_nodes.row(1).array()) + (y3 - y1) * (reference_gauss_weight_nodes.row(2).array());

		//cout << "local gauss nodes" << endl;
		//cout << this->gauss_weight_nodes << endl;
	
}

double FE_solver_2D::FE_basis_local_fun_test(double x, double y,int basis_index, int basis_der_x,int basis_der_y)
{
	return 0;
}


double FE_solver_2D::FE_basis_local_fun_trial(double x, double y, int basis_index, int basis_der_x, int basis_der_y)
{
	//
	double jacobi;
	double local_result;

	//Caculate_vertices(n);

	/*double xn1 = this->gauss_weight_nodes(1, 0);
	double xn2 = this->gauss_weight_nodes(1, 1);
	double xn3 = this->gauss_weight_nodes(1, 2);
	double yn1 = this->gauss_weight_nodes(2, 0);
	double yn2 = this->gauss_weight_nodes(2, 1);
	double yn3 = this->gauss_weight_nodes(2, 2);*/

	double xn1 = this->vertices_(0, 0);
	double xn2 = this->vertices_(0, 1);
	double xn3 = this->vertices_(0, 2);
	double yn1 = this->vertices_(1, 0);
	double yn2 = this->vertices_(1, 1);
	double yn3 = this->vertices_(1, 2);


	jacobi = (xn2 - xn1) * (yn3 - yn1) - (xn3 - xn1) * (yn2 - yn1);

	double xh = ((yn3 - yn1) * (x - xn1) - (xn3 - xn1) * (y - yn1)) / jacobi;
	double yh = (-(yn2 - yn1) * (x - xn1) + (xn2 - xn1) * (y - yn1)) / jacobi;

	//cout << "\txh:" << xh << "\tyh:" << yh << endl;

	if (basis_der_x == 0 && basis_der_y == 0)
	{
		local_result = reference_basis_2D(xh, yh, basis_index, basis_der_x, basis_der_y);
	}
	else if (basis_der_x == 1 && basis_der_y == 0)
	{
		local_result = ((yn3 - yn1) / jacobi) * reference_basis_2D(xh, yh, basis_index, 1, 0) + \
			((yn1 - yn2) / jacobi) * reference_basis_2D(xh, yh, basis_index, 0, 1);
	}
	else if (basis_der_x == 0 && basis_der_y == 1)
	{
		local_result = ((xn1 - xn3) * reference_basis_2D(xh, yh, basis_index, 1, 0) + reference_basis_2D(xh, yh, basis_index, 0, 1) * (xn2 - xn1)) / jacobi;
	}
	else if (basis_der_x == 1 && basis_der_y == 1)
	{
		local_result = \
			(reference_basis_2D(xh, yh, basis_index, 2, 0) * (xn1 - xn3) * (yn3 - yn1) + \
				reference_basis_2D(xh, yh, basis_index, 1, 1) * (xn1 - xn3) * (yn1 - yn2) + \
				reference_basis_2D(xh, yh, basis_index, 1, 1) * (xn2 - xn1) * (yn3 - yn1) + \
				reference_basis_2D(xh, yh, basis_index, 0, 2) * (xn2 - xn1) * (yn1 - yn2)\
				) / (pow(jacobi, 2));
	}

	else if (basis_der_x == 2 && basis_der_y == 0)
	{
		local_result = \
			(reference_basis_2D(xh, yh, basis_index, 2, 0) * (pow((yn3 - yn1), 2)) + \
				3 * reference_basis_2D(xh, yh, basis_index, 1, 1) * (yn1 - yn2) * (yn3 - yn1) + \
				reference_basis_2D(xh, yh, basis_index, 0, 2) * (pow((yn1 - yn2), 2)) \
				) / (pow(jacobi, 2));

	}
	else if (basis_der_x == 0 && basis_der_y == 2)
	{
		local_result = \
			(reference_basis_2D(xh, yh, basis_index, 2, 0) * (pow((xn3 - xn1), 2)) + \
				3 * reference_basis_2D(xh, yh, basis_index, 1, 1) * (xn1 - xn3) * (xn2 - xn1) + \
				reference_basis_2D(xh, yh, basis_index, 0, 2) * (pow((xn1 - xn2), 2)) \
				) / (pow(jacobi, 2));
	}

	return local_result;

	//

 }



void FE_solver_2D::Print_message_normal()
{

 }


double FE_solver_2D::reference_basis_2D(double xh, double yh ,int basis_index, int basis_der_x,int basis_der_y)  //参考局部奇函数
{

	if (this->basis_type_test_ == 201)        //二维线性函数  ppt  32/103
	{
		if(basis_index==0)
		{ 
			if (basis_der_x == 0 && basis_der_y == 0)
			{
				return -xh - yh + 1;
			}
			else if (basis_der_x == 1 && basis_der_y == 0)
			{
				return -1;
			}
			else if (basis_der_x == 0 && basis_der_y == 1)
			{
				return -1;
			}
			else
			{
				return 0;
			}
		}
		else if(basis_index==1)
		{
			if (basis_der_x == 0 && basis_der_y == 0)
			{
				return xh;
			}
			else if (basis_der_x == 1 && basis_der_y == 0)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
		else if (basis_index == 2)
		{
			if (basis_der_x == 0 && basis_der_y == 0)
			{
				return yh;
			}
			else if (basis_der_x == 0 && basis_der_y == 1)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}

	}
	else if (this->basis_type_test_ == 202)      //二维二次函数   ppt53/103
	{

		if (basis_index == 0)
		{
			if (basis_der_x == 0 && basis_der_y == 0)
			{
				return 2 * pow(xh , 2) + 2 * pow(yh , 2)+4*xh*yh-3*yh-3*xh+1;
			}
			else if (basis_der_x == 1 && basis_der_y == 0)
			{
				return 4 * (xh + yh) - 3;
			}
			else if (basis_der_x == 0 && basis_der_y == 1)
			{
				return 4 * (xh + yh) - 3;
			}
			else if (basis_der_x == 1 && basis_der_y == 1)
			{
				return 4;
			}
			else if (basis_der_x == 2 && basis_der_y == 0)
			{
				return 4;
			}
			else if (basis_der_x == 0 && basis_der_y == 2)
			{
				return 4;
			}
			else
			{
				return 0;
			}

		}
		else if (basis_index == 1)
		{
			if (basis_der_x == 0 && basis_der_y == 0)
			{
				return 2 * pow(xh, 2) -  xh ;
			}
			else if (basis_der_x == 1 && basis_der_y == 0)
			{
				return 4 * xh - 1;
			}

			else if (basis_der_x == 2 && basis_der_y == 0)
			{
				return 4;
			}
			else
			{
				return 0;
			}

		}
		else if (basis_index == 2)
		{
			if (basis_der_x == 0 && basis_der_y == 0)
			{
				return 2 * pow(yh, 2) - yh;
			}
			else if (basis_der_x == 0 && basis_der_y == 1)
			{
				return 4 * yh - 1;
			}

			else if (basis_der_x == 0 && basis_der_y == 2)
			{
				return 4;
			}
			else
			{
				return 0;
			}
		}
		else if (basis_index == 3)
		{
			if (basis_der_x == 0 && basis_der_y == 0)
			{
				return -4 * pow(xh, 2) - 4*xh*yh+4*xh;
			}
			else if (basis_der_x == 1 && basis_der_y == 0)
			{
				return -8 * xh - -4*yh+4;
			}
			else if (basis_der_x == 0 && basis_der_y == 1)
			{
				return -4 * xh;
			}

			else if (basis_der_x == 2 && basis_der_y == 0)
			{
				return -8;
			}
			else
			{
				return 0;
			}
		}
		else if (basis_index == 4)
		{
			if (basis_der_x == 0 && basis_der_y == 0)
			{
				return 4*xh*yh;
			}
			else if (basis_der_x == 1 && basis_der_y == 0)
			{
				return 4 * yh;
			}

			else if (basis_der_x == 0 && basis_der_y == 1)
			{
				return 4 * xh;
			}
			else
			{
				return 0;
			}
		}
		else if (basis_index == 5)
		{
		if (basis_der_x == 0 && basis_der_y == 0)
		{
			return -4 * pow(yh, 2) - 4 * xh * yh + 4 * yh;
		}
		else if (basis_der_x == 0 && basis_der_y == 1)
		{
			return -8 * yh - -4 * xh + 4;
		}
		else if (basis_der_x == 1 && basis_der_y == 0)
		{
			return -4 * yh;
		}

		else if (basis_der_x == 0 && basis_der_y == 2)
		{
			return -8;
		}
		else
		{
			return 0;
		}
		}
	}
}


void FE_solver_2D::Caculate_vertices(int n)
{
	this->vertices_.col(0) = this->p_.col(this->t_(0, n) - 1);
	this->vertices_.col(1) = this->p_.col(this->t_(1, n) - 1);
	this->vertices_ .col(2)= this->p_.col(this->t_(2, n) - 1);

}


double FE_solver_2D::Gauss_qual_trial_test_2D(int alpha, int belta,int r,int s,int p,int q,bool T)
{
	double int_value = 0;
	for (int k = 0; k < this->number_of_gauss_points; k++)   //三点高斯插值
	{
		double cx = 0;
		if (T = false)	//稳态
		{
			cx = this->Cx(this->gauss_weight_nodes(1, k), this->gauss_weight_nodes(2, k)); //计算Cx（xi,yi）
		}
		else  //非稳态
		{
			cx= this->Cx(this->gauss_weight_nodes(1, k), this->gauss_weight_nodes(2, k),1);
		}
		

		double fai_trial_x = 0;
		fai_trial_x=this->FE_basis_local_fun_trial(this->gauss_weight_nodes(1, k), this->gauss_weight_nodes(2, k), alpha, r, s); //计算∂ψnα/∂x(xi,yi)
		
		double fai_test_x = 0;
		fai_test_x = this->FE_basis_local_fun_trial(this->gauss_weight_nodes(1, k), this->gauss_weight_nodes(2, k), belta, p, q);//计算∂ψβ/∂x(xi,yi)

		int_value += this->gauss_weight_nodes(0, k) * cx * fai_trial_x * fai_test_x;
	}
	
	return int_value;
}


double FE_solver_2D::Gauss_qual_fx_test_2D(int belta, int r, int s, int p, int q,bool T)
{
	double int_value = 0;
	for (int k = 0; k < this->number_of_gauss_points; k++)   //三点高斯插值
	{
		double fx_value = 0;
		if (T == false)
		{
			fx_value = this->fx(this->gauss_weight_nodes(1, k), this->gauss_weight_nodes(2, k)); //计算Cx（xi,yi）
		}
		else
		{
			fx_value = this->fx(this->gauss_weight_nodes(1, k), this->gauss_weight_nodes(2, k),1); //计算Cx（xi,yi）
		}
		

		double fai_test_x = 0;
		fai_test_x = this->FE_basis_local_fun_trial(this->gauss_weight_nodes(1, k), this->gauss_weight_nodes(2, k), belta, p, q);//计算ψβ(xi,yi)

		int_value += this->gauss_weight_nodes(0, k) * fx_value * fai_test_x;
	}

	return int_value;
}

double  FE_solver_2D::Gauss_qual_neumann_test_2D(int belta, int r, int s, int p, int q)
{
	double int_value = 0;
	for (int k = 0; k < this->number_of_gauss_points; k++)   //三点高斯插值
	{
		//niumann边上的高斯积分还要另写
		double cp_value = 0;
		cp_value = this->cp(this->gauss_weight_nodes(1, k), this->gauss_weight_nodes(2, k)); //计算Cp（xi,yi）

		double fai_test_x = 0;
		fai_test_x = this->FE_basis_local_fun_trial(this->gauss_weight_nodes(1, k), this->gauss_weight_nodes(2, k), belta, p, q);//计算ψβ(xi,yi)




		int_value += this->gauss_weight_nodes(0, k) * cp_value * fai_test_x;
	}

	return int_value;
}




double FE_solver_2D::Cx(double x,double y)
{
	return 1;
}

double FE_solver_2D::Cx(double x, double y,double t)
{
	return 1;
}


double  FE_solver_2D::fx(double x, double y)
{

	return -y * (1 - y) * (1 - x - pow(x , 2) / 2.0)*exp(x+y)+x*(1-x/2.0)*(3*y+pow(y,2))*exp(x+y);
}

double  FE_solver_2D::fx(double x, double y,double t)
{

	return 1;
}

double  FE_solver_2D::cp(double x, double y)
{

	return 0;//？还没写
}

//计算边界值
double  FE_solver_2D::g_boundary(double x, double y)
{
	if (x == -1)
	{
		return -1.5 * y * (1 - y) * exp(y - 1);
	}
	else if (x == 1)
	{
		return 0.5 * y * (1 - y) * exp(y + 1);
	}
	else if (y == -1)
	{
		return -2 * x* (1 - x/2.0) * exp(x- 1);
	}
	else if (y == 1)
	{
		return 0;
	}
}




double FE_solver_2D::Real_Ux(double x, double y)
{
	return x * y * (1 - x / 2.0) * (1 - y) * exp(x + y);
}