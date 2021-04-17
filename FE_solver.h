/*************************************************
Copyright:杨昌鹄
Author:杨昌鹄
Date:2021-01-15
Description:有限元求解器的基类
**************************************************/

#pragma once
#include<iostream>
#include<string>
#include<math.h>
#include <Eigen/Dense>
using namespace Eigen;     // 改成这样亦可 using Eigen::MatrixXd; 
#include<iostream>

using namespace std;



class FE_solver
{
public:

	virtual void autoRun() = 0;

	//计算P、Pb_trial、Pb_test
	virtual void Generate_PT()=0;

	//重载虚函数，因为二维网格类型有三角形和四边形
	virtual void Generate_PT(int mesh_type) = 0;

	//设定边界条件
	virtual void Generate_BoundaryNodes()=0;

	//组装A矩阵
	virtual void Assemble_matrix_A()=0;

	//组装b向量
	virtual void Assemble_b()=0;

	//处理边界条件
	virtual void Treat_Boundary()=0;


	//计算最大误差
	virtual void Compute_Error()=0;


	//计算高斯点
	virtual void Compute_Gauss(int n) = 0;

	//显示基本信息
	virtual void Print_message_normal() = 0;


	void Solution()
	{

		this->solution_ = MatrixXd::Zero(this->nb_test_, 1);
		this->solution_ = this->a_matrix_.inverse() * (this->b_vector_);
		cout << "solution_:" << endl;
		cout << this->solution_ << endl;

	}
	//1D:左边界值g(a)
	//2D:左下边界值
	double ga_;

	//1D:右边界值g(b)
	//2D:右上边界值
	double gb_;




	//试验函数基函数类型
	int basis_type_trial_;

	//测试函数基函数类型
	int basis_type_test_;
	 
	//边界类型  11:Dirichlet-Dirichlet，12：Dirichlet-Neumann，13：Dirichlet-Robin
	int boundary_;    

	//高斯类型
	int gauss_type_;

	//网格分成n分
	int n_;

	//边界节点数
	int nbn_;
	//ux
	//MatrixXd solution_;

	//边界边数
	int nbe_;


	//1D：左边界a
	//2D:左下边界
	ArrayXd a_;

	//1D：右边界
	//2D:右上边界
	ArrayXd b_;



	//网格矩阵T_,t_的第i列表示第i个网格单元中节点的全局索引;
	MatrixXi t_;

	//单元矩阵Tb_test
	MatrixXi tb_test_;

	//trial单元矩阵Tb_trial
	MatrixXi tb_trial_;

	//边界节点
	MatrixXi boundary_nodes_;

	//A矩阵
	MatrixXd a_matrix_;

	//b向量
	MatrixXd b_vector_;

	//第i个网格单元的3/4个顶点的坐标
	//第i列表示第i个节点的坐标 
	MatrixXd vertices_;

	//局部基trial函数数量 节点数
	//一维线性：2
	//一维二次：3
	//二维线性：3
	//二维二次：6
	int number_of_local_basis_trial_;

	//局部基test函数数量 节点数
	//一维线性：2
	//一维二次：3
	//二维线性：3
	//二维二次：6
	int number_of_local_basis_test_;

	//基函数类别,就是1：区间内单调递减，2：区间内    单调递增
	int basis_index_;

	//基函数对x的几阶导
	int basis_der_x;

	//基函数对y的几阶导
	int basis_der_y;

	//高斯权值与节点矩阵，四点高斯积分的权重和节点 第一行是权重，第二行是积分点(对于二维，第二行是积分点的x坐标，第三行是积分点的y坐标)
	MatrixXd gauss_weight_nodes;

	//高斯插值点个数
	int number_of_gauss_points;

	//test函数有限元节点数，用于设定A的行数，线性基函数时，nb=n_m_，二次元时不一样
	int nb_test_; 

	//trial函数单元节点数,用于设定A的列数，线性基函数时，nb=n_m_，二次元时不一样
	int nb_trial_;

	//网格节点数
	int n_m_;

	//解
	MatrixXd solution_;

	//Robin q(b)*u(b)
	double qbub_;
};
