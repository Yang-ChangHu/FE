﻿/*************************************************
Copyright:杨昌鹄
Author:杨昌鹄
Date:2021-01-15
Description:二维有限元求解器 派生类
**************************************************/

#pragma once
#include<iostream>
#include<string>
#include <Eigen/Dense>
using namespace Eigen;
#include<iostream>
using namespace std;
//using namespace cv;
#include"FE_solver.h"

//二维求解器类
class FE_solver_2D :public FE_solver
{
public:

	/*FE_solver_1D构造函数
	输入：
			N1_:水平方向网格份数
			N2_:竖直方向网格分数
			gauss_type_：
			a_x,a_y:区域左下角坐标
			b_x,b_y：区域右上角坐标
			basis_type_trial,basis_type_test_:基函数类型 201 线性基函数 202 二次基函数*/
	FE_solver_2D(int N1_, int N2_, int mesh_type, double a_x, double a_y, double b_x, double b_y, int basis_type_trial_, int basis_type_test_);

	FE_solver_2D() {}
	//自动调用各个成员函数，实现计算
	virtual void autoRun();

	/*计算P、Pb_trial、Pb_test:
	输入：	mesh_type ：3：三角形网格:， 4：四边形网格
	调用：
	输出：	this->p_ 、this->t_、this->pb_trial_、this->tb_trial_ 、this->pb_test_、this->tb_test_ ;
	 void Generate_PT(int mesh_type);

	 计算P、Pb_trial、Pb_test:
	 输入：	mesh_type ：3：三角形网格:， 4：四边形网格
	 调用：
	 输出：	this->p_ 、this->t_、this->pb_trial_、this->tb_trial_ 、this->pb_test_、this->tb_test_ ;*/
	virtual void Generate_PT();

	/*生成边界边和边界点矩阵
	输入：
			this->N1_  ,this->N2_
	调用：
			Generate_boundary_edge()
			Generate_boundary_nodes()
	输出：
			this->boundary_nodes_
			this->boundary_edges_*/
	void Generate_BoundaryNodes();

	//组装A矩阵
	//输入：nb_test_,nb_trial,tb_test,tb_trial_,number_of_local_basis_trial_,number_of_local_basis_test_
	//调用:
	//		Gauss_qual_trial_test_2D();//用高斯积分法计算积分
	//输出：a_matrix_
	virtual void Assemble_matrix_A();

	/*组装b向量
	输入：
			n_，nb_test_,number_of_local_basis_test_
	调用：
			Gauss_qual_fx_test_2D()
	输出：
			b_vector*/
	virtual void Assemble_b();

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

	//计算无穷范误差
	virtual void Compute_Error();

	/*test基函数(目前跟trial一样，所以直接用FE_basis_local_fun_trial)
	输入：
			x:转换后的高斯积分点x坐标
			y:转换后的高斯积分点y坐标
			basis_index:test基函数类型(线性基函数：3条,二次基函数：6条)
			basis_der_x/basis_der_y:对x、y的偏导次数
	调用：
	输出：
			ψ对x和y的复合偏导.∂(basis_der_x, int basis_der_y)ψ/(∂(basis_der_x)x*∂(basis_der_y)y)*/
	double FE_basis_local_fun_test(double x, double y, int n, int basis_index, int basis_der_x, int basis_der_y);

	/*计算第n个单元的trial基函数:fai_trial(xi,yi)
	输入：
			x :x坐标
			y：y坐标
			n：第n个单元
			basis_index:基函数类型，对于三角形单元线性基函数 basis_index：0,1,2 对于三角形单元二次基函数basis_index:0，1,2,3,4,5
			basis_der_x:对x的偏导数次数
			basis_der_y:对y的偏导数次数
	调用：
			Caculate_vertices(n)
			reference_basis_2D()
	输出：
		局部基函数ψ_trial(xi,yi)及其偏导数*/
	double FE_basis_local_fun_trial(double x, double y, int n, int basis_index, int basis_der_x, int basis_der_y);

	//计算高斯积分的权重和节点(local高斯公式) 四节点（废了）
	virtual void Compute_Gauss(int n);

	/*计算第i个网格单元的n个高斯插值点权重
	输入 ：
			n:插值点个数
			i:第i个网格单元
	调用：
			 Caculate_vertices():计算出第i个网格单元的几个顶点的坐标
	输出：
			local_gauss_weight_nodes:第i个网格单元的高斯积分点以及权重矩阵*/
	MatrixXd Compute_Gauss(int n, int i);

	//显示基本信息
	virtual void Print_message_normal();

	//计算边界边矩阵 （有限元概念）
	void Generate_boundary_edge();

	//计算边界点矩阵（有限元概念）
	void Generate_boundary_nodes();

	//c(x,y)
	double  Cx(double x, double y);

	double Qx(double x, double y);

	//r(x,y) Robin边界的r(x,y)
	double rxy(double x, double y);

	//f(x)
	double fx(double x, double y);

	//Neumann边界的cp（x，y）
	double  cp(double x, double y);

	//Robin边界的cq（x，y）
	double  cq(double x, double y);

	double cr(double x, double y);

	/*计算第n个单元的trial_test高斯积分
	输入：
			r,s,p,q:trial基函数和test基函数分别对x、y的偏导数
			n：第n个单元
	调用:
		   Compute_Gauss（）
			Cx(x,y)
			FE_basis_local_fun_trial：∂ψ(r+s)/∂x(r)∂y(s)
		   FE_basis_local_fun_test：∂ψ(p+q)/∂x(p)∂y(q)  目前test基函数==trial基函数
	输出:
			integral（c（xn,yn）*∂ψtrial(r+s)/∂x(r)∂y(s)*∂ψtest(p+q)/∂x(p)∂y(q)）*/
	double Gauss_qual_trial_test_2D(int alpha, int belta, int n, int r, int s, int p, int q);

	/*计算第n个网格单元的fx_test高斯积分
	输入：
			p,q:test基函数分别对x、y的偏导数
	调用:
		   Compute_Gauss（）
			fx(xi,yi)
			FE_basis_local_fun_test：ψtest(xi,yi)
	输出:
			integral（fx（xi,yi）*ψ_nβ(xi,yi)）*/
	double Gauss_qual_fx_test_2D(int belta, int n, int p, int q);

	/*计算第k条（在第nk个网格单元上）纽曼边界上c(x,y)*p(x,y)*basis_function_test(x,y)的积分  ： P96/138中的r
	输入：boundary_edges_，pb_test_
	调用：
			cp(x,y):c(x,y)*p(x,y)
			 Compute_neumann_line_Gauss(x,y):计算边界边上的高斯积分点和权重
			FE_basis_local_fun_test：basis_function_test(x,y)*/
	double Gauss_qual_neumann_test_2D(int belta, int nk, int k, int p, int q);

	double Gauss_qual_Robin_cr_2D(int alpha, int belta, int nk, int k, int r, int s, int p, int q);

	double Gauss_qual_Robin_cq_2D(int belta, int nk, int k, int p, int q);

	/*输入：纽曼边的（a，b）
	输出： 高斯4个插值点和权重*/
	MatrixXd Compute_neumann_line_Gauss(double a, double b);

	//真实u（x）值
	double Real_Ux(double x, double y);

	/*参考局部基函数
	输入：
			xh:x_head
			yh:y_head
			basis_index:基函数类型，对于三角形单元线性基函数 basis_index：0,1,2 对于三角形单元二次基函数basis_index:0，1,2,3,4,5
						四边形单元 basis_index：0,1,2,3
			basis_der_x:对x的偏导数次数
			basis_der_y:对y的偏导数次数
	输出：
		参考基函数ψ_hat(x_hat,y_hat)或其偏导数*/
	double reference_basis_2D(double xh, double yh, int basis_index, int basis_der_x, int basis_der_y);

	/*计算第n个网格各个节点的坐标
	输入：
			n:第n个网格单元
			this->basis_type_trial_：单元类型
			this->mesh_type：网格类型
	调用：

	输出：
		:MatrixXd::Zero(2, n) n表示第n个网格节点，第一行是x坐标，第二行是y坐标*/
	MatrixXd Caculate_vertices(int n);

	//Dirichlet边界有限元节点的值
	double g_boundary(double x, double y);

	//水平方向分成N1_份网格单元
	int N1_;

	//竖直方向分成N1_份网格单元
	int N2_;

	//对y的导数阶
	int basis_der_y;

	//网格数组p_，p的第i列表示的是第i个网格点的真实坐标;
	MatrixXd p_;

	//单元数组pb_test，pb_test的第i列表示的是第i个有限元节点的真实坐标，线性函数时，p_==pb;
	MatrixXd pb_test_;

	//trial单元数组 pb_trial的第i列表示的是第i个有限元节点的真实坐标，线性函数时，p_==pb;
	MatrixXd pb_trial_;

	//边界边，第一行是边界条件类型（1：Dirichelet，2：Neumann，3：Robin），第二行是边界边属于的网格单元编号，第三四行分别是第k个边界的两个网格
	//定点编号
	MatrixXi boundary_edges_;

	//边界有限元节点
	MatrixXi boundary_nodes_;
};
