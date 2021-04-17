﻿#pragma once
#include<iostream>
#include<string>
#include <Eigen/Dense>
using namespace Eigen;
#include"FE_solver.h"
using namespace std;
class  FE_solver_2D_heat :public FE_solver
{
public:
	//非稳态c(x,y)，f（x,y,t），u（x,y,t)
	FE_solver_2D_heat(double start, double end, double dt, int N1_, int N2_, int gauss_type_, double a_x, double a_y, double b_x, double b_y, int basis_type_trial_, int basis_type_test_);
	

	virtual void autoRun();



	

	//计算P、Pb_trial、Pb_test:
	//输入：	mesh_type 三角形网格/四边形网格
	//调用：
	//输出：	this->p_ 、this->t_、this->pb_trial_、this->tb_trial_ 、this->pb_test_、this->tb_test_ ;
	void Generate_PT(int mesh_type);            //二维的PT矩阵要根据网格单元类型来生成

	virtual void Generate_PT();            //空实现

	//设定边界条件
	//输入：
	//		this->N1_  ,this->N2_
	//调用：
	//		Generate_boundary_edge()
	//		Generate_boundary_nodes()
	//输出：
	//		this->boundary_nodes_
	//		this->boundary_edges_		
	virtual void Generate_BoundaryNodes();

	//组装A矩阵
	//调用:
	//		Caculate_vertices(n);网格单元节点坐标
	//		Compute_Gauss(); //计算高斯点权重(local)
	//		this->Gauss_qual_trial_test_2D;//用高斯积分法计算积分
	//				Cx(x,y)
	//				this->FE_basis_local_fun_trial(x,y, alpha, r, s) 计算偏导
	virtual void Assemble_matrix_A();


	void IterateInTime(double theta);

	//组装b向量,在这是空实现
	virtual void Assemble_b();

	//组装b（tm）.b(tm+1)
	void Assemble_b(double tm, double tmp1);


	//处理边界条件
	virtual void Treat_Boundary();


	void Solution_heat();


	//计算最大误差
	virtual void Compute_Error();

	//trial local基函数
	//输入：
	//		x
	//		y
	//		basis_index:基函数类型，对于三角形单元线性基函数 basis_index：1,2,3 对于三角形单元二次基函数basis_index:1,2,3,4,5,6
	//		basis_der_x:对x的偏导数次数
	//		basis_der_y:对y的偏导数次数
	//		n:第n个网格单元
	//调用：
	//	reference_basis_2D(double xh, double yh, int basis_index, int basis_der_x, int basis_der_y)
	//	Caculate_vertices(n)
	//输出：
	//	局部基函数ψ(x,y)及其偏导数
	double FE_basis_local_fun_trial(double x, double y, int basis_index, int basis_der_x, int basis_der_y);


	//test基函数(目前跟test一样，所以直接用FE_basis_local_fun_trial)
	//输入：
	//		x:转换后的高斯积分点x坐标
	//		y:转换后的高斯积分点y坐标
	//		basis_index:test基函数类型(线性基函数：3条,二次基函数：6条)
	//		basis_der_x/basis_der_y:对x、y的偏导次数
	//调用：
	//输出：
	//		ψ对x和y的复合偏导.∂(basis_der_x, int basis_der_y)ψ/(∂(basis_der_x)x*∂(basis_der_y)y)
	virtual double FE_basis_local_fun_test(double x, double y, int basis_index, int basis_der_x, int basis_der_y);

	//计算高斯积分的权重和节点(local高斯公式) 四节点
	virtual void Compute_Gauss(int n);



	//显示基本信息
	virtual void Print_message_normal();

	//边界边矩阵 （有限元概念）
	void Generate_boundary_edge();

	//边界点矩阵（有限元概念）
	void Generate_boundary_nodes();



	//c(x,y,t)
	double  Cx(double x, double y);

	//f(x,y,t)
	double fx(double x, double y, double t);

	//Neumann边界的cp（x，y）
	double  cp(double x, double y);

	//计算trial_test高斯积分
	//输入：
	//		r,s,p,q:trial基函数和test基函数分别对x、y的偏导数
	//		gauss_weight_nodes：修正后的高斯积分节点及其权重
	//调用:
	//		Cx(x,y)
	//		FE_basis_local_fun_trial：∂ψ/∂x∂y
	//输出:
	//		integral（c（x,y）*▽ψ_nα*▽ψ_nβ）
	//double Gauss_qual_trial_test(int alpha, int belta);
	double Gauss_qual_trial_test_2D(int alpha, int belta, int r, int s, int p, int q,bool M);

	//计算fx_test高斯积分(非稳态)
	double Gauss_qual_fx_test_2D(int belta, int r, int s, int p, int q, double t);

	//纽曼边界积分
	double Gauss_qual_neumann_test_2D(int belta, int r, int s, int p, int q);


	//非稳态的u（x，y，t）
	double Real_Ux(double x, double y, double t);

	//参考基函数
	//输入：
	//		xh:x_head
	//		yh:y_head
	//		basis_index:基函数类型，对于三角形单元线性基函数 basis_index：1,2,3 对于三角形单元二次基函数basis_index:1,2,3,4,5,6
	//		basis_der_x:对x的偏导数次数
	//		basis_der_y:对y的偏导数次数
	//输出：
	//	参考基函数ψ_hat(x_hat,y_hat)及其偏导数
	double reference_basis_2D(double xh, double yh, int basis_index, int basis_der_x, int basis_der_y);




	//计算第n个网格三个节点的坐标
	//输入：
	//		n:第n个网格单元
	//调用：
	//	
	//输出：
	//	this->vertices_  :MatrixXd::Zero(2, 3)  每列表示网格节点，第一行是x坐标，第二行是y坐标
	void Caculate_vertices(int n);

	//边界有限元节点的值
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