/*************************************************
Copyright:�����
Author:�����
Date:2021-01-15
Description:һά����Ԫ����� ������
**************************************************/

#pragma once
#include<iostream>
#include<string>
#include <Eigen/Dense>
using namespace Eigen;     // �ĳ�������� using Eigen::MatrixXd; 
#include<iostream>
using namespace std;
//using namespace cv;
#include"FE_solver.h"

//һά�������
class FE_solver_1D:public FE_solver
{
public:

	//FE_solver_1D���캯��
	FE_solver_1D(int a_,int b_,int n_,int gauss_type_,double ga_,double gb_,int basis_type_trial_,int basis_type_test_,int boundary_,double qbub_=0);

	//����P��Pb_trial��Pb_test:
	virtual void Generate_PT();            //������д������麯�����ߴ��麯����virtual�ؼ��ֿ�ɾ��Ҳ�ɲ�ɾ��

	virtual void Generate_PT(int mesh_type);  //��ʵ�֣���һά��������

	//�趨�߽�����
	virtual void Generate_BoundaryNodes();

	//��װA����
	virtual void Assemble_matrix_A(bool T=false);

	//��װb����
	virtual void Assemble_b(bool T=false);

	//����߽�����
	virtual void Treat_Boundary();

	//���Uh
	virtual void Solution();

	//����������
	virtual void Compute_Error();

	//trial������
	 double FE_basis_local_fun_trial(double x, int basis_index, int basis_der_x);

	//virtual double FE_basis_local_fun_trial(double x, double y, int basis_index, int basis_der_x, int basis_der_y) ;

	//test������
	 double FE_basis_local_fun_test(double x, int basis_index, int basis_der_x);

	//virtual double FE_basis_local_fun_test(double x, double y, int basis_index, int basis_der_x, int basis_der_y) ;

	//�����˹���ֵ�Ȩ�غͽڵ�
	virtual void Compute_Gauss(int n);

	//��ʾ������Ϣ
	virtual void Print_message_normal();



	//c(x)
	double  Cx(double x);

	//f(x)
	double fx(double x);

	//����trial_test��˹����
	double Gauss_qual_trial_test( int alpha, int belta) ;

	//����fx_test��˹����
	double Gauss_qual_fx_test( int belta);

	//��ʵu��x��ֵ
	double Real_Ux(double x);

	//��������p_;
	RowVectorXd p_;

	//��Ԫ����pb_test
	RowVectorXd pb_test_;

	//trial��Ԫ���� pb_trial
	RowVectorXd pb_trial_;



};
