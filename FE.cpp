#include<iostream>
using namespace std;
#include "FE_solver.h"
#include "FE_solver_1D.h"
#include "FE_solver_2D.h"
#include"FE_solver_2D_heat.h"
#include<string>
#include<math.h>
#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;     
using namespace std;


//һά����Ԫ�����
void test_1D()
{
	FE_solver* fe = NULL;

	//һά����Ԫ
	fe = new FE_solver_1D(0, 1, 4, 1, 0.0, cos(1), 101, 101, 11);
	fe->Generate_PT();
	fe->Generate_BoundaryNodes();
	fe->Print_message_normal();
	fe->Assemble_matrix_A();
	fe->Assemble_b();
	fe->Treat_Boundary_Neumann();
	fe->Treat_Boundary_Robin();
	fe->Treat_Boundary_Dirichlet();
	fe->Solution();
	fe->Compute_Error();
}

//��ά����Ԫ�����
void test_2D()
{
	FE_solver* fe = NULL;
	fe = new FE_solver_2D(16,16,1,-1,-1,1,1,201,201);
	//int N1_, int N2_, int gauss_type_, double a_x, double a_y, double b_x, double b_y, int basis_type_trial_, int basis_type_test_,
	fe->autoRun();


	
}

//��ά����Ԫ���������̬
void test_2D_heat()
{
	FE_solver* fe = NULL;
	fe = new FE_solver_2D_heat(0,1,1/4.0,8,4,1,0,0,2,1,201,201);
	//double start, double end, double dt, int N1_, int N2_, int gauss_type_, double a_x, double a_y, double b_x, double b_y, int basis_type_trial_, int basis_type_test_
	fe->autoRun();


}

int main()
{
	//һά����Ԫ�����
	//test_1D();

	//��ά����Ԫ
	test_2D();

	//��ά����̬
	//test_2D_heat();

	system("pause");
	return 0;
}


