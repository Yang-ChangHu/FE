#include<iostream>
using namespace std;
#include "FE_solver.h"
#include "FE_solver_1D.h"
#include "FE_solver_2D.h"
#include<string>
#include<math.h>
#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;     
using namespace std;


//一维有限元求解器
void test_1D()
{
	FE_solver* fe = NULL;

	//一维有限元
	fe = new FE_solver_1D(0, 1, 4, 1, 0.0, cos(1), 101, 101, 11);
	fe->Generate_PT();
	fe->Generate_BoundaryNodes();
	fe->Print_message_normal();
	fe->Assemble_matrix_A();
	fe->Assemble_b();
	fe->Treat_Boundary();
	fe->Solution();
	fe->Compute_Error();
}

//二维有限元求解器
void test_2D()
{
	FE_solver* fe = NULL;
	fe = new FE_solver_2D(16,16,1,-1,-1,1,1,201,201);
	//int N1_, int N2_, int gauss_type_, double a_x, double a_y, double b_x, double b_y, int basis_type_trial_, int basis_type_test_,
	fe->Generate_PT(3);
	fe->Generate_BoundaryNodes();
	fe->Assemble_matrix_A();
	fe->Assemble_b();
	fe->Treat_Boundary();
	fe->Solution();
	fe->Compute_Error();

	
}

int main()
{
	//一维有限元求解器
	//test_1D();

	//二维有限元
	test_2D();

	system("pause");
	return 0;
}


