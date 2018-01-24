#include "Constitutive equation.h"
#include <fstream>


using namespace std;
ofstream Material_Mechanical_Properties("Material_Mechanical_Properties.txt");

void Orthotropic_constitutive_equation::Initialize()
{
	Material_Mechanical_Properties << "Elastic_Modulus_1" << "\t" << "Elastic_Modulus_2" << "\t" << "Elastic_Modulus_3" << "Poisson_ratio_12" << "\t" << "Poisson_ratio_13" << "\t" << "Poisson_ratio_23" << endl;

}
void Orthotropic_constitutive_equation::Incremental_Constitutive_Equation(double Elastic_modulus_1, double Elastic_modulus_2, double Elastic_modulus_3, double Poisson_ratio_12, double Poisson_ratio_13, double Poisson_ratio_21, double Poisson_ratio_23, double Poisson_ratio_31, double Poisson_ratio_32)
{
	/*

	Stress*D = Strain

	Stress = C *(Strain) = H*E*(Strain) = E*(Equivelaent-Uniaxial-Strain)

		/ E1  0   0 \
	E=|  0   E2  0   |
		\ 0   0   E3/

		/   V11     V12   V13   \
	V=|    V21     V22   V23     |
		\   V31     V32   V33   /

		/	   1/E1		    v12/E2		 v13/E3	 \
	D=|    v21/E1		     1/E2		 v23/E3	  |
		\    v31/E1		    v32/E2		   1/E3	 /

					/	E1(-v23 v32 + v22 v33)	  E1(v13 v32 +v12 v33)	    E1 (v13 v22 + v12 v23)	   \
	C=(1/omega)		|   E2(v23 v31 + v21 v33)	  E2(-v13 v31 + v11 v33)	E2 (v13 v21 + v11 v23)	    |
					\   E3(v22 v31 + v21 v32)	  E3(v12 v31 + v11 v32)	    E3(-v12 v21 + v11 v22)	   /

				/ (-v23 v32 + v22 v33)	 (v13 v32 -v12 v33)    (-v13 v22 + v12 v23)  \
	H=(1/omega)|  (v23 v31 - v21 v33)	 (-v13 v31 + v11 v33)  (v13 v21 - v11 v23)     |
				\ (-v22 v31 + v21 v32)	 (v12 v31 - v11 v32)   (-v12 v21 + v11 v22)	 /

	omega=-v13 v22 v31 + v12 v23 v31 + v13 v21 v32 - v11 v23 v32 - v12 v21 v33 +	v11 v22 v33

	L = INVERSE(H)

					/ (H22H33-H23H32) (H13H32-H12H33) (H12H23-H13H22)	\
	L=(1/omegaH)	| (H23H31-H21H33) (H11H33-H13H31) (H21H13-H23H11)    |
					\ (H32H21-H31H22) (H31H12-H32H11) (H11H22-H12H21)	/

	omegaH = H11H22H33 - H11H23H32 - H22H13H31 - H33H12H21 + H12H23H31 + H13H32H21
	*/


	/* E */

	Matrix_E[0][0] = Elastic_modulus_1;
	Matrix_E[0][1] = 0.0;
	Matrix_E[0][2] = 0.0;
	Matrix_E[1][0] = 0.0;
	Matrix_E[1][1] = Elastic_modulus_2;
	Matrix_E[1][2] = 0.0;
	Matrix_E[2][0] = 0.0;
	Matrix_E[2][1] = 0.0;
	Matrix_E[2][2] = Elastic_modulus_3;

	/* v */

	Matrix_v[0][0] = 1.0;
	Matrix_v[0][1] = Poisson_ratio_12;
	Matrix_v[0][2] = Poisson_ratio_13;

	Matrix_v[1][0] = Poisson_ratio_21;
	Matrix_v[1][1] = 1.0;
	Matrix_v[1][2] = Poisson_ratio_23;

	Matrix_v[2][0] = Poisson_ratio_31;
	Matrix_v[2][1] = Poisson_ratio_32;
	Matrix_v[2][2] = 1.0;

	/*H*/

	Omega = (-Matrix_v[0][2] * Matrix_v[1][1] * Matrix_v[2][0]) - (Matrix_v[0][1] * Matrix_v[1][2] * Matrix_v[2][0]) - (Matrix_v[0][2] * Matrix_v[1][0] * Matrix_v[2][1]) - (Matrix_v[0][0] * Matrix_v[1][2] * Matrix_v[2][1]) - (Matrix_v[0][1] * Matrix_v[1][0] * Matrix_v[2][2]) + (Matrix_v[0][0] * Matrix_v[1][1] * Matrix_v[2][2]);

	Matrix_H[0][0] = (-Matrix_v[1][2] * Matrix_v[2][1] + Matrix_v[1][1] * Matrix_v[2][2]) / Omega;
	Matrix_H[0][1] = (Matrix_v[0][2] * Matrix_v[2][1] + Matrix_v[0][1] * Matrix_v[2][2]) / Omega;
	Matrix_H[0][2] = (Matrix_v[0][2] * Matrix_v[1][1] + Matrix_v[0][1] * Matrix_v[1][2]) / Omega;
	Matrix_H[1][0] = (Matrix_v[1][2] * Matrix_v[2][0] + Matrix_v[1][0] * Matrix_v[2][2]) / Omega;
	Matrix_H[1][1] = (-Matrix_v[0][2] * Matrix_v[2][0] + Matrix_v[0][0] * Matrix_v[2][2]) / Omega;
	Matrix_H[1][2] = (Matrix_v[0][2] * Matrix_v[1][0] + Matrix_v[0][0] * Matrix_v[1][2]) / Omega;
	Matrix_H[2][0] = (Matrix_v[1][1] * Matrix_v[2][0] + Matrix_v[1][0] * Matrix_v[2][1]) / Omega;
	Matrix_H[2][1] = (Matrix_v[0][1] * Matrix_v[2][0] + Matrix_v[0][0] * Matrix_v[2][1]) / Omega;
	Matrix_H[2][2] = (-Matrix_v[0][1] * Matrix_v[1][0] + Matrix_v[0][0] * Matrix_v[1][1]) / Omega;

	/*D*/

	Matrix_D[0][0] = Matrix_v[0][0] / Matrix_E[0][0];
	Matrix_D[0][1] = -Matrix_v[0][1] / Matrix_E[1][1];
	Matrix_D[0][2] = -Matrix_v[0][2] / Matrix_E[2][2];
	Matrix_D[1][0] = -Matrix_v[1][0] / Matrix_E[0][0];
	Matrix_D[1][1] = Matrix_v[1][1] / Matrix_E[1][1];
	Matrix_D[1][2] = -Matrix_v[1][2] / Matrix_E[2][2];
	Matrix_D[2][0] = -Matrix_v[2][0] / Matrix_E[0][0];
	Matrix_D[2][1] = -Matrix_v[2][1] / Matrix_E[1][1];
	Matrix_D[2][2] = Matrix_v[2][2] / Matrix_E[2][2];


	/*C*/
	Matrix_C[0][0] = Matrix_E[0][0] * (-Matrix_v[1][2] * Matrix_v[2][1] + Matrix_v[1][1] * Matrix_v[2][2]) / Omega;
	Matrix_C[0][1] = Matrix_E[0][0] * (Matrix_v[0][2] * Matrix_v[2][1] + Matrix_v[0][1] * Matrix_v[2][2]) / Omega;
	Matrix_C[0][2] = Matrix_E[0][0] * (Matrix_v[0][2] * Matrix_v[1][1] + Matrix_v[0][1] * Matrix_v[1][2]) / Omega;
	Matrix_C[1][0] = Matrix_E[1][1] * (Matrix_v[1][2] * Matrix_v[2][0] + Matrix_v[1][0] * Matrix_v[2][2]) / Omega;
	Matrix_C[1][1] = Matrix_E[1][1] * (-Matrix_v[0][2] * Matrix_v[2][0] + Matrix_v[0][0] * Matrix_v[2][2]) / Omega;
	Matrix_C[1][2] = Matrix_E[1][1] * (Matrix_v[0][2] * Matrix_v[1][0] + Matrix_v[0][0] * Matrix_v[1][2]) / Omega;
	Matrix_C[2][0] = Matrix_E[2][2] * (Matrix_v[1][1] * Matrix_v[2][0] + Matrix_v[1][0] * Matrix_v[2][1]) / Omega;
	Matrix_C[2][1] = Matrix_E[2][2] * (Matrix_v[0][1] * Matrix_v[2][0] + Matrix_v[0][0] * Matrix_v[2][1]) / Omega;
	Matrix_C[2][2] = Matrix_E[2][2] * (-Matrix_v[0][1] * Matrix_v[1][0] + Matrix_v[0][0] * Matrix_v[1][1]) / Omega;


	/*L*/
	OmegaH = (-Matrix_H[0][2] * Matrix_H[1][1] * Matrix_H[2][0]) + (Matrix_H[0][1] * Matrix_H[1][2] * Matrix_H[2][0]) + (Matrix_H[0][2] * Matrix_H[1][0] * Matrix_H[2][1]) - (Matrix_H[0][0] * Matrix_H[1][2] * Matrix_H[2][1]) - (Matrix_H[0][1] * Matrix_H[1][0] * Matrix_H[2][2]) + (Matrix_H[0][0] * Matrix_H[1][1] * Matrix_H[2][2]);

	Matrix_L[0][0] = (-Matrix_H[1][2] * Matrix_H[2][1] + Matrix_H[1][1] * Matrix_H[2][2]) / OmegaH;
	Matrix_L[0][1] = (Matrix_H[0][2] * Matrix_H[2][1] - Matrix_H[0][1] * Matrix_H[2][2]) / OmegaH;
	Matrix_L[0][2] = (-Matrix_H[0][2] * Matrix_H[1][1] + Matrix_H[0][1] * Matrix_H[1][2]) / OmegaH;
	Matrix_L[1][0] = (Matrix_H[1][2] * Matrix_H[2][0] - Matrix_H[1][0] * Matrix_H[2][2]) / OmegaH;
	Matrix_L[1][1] = (-Matrix_H[0][2] * Matrix_H[2][0] + Matrix_H[0][0] * Matrix_H[2][2]) / OmegaH;
	Matrix_L[1][2] = (Matrix_H[0][2] * Matrix_H[1][0] - Matrix_H[0][0] * Matrix_H[1][2]) / OmegaH;
	Matrix_L[2][0] = (-Matrix_H[1][1] * Matrix_H[2][0] + Matrix_H[1][0] * Matrix_H[2][1]) / OmegaH;
	Matrix_L[2][1] = (Matrix_H[0][1] * Matrix_H[2][0] - Matrix_H[0][0] * Matrix_H[2][1]) / OmegaH;
	Matrix_L[2][2] = (-Matrix_H[0][1] * Matrix_H[1][0] + Matrix_H[0][0] * Matrix_H[1][1]) / OmegaH;

	/*
	cout<<"V:"<<endl;
	for (int i=0.0;i<=2.0;i++)
	{
	for (int j=0.0;j<=2.0;j++)
	{
	cout<<Matrix_v[i][j]<<"\t";
	}
	cout<<endl;
	}
	cout<<"---------------------------------"<<endl;


	cout<<"E:"<<endl;
	for (int i=0.0;i<=2.0;i++)
	{
	for (int j=0.0;j<=2.0;j++)
	{

	cout<<Matrix_E[i][j]<<"\t";
	}
	cout<<endl;
	}
	cout<<"---------------------------------"<<endl;

	cout<<"H:"<<endl;
	for (int i=0.0;i<=2.0;i++)
	{
	for (int j=0.0;j<=2.0;j++)
	{

	cout<<Matrix_H[i][j]<<"\t";
	}
	cout<<endl;
	}
	cout<<"---------------------------------"<<endl;

	cout<<"C:"<<endl;
	for (int i=0.0;i<=2.0;i++)
	{
	for (int j=0.0;j<=2.0;j++)
	{

	cout<<Matrix_C[i][j]<<"\t";
	}
	cout<<endl;
	}
	cout<<"---------------------------------"<<endl;

	cout<<"D:"<<endl;
	for (int i=0.0;i<=2.0;i++)
	{
	for (int j=0.0;j<=2.0;j++)
	{

	cout<<Matrix_D[i][j]<<"\t";
	}
	cout<<endl;
	}
	cout<<"---------------------------------"<<endl;

	cout<<"L:"<<endl;
	for (int i=0.0;i<=2.0;i++)
	{
	for (int j=0.0;j<=2.0;j++)
	{

	cout<<Matrix_L[i][j]<<"\t";
	}
	cout<<endl;
	}
	cout<<"---------------------------------"<<endl;*/

	Material_Mechanical_Properties << Elastic_modulus_1 << "\t" << Elastic_modulus_2 << "\t" << Elastic_modulus_3 << "\t" << Poisson_ratio_12 << "\t" << Poisson_ratio_13 << "\t" << Poisson_ratio_23 << endl;

}

void Orthotropic_constitutive_equation::Normal_stress_to_normal_strain(double dStress_1, double dStress_2, double dStress_3)
{

	dStrain_1 = Matrix_D[0][0] * dStress_1 + Matrix_D[0][1] * dStress_2 + Matrix_D[0][2] * dStress_3;
	dStrain_2 = Matrix_D[1][0] * dStress_1 + Matrix_D[1][1] * dStress_2 + Matrix_D[1][2] * dStress_3;
	dStrain_3 = Matrix_D[2][0] * dStress_1 + Matrix_D[2][1] * dStress_2 + Matrix_D[2][2] * dStress_3;

}

void Orthotropic_constitutive_equation::Normal_strain_to_normal_stress(double dStrain_1, double dStrain_2, double dStrain_3)
{
	dStress_1 = Matrix_C[0][0] * dStrain_1 + Matrix_C[0][1] * dStrain_2 + Matrix_C[0][2] * dStrain_3;
	dStress_2 = Matrix_C[1][0] * dStrain_1 + Matrix_C[1][1] * dStrain_2 + Matrix_C[1][2] * dStrain_3;
	dStress_3 = Matrix_C[2][0] * dStrain_1 + Matrix_C[2][1] * dStrain_2 + Matrix_C[2][2] * dStrain_3;

}

void Orthotropic_constitutive_equation::Normal_strain_to_equiv_strain(double dStrin_1, double dStrin_2, double dStrin_3)
{

	dEquiv_strain_1 = Matrix_H[0][0] * dStrin_1 + Matrix_H[0][1] * dStrin_2 + Matrix_H[0][2] * dStrin_3;
	dEquiv_strain_2 = Matrix_H[1][0] * dStrin_1 + Matrix_H[1][1] * dStrin_2 + Matrix_H[1][2] * dStrin_3;
	dEquiv_strain_3 = Matrix_H[2][0] * dStrin_1 + Matrix_H[2][1] * dStrin_2 + Matrix_H[2][2] * dStrin_3;

}

void Orthotropic_constitutive_equation::Equiv_strain_to_normal_strain(double dEqui_Uni_strain_1, double dEqui_Uni_strain_2, double dEqui_Uni_strain_3)
{
	dStrain_1 = Matrix_L[0][0] * dEqui_Uni_strain_1 + Matrix_L[0][1] * dEqui_Uni_strain_2 + Matrix_L[0][2] * dEqui_Uni_strain_3;
	dStrain_2 = Matrix_L[1][0] * dEqui_Uni_strain_1 + Matrix_L[1][1] * dEqui_Uni_strain_2 + Matrix_L[1][2] * dEqui_Uni_strain_3;
	dStrain_3 = Matrix_L[2][0] * dEqui_Uni_strain_1 + Matrix_L[2][1] * dEqui_Uni_strain_2 + Matrix_L[2][2] * dEqui_Uni_strain_3;

}