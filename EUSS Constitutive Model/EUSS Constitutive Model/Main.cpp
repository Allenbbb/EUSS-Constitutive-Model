/*
EUSS Code remastered

This code is algorithm of "Equivalent uniaxial strain-stress constitutive model"
and the model contained "CONCRETE", "ROCK", "METAL"

2017/11/23 Allen

*/

#include "Main.h"
#include "Constitutive equation.h"
#include "Equivalent Uniaxial Strain Stress Model.h"

using namespace std;


int main()
{
	ofstream Stress_Strain_Output("Stress_Strain_Output.txt");
	ofstream Ultimate_Hydro_Devia_Output("Ultimate_Hydro_Devia_Output.txt");

	Ultimate_Hydro_Devia_Output << "ult_hydro_stress" << "\t" << "ult_devia_stress" << "\t" << "ult_theta" << "\t" << "ult_hydro_strain" << "\t" << "ult_devia_strain" << "\t" << "ult_theta_s" << endl;
	Stress_Strain_Output << "Loading_State_Index" << "\t" << "Ultimate_stress_1" << "\t" << "Ultimate_stress_2" << "\t" << "Ultimate_stress_3" << "\t" << "Stress_1" << "\t" << "Stress_2" << "\t" << "Stress_3" << "\t" << "Ultimate_strain_1" << "\t" << "Ultimate_strain_2" << "\t" << "Ultimate_strain_3" << "\t" << "Strain_1" << "\t" << "Strain_2" << "\t" << "Strain_3" << endl;


	/*Input essential parameters*/
	Target_stress = -0.;
	Elastic_Modulus_1 = Elastic_Modulus_2 = Elastic_Modulus_3 = 40700;
	Poisson_ratio_12 = Poisson_ratio_13 = Poisson_ratio_21 = Poisson_ratio_23 = Poisson_ratio_31 = Poisson_ratio_32 = 0.3;
	Compressive_stress = -68.;
	Tensile_stress = -.1*Compressive_stress;
	Compressive_strain = -0.0028;
	Tensile_strain = -.1*Compressive_strain;

	/**/
	General_functions gf;
	Orthotropic_constitutive_equation oce;
	Closed_Menetrey_Willam_failure_surface cmwfs;
	Uniaxial_stress_strain_model ussm;


	/*initialize*/
	Loading_State_Index = 0;
	Material_type = 1;
	Criterion_type = 1;
	dt = 1.;

	gf.Initialize();
	gf.Modification_of_failure_surface(Criterion_type);
	oce.Initialize();
	oce.Incremental_Constitutive_Equation(Elastic_Modulus_1, Elastic_Modulus_2, Elastic_Modulus_3, Poisson_ratio_12, Poisson_ratio_13, Poisson_ratio_21, Poisson_ratio_23, Poisson_ratio_31, Poisson_ratio_32);
	ussm.Initiallize();
	cmwfs.Initiallize();

	for (int i = 0; i < 150; i++)
	{
		Ultimate_Hydro_Devia_Output << gf.Get_ult_hydro_stress() << "\t" << gf.Get_ult_devia_stress() << "\t" << gf.Get_ult_theta() << "\t" << gf.Get_ult_hydro_strain() << "\t" << gf.Get_ult_devia_strain() << "\t" << gf.Get_ult_theta_s() << endl;
		Stress_Strain_Output << Loading_State_Index << "\t" << Ultimate_Stress_1 << "\t" << Ultimate_Stress_2 << "\t" << Ultimate_Stress_3 << "\t" << abs(Stress_1) << "\t" << abs(Stress_2) << "\t" << abs(Stress_3) << "\t" << Ultimate_Strain_1 << "\t" << Ultimate_Strain_2 << "\t" << Ultimate_Strain_3 << "\t" << Strain_1 << "\t" << Strain_2 << "\t" << Strain_3 << endl;

		/*stress control*/
		if (abs(Stress_2) < abs(Target_stress))
		{

			dStress_1 = -1.;
			dStress_2 = -0.;
			dStress_3 = -0.;

			oce.Normal_stress_to_normal_strain(dStress_1, dStress_2, dStress_3);
			dStrain_1 = oce.dStrain_1;
			dStrain_2 = oce.dStrain_2;
			dStrain_3 = oce.dStrain_3;

			cout << "stress control" << endl;
		}

		/*strain control*/
		else if (abs(Stress_2) >= abs(Target_stress))
		{
			dStress_1 = -1.;
			dStress_2 = -0.;
			dStress_3 = -0.;

			oce.Normal_stress_to_normal_strain(dStress_1, dStress_2, dStress_3);
			dStrain_1 = oce.dStrain_1;
			dStrain_2 = oce.dStrain_2;
			dStrain_3 = oce.dStrain_3;

			// 			dStrain_1 = -0.00005;
			// 			dStrain_2 = 0.;
			// 			dStrain_3 = 0.;

			oce.Normal_strain_to_normal_stress(dStrain_1, dStrain_2, dStrain_3);
			dStress_1 = oce.dStress_1;
			dStress_2 = oce.dStress_2;
			dStress_3 = oce.dStress_3;

			
			cout << "strain control" << endl;
		}



		Strain_1 += dStrain_1;
		Strain_2 += dStrain_2;
		Strain_3 += dStrain_3;

		oce.Normal_strain_to_equiv_strain(dStrain_1, dStrain_2, dStrain_3);
		dEquiv_strain_1 = oce.dEquiv_strain_1;
		dEquiv_strain_2 = oce.dEquiv_strain_2;
		dEquiv_strain_3 = oce.dEquiv_strain_3;

		gf.Stress_invariants_to_HW(dStress_1, dStress_2, dStress_3);
		gf.Strain_invariants_to_HW(dEquiv_strain_1, dEquiv_strain_2, dEquiv_strain_3);


		if (dt >= 1.)		/*Before apex*/
		{
			Loading_State_Index = 0;

			Stress_1 += dStress_1;
			Stress_2 += dStress_2;
			Stress_3 += dStress_3;

			switch (Criterion_type)		/*Choose failure surface*/
			{
			case 1:
				cmwfs.Stress_failure_surface(Compressive_stress, Tensile_stress);
				cmwfs.Strain_failure_surface(Compressive_strain, Tensile_strain);
				break;

			default:
				break;
			}

			gf.Determination_of_ultimate_stresses();
			gf.Determination_of_ultimate_strains();
			gf.Coordinate_transform_to_principal_stresses();
			gf.Coordinate_transform_to_principal_strains();
			gf.Ultimate_stress_strain_correction();

			Ultimate_Stress_1 = gf.Get_Ulti_stress_1();		/*Calculate ultimate stresses*/
			Ultimate_Stress_2 = gf.Get_Ulti_stress_2();
			Ultimate_Stress_3 = gf.Get_Ulti_stress_3();
			Ultimate_Strain_1 = gf.Get_Ulti_strain_1();		/*Calculate ultimate strains*/
			Ultimate_Strain_2 = gf.Get_Ulti_strain_2();
			Ultimate_Strain_3 = gf.Get_Ulti_strain_3();

		}

		else if (dt < 1.)	/*After apex*/
		{
			Loading_State_Index = 1;

			Stress_1 -= dStress_1;
			Stress_2 -= dStress_2;
			Stress_3 -= dStress_3;

			switch (Criterion_type)		/*Choose softening surface*/
			{
			case 1:
				cmwfs.Stress_softening_surface(Compressive_stress, Tensile_stress);
				cmwfs.Strain_softening_surface(Compressive_strain, Tensile_strain);
				break;

			default:
				break;
			}


		}


		/*Determined loading state*/
		//Loading_State_Index = gf.Current_loading_state();

		/**/



		//ussm.Popovics_function(1, dEquiv_strain_1, dEquiv_strain_2, dEquiv_strain_3, Elastic_Modulus_1, Elastic_Modulus_2, Elastic_Modulus_3, Ultimate_Stress_1, Ultimate_Stress_2, Ultimate_Stress_3, Ultimate_Strain_1, Ultimate_Strain_2, Ultimate_Strain_3);
		ussm.Saenz_cubic_function(1, dEquiv_strain_1, dEquiv_strain_2, dEquiv_strain_3, Elastic_Modulus_1, Elastic_Modulus_2, Elastic_Modulus_3, Ultimate_Stress_1, Ultimate_Stress_2, Ultimate_Stress_3, Ultimate_Strain_1, Ultimate_Strain_2, Ultimate_Strain_3);
		Elastic_Modulus_1 = ussm.Get_Elas_modu_1();
		Elastic_Modulus_2 = ussm.Get_Elas_modu_2();
		Elastic_Modulus_3 = ussm.Get_Elas_modu_3();
		/*

				ussm.Kupfer_function();
				Poisson_ratio_12 = ussm.Get_Poisson_ratio_12();
				Poisson_ratio_13 = ussm.Get_Poisson_ratio_13();
				Poisson_ratio_21 = ussm.Get_Poisson_ratio_21();
				Poisson_ratio_23 = ussm.Get_Poisson_ratio_23();
				Poisson_ratio_31 = ussm.Get_Poisson_ratio_31();
				Poisson_ratio_32 = ussm.Get_Poisson_ratio_32();
		*/

		/*Update mechanical properties*/
		oce.Incremental_Constitutive_Equation(Elastic_Modulus_1, Elastic_Modulus_2, Elastic_Modulus_3, Poisson_ratio_12, Poisson_ratio_13, Poisson_ratio_21, Poisson_ratio_23, Poisson_ratio_31, Poisson_ratio_32);

		dt = gf.General_functions_output();
	}


	cmwfs.Plot_failure_surface(Compressive_stress, Compressive_strain);


	system("pause");
	return 0;

}