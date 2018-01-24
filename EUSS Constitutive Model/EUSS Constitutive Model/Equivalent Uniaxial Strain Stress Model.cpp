#include "Equivalent Uniaxial Strain Stress Model.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

ofstream Haigh_Westergaard_Output("Haigh_Westergaard_Output.txt");
ofstream Uniaxial_Model_Output("Uniaxial_model_output.txt");


/*Index*/
int Stress_state_type, Modulus_type;
int Ultimate_stress_location, Ultimate_strain_location;
int Loading_state_index;
int i, j, k;

/*Haigh-Westergaard*/
double dI1, dJ2, dJ3;
double dI1_s, dJ2_s, dJ3_s;
double dHydro_stress, dDevia_stress, _Theta, Cosine3Theta;
double dHydro_strain, dDevia_strain, _Theta_s, Cosine3Theta_s;

/*Point 0 a b*/
double Point_0, Point_a, Point_b;
double Point_0_s, Point_a_s, Point_b_s;

/**/
double Hydro_stress, Devia_stress, Theta, Hydro_strain, Devia_strain, Theta_s;
double Devia_stress_a, dt_Meridian, dt_Cap, Hydro_stress_M, Hydro_stress_C, Devia_stress_M, Devia_stress_C;
double Devia_strain_a, dt_Meridian_s, dt_Cap_s, Hydro_strain_M, Hydro_strain_C, Devia_strain_M, Devia_strain_C;

double w0, w1, w2, w3, w4;
double w0_s, w1_s, w2_s, w3_s, w4_s;
double Equiv_strain_1, Equiv_strain_2, Equiv_strain_3;
double Elas_modu_1, Elas_modu_2, Elas_modu_3;

double Plot_devia_stress, Plot_devia_strain;	//Plot

/*General functions*/


void General_functions::Initialize()
{
	i = j = k = 1;
	Temp_1 = Temp_2 = Temp_3 = 0.;
	Ult_hydro_stress = Ult_devia_stress = Ult_theta = 0.;
	Ult_hydro_strain = Ult_devia_strain = Ult_theta_s = 0.;
	Hydro_stress = Devia_stress = Theta = 0;
	dt_Meridian = dt_Cap = dt_Meridian_s = dt_Cap_s = 0.;

	Haigh_Westergaard_Output << "dt" << "\t" << "Hydro_stress" << "\t" << "Devia_stress" << "\t" << "Theta" << "\t" << "Hydro_strain" << "\t" << "Devia_strain" << "\t" << "Theta_s" << endl;

	Haigh_Westergaard_Output << dt_Meridian << "\t" << Hydro_stress << "\t" << Devia_stress << "\t" << Theta << "\t" << Hydro_strain << "\t" << Devia_strain << "\t" << Theta_s << endl;

}

void General_functions::Stress_invariants_to_HW(double dSigma_1, double dSigma_2, double dSigma_3)
{
	dI1 = dSigma_1 + dSigma_2 + dSigma_3;
	dJ2 = (pow(dSigma_1 - dSigma_2, 2) + pow(dSigma_2 - dSigma_3, 2) + pow(dSigma_3 - dSigma_1, 2)) / 6;
	dJ3 = ((2 * dSigma_1 - dSigma_2 - dSigma_3)*(2 * dSigma_2 - dSigma_1 - dSigma_3)*(2 * dSigma_3 - dSigma_1 - dSigma_2)) / 27;

	dHydro_stress = dI1 / sqrt(3.0);
	dDevia_stress = sqrt(2 * dJ2);
	Cosine3Theta = (3 * sqrt(3.0) / 2) * dJ3 / pow(dJ2, 1.5);

	if (Cosine3Theta == -0) Cosine3Theta = 0;
	if (dI1 == 0 && dJ2 == 0) Cosine3Theta = 0;

	Hydro_stress += dHydro_stress;
	Devia_stress += dDevia_stress;

	if (Cosine3Theta > 0) {

		if (Cosine3Theta >= 1 && Cosine3Theta <= 1.00005) {//注意軸壓力-1:-0.52:0與-1:-1:0

			Cosine3Theta = 1.0;
			Theta = (1.0 / 3.0) * acos(Cosine3Theta) * 180. / 3.14159265;
			//cout<<"dTheta(Cos 3theta)= "<<dTheta<<endl;
		}

		else {

			Theta = (1.0 / 3.0) * acos(Cosine3Theta) * 180. / 3.14159265;
			//cout<<"dTheta(Cos 3theta)= "<<dTheta<<endl;
		}
	}

	//else if (dJ2 == 0) 

	else if (Cosine3Theta < 0) {

		if (Cosine3Theta <= -1 && Cosine3Theta >= -1.00005) {

			Cosine3Theta = -1.0;
			Theta = (1.0 / 3.0) * acos(Cosine3Theta) * 180. / 3.14159265;
			//cout<<"dTheta(Cos 3theta)= "<<dTheta<<endl;
		}
		else {

			Theta = (1.0 / 3.0) * acos(Cosine3Theta) * 180. / 3.14159265;
			//cout<<"dTheta(Cos 3theta)= "<<dTheta<<endl;
		}
	}
	else {

		Theta = (1.0 / 3.0) * acos(Cosine3Theta) * 180. / 3.14159265;
		//cout<<"dTheta(Cos 3theta)= "<<dTheta<<endl;
	}

	if (dSigma_2 > dSigma_1 && dSigma_1 > dSigma_3)
	{
		_Theta = 120 - Theta;
	}
	else	if (dSigma_2 > dSigma_1 && dSigma_1 == dSigma_3)
	{
		_Theta = Theta + 120.0;
	}
	else	if (dSigma_2 > dSigma_3 && dSigma_3 > dSigma_1)
	{
		_Theta = Theta + 120.0;
	}
	else if (dSigma_2 == dSigma_3 && dSigma_3 > dSigma_1)
	{
		_Theta = Theta + 120.0;
	}
	else if (dSigma_3 > dSigma_2 && dSigma_2 > dSigma_1)
	{
		_Theta = 240.0 - Theta;
	}
	else if (dSigma_3 > dSigma_2 && dSigma_2 == dSigma_1)
	{
		_Theta = Theta + 240;
	}
	else if (dSigma_3 > dSigma_1 && dSigma_1 > dSigma_2)
	{
		_Theta = Theta + 240.0;
	}
	else if (dSigma_3 == dSigma_1 && dSigma_1 > dSigma_2)
	{
		_Theta = Theta + 240.0;
	}
	else if (dSigma_3 > dSigma_1 && dSigma_1 > dSigma_2)
	{
		_Theta = 360.0 - Theta;
	}


	// 	cout << "Hydro_stress=" << Hydro_stress << endl;
	// 	cout << "Devia_stress=" << Devia_stress << endl;
}

void General_functions::Strain_invariants_to_HW(double dEqui_Uni_strain_1, double dEqui_Uni_strain_2, double dEqui_Uni_strain_3)
{
	dI1_s = dEqui_Uni_strain_1 + dEqui_Uni_strain_2 + dEqui_Uni_strain_3;
	dJ2_s = (pow(dEqui_Uni_strain_1 - dEqui_Uni_strain_2, 2) + pow(dEqui_Uni_strain_2 - dEqui_Uni_strain_3, 2) + pow(dEqui_Uni_strain_3 - dEqui_Uni_strain_1, 2)) / 6;
	dJ3_s = ((2 * dEqui_Uni_strain_1 - dEqui_Uni_strain_2 - dEqui_Uni_strain_3) * (2 * dEqui_Uni_strain_2 - dEqui_Uni_strain_1 - dEqui_Uni_strain_3) * (2 * dEqui_Uni_strain_3 - dEqui_Uni_strain_1 - dEqui_Uni_strain_2)) / 27;

	dHydro_strain = dI1_s / sqrt(3);
	dDevia_strain = sqrt(2.* dJ2_s);
	Cosine3Theta_s = (3.* sqrt(3) / 2.) * dJ3_s / pow(dJ2_s, 1.5);

	if (Cosine3Theta_s == -0)  Cosine3Theta_s = 0;
	if (dI1_s == 0 && dJ2_s == 0)  Cosine3Theta_s = 0;

	Hydro_strain += dHydro_strain;
	Devia_strain += dDevia_strain;

	if (Cosine3Theta_s > 0)
	{

		if (Cosine3Theta_s >= 1 && Cosine3Theta_s <= 1.00005)
		{
			Cosine3Theta_s = 1.0;
			Theta_s = (1.0 / 3.0)*acos(Cosine3Theta_s)*180. / 3.14159265;
			//cout<<"Theta_s= "<<Theta_s<<endl;
		}

		else
		{
			Theta_s = (1.0 / 3.0)*acos(Cosine3Theta_s)*180. / 3.14159265;
			//cout<<"Theta_s= "<<Theta_s<<endl;
		}
	}

	else if (dJ2_s == 0)
	{
		//cout<<"There no Lord angle about this stress state"<<endl;
	}

	else if (Cosine3Theta_s < 0)
	{
		if (Cosine3Theta_s <= -1 && Cosine3Theta_s >= -1.00005)
		{

			Cosine3Theta_s = -1.0;
			Theta_s = (1.0 / 3.0) * acos(Cosine3Theta_s) * 180. / 3.14159265;
			//cout<<"Theta_s= "<<Theta_s<<endl;
		}
		else
		{

			Theta_s = (1.0 / 3.0) * acos(Cosine3Theta_s) * 180. / 3.14159265;
			//cout<<"Theta_s= "<<Theta_s<<endl;
		}
	}

	else
	{
		Theta_s = (1.0 / 3.0) * acos(Cosine3Theta_s) * 180. / 3.14159265;
		//cout<<"dTheta= "<<Theta_s<<endl;       
	}

	// 	cout << "Hydro_strain=" << Hydro_strain << endl;
	// 	cout << "Devia_strain=" << Devia_strain << endl;
}

void General_functions::Modification_of_failure_surface(int Criterion_type)
{
	switch (Criterion_type)
	{
	case 1:
		w0 = 1.;
		w1 = 1.;
		w2 = 1.;
		w3 = 1.;
		w4 = 1.;

		w0_s = 1.;
		w1_s = 1.;
		w2_s = 1.;
		w3_s = 1.;
		w4_s = 1.;

		break;

	case 2:


	default:
		break;
	}
}

int General_functions::Current_loading_state()
{
	switch (Ultimate_stress_location)
	{

	case 1:
		if (dt_Meridian <= 1.)
		{
			Loading_state_index = 1;
			//cout << "dt_Meridian=" << dt_Meridian << endl;
		}
		break;
	case 2:
		if (dt_Meridian <= 1.)
		{
			Loading_state_index = 1;
		}
		break;
	case 3:
		if (dt_Cap <= 1.)
		{
			Loading_state_index = 1;
			//cout << "dt_Cap=" << dt_Cap << endl;

		}
		break;
	default:
		Loading_state_index = 0;
		break;
	}

	return Loading_state_index;
}

void General_functions::Determination_of_ultimate_stresses()
{
	if (dJ2 > 0)
	{
		if (Point_a < Hydro_stress_M && Hydro_stress_M < Point_0)
		{
			/*cout<<"極限狀態在子午線上"<<endl;*/
			Ultimate_stress_location = 1;
			cout << "dt_Meridian=" << dt_Meridian << endl;
		}

		else if (abs(Hydro_stress_C - Hydro_stress_M) < 0.01)
		{
			/*cout<<"極限狀態在帽蓋與子午線交會處"<<endl;*/
			Ultimate_stress_location = 2;

		}

		else if (Point_b < Hydro_stress_C && Hydro_stress_C < Point_a)
		{
			/*cout<<"極限狀態在帽蓋上"<<endl;*/
			Ultimate_stress_location = 3;
			cout << "dt_Cap=" << dt_Cap << endl;

		}

	}

	else if (dJ2 == 0)
	{
		if (dHydro_stress < 0.0)
		{
			/*cout<<"極限狀態在帽蓋與靜水壓力軸交點上"<<endl;*/
			Ultimate_stress_location = 4;
		}

		else if (dHydro_stress > 0.0)
		{
			/*cout<<"極限狀態在子午線與靜水壓力軸交點上"<<endl;*/
			Ultimate_stress_location = 5;
		}
	}

	else
	{
		/*cout<<"無法判斷位置"<<endl;*/
		Ultimate_stress_location = 6;
	}

	switch (Ultimate_stress_location)
	{
	case 1:
		Ult_hydro_stress = Hydro_stress_M;
		Ult_devia_stress = Devia_stress_M;
		Ult_theta = Theta;
		break;

	case 2:
		Ult_hydro_stress = Hydro_stress_M;
		Ult_devia_stress = Devia_stress_M;
		Ult_theta = Theta;
		break;

	case 3:
		Ult_hydro_stress = Hydro_stress_C;
		Ult_devia_stress = Devia_stress_C;
		Ult_theta = Theta;
		break;

	case 4:
		Ult_hydro_stress = Point_b;
		Ult_devia_stress = 0.;
		Ult_theta = Theta;
		break;

	case 5:
		Ult_hydro_stress = Point_0;
		Ult_devia_stress = 0.;
		Ult_theta = Theta;
		break;

	case 6:

	default:
		break;
	}



	cout << "Ultimate_Stress_Location=" << Ultimate_stress_location << endl;
	cout << "Ult_hydro_stress=" << Ult_hydro_stress << endl;
	cout << "Ult_devia_stress=" << Ult_devia_stress << endl;
	cout << endl;

}

void General_functions::Determination_of_ultimate_strains()
{
	if (dJ2_s > 0.0)
	{
		if (Point_a_s <= Hydro_strain_M && Hydro_strain_M < Point_0_s)
		{
			//cout << "極限狀態在子午線上" << endl;
			Ultimate_strain_location = 1;
		}

		else if (abs(Hydro_strain_C - Hydro_strain_M) < 1e-10)
		{
			/*cout<<"極限狀態在帽蓋與子午線交會處"<<endl;*/
			Ultimate_strain_location = 2;
		}

		else if (Point_b_s <= Hydro_strain_C && Hydro_strain_C < Point_a_s)
		{
			/*cout<<"極限狀態在帽蓋上"<<endl;*/
			Ultimate_strain_location = 3;
		}
		else
		{
			Ultimate_strain_location = 0;
		}

	}

	else if (dJ2_s <= 0)
	{
		if (dHydro_strain < 0.)
		{
			/*cout<<"極限狀態在帽蓋與靜水壓力軸交點上"<<endl;*/
			Ultimate_strain_location = 4;
		}

		else if (dHydro_strain > 0.0)
		{
			/*cout<<"極限狀態在子午線與靜水壓力軸交點上"<<endl;*/
			Ultimate_strain_location = 5;
		}
		else
		{
			Ultimate_strain_location = 0;
		}

	}

	switch (Ultimate_strain_location)
	{
	case 1:
		Ult_hydro_strain = Hydro_strain_M;
		Ult_devia_strain = Devia_strain_M;
		Ult_theta = Theta_s;
		break;

	case 2:
		Ult_hydro_strain = Hydro_strain_M;
		Ult_devia_strain = Devia_strain_M;
		Ult_theta = Theta_s;
		break;

	case 3:
		Ult_hydro_strain = Hydro_strain_C;
		Ult_devia_strain = Devia_strain_C;
		Ult_theta = Theta_s;
		break;

	case 4:
		Ult_hydro_strain = Point_b_s;
		Ult_devia_strain = 0.;
		Ult_theta = Theta_s;
		break;

	case 5:
		Ult_hydro_strain = Point_0_s;
		Ult_devia_strain = 0.;
		Ult_theta = Theta_s;
		break;

	case 6:

	default:
		break;
	}
	/*
		cout << "Ultimate_Strain_Location=" << Ultimate_Strain_Location << endl;
		cout << "Ult_hydro_strain=" << Ult_hydro_strain << endl;
		cout << "Ult_devia_strain=" << Ult_devia_strain << endl;
		cout << endl;*/
}

void General_functions::Coordinate_transform_to_principal_stresses()
{
	double a = Ult_theta / 180. *3.14159265;
	double b = (Ult_theta - 120.) / 180. *3.14159265;
	double c = (Ult_theta + 120.) / 180. *3.14159265;

	Ulti_stress_1 = Ult_hydro_stress / sqrt(3.) + sqrt(2.) / sqrt(3.) * Ult_devia_stress * cos(a);
	Ulti_stress_2 = Ult_hydro_stress / sqrt(3.) + sqrt(2.) / sqrt(3.) * Ult_devia_stress * cos(b);
	Ulti_stress_3 = Ult_hydro_stress / sqrt(3.) + sqrt(2.) / sqrt(3.) * Ult_devia_stress * cos(c);

	if (Ulti_stress_1 >= Ulti_stress_2 && Ulti_stress_2 >= Ulti_stress_3) {

		Temp_1 = Ulti_stress_1;
		Temp_2 = Ulti_stress_2;
		Temp_3 = Ulti_stress_3;

	}
	else if (Ulti_stress_1 >= Ulti_stress_3 && Ulti_stress_3 >= Ulti_stress_2) {

		Temp_1 = Ulti_stress_1;
		Temp_2 = Ulti_stress_3;
		Temp_3 = Ulti_stress_2;

	}
	else if (Ulti_stress_2 >= Ulti_stress_1 && Ulti_stress_1 >= Ulti_stress_3) {

		Temp_1 = Ulti_stress_2;
		Temp_2 = Ulti_stress_1;
		Temp_3 = Ulti_stress_3;

	}
	else if (Ulti_stress_2 >= Ulti_stress_3 && Ulti_stress_3 >= Ulti_stress_1) {

		Temp_1 = Ulti_stress_2;
		Temp_2 = Ulti_stress_3;
		Temp_3 = Ulti_stress_1;

	}
	else if (Ulti_stress_3 >= Ulti_stress_1 && Ulti_stress_1 >= Ulti_stress_2) {

		Temp_1 = Ulti_stress_3;
		Temp_2 = Ulti_stress_1;
		Temp_3 = Ulti_stress_2;

	}
	else if (Ulti_stress_3 >= Ulti_stress_2 && Ulti_stress_2 >= Ulti_stress_1) {

		Temp_1 = Ulti_stress_3;
		Temp_2 = Ulti_stress_2;
		Temp_3 = Ulti_stress_1;

	}

	if (Ulti_stress_1 >= 0 && Ulti_stress_2 >= 0 && Ulti_stress_3 >= 0) Stress_state_type = 1;
	else if (Ulti_stress_1 <= 0 && Ulti_stress_2 <= 0 && Ulti_stress_3 <= 0) Stress_state_type = 2;
	else if (Ulti_stress_1*Ulti_stress_2*Ulti_stress_3 <= 0) Stress_state_type = 3;
	else if (Ulti_stress_1*Ulti_stress_2*Ulti_stress_3 >= 0) Stress_state_type = 4;

	switch (Stress_state_type)
	{
	case 1:
		Ulti_stress_1 = Temp_1;
		Ulti_stress_2 = Temp_2;
		Ulti_stress_3 = Temp_3;
		break;

	case 2:
		Ulti_stress_1 = Temp_3;
		Ulti_stress_2 = Temp_2;
		Ulti_stress_3 = Temp_1;
		break;

	case 3:
		Ulti_stress_1 = Temp_3;
		Ulti_stress_2 = Temp_2;
		Ulti_stress_3 = Temp_1;
		break;

	case 4:
		Ulti_stress_1 = Temp_3;
		Ulti_stress_2 = Temp_2;
		Ulti_stress_3 = Temp_1;
		break;

	default:
		break;
	}
}

void General_functions::Coordinate_transform_to_principal_strains()
{
	double a = Ult_theta / 180. *3.14159265;
	double b = (Ult_theta - 120.) / 180. *3.14159265;
	double c = (Ult_theta + 120.) / 180. *3.14159265;

	Ulti_strain_1 = Ult_hydro_strain / sqrt(3.) + sqrt(2.) / sqrt(3.) * Ult_devia_strain * cos(a);
	Ulti_strain_2 = Ult_hydro_strain / sqrt(3.) + sqrt(2.) / sqrt(3.) * Ult_devia_strain * cos(b);
	Ulti_strain_3 = Ult_hydro_strain / sqrt(3.) + sqrt(2.) / sqrt(3.) * Ult_devia_strain * cos(c);

	if (Ulti_strain_1 >= Ulti_strain_2 && Ulti_strain_2 >= Ulti_strain_3) {

		Temp_1_ = Ulti_strain_1;
		Temp_2_ = Ulti_strain_2;
		Temp_3_ = Ulti_strain_3;

	}
	else if (Ulti_strain_1 >= Ulti_strain_3 && Ulti_strain_3 >= Ulti_strain_2) {

		Temp_1_ = Ulti_strain_1;
		Temp_2_ = Ulti_strain_3;
		Temp_3_ = Ulti_strain_2;

	}
	else if (Ulti_strain_2 >= Ulti_strain_1 && Ulti_strain_1 >= Ulti_strain_3) {

		Temp_1_ = Ulti_strain_2;
		Temp_2_ = Ulti_strain_1;
		Temp_3_ = Ulti_strain_3;

	}
	else if (Ulti_strain_2 >= Ulti_strain_3 && Ulti_strain_3 >= Ulti_strain_1) {

		Temp_1_ = Ulti_strain_2;
		Temp_2_ = Ulti_strain_3;
		Temp_3_ = Ulti_strain_1;

	}
	else if (Ulti_strain_3 >= Ulti_strain_1 && Ulti_strain_1 >= Ulti_strain_2) {

		Temp_1_ = Ulti_strain_3;
		Temp_2_ = Ulti_strain_1;
		Temp_3_ = Ulti_strain_2;

	}
	else if (Ulti_strain_3 >= Ulti_strain_2 && Ulti_strain_2 >= Ulti_strain_1) {

		Temp_1_ = Ulti_strain_3;
		Temp_2_ = Ulti_strain_2;
		Temp_3_ = Ulti_strain_1;

	}

	if (Ulti_strain_1 >= 0 && Ulti_strain_2 >= 0 && Ulti_strain_3 >= 0) Stress_state_type = 1;
	else if (Ulti_strain_1 <= 0 && Ulti_strain_2 <= 0 && Ulti_strain_3 <= 0) Stress_state_type = 2;
	else if (Ulti_strain_1*Ulti_strain_2*Ulti_strain_3 <= 0) Stress_state_type = 3;
	else if (Ulti_strain_1*Ulti_strain_2*Ulti_strain_3 >= 0) Stress_state_type = 4;

	switch (Stress_state_type)
	{
	case 1:
		Ulti_strain_1 = Temp_1_;
		Ulti_strain_2 = Temp_2_;
		Ulti_strain_3 = Temp_3_;
		break;

	case 2:
		Ulti_strain_1 = Temp_3_;
		Ulti_strain_2 = Temp_2_;
		Ulti_strain_3 = Temp_1_;
		break;

	case 3:
		Ulti_strain_1 = Temp_3_;
		Ulti_strain_2 = Temp_2_;
		Ulti_strain_3 = Temp_1_;
		break;

	case 4:
		Ulti_strain_1 = Temp_3_;
		Ulti_strain_2 = Temp_2_;
		Ulti_strain_3 = Temp_1_;
		break;

	default:
		break;
	}
}

void General_functions::Ultimate_stress_strain_correction()
{
	if (abs(Ulti_stress_1) < 1e-7)
	{
		Ulti_stress_1 = 0.;
		Ulti_strain_1 = 0.;
		i = 0;
	}

	if (abs(Ulti_stress_2) < 1e-7)
	{
		Ulti_stress_2 = 0.;
		Ulti_strain_2 = 0.;
		j = 0;
	}

	if (abs(Ulti_stress_3) < 1e-7)
	{
		Ulti_stress_3 = 0.;
		Ulti_strain_3 = 0.;
		k = 0;
	}

	cout << "Loading_state=" << i + j + k << endl;

}

double General_functions::General_functions_output()
{
	double dt=0;

	if (Ultimate_stress_location == 1)
	{
		dt = dt_Meridian;
		Haigh_Westergaard_Output << dt << "\t" << Hydro_stress << "\t" << Devia_stress << "\t" << Theta << "\t" << Hydro_strain << "\t" << Devia_strain << "\t" << Theta_s << endl;

	}
	else if (Ultimate_stress_location == 3)
	{
		dt = dt_Cap;
		Haigh_Westergaard_Output << dt_Cap << "\t" << Hydro_stress << "\t" << Devia_stress << "\t" << Theta << "\t" << Hydro_strain << "\t" << Devia_strain << "\t" << Theta_s << endl;

	}

	return dt;
}



/*CMW failrue surface*/

ofstream Stress_Meridian_Output("Stress_Meridian_Output.txt");
ofstream Strain_Meridian_Output("Strain_Meridian_Output.txt");
ofstream Stress_Cap_Output("Stress_Cap_Output.txt");
ofstream Strain_Cap_Output("Strain_Cap_Output.txt");
ofstream Ellipse_Function_Output("Ellipse_Function_Output.txt");

void Closed_Menetrey_Willam_failure_surface::Initiallize()
{
	Ellipse_Function_Output << "Radious_of_Ellipse" << "\t" << "Radious_of_Ellipse_s" << "\t" << "Lode_angle" << "\t" << "Lode_angle + 3.14159265 * 2 / 3" << "\t" << "Lode_angle + 3.14159265 * 4 / 3" << endl;
	Stress_Meridian_Output << "dHydro_stress" << "\t" << "Plot_devia_stress" << "\t" << "X1" << "\t" << "Y1" << "\t" << "X2" << "\t" << "Y2" << "\t" << "X3" << "\t" << "Y3" << endl;
	Stress_Cap_Output << "dHydro_stress" << "\t" << "Plot_devia_stress" << "\t" << "X1" << "\t" << "Y1" << "\t" << "X2" << "\t" << "Y2" << "\t" << "X3" << "\t" << "Y3" << endl;
	Strain_Meridian_Output << "dHydro_strain" << "\t" << "Plot_devia_strain" << "\t" << "X1" << "\t" << "Y1" << "\t" << "X2" << "\t" << "Y2" << "\t" << "X3" << "\t" << "Y3" << endl;
	Strain_Cap_Output << "dHydro_strain" << "\t" << "Plot_devia_strain" << "\t" << "X1" << "\t" << "Y1" << "\t" << "X2" << "\t" << "Y2" << "\t" << "X3" << "\t" << "Y3" << endl;

}

void Closed_Menetrey_Willam_failure_surface::Stress_failure_surface(double Compressive_stress, double Tensile_stress)
{
	Alpha = Tensile_stress / abs(Compressive_stress);

	e = (2.0 + Alpha) / (4.0 - Alpha);

	m = sqrt(3.0)* (pow(Compressive_stress, 2) - pow(Tensile_stress, 2))*e / (Compressive_stress * Tensile_stress * (e + 1));

	/*Elliptical function*/
	if (dJ2 == 0)
	{
		Radious_of_Ellipse = 0;
		//cout << "The state is on the hydrostatic axis" << endl;
	}
	else
	{

		Y = Theta / 180.0 * 3.14159265;

		H = 4.0 * (1.0 - pow(e, 2.0)) * pow(cos(Y), 2.0) + pow(2.0 * e - 1.0, 2.0);

		N = 2.0 * (1.0 - pow(e, 2.0)) * cos(Y) + (2.0 * e - 1.0) * pow(4 * (1.0 - pow(e, 2.0)) * pow(cos(Y), 2.0) + 5 * pow(e, 2.0) - 4 * e, 0.5);

		Radious_of_Ellipse = H / N;		// In this criterion

		/*cout<<"Radious_of_Ellipse= "<<Radious_of_Ellipse<<endl;*/

	}

	Point_0 = sqrt(3.0)*w3*Compressive_stress / (m*w2);
	Point_a = 1.*Compressive_stress*w4;
	Point_b = 2.*Compressive_stress*w4;

	/*Radius at point a*/
	w = w0*1.5 / pow(Compressive_stress, 2);
	s = w1*m*Radious_of_Ellipse / (sqrt(6.0)*Compressive_stress);
	x = w2*m*Point_a / (sqrt(3.0)*Compressive_stress) - w3;

	double temp1 = (-s + sqrt(pow(s, 2) - 4. * w * x)) / (2 * w);
	double temp2 = (-s - sqrt(pow(s, 2) - 4. * w * x)) / (2 * w);

	if (temp1 > temp2)
	{
		Devia_stress_a = temp1;
	}
	else Devia_stress_a = temp2;

	/*Meridian*/
	a = 1.5*w0 * pow(Devia_stress, 2) / pow(Compressive_stress, 2);
	b = w1*m * Radious_of_Ellipse * Devia_stress / (sqrt(6.0) * Compressive_stress) + w2*m*Hydro_stress / (sqrt(3.0)*Compressive_stress);
	c = -w3;

	double temp3 = (-b - sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
	double temp4 = (-b + sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);

	if (temp3 > temp4)
	{
		dt_Meridian = temp3;
	}
	else dt_Meridian = temp4;

	Hydro_stress_M = dt_Meridian*Hydro_stress;
	Devia_stress_M = dt_Meridian*Devia_stress;

	/*Cap*/
	double J = -Devia_stress_a / pow(Point_a - Point_b, 2.0);
	double d = J * pow(Hydro_stress, 2.0);
	double e = -2.0 * J * Point_a * Hydro_stress - Devia_stress;
	double f = 2.*J*Point_a*Point_b - J * pow(Point_b, 2);

	double temp5 = (-e - sqrt(pow(e, 2.) - 4. * d * f)) / (2. * d);
	double temp6 = (-e + sqrt(pow(e, 2.) - 4. * d * f)) / (2. * d);

	if (temp5 > temp6)
	{
		dt_Cap = temp5;
	}
	else dt_Cap = temp6;

	Hydro_stress_C = dt_Cap*Hydro_stress;
	Devia_stress_C = dt_Cap*Devia_stress;

	/*
		cout << "Point_0=" << Point_0 << endl;
		cout << "Point_a=" << Point_a << endl;
		cout << "Point_b=" << Point_b << endl;
		cout << "Hydro_stress_M=" << Hydro_stress_M << endl;
		cout << "Devia_stress_M=" << Devia_stress_M << endl;
		cout << "Devia_stress_a=" << Devia_stress_a << endl;
		cout << "temp1=" << temp1 << endl;
		cout << "temp2=" << temp2 << endl;
		cout << "temp3=" << temp3 << endl;
		cout << "temp4=" << temp4 << endl;
		cout << "temp5=" << temp5 << endl;
		cout << "temp6=" << temp6 << endl;
		*/

}

void Closed_Menetrey_Willam_failure_surface::Strain_failure_surface(double Compressive_strain, double Tensile_strain)
{
	Alpha_s = Tensile_strain / abs(Compressive_strain);
	//cout<<"Alpha_s= "<<Alpha_s<<endl;

	e_s = (2.0 + Alpha_s) / (4.0 - Alpha_s);
	//cout<<"e_s= "<<e_s<<endl;

	double a = sqrt(3.0)* (pow(Compressive_strain, 2) - pow(Tensile_strain, 2))* e_s;

	double b = (Compressive_strain * Tensile_strain * (e_s + 1));

	m_s = a / b;

	/*Elliptical function*/
	if (dJ2 == 0)
	{
		Radious_of_Ellipse = 0;
		//cout << "The state is on the hydrostatic axis" << endl;
	}
	else
	{

		Y_s = Theta_s / 180.0 * 3.14159265;

		H_s = 4.0 * (1.0 - pow(e_s, 2.0)) * pow(cos(Y_s), 2.0) + pow(2.0 * e_s - 1.0, 2.0);

		N_s = 2.0 * (1.0 - pow(e_s, 2.0)) * cos(Y_s) + (2.0 * e_s - 1.0) * pow(4 * (1.0 - pow(e_s, 2.0)) * pow(cos(Y_s), 2.0) + 5 * pow(e_s, 2.0) - 4 * e_s, 0.5);

		Radious_of_Ellipse_s = H_s / N_s;		// In this criterion

		/*cout<<"Radious_of_Ellipse= "<<Radious_of_Ellipse<<endl;*/

	}

	Point_0_s = sqrt(3.0)*Compressive_strain / m_s;
	Point_a_s = 1.0*Compressive_strain*w4_s;
	Point_b_s = 2.0*Compressive_strain*w4_s;

	/*Radius at point a*/
	w_s = w0_s*1.5 / pow(Compressive_strain, 2);
	s_s = w1_s*m_s*Radious_of_Ellipse_s / (sqrt(6.0)*Compressive_strain);
	x_s = w2_s*m_s*Point_a_s / (sqrt(3.0)*Compressive_strain) - w3_s;

	double temp7 = (-s_s + sqrt(pow(s_s, 2) - 4. * w_s * x_s)) / (2 * w_s);
	double temp8 = (-s_s - sqrt(pow(s_s, 2) - 4. * w_s * x_s)) / (2 * w_s);

	if (temp7 > temp8)
	{
		Devia_strain_a = temp7;
	}
	else Devia_strain_a = temp8;

	/**/
	a_s = w0_s*1.5*pow(Devia_strain, 2) / pow(Compressive_strain, 2);
	b_s = w1_s*m_s * Radious_of_Ellipse_s * Devia_strain / (sqrt(6.0)*Compressive_strain) + w2_s*m_s*Hydro_strain / (sqrt(3.0)*Compressive_strain);
	c_s = -w3_s;

	double temp9 = (-b_s + sqrt(pow(b_s, 2) - 4 * a_s * c_s)) / (2 * a_s);
	double temp10 = (-b_s - sqrt(pow(b, 2) - 4 * a_s * c_s)) / (2 * a_s);

	if (temp9 > temp10)
	{
		dt_Meridian_s = temp9;
	}
	else dt_Meridian_s = temp10;

	Hydro_strain_M = dt_Meridian_s*Hydro_strain;
	Devia_strain_M = dt_Meridian_s*Devia_strain;

	/**/
	double J_s = -Devia_strain_a / pow(Point_a_s - Point_b_s, 2.0);
	double d_s = J_s * pow(Hydro_strain, 2.0);
	double e_s = -2.0 * J_s * Point_a_s * Hydro_strain - Devia_strain;
	double f_s = 2.*J_s * Point_a_s*Point_b_s - J_s*pow(Point_b_s, 2);

	double temp11 = (-e_s - sqrt(pow(e_s, 2) - 4 * d_s * f_s)) / (2 * d_s);
	double temp12 = (-e_s + sqrt(pow(e_s, 2) - 4 * d_s * f_s)) / (2 * d_s);

	if (temp11 > temp12)
	{
		dt_Cap_s = temp11;
	}
	else dt_Cap_s = temp12;

	Hydro_strain_C = dt_Cap_s*Hydro_strain;
	Devia_strain_C = dt_Cap_s*Devia_strain;
	/*

		cout << "Point_0_s=" << Point_0_s << endl;
		cout << "Point_a_s=" << Point_a_s << endl;
		cout << "Point_b_s=" << Point_b_s << endl;
		cout << "Hydro_stress_M=" << Hydro_strain_M << endl;
		cout << "Devia_stress_M=" << Devia_strain_M << endl;
		cout << "Devia_strain_a=" << Devia_strain_a << endl;
		cout << "temp7=" << temp7 << endl;
		cout << "temp8=" << temp8 << endl;
		cout << "temp9=" << temp9 << endl;
		cout << "temp10=" << temp10 << endl;
		cout << "temp11=" << temp11 << endl;
		cout << "temp12=" << temp12 << endl;
		*/
}

void Closed_Menetrey_Willam_failure_surface::Plot_failure_surface(double Compressive_stress, double Compressive_strain)
{
	for (double Lode_angle = 3.14159265 / 3; Lode_angle <= 3.14159265 / 3; Lode_angle += 0.05)
	{
		double H = 4.0 * (1.0 - pow(e, 2.0)) * pow(cos(Lode_angle), 2.0) + pow(2.0 * e - 1.0, 2.0);

		double N = 2.0 * (1.0 - pow(e, 2.0)) * cos(Lode_angle) + (2.0 * e - 1.0) * pow(4 * (1.0 - pow(e, 2.0)) * pow(cos(Lode_angle), 2.0) + 5 * pow(e, 2.0) - 4 * e, 0.5);

		Radious_of_Ellipse = H / N; /*plot is different*/

		double H_s = 4.0 * (1.0 - pow(e_s, 2.0)) * pow(cos(Lode_angle), 2.0) + pow(2.0 * e_s - 1.0, 2.0);

		double N_s = 2.0 * (1.0 - pow(e_s, 2.0)) * cos(Lode_angle) + (2.0 * e_s - 1.0) * pow(4 * (1.0 - pow(e_s, 2.0)) * pow(cos(Lode_angle), 2.0) + 5 * pow(e_s, 2.0) - 4 * e_s, 0.5);

		Radious_of_Ellipse_s = H_s / N_s;		// In this criterion

		Ellipse_Function_Output << Radious_of_Ellipse << "\t" << Radious_of_Ellipse_s << "\t" << Lode_angle << "\t" << Lode_angle + 3.14159265 * 2 / 3 << "\t" << Lode_angle + 3.14159265 * 4 / 3 << endl;

		/*stress surface*/
		for (dHydro_stress = Point_0; dHydro_stress >= Point_a; dHydro_stress -= 0.5)
		{
			/**/
			double a = w0*1.5 / pow(Compressive_stress, 2);
			double b = w1*m * Radious_of_Ellipse / (sqrt(6.0)*Compressive_stress);
			double c = w2*m*dHydro_stress / (sqrt(3.0)*Compressive_stress) - w3;

			double temp1 = (-b - sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
			double temp2 = (-b + sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);

			if (temp1 > temp2)
			{
				Plot_devia_stress = temp1;
			}
			else Plot_devia_stress = temp2;

			Stress_Meridian_Output << dHydro_stress << "\t" << Plot_devia_stress << "\t" << Plot_devia_stress*cos(Lode_angle) << "\t" << Plot_devia_stress*sin(Lode_angle) << "\t" << Plot_devia_stress*cos(Lode_angle + 3.14159 / 3 * 2) << "\t" << Plot_devia_stress*sin(Lode_angle + 3.14159 / 3 * 2) << "\t" << Plot_devia_stress*cos(Lode_angle + 3.14159 / 3 * 4) << "\t" << Plot_devia_stress*sin(Lode_angle + 3.14159 / 3 * 4) << endl;
		}

		/*Radius at point a*/
		double w = w0*1.5 / pow(Compressive_stress, 2);
		double s = w1*m*Radious_of_Ellipse / (sqrt(6.0)*Compressive_stress);
		double x = w2*m*Point_a / (sqrt(3.0)*Compressive_stress) - w3;

		double temp1 = (-s + sqrt(pow(s, 2) - 4. * w * x)) / (2 * w);
		double temp2 = (-s - sqrt(pow(s, 2) - 4. * w * x)) / (2 * w);

		if (temp1 > temp2)
		{
			Devia_stress_a = temp1;
		}
		else Devia_stress_a = temp2;

		/*Cap*/
		for (dHydro_stress = Point_a; dHydro_stress >= Point_b; dHydro_stress -= 0.5)
		{
			Plot_devia_stress = -Devia_stress_a / pow(Point_a - Point_b, 2)*(dHydro_stress*dHydro_stress - 2 * Point_a*(dHydro_stress - Point_b) - Point_b*Point_b);

			Stress_Cap_Output << dHydro_stress << "\t" << Plot_devia_stress << "\t" << Plot_devia_stress*cos(Lode_angle) << "\t" << Plot_devia_stress*sin(Lode_angle) << "\t" << Plot_devia_stress*cos(Lode_angle + 3.14159 / 3 * 2) << "\t" << Plot_devia_stress*sin(Lode_angle + 3.14159 / 3 * 2) << "\t" << Plot_devia_stress*cos(Lode_angle + 3.14159 / 3 * 4) << "\t" << Plot_devia_stress*sin(Lode_angle + 3.14159 / 3 * 4) << endl;
		}

		/*strain surface*/
		for (dHydro_strain = Point_0_s; dHydro_strain >= Point_a_s; dHydro_strain -= 0.00001)
		{
			a_s = w0_s*1.5*pow(Devia_strain, 2) / pow(Compressive_strain, 2);
			b_s = w1_s*m_s * Radious_of_Ellipse_s * dDevia_strain / (sqrt(6.0)*Compressive_strain);
			c_s = w2_s*m_s*dHydro_strain / (sqrt(3.0)*Compressive_strain) - w3_s;

			double temp3 = (-b_s + sqrt(pow(b_s, 2) - 4 * a_s * c_s)) / (2 * a_s);
			double temp4 = (-b_s - sqrt(pow(b, 2) - 4 * a_s * c_s)) / (2 * a_s);

			if (temp3 > temp4)
			{
				Plot_devia_strain = temp3;
			}
			else Plot_devia_strain = temp4;

			Strain_Meridian_Output << dHydro_strain << "\t" << Plot_devia_strain << "\t" << Plot_devia_strain*cos(Lode_angle) << "\t" << Plot_devia_strain*sin(Lode_angle) << "\t" << Plot_devia_strain*cos(Lode_angle + 3.14159 / 3 * 2) << "\t" << Plot_devia_strain*sin(Lode_angle + 3.14159 / 3 * 2) << "\t" << Plot_devia_strain*cos(Lode_angle + 3.14159 / 3 * 4) << "\t" << Plot_devia_strain*sin(Lode_angle + 3.14159 / 3 * 4) << endl;
		}

		for (dHydro_strain = Point_a_s; dHydro_strain >= Point_b_s; dHydro_strain -= 0.00001)
		{
			Plot_devia_strain = -Devia_strain_a / pow(Point_a_s - Point_b_s, 2)*(dHydro_strain*dHydro_strain - 2 * Point_a_s*(dHydro_strain - Point_b_s) - Point_b_s*Point_b_s);

			Strain_Cap_Output << dHydro_strain << "\t" << Plot_devia_strain << "\t" << Plot_devia_strain*cos(Lode_angle) << "\t" << Plot_devia_strain*sin(Lode_angle) << "\t" << Plot_devia_strain*cos(Lode_angle + 3.14159 / 3 * 2) << "\t" << Plot_devia_strain*sin(Lode_angle + 3.14159 / 3 * 2) << "\t" << Plot_devia_strain*cos(Lode_angle + 3.14159 / 3 * 4) << "\t" << Plot_devia_strain*sin(Lode_angle + 3.14159 / 3 * 4) << endl;
		}
	}

}

void Closed_Menetrey_Willam_failure_surface::Stress_softening_surface(double Compressive_stress, double Tensile_stress)
{
}

void Closed_Menetrey_Willam_failure_surface::Strain_softening_surface(double Compressive_strain, double Tensile_strain)
{
}

/**/

void Uniaxial_stress_strain_model::Initiallize()
{
	Previous_equiv_strain_1 = Previous_equiv_strain_2 = Previous_equiv_strain_3 = 0.;
	Current_equiv_strain_1 = Current_equiv_strain_2 = Current_equiv_strain_3 = 0.;

	Uniaxial_Model_Output << "Equiv_Stress_1" << "\t" << "Equiv_Stress_2" << "\t" << "Equiv_Stress_3" << "\t" << "Equiv_strain_1" << "\t" << "Equiv_strain_2" << "\t" << "Equiv_strain_3" << endl;
}

void Uniaxial_stress_strain_model::Saenz_cubic_function(int Modulus_type, double dEquiv_strain_1, double dEquiv_strain_2, double dEquiv_strain_3, double Ini_elastic_1, double Ini_elastic_2, double Ini_elastic_3, double Ulti_stress_1, double Ulti_stress_2, double Ulti_stress_3, double Ulti_strain_1, double Ulti_strain_2, double Ulti_strain_3)
{
	Equiv_strain_1 += dEquiv_strain_1;
	Equiv_strain_2 += dEquiv_strain_2;
	Equiv_strain_3 += dEquiv_strain_3;

	Current_equiv_strain_1 = Equiv_strain_1;
	Current_equiv_strain_2 = Equiv_strain_2;
	Current_equiv_strain_3 = Equiv_strain_3;

	double Strain_ratio_1_ = Previous_equiv_strain_1 / Ulti_strain_1;
	double Strain_ratio_2_ = Previous_equiv_strain_2 / Ulti_strain_2;
	double Strain_ratio_3_ = Previous_equiv_strain_3 / Ulti_strain_3;

	double Strain_ratio_1 = Equiv_strain_1 / Ulti_strain_1;
	double Strain_ratio_2 = Equiv_strain_2 / Ulti_strain_2;
	double Strain_ratio_3 = Equiv_strain_3 / Ulti_strain_3;

	double Secant_modulus_1 = Ulti_stress_1 / Ulti_strain_1;
	double Secant_modulus_2 = Ulti_stress_2 / Ulti_strain_2;
	double Secant_modulus_3 = Ulti_stress_3 / Ulti_strain_3;

	double Softening_Stress_1 = Ulti_stress_1 * 0.85;
	double Softening_Stress_2 = Ulti_stress_2 * 0.85;
	double Softening_Stress_3 = Ulti_stress_3 * 0.85;

	double Softening_Strain_1 = Ulti_strain_1 * 1.41;
	double Softening_Strain_2 = Ulti_strain_2 * 1.41;
	double Softening_Strain_3 = Ulti_strain_3 * 1.41;

	double Ks_1 = Ulti_stress_1 / Softening_Stress_1;
	double Ks_2 = Ulti_stress_2 / Softening_Stress_2;
	double Ks_3 = Ulti_stress_3 / Softening_Stress_3;

	/*

		cout << "Ks_1= " << Ks_1 << endl;
		cout << "Ks_2= " << Ks_2 << endl;
		cout << "Ks_3= " << Ks_3 << endl;
	*/

	double Ke_1 = Softening_Strain_1 / Ulti_strain_1;
	double Ke_2 = Softening_Strain_2 / Ulti_strain_2;
	double Ke_3 = Softening_Strain_3 / Ulti_strain_3;

	/*
		cout << "Ke_1= " << Ke_1 << endl;
		cout << "Ke_2= " << Ke_2 << endl;
		cout << "Ke_3= " << Ke_3 << endl;
	*/

	double K1 = Ini_elastic_1*Ulti_strain_1 / Ulti_stress_1;
	double K2 = Ini_elastic_2*Ulti_strain_2 / Ulti_stress_2;
	double K3 = Ini_elastic_3*Ulti_strain_3 / Ulti_stress_3;

	/*

		cout << "K1= " << K1 << endl;
		cout << "K2= " << K2 << endl;
		cout << "K3= " << K3 << endl;
	*/


	double C1 = K1 * (Ks_1 - 1.0) / pow(Ke_1 - 1.0, 2) - 1.0 / Ke_1;
	double C2 = K2 * (Ks_2 - 1.0) / pow(Ke_2 - 1.0, 2) - 1.0 / Ke_2;
	double C3 = K3 * (Ks_3 - 1.0) / pow(Ke_3 - 1.0, 2) - 1.0 / Ke_3;


	/*
		cout << "C1= " << C1 << endl;
		cout << "C2= " << C2 << endl;
		cout << "C3= " << C3 << endl;
	*/


	double B1 = 1.0 - 2.0*C1;
	double B2 = 1.0 - 2.0*C2;
	double B3 = 1.0 - 2.0*C3;

	/*

		cout << "B1= " << B1 << endl;
		cout << "B2= " << B2 << endl;
		cout << "B3= " << B3 << endl;
	*/


	double A1 = C1 + K1 - 2.0;
	double A2 = C2 + K2 - 2.0;
	double A3 = C3 + K3 - 2.0;

	/*

		cout << "A1= " << A1 << endl;
		cout << "A2= " << A2 << endl;
		cout << "A3= " << A3 << endl;

		cout << "Strain_ratio_1= " << Strain_ratio_1 << endl;
		cout << "Strain_ratio_2= " << Strain_ratio_2 << endl;
		cout << "Strain_ratio_3= " << Strain_ratio_3 << endl;

		cout << "Equivalent_Uniaxial_strain_1= " << Equiv_strain_1 << endl;
		cout << "Equivalent_Uniaxial_strain_2= " << Equiv_strain_2 << endl;
		cout << "Equivalent_Uniaxial_strain_3= " << Equiv_strain_3 << endl;
	*/

	D1_ = (1. + A1 * Strain_ratio_1_ + B1 * pow(Strain_ratio_1_, 2) + C1 * pow(Strain_ratio_1_, 3));
	D2_ = (1. + A2 * Strain_ratio_2_ + B2 * pow(Strain_ratio_2_, 2) + C2 * pow(Strain_ratio_2_, 3));
	D3_ = (1. + A3 * Strain_ratio_3_ + B3 * pow(Strain_ratio_3_, 2) + C3 * pow(Strain_ratio_3_, 3));

	D1 = (1. + A1 * Strain_ratio_1 + B1 * pow(Strain_ratio_1, 2) + C1 * pow(Strain_ratio_1, 3));
	D2 = (1. + A2 * Strain_ratio_2 + B2 * pow(Strain_ratio_2, 2) + C2 * pow(Strain_ratio_2, 3));
	D3 = (1. + A3 * Strain_ratio_3 + B3 * pow(Strain_ratio_3, 2) + C3 * pow(Strain_ratio_3, 3));

	/*Modulus*/
	switch (Modulus_type)
	{
	case 1:		/*Tangent*/

		Previous_equiv_stress_1 = 1.*Ulti_stress_1 * K1 * Strain_ratio_1_ / D1_;
		Previous_equiv_stress_2 = 1.*Ulti_stress_2 * K2 * Strain_ratio_2_ / D2_;
		Previous_equiv_stress_3 = 1.*Ulti_stress_3 * K3 * Strain_ratio_3_ / D3_;

		Current_equiv_stress_1 = 1.*Ulti_stress_1 * K1 * Strain_ratio_1 / D1;
		Current_equiv_stress_2 = 1.*Ulti_stress_2 * K2 * Strain_ratio_2 / D2;
		Current_equiv_stress_3 = 1.*Ulti_stress_3 * K3 * Strain_ratio_3 / D3;

		Elas_modu_1 = (Current_equiv_stress_1 - Previous_equiv_stress_1) / dEquiv_strain_1;
		Elas_modu_2 = (Current_equiv_stress_2 - Previous_equiv_stress_2) / dEquiv_strain_2;
		Elas_modu_3 = (Current_equiv_stress_3 - Previous_equiv_stress_3) / dEquiv_strain_3;

		Uniaxial_Model_Output << A1 << "\t" << B1 << "\t" << C1 << "\t" << Strain_ratio_1 << "\t" << D1 << "\t" << Equiv_strain_3 << endl;

		break;

	case 2:		/*Secant*/

		Equiv_stress_1 = 1.*Ulti_stress_1 * K1 * Strain_ratio_1 / D1;
		Equiv_stress_2 = 1.*Ulti_stress_2 * K2 * Strain_ratio_2 / D2;
		Equiv_stress_3 = 1.*Ulti_stress_3 * K3 * Strain_ratio_3 / D3;

		Elas_modu_1 = Equiv_stress_1 / Equiv_strain_1;
		Elas_modu_2 = Equiv_stress_2 / Equiv_strain_2;
		Elas_modu_3 = Equiv_stress_3 / Equiv_strain_3;

		Uniaxial_Model_Output << Equiv_stress_1 << "\t" << Equiv_stress_2 << "\t" << Equiv_stress_3 << "\t" << Equiv_strain_1 << "\t" << Equiv_strain_2 << "\t" << Equiv_strain_3 << endl;

		break;

	default:
		break;
	}


	cout << "Ulti_stress_1=" << Ulti_stress_1 << endl;
	cout << "K1=" << K1 << endl;
	cout << "Equiv_strain_1=" << Equiv_strain_1 << endl;
	cout << "Ulti_strain_1=" << Ulti_strain_1 << endl;
	cout << "Strain_ratio_1=" << Strain_ratio_1 << endl;
	cout << "D1=" << D1 << endl;
	// 	cout << "Previous_equiv_stress_1=" << Previous_equiv_stress_1 << endl;

	if (i == 0)	Elas_modu_1 = Ini_elastic_1;
	if (j == 0)	Elas_modu_2 = Ini_elastic_2;
	if (k == 0)	Elas_modu_3 = Ini_elastic_3;

	Elastic_Modulus_1 = Elas_modu_1;
	Elastic_Modulus_2 = Elas_modu_2;
	Elastic_Modulus_3 = Elas_modu_3;

	Previous_equiv_strain_1 = Equiv_strain_1;
	Previous_equiv_strain_2 = Equiv_strain_2;
	Previous_equiv_strain_3 = Equiv_strain_3;


}

void Uniaxial_stress_strain_model::Popovics_function(int Modulus_type, double dEquiv_strain_1, double dEquiv_strain_2, double dEquiv_strain_3, double Ini_elastic_1, double Ini_elastic_2, double Ini_elastic_3, double Ulti_stress_1, double Ulti_stress_2, double Ulti_stress_3, double Ulti_strain_1, double Ulti_strain_2, double Ulti_strain_3)
{
	Equiv_strain_1 += dEquiv_strain_1;
	Equiv_strain_2 += dEquiv_strain_2;
	Equiv_strain_3 += dEquiv_strain_3;

	double K1 = Ini_elastic_1*Ulti_strain_1 / Ulti_stress_1;
	double K2 = Ini_elastic_2*Ulti_strain_2 / Ulti_stress_2;
	double K3 = Ini_elastic_3*Ulti_strain_3 / Ulti_stress_3;

	double Strain_ratio_1 = Equiv_strain_1 / Ulti_strain_1;
	double Strain_ratio_2 = Equiv_strain_2 / Ulti_strain_2;
	double Strain_ratio_3 = Equiv_strain_3 / Ulti_strain_3;

	double Strain_ratio_1_ = Previous_equiv_strain_1 / Ulti_strain_1;
	double Strain_ratio_2_ = Previous_equiv_strain_2 / Ulti_strain_2;
	double Strain_ratio_3_ = Previous_equiv_strain_3 / Ulti_strain_3;

	double D1 = 1. + (K1 - 1.)*pow(Strain_ratio_1, (K1 / (K1 - 1.)));
	double D2 = 1. + (K2 - 1.)*pow(Strain_ratio_2, (K2 / (K2 - 1.)));
	double D3 = 1. + (K3 - 1.)*pow(Strain_ratio_3, (K3 / (K3 - 1.)));

	double D1_ = 1. + (K1 - 1.)*pow(Strain_ratio_1_, (K1 / (K1 - 1.)));
	double D2_ = 1. + (K2 - 1.)*pow(Strain_ratio_2_, (K2 / (K2 - 1.)));
	double D3_ = 1. + (K3 - 1.)*pow(Strain_ratio_3_, (K3 / (K3 - 1.)));

	Previous_equiv_stress_1 = 1.*Ulti_stress_1 * K1 * Strain_ratio_1_ / D1_;
	Previous_equiv_stress_2 = 1.*Ulti_stress_2 * K2 * Strain_ratio_2_ / D2_;
	Previous_equiv_stress_3 = 1.*Ulti_stress_3 * K3 * Strain_ratio_3_ / D3_;

	Equiv_stress_1 = 1.*Ulti_stress_1 * K1 * Strain_ratio_1 / D1;
	Equiv_stress_2 = 1.*Ulti_stress_2 * K2 * Strain_ratio_2 / D2;
	Equiv_stress_3 = 1.*Ulti_stress_3 * K3 * Strain_ratio_3 / D3;



	switch (Modulus_type)
	{
	case 1:
		Previous_equiv_stress_1 = 1.*Ulti_stress_1 * K1 * Strain_ratio_1_ / D1_;
		Previous_equiv_stress_2 = 1.*Ulti_stress_2 * K2 * Strain_ratio_2_ / D2_;
		Previous_equiv_stress_3 = 1.*Ulti_stress_3 * K3 * Strain_ratio_3_ / D3_;

		Current_equiv_stress_1 = 1.*Ulti_stress_1 * K1 * Strain_ratio_1 / D1;
		Current_equiv_stress_2 = 1.*Ulti_stress_2 * K2 * Strain_ratio_2 / D2;
		Current_equiv_stress_3 = 1.*Ulti_stress_3 * K3 * Strain_ratio_3 / D3;

		Elas_modu_1 = (Current_equiv_stress_1 - Previous_equiv_stress_1) / dEquiv_strain_1;
		Elas_modu_2 = (Current_equiv_stress_2 - Previous_equiv_stress_2) / dEquiv_strain_2;
		Elas_modu_3 = (Current_equiv_stress_3 - Previous_equiv_stress_3) / dEquiv_strain_3;
		break;

	case 2:
		Equiv_stress_1 = 1.*Ulti_stress_1 * K1 * Strain_ratio_1 / D1;
		Equiv_stress_2 = 1.*Ulti_stress_2 * K2 * Strain_ratio_2 / D2;
		Equiv_stress_3 = 1.*Ulti_stress_3 * K3 * Strain_ratio_3 / D3;

		Elas_modu_1 = Equiv_stress_1 / Equiv_strain_1;
		Elas_modu_2 = Equiv_stress_2 / Equiv_strain_2;
		Elas_modu_3 = Equiv_stress_3 / Equiv_strain_3;
		break;

	default:
		break;
	}

	Previous_equiv_strain_1 = Equiv_strain_1;
	Previous_equiv_strain_2 = Equiv_strain_2;
	Previous_equiv_strain_3 = Equiv_strain_3;
}

void Uniaxial_stress_strain_model::Kupfer_function()
{

	Eqiv_Poisson_1 = 0.3 * D1;
	Eqiv_Poisson_2 = 0.3 * D2;
	Eqiv_Poisson_3 = 0.3 * D3;

	Poisson_ratio_12 = sqrt(Eqiv_Poisson_1*Eqiv_Poisson_2*abs(Elas_modu_1 / Elas_modu_2));
	Poisson_ratio_13 = sqrt(Eqiv_Poisson_1*Eqiv_Poisson_3*abs(Elas_modu_1 / Elas_modu_3));
	Poisson_ratio_21 = sqrt(Eqiv_Poisson_2*Eqiv_Poisson_1*abs(Elas_modu_2 / Elas_modu_1));
	Poisson_ratio_23 = sqrt(Eqiv_Poisson_2*Eqiv_Poisson_3*abs(Elas_modu_2 / Elas_modu_3));
	Poisson_ratio_31 = sqrt(Eqiv_Poisson_3*Eqiv_Poisson_1*abs(Elas_modu_3 / Elas_modu_1));
	Poisson_ratio_32 = sqrt(Eqiv_Poisson_3*Eqiv_Poisson_2*abs(Elas_modu_3 / Elas_modu_2));

	if (i == 0)	Poisson_ratio_21 = Poisson_ratio_31 = 0.3;
	if (j == 0)	Poisson_ratio_12 = Poisson_ratio_32 = 0.3;
	if (k == 0)	Poisson_ratio_13 = Poisson_ratio_23 = 0.3;

}



/*
void Uniaxial_stress_strain_model::Balan_function()
{
	Kv_1 = Initial_Elastic_Modulus*Ultimate_strain_1 / (2 * Ultimate_stress_1);
	Kv_2 = Initial_Elastic_Modulus*Ultimate_strain_2 / (2 * Ultimate_stress_2);
	Kv_3 = Initial_Elastic_Modulus*Ultimate_strain_3 / (2 * Ultimate_stress_3);

	/ *
	b1=B1;
	b2=B2;
	b3=B3;

	c1=C1;
	c2=C2;
	c3=C3;

	a1=c1+0.5*Initial_Poisson_Ratio;
	a2=c2+0.5*Initial_Poisson_Ratio;
	a3=c3+0.5*Initial_Poisson_Ratio;


	d1 = (1.0 + Kv_1*(A1 * StrainRatio_1 +B1 * pow(StrainRatio_1, 2) +C1 * pow(StrainRatio_1, 3)));
	d2 = (1.0 + Kv_2*(A2 * StrainRatio_2 +B2 * pow(StrainRatio_2, 2) +C2 * pow(StrainRatio_2, 3)));
	d3 = (1.0 + Kv_3*(A3 * StrainRatio_3 +B3 * pow(StrainRatio_3, 2) +C3 * pow(StrainRatio_3, 3)));* /




}*/

