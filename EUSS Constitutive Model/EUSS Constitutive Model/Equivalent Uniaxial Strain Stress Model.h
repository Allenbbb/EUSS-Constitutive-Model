class General_functions
{
public:
	void Initialize();
	void Stress_invariants_to_HW(double dStress_1, double dStress_2, double dStress_3);
	void Strain_invariants_to_HW(double dEquiv_strain_1, double dEquiv_strain_2, double dEquiv_strain_3);

	void Modification_of_failure_surface(int Criterion_type);

	int Current_loading_state();
	
	void Determination_of_ultimate_stresses();
	double Get_ult_hydro_stress() { return Ult_hydro_stress; }
	double Get_ult_devia_stress() { return Ult_devia_stress; }
	double Get_ult_theta() { return Ult_theta; }

	void Determination_of_ultimate_strains();
	double Get_ult_hydro_strain() { return Ult_hydro_strain; }
	double Get_ult_devia_strain() { return Ult_devia_strain; }
	double Get_ult_theta_s() { return Ult_theta_s; }

	void Coordinate_transform_to_principal_stresses();
	double Get_Ulti_stress_1() { return Ulti_stress_1; }
	double Get_Ulti_stress_2() { return Ulti_stress_2; }
	double Get_Ulti_stress_3() { return Ulti_stress_3; }

	void Coordinate_transform_to_principal_strains();
	double Get_Ulti_strain_1() { return Ulti_strain_1; }
	double Get_Ulti_strain_2() { return Ulti_strain_2; }
	double Get_Ulti_strain_3() { return Ulti_strain_3; }
	
	void Ultimate_stress_strain_correction();

	double General_functions_output();
private:
	double Temp_1, Temp_2, Temp_3;
	double Temp_1_, Temp_2_, Temp_3_;
	double Ulti_stress_1, Ulti_stress_2, Ulti_stress_3;
	double Ulti_strain_1, Ulti_strain_2, Ulti_strain_3;
	double Ult_hydro_stress, Ult_devia_stress, Ult_theta, Ult_hydro_strain, Ult_devia_strain, Ult_theta_s;
	

};

class Closed_Menetrey_Willam_failure_surface
{
public:
	void Initiallize();
	void Stress_failure_surface(double Compressive_stress, double Tensile_stress);
	void Strain_failure_surface(double Compressive_strain, double Tensile_strain);
	void Plot_failure_surface(double Compressive_stress, double Compressive_strain);


	void Stress_softening_surface(double Compressive_stress, double Tensile_stress);
	void Strain_softening_surface(double Compressive_strain, double Tensile_strain);

private:

	double Alpha, e, m;
	double Alpha_s, e_s, m_s;
	double Radious_of_Ellipse;
	double Radious_of_Ellipse_s;

	double Y, H, N;
	double w, s, x;
	double a, b, c;

	double Y_s, H_s, N_s;
	double w_s, s_s, x_s;
	double a_s, b_s, c_s;

};

class Uniaxial_stress_strain_model
{
public:
	void Initiallize();
	void Saenz_cubic_function(int Modulus_type, double dEquiv_strain_1, double dEquiv_strain_2, double dEquiv_strain_3, double Ini_elastic_1, double Ini_elastic_2, double Ini_elastic_3, double Ulti_stress_1, double Ulti_stress_2, double Ulti_stress_3, double Ulti_strain_1, double Ulti_strain_2, double Ulti_strain_3);
	void Popovics_function(int Modulus_type, double dEquiv_strain_1, double dEquiv_strain_2, double dEquiv_strain_3, double Ini_elastic_1, double Ini_elastic_2, double Ini_elastic_3, double Ulti_stress_1, double Ulti_stress_2, double Ulti_stress_3, double Ulti_strain_1, double Ulti_strain_2, double Ulti_strain_3);
	double Get_Elas_modu_1() { return Elastic_Modulus_1; };
	double Get_Elas_modu_2() { return Elastic_Modulus_2; };
	double Get_Elas_modu_3() { return Elastic_Modulus_3; };

	void Kupfer_function();
	double Get_Poisson_ratio_12() { return Poisson_ratio_12; };
	double Get_Poisson_ratio_13() { return Poisson_ratio_13; };
	double Get_Poisson_ratio_21() { return Poisson_ratio_21; };
	double Get_Poisson_ratio_23() { return Poisson_ratio_23; };
	double Get_Poisson_ratio_31() { return Poisson_ratio_31; };
	double Get_Poisson_ratio_32() { return Poisson_ratio_32; };

	void Balan_function();
		
private:

	double D1, D2, D3;
	double D1_, D2_, D3_;

	double Equiv_stress_1, Equiv_stress_2, Equiv_stress_3;
	double Previous_equiv_stress_1, Previous_equiv_stress_2, Previous_equiv_stress_3;
	double Current_equiv_stress_1, Current_equiv_stress_2, Current_equiv_stress_3;
	double Previous_equiv_strain_1, Previous_equiv_strain_2, Previous_equiv_strain_3;
	double Current_equiv_strain_1, Current_equiv_strain_2, Current_equiv_strain_3;

	double Elastic_Modulus_1, Elastic_Modulus_2, Elastic_Modulus_3;

	double Eqiv_Poisson_1, Eqiv_Poisson_2, Eqiv_Poisson_3;
	double Poisson_ratio_12, Poisson_ratio_13, Poisson_ratio_21, Poisson_ratio_23, Poisson_ratio_31, Poisson_ratio_32;

};

