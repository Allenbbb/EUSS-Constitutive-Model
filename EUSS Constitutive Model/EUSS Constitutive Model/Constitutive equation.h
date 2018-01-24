class Orthotropic_constitutive_equation
{
public:
	void Initialize();
	void Incremental_Constitutive_Equation(double Elastic_modulus_1, double Elastic_modulus_2, double Elastic_modulus_3, double Poisson_ratio_12, double Poisson_ratio_13, double Poisson_ratio_21, double Poisson_ratio_23, double Poisson_ratio_31, double Poisson_ratio_32);
	void Normal_stress_to_normal_strain(double dStress_1, double dStress_2, double dStress_3);
	void Normal_strain_to_normal_stress(double dStrain_1, double dStrain_2, double dStrain_3);
	void Normal_strain_to_equiv_strain(double dStrain_1, double dStrain_2, double dStrain_3);
	void Equiv_strain_to_normal_strain(double dEqui_Uni_strain_1, double dEqui_Uni_strain_2, double dEqui_Uni_strain_3);

	/*Incremental*/
	double dStrain_1, dStrain_2, dStrain_3;
	double dEquiv_strain_1, dEquiv_strain_2, dEquiv_strain_3;
	double dStress_1, dStress_2, dStress_3;

private:
	double Omega, OmegaH;
	double Matrix_E[3][3], Matrix_v[3][3], Matrix_D[3][3], Matrix_H[3][3], Matrix_C[3][3], Matrix_L[3][3];



};






