#include <iostream>
#include <cmath>
#include <fstream>

int Material_type, Criterion_type;
int Loading_State_Index;
double Target_stress;
double dt;

/*Material properties*/
double Compressive_stress, Tensile_stress;
double Compressive_strain, Tensile_strain;
double Elastic_Modulus_1, Elastic_Modulus_2, Elastic_Modulus_3;
double Poisson_ratio_12, Poisson_ratio_13, Poisson_ratio_21, Poisson_ratio_23, Poisson_ratio_31, Poisson_ratio_32;

/*Incremental*/
double dStrain_1, dStrain_2, dStrain_3;
double dEquiv_strain_1, dEquiv_strain_2, dEquiv_strain_3;
double dStress_1, dStress_2, dStress_3;

/*Total*/
double Strain_1, Strain_2, Strain_3;
double Stress_1, Stress_2, Stress_3;

/**/
double Hydrostatic_Stress, Deviatoric_Stress, Lode_angle;
double Hydrostatic_Strain, Deviatoric_Strain, Lode_angle_s;

/**/
double Ultimate_Stress_1, Ultimate_Stress_2, Ultimate_Stress_3;
double Ultimate_Strain_1, Ultimate_Strain_2, Ultimate_Strain_3;