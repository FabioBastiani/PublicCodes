#include "udf.h"
#include "mem.h"
#include "prop.h"
#include "math.h"
#include "sg_mem.h"
#include "materials.h"

/* Constant declaration */
    real C_b = 1.0;
    real pi = 3.14159265359;
    real rho_s = 1900.0;        /* kg/m^3 */
    real C_a = 54.0;            /* 1/s */
    real N_a = 6.022137e+26;    /* 1/kmol */
    real T_a = 21000.0;         /* K */
    real M_F = 26.03788;               /* kg/kmol */
    real R = 8314.46261815324; /* J/kmol K */
    real M_p = 144.0;           /* kg/kmol */
    real C_y = 11700.0;         /* kg m / kmol s */ 
    real T_y = 12100.0;         /* K */
/* Defining soot nuclei source */
DEFINE_SOURCE(nuclei_source,c,t,dS,eqn)
{
    real Y_s, b, N;

    /* Declaration of needed variables */
    real N_norm = 1.0e+15;
    b = C_UDSI(c, t, 0);
    Y_s = C_UDSI(c, t, 1);
    /* Calculating particle number density */
    N = C_R(c, t) * b * N_norm;

    /* Preventing negative values */
    if (b < 0.0) {b = 0.0;}
    if (Y_s < 0.0) {Y_s = 0.0;}
    if (N < 0.0) {N = 0.0;}

    /* Calculating nucleation of particles rate */
    real rho = C_R(c, t);       /* kg/m^3 */
    real T = C_T(c, t);         /* K */
    real P = C_P(c, t); ; /* Pressure */
    real X_F;                   /* Molar fraction of fuel */
    real nuc_rate;              /* 1/m^3/s */

    /* Calculate molar fraction of fuel */
    X_F = C_YI(c, t, 22) * rho * R * T / (P * M_F);

    /* Calculating nucleation rate */
    nuc_rate = C_a * N_a * (X_F * P / (R * T)) * exp(- 1.0 * T_a / T);
    /* Saving nucleation rate to memory 0 */
    C_UDMI(c, t, 0) = nuc_rate / N_norm;

    /* Calculating coagulation rate */
    real dp;
    real M;

    /* Mass concentration of soot */
    M = Y_s * rho;

    /* Calculating particle mean diameter */
    if (Y_s > 0 && N > 0) {
        dp = pow((6.0 * M / (pi * N * rho_s)), (1.0 / 3.0));
    } else {
        dp = 0.0;
    }
    C_UDMI(c,t,10) = dp * 1.0e+09; /* Particle mean diameter in nm */

    /* Calculating coagulation rate */
    real coag_rate;
    coag_rate = C_b * pow((24.0 * R * T / (rho_s * N_a)), 1.0 / 2.0) * pow(6 / ( pi * rho_s ), 1.0 / 6.0) * pow(M, 1.0 / 6.0) * pow(N, 11.0 / 6.0);

    /* Saving coagulation rate to memory 1 */
    C_UDMI(c, t, 1) = coag_rate / N_norm;

    /* Returning source term */
    real dNdt;
    dNdt = (nuc_rate - coag_rate) / N_norm;

    /* Saving source term to memory 2 */
    C_UDMI(c, t, 2) = dNdt;
    dS[eqn] = 0.0;
    return dNdt;
}

/* Defining soot mass fraction source term */
DEFINE_SOURCE(massfraction_source,c,t,dS,eqn)
{
    real Y_s, b, N;

    /* Declaration of needed variables */
    real N_norm = 1.0e+15;
    b = C_UDSI(c, t, 0);
    Y_s = C_UDSI(c, t, 1);
    /* Calculating particle number density */
    N = C_R(c, t) * b * N_norm;

    /* Preventing negative values */
    if (b < 0.0) {b = 0.0;}
    if (Y_s < 0.0) {Y_s = 0.0;}
    if (N < 0.0) {N = 0.0;}

    /* Calculating nucleation of particles rate */
    real P = C_P(c,t);          /* Pressure */
    real rho = C_R(c, t);       /* kg/m^3 */
    real T = C_T(c, t);         /* K */
    real X_F;                   /* Molar fraction of fuel */

    /* Calculate the molar fraction of fuel */
    X_F = C_YI(c, t, 22) * rho * R * T / (P * M_F);

    /* Calculating nucleation rate */
    real nuc_rate;              /* 1/m^3/s */
    nuc_rate = M_p * C_a * (X_F * P / (R * T)) * exp(- 1.0 * T_a / T);
    /* Saving nucleation rate of soot to memory 3 */
    C_UDMI(c, t, 3) = nuc_rate;

    /* Calculating Surface Growth */
    real M;                 /* Mass concentration of soot */
    M = Y_s * rho;

    /* Calculating surface growth rate */
    real surf_growth;

    if (Y_s > 0.0 && N > 0.0) {
        surf_growth = C_y * (X_F * P / (R * T)) * exp(- 1.0 * T_y / T) * pow(N * pi, 1.0 / 3.0) * pow(6 * M / rho_s, 2.0 / 3.0);
    }
    else {
        surf_growth = 0.0;
    }
    /* Saving surface growth rate to memory 4 */
    C_UDMI(c, t, 4) = surf_growth;
    
    /* Calculating Oxidation */
    real C_w = 105.8125;             /* kg m / kmol K^0.5 s */  
    real T_O = 19778.0;              /* K */
    real C_wo = 8903.51;            /* kg m / kmol K^0.5 s */
    real M_OH = 17.0;                /* kg/kmol */
    real M_O = 32.0;                 /* kg/kmol */
    real n_coll;
    real C_oxid;
    n_coll = 0.04;
    C_oxid = 1.0;
    /* Calculating oxidation rate for OH */
    real oxid_rate_oh;
    oxid_rate_oh = C_oxid * C_w * n_coll * (C_YI(c, t, 4) * rho / M_OH) * pow(T, 1.0 / 2.0) * pow(pi * N, 1.0 / 3.0) * pow(6 * M / rho_s, 2.0 / 3.0); 
    /* Calculating oxidation rate for O2 */
    real oxid_rate_o2;
    oxid_rate_o2 = C_oxid * C_wo * (C_YI(c, t, 3) * rho / M_O) * exp(- T_O / T) * pow(T, 1.0 / 2.0) * pow(pi * N, 1.0 / 3.0) * pow(6 * M / rho_s, 2.0 / 3.0); 
    /* When not considering Lee Oxidation model (O2) */
    oxid_rate_o2 = 0.0;
    
    real oxid_rate;
    oxid_rate = oxid_rate_oh + oxid_rate_o2;
    /* Saving oxidation rate to memory 5 */
    C_UDMI(c, t, 5) = oxid_rate;

    /* Calculating source term */
    real dYsdt;
    dYsdt = (nuc_rate + surf_growth - oxid_rate);

    /* Saving source term to memory 6 */
    C_UDMI(c, t, 6) = dYsdt;
    dS[eqn] = 0.0;
    
    /* Volume Fraction Calculation */
    C_UDMI(c, t, 7) = Y_s * C_R(c, t) / rho_s;
    /* Saving N in memory */
    C_UDMI(c, t, 8) = N;
    return dYsdt;
}

/* Widmann (2003) model */
DEFINE_WSGGM_ABS_COEFF(user_wsggm_abs_coeff,c,t,xi,p_t,s,soot_conc,Tcell,nb,ab_wsggm,ab_soot)
{
 real a_s;        
 real a_g;	   

 a_g = *ab_wsggm;
 a_s = 2370.0 * C_T(c,t) * C_R(c,t) * C_UDSI(c,t,1) / rho_s; 

 *ab_wsggm = a_s+a_g;
 C_UDMI(c,t,11) = a_s; 
 C_UDMI(c,t,12) = a_g; 
 C_UDMI(c,t,13) = *ab_wsggm; 
}

DEFINE_DIFFUSIVITY(scalar_diff,c,t,i)
{
 real D;
 D = 1.0e-5 + C_MU_T(c,t) / 0.7;
 return D;
}

