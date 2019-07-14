/* The main header for ESP. Yes, I know, I should arrange it better :) */

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>

#define ROUND(x) floor(x+0.5)

/* Max No of species */
#define MAXSPEC 20 // only used to initialize various arrays, so can be removed at one point
/* Max No of chemical reactions */
#define MAXCHIM 10 // only used to initialize various arrays, so can be removed at one point
/* Max No of redox */
#define MAXREDOX 10 // only used to initialize various arrays, so can be removed at one point
/* Max No of space grid */
#define MAXGRID 1000 // is used heavily in the simul algorithm
/* int constant for "pass" value */
#define PASS -32760
/* double constant for "pass" value */
#define PASSF -99999.
/* Max No of points that can be stored */
#define MAXNP 40000 // only used to initialize s.Cur and s.Pot arrays, so can be removed at one point
/* Default maximum Diffusion coefficient */
#define DSTAR 0.45 // used for normalization and space grid calculation

/* Known techniques */
enum technique {mode_CV, mode_SWV, mode_CA, mode_SDC};
/* IR drop calculation: None, Resistance (RU) , Resistance+Capacitance (DL) */
enum ir_mode {ir_none, ir_ru, ir_dl};
/* Electrode type */
enum elec_type {WE_Solid, WE_HMDE, WE_DME};

struct electro {        /* Define an electrode redox reaction */
    int Ox;	/* Ox + n e- <==> Red */
    int Red;
    int n;        /* number of electron */
    int E;	/* potential of reduction */
    double Ke;	/* heterogeneous rate constant */
    double a;	/* alpha coefficient */
};
typedef struct electro ELECTRO;

struct chim {	/* Define a homogeneous chemistry reaction */
    int a;
    int b;        /* the reaction is : */
    int c;        /* a + b <==> c + d */
    int d;
    double Kf;    /* forward homogenous constant */
    double Kb;	/* backward homogeneous constant */
};
typedef struct chim CHIM;

struct var_exp {/* experimental variables (used also by simulation) */
    char CO[67];	/* Comment */
    char DR[10];	/* Date Recording */
    char TR[10];	/* Time Recording */
    int CP;	/* Condition Potential (mV) */
    double CT;	/* Contidion Time (s) */
    double ET;	/* Equilibration Time (s) */
    int SI;	/* Scan Increment (mV) */
    double SR;	/* CV Scan Rate (V/s) */
    double ST;	/* Step Time/Drop Time (s), also Time per Points (TP) in CA */
    int IP;	/* Initial Potential (mV) */
    int V1;	/* Vertex 1 Potential (mV), also E1 in CA */
    int V2;	/* Vertex 2 Potential (mV), also E2 in CA */
    double VD;	/* Vertex Delay (s), also T1 in CA */
    double T2;	/* T2 in CA (s) */
    int FP;	/* Final Potential (mV) */
    int NC;	/* Number of Cycle */
    int SC;	/* Store Cycle */
    int AM;	/* Acquisition Mode {1,2,3,4,All=0,Ramp=5} */
    int NP;	/* Number of point (max MAXNP); calculated by program */
    double FR;	/* Frequencies (Hz) for SWV */
    int PH;	/* Pulse Height (mV) for SWV */
    enum technique Mode;	/* mode_CV, mode_SWV, mode_CA, mode_SDC */
    enum ir_mode IR;  /* {ir_none, ir_dl, ir_ru} -> {none, Double layer, RU */
    double DL;	/* Double Layer Capacitance (F) */
    double RU;	/* uncompensated resistance m� (mohm) */
    int AP;	/* Approx. Chem. 1 (Fast) or 0 (Exact: Chem. in all boxes) */
    int ncycle;	/* Set half, one, one and half cycle:
    	   {1,2,3} --> {IP-FP, IP-V1-FP, IP-V1-V2-FP} */
    double TE;	/* Temperature (K) (298.15=25 �C) */
    double AR;	/* Electrode Area (cm2), also mercury flow MF (mg/s) for DME */
    enum elec_type WE;    /* Electrode type */
};
typedef struct var_exp VAR_EXP;

struct mec_exp {/* More about mechanism */
    int nspec;	/* Number of chemical species in solution */
    int nchim;	/* Number of homogeneous chemical reactions */
    int nredox;	/* Number of heterogeneous redox */
    double *C;	/* Concentrations */
    double D[MAXSPEC];		/* Diffusion coefficients */
    ELECTRO Re[MAXREDOX];		/* Electrode reactions */
    CHIM Ch[MAXCHIM];		/* Chemical reactions */
};
typedef struct mec_exp MEC_EXP;

struct var_sim {	/* Simulation variables */
    double FRT;			/* store F/RT */
    double SkipChim;
    int Pot[MAXNP];		/* Simulated Potential (mV) */
    double Cur[MAXNP];		/* Simulated Current */
    double FluxJ[MAXSPEC];	/* The fluxs at electrode surface */
    double Deltat;		/* Time Grid */
    double Deltax;		/* Space Grid */
    int NS;			/* Number of Space Grid */
    int NI;	/* Number of time-division of SI. It set the approximation.
    	   Need to be a multiple of 4 (due to AM=Ramp) */
    int NIDIV4;		/* NI/4 */
    int NIDIV2;		/* NI/2 */
    int c_idx[MAXSPEC][MAXCHIM];	/* The chemical index (for homogeneous
    			   chemical reactions) */
};
typedef struct var_sim VAR_SIM;


// interfaced functions:
uint64_t get_posix_clock_time(void);
int _setup(void);
int _destroy(void);
int _set_params(int, double, double, int, int, int, double, int, int, double, int, int, int, int, int, int, double, int, int, double, double, int, double, double, int);
int _simulate(void);
int _add_species(double, double);
int _add_redox(int, int, int, double, double, double);
int _add_chemical(int, int, int, int, double, double);
// internal functions:
void Do_Simul(MEC_EXP *);
void Scan(int, int, int *, MEC_EXP *);
void Aspetta(double, int, MEC_EXP *);
void Step_CV_All(int, int *, int *, MEC_EXP *);
void Step_CV_1(int, int *, int *, MEC_EXP *);
void Step_CV_2(int, int *, int *, MEC_EXP *);
void Step_CV_3(int, int *, int *, MEC_EXP *);
void Step_CV_4(int, int *, int *, MEC_EXP *);
void Step_SDC(int, int *, int *, MEC_EXP *);
void Step_CV_Ramp(int, int *, int *, MEC_EXP *);
void Step_SWV(int, int *, int *, MEC_EXP *);
double RDCstep(int, int, MEC_EXP *);
double RDCstep_reduced(int, int , MEC_EXP *);
int correct_IR(int);
double Redox(double, MEC_EXP *);
void Chimica(MEC_EXP *);
void Diffusione(MEC_EXP *);
void Set_Simul(MEC_EXP *);
void Calcola_NI(MEC_EXP *);
int CalcolaNP(void);
