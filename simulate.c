
#include "simulate.h"

double slope=1.f, intcp=0.f;     /* for best-fit curve */

VAR_EXP v;
MEC_EXP m;
VAR_SIM s;
double C0[MAXSPEC];

/* The variable Step_ contain the address of the function for the correct
   electrochemical technique; can be SWV or CV (wich depend by AM) */
static void (*Step_)(int, int *, int *, MEC_EXP *);
static double current_norm, Capacity, icap0=0.0;
static int oldpot;
static uint64_t start_time, elapsed_redox, elapsed_chim, elapsed_diff;


uint64_t get_posix_clock_time(void)
{
    struct timespec ts;

    if (clock_gettime (CLOCK_MONOTONIC, &ts) == 0)
        return (uint64_t) (ts.tv_sec * 1000000000 + ts.tv_nsec);
    else
        return 0;
}

int _setup(void)
{
	// ALL MECHANISTIC SETTINGS
	// allocate memory for concentration matrix:
	m.C= (double *) malloc(sizeof(double) * MAXSPEC * MAXGRID);
	if (m.C==NULL) return(-1);
	memset(m.C, 0, sizeof(double)*MAXSPEC*MAXGRID); // init to 0
	// start with nothing happening:
	m.nspec = m.nredox = m.nchim = 0;

	// ALL EXPERIMENTAL VARIABLES
	v.CP = 0;		/* Condition Potential (mV), only used if v.CT > 0 */
	v.CT = -1.0;	/* Condition Time (s) */
	v.ET = -1.0;	/* Equilibration Time (s) */
	v.IP = 500;		/* Initial Potential (mV) */
	v.Mode = mode_CV;  /* mode_CV, mode_SWV, mode_CA, mode_SDC */
	v.SI = 2;		/* Scan Increment (mV) */
	v.ST = 0.01;	/* Step Time/Drop Time (s), also Time per Points (TP) in CA */
	v.V1 = -500;	/* Vertex 1 Potential (mV), also E1 in CA */
	v.V2 = 0;		/* Vertex 2 Potential (mV), also E2 in CA */
	v.VD = -1.0;	/* Vertex Delay (s), also T1 in CA */
	v.FP = 500;		/* Final Potential (mV) */
	v.NC = 1;		/* Number of Cycle */
	v.SC = 1;		/* Store Cycle */
	v.AM = 4;		/* Acquisition Mode {1,2,3,4,All=0,Ramp=5} */
	v.ncycle = 2;	/* Set half, one, one and half cycle:
			 		   {1,2,3} --> {IP-FP, IP-V1-FP, IP-V1-V2-FP} */
	v.T2 = 0;		/* T2 in CA (s) */
	v.FR = 60.0;	/* Frequencies (Hz) for SWV */
	v.PH = 30;		/* Pulse Height (mV) for SWV */
	v.IR = ir_none; /* {ir_none, ir_dl, ir_ru} -> {none, Double layer, RU */
	v.DL = 0.0;		/* Double Layer Capacitance (F) */
	v.RU = 0.0;		/* uncompensated resistance mohm */
	v.AP = 1;		/* Approx. Chem. 1 (Fast) or 0 (Exact: Chem. in all boxes) */
	v.TE = 298.15f;	/* Temperature (K) (298.15=25 C) */
	v.AR = 0.01;	/* Electrode Area (cm2), also mercury flow MF (mg/s) for DME */
	v.WE = WE_Solid;    /* Electrode type */

	v.NP = CalcolaNP();	/* Number of point (max MAXNP); calculated by program */

	// ALL SIMULATION SETTINGS:
	// (all other s.params are set in simul.c)
	s.NI = 1; /* a value lower than 4 enable the automatic setting of NI */

	// initialize all elapsed times:
	start_time = get_posix_clock_time();
	elapsed_redox = start_time;
	elapsed_chim = start_time;
	elapsed_diff = start_time;

	return(0);
}

int _set_params(int CP, double CT, double ET, int IP, int Mode, int SI, double ST, int V1, int V2,
		   double VD, int FP, int NC, int SC, int AM, int ncycle, int T2, double FR, int PH,
		   int IR, double DL, double RU, int AP, double TE, double AR, int WE)
{
	v.CP=CP; v.CT=CT; v.ET=ET; v.IP=IP; v.Mode=Mode; v.SI=SI; v.ST=ST; v.V1=V1; v.V2=V2; v.VD=VD;
	v.FP=FP; v.NC=NC; v.SC=SC; v.AM=AM; v.ncycle=ncycle; v.T2=T2; v.FR=FR; v.PH=PH; v.IR=IR;
	v.DL=DL; v.RU=RU; v.AP=AP; v.TE=TE; v.AR=AR; v.WE=WE;

	v.SR=(double)v.SI/(v.ST*1000.);
	v.NP = CalcolaNP();

	return(0);
}

int _destroy(void)
{
	free(m.C);

	return(0);
}

int _simulate(void)
{
	MEC_EXP* m2;

	m2 = (MEC_EXP*) malloc(sizeof(MEC_EXP));
	if (m2 == NULL) return(-1);
	*m2 = m;
	Do_Simul(m2);
	free(m2);

	printf("Redox time elapsed: %lf sec\n", (double) (elapsed_redox - start_time)/(double) 1e9);
	printf("Chimica time elapsed: %lf sec\n", (double) (elapsed_chim - start_time)/(double) 1e9);
	printf("Diffusion time elapsed: %lf sec\n", (double) (elapsed_diff - start_time)/(double) 1e9);

    return(v.NP);
}

int _add_species(double conc, double diff_const)
{
	int new_idx;

	if (m.nspec < MAXSPEC) {
		new_idx = m.nspec++;
		C0[new_idx] = conc;
		m.D[new_idx] = diff_const;
		return new_idx+1;
	} else {
		// cannot add more than MAXSPEC species to mechanism
		return(0);
	}
}

int _add_redox(int Ox, int Red, int n, double E, double Ke, double a)
{
	int new_idx;

	if (m.nredox < MAXREDOX) {
		new_idx = m.nredox++;
		m.Re[new_idx].Ox = Ox;
		m.Re[new_idx].Red = Red;
		m.Re[new_idx].n = n;
		m.Re[new_idx].E = E;
		m.Re[new_idx].Ke = Ke;
		m.Re[new_idx].a = a;
		return(new_idx+1);
	} else {
		// cannot add more than MAXREDOX redox steps to mechanism
		return(0);
	}
}

int _add_chemical(int a, int b, int c, int d, double Kf, double Kb)
{
	int new_idx;

	if (m.nchim < MAXCHIM) {
		new_idx = m.nchim++;
		m.Ch[new_idx].a = a;
		m.Ch[new_idx].b = b;
		m.Ch[new_idx].c = c;
		m.Ch[new_idx].d = d;
		m.Ch[new_idx].Kf = Kf;
		m.Ch[new_idx].Kb = Kb;
		return(new_idx+1);
	} else {
		// cannot add more than MAXCHIM chemical reactions to mechanism
		return(0);
	}
}

void Do_Simul(MEC_EXP *m)
{

    int i = 0, oldi = 0; // i is the current data point
	int cycle_num; // cycle_num is the current cycle
	int k, np1;
    double tmp[4];

    Set_Simul(m); /* warning: m is modified */

    if (v.CT>0.0) Aspetta(v.CT, v.CP, m);
    if (v.ET>0.0) Aspetta(v.ET, v.IP, m);

	// run as many cycles as necessary (e.g. cyclic voltammetry total cycles)
	for (cycle_num = 1; cycle_num <= v.NC; cycle_num++) {
        /* Select Step_ and take the first point according the technique selected */
        switch (v.Mode) {
            case mode_SWV:
                Step_ = Step_SWV;


                s.Pot[i] = v.IP;
				oldpot = v.IP;
                s.Cur[i++] = RDCstep(s.NI, v.IP, m);

				// set correct sign of pulse height:
                if ( ((v.ncycle == 1) && (v.IP < v.FP)) || ((v.ncycle != 1) && (v.IP < v.V1)) ) v.PH = -v.PH;

                s.Pot[i] = v.IP - v.PH;
                s.Cur[i] = RDCstep(s.NIDIV2, s.Pot[i], m);
				i++;

                s.Pot[i] = v.IP + v.PH;
                s.Cur[i] = RDCstep(s.NIDIV2, s.Pot[i], m);
				i++;

                v.PH = abs(v.PH);
                break;

            case mode_SDC:
                Step_ = Step_SDC;
                s.Pot[i] = v.IP;
                s.Cur[i++] = RDCstep(s.NI, v.IP, m);
                break;

            case mode_CV:
                for (k=0; k<4; k++) tmp[k] = RDCstep(s.NIDIV4, v.IP, m);
                s.Pot[i] = v.IP;
                switch (v.AM) {
                    case 0:
                        Step_=Step_CV_All;
                        for (k=0; k<4; k++) {
                        	s.Cur[i] = tmp[k];
                        	s.Pot[i++] = v.IP;
                        } break;
                    case 1:
                        Step_=Step_CV_1;
                        s.Cur[i++] = tmp[0];
                        break;
                    case 2:
                        Step_=Step_CV_2;
                        s.Cur[i++] = tmp[1];
                        break;
                    case 3:
                        Step_=Step_CV_3;
                        s.Cur[i++] = tmp[2];
                        break;
                    case 4:
                        Step_=Step_CV_4;
                        s.Cur[i++] = tmp[3];
                        break;
                    case 5:
                        Step_=Step_CV_Ramp;
                        s.Cur[i] = (tmp[0] + tmp[1] + tmp [2] + tmp [3]) / 4.0f;
                        s.Pot[i++] = v.IP;
                        break;
                }
                break;

            case mode_CA:
                for (k=0; k<10; k++) {
                    s.Cur[i]=0;
                    s.Pot[i++]=v.IP;
                }
                np1=ROUND(v.VD/v.ST)+10;
                for (; i<np1; i++) {
                    s.Pot[i]=v.V1;
                    s.Cur[i]=RDCstep(s.NI, v.V1, m);
                    //Plot(&i, COL_RUN);
                }
                if (v.T2>0.0) {
                    np1+=ROUND(v.T2/v.ST);
                    for (; i<np1; i++) {
                        s.Pot[i]=v.V2;
                        s.Cur[i]=RDCstep(s.NI, v.V2, m);
                        //Plot(&i, COL_RUN);
                    }
                }
                goto Noscan;
        }

        /* Begin the simulation */
        switch (v.ncycle) {
            case 1:
                Scan(v.IP, v.FP, &i, m);
                break;
            case 2:
                Scan(v.IP, v.V1, &i, m);
                if (v.VD>0.0) Aspetta(v.VD, v.V1, m);
                Scan(v.V1, v.FP, &i, m);
                break;
            case 3:
                Scan(v.IP, v.V1, &i, m);
                if (v.VD>0.0) Aspetta(v.VD, v.V1, m);
                Scan(v.V1, v.V2, &i, m);
                if (v.VD>0.0) Aspetta(v.VD, v.V2, m);
                Scan(v.V2, v.FP, &i, m);
                break;
        } // end switch

Noscan:
        if (cycle_num == v.SC) oldi=i;
        if (cycle_num < v.NC) {
			// do this between cycles
            if (v.VD>0.0) Aspetta(v.VD, v.FP, m);
            if (cycle_num < v.SC) i=1; else i=oldi;
        }
    } // end while
}

void Scan(int pot, int fine, int *idx, MEC_EXP *m)
{
    int xp, i;

    xp = abs(pot - fine) / v.SI;
    if (fine < pot) v.SI=-v.SI; else v.PH=-v.PH;

    for (i=0; i<xp; i++) {
        Step_(1, &pot, idx, m);
        //Plot(idx, COL_RUN);
    }

    v.SI=abs(v.SI);
    v.PH=abs(v.PH);
}

void Aspetta(double wtime, int pot, MEC_EXP *m)
{
    double t, st;
    int i;

    i = ROUND(0.1 / s.Deltat);
    st = s.Deltat * i;
    for (t=0; t<wtime; t+=st) {
        RDCstep(i, pot, m);
    }
}

/* Here, with Step_ functions, I have increased the code size; however,
   the speed is increased, since there are few function calls and no if

   About RU/DL calculations: if there is no RU, v.RU=0, so calculations
   are done but have no effects.
   The calculation follows different paths depending if we have fitting
   or not. Whitout fit it is necessary after the RU/DL calculation to
   cycle if the "new" corrected potential (pot_2) is different from the
   old (pot_1).
*/

void Step_CV_All(int xp, int *pot, int *idx, MEC_EXP *m)
{
    double pot_1;
	int i;

    do {
        *pot += v.SI;

        for (i=0; i<4; i++) {
			pot_1 = correct_IR(*pot);

            s.Cur[*idx] = RDCstep_reduced(s.NIDIV4, pot_1, m);
            s.Pot[(*idx)++] = *pot;
        }
    } while (--xp);
}

void Step_CV_1(int xp, int *pot, int *idx, MEC_EXP *m)
{
    double pot_1;

    do {
        *pot += v.SI;
		pot_1 = correct_IR(*pot);

		s.Cur[*idx] = RDCstep_reduced(s.NIDIV4, pot_1, m);
        s.Pot[(*idx)++] = *pot;
		RDCstep_reduced(s.NIDIV4*3, pot_1, m);
    } while (--xp);
}

void Step_CV_2(int xp, int *pot, int *idx, MEC_EXP *m)
{
    double pot_1;

    do {
        *pot += v.SI;
		pot_1 = correct_IR(*pot);

        s.Cur[*idx] = RDCstep_reduced(s.NIDIV2, pot_1, m);
        s.Pot[(*idx)++] = *pot;
        RDCstep_reduced(s.NIDIV2, pot_1, m);
    } while (--xp);
}

void Step_CV_3(int xp, int *pot, int *idx, MEC_EXP *m)
{
    double pot_1;

    do {
        *pot += v.SI;
		pot_1 = correct_IR(*pot);

    	s.Cur[*idx] = RDCstep_reduced(s.NIDIV4*3, pot_1, m);
        RDCstep_reduced(s.NIDIV4, pot_1, m);
        s.Pot[(*idx)++] = *pot;
    } while (--xp);
}

void Step_CV_4(int xp, int *pot, int *idx, MEC_EXP *m)
{
    double pot_1;

    do {
        *pot += v.SI;
		pot_1 = correct_IR(*pot);

        s.Cur[*idx] = RDCstep_reduced(s.NI, pot_1, m);
        s.Pot[(*idx)++] = *pot;
    } while (--xp);
}

void Step_SDC(int xp, int *pot, int *idx, MEC_EXP *m)
{
    int i, j;
    double pot_1;

    do {
        for (i=0; i < m->nspec; i++)
            for (j=0; j<=s.NS; j++) *(m->C+i*MAXGRID+j) = *(m->C+i*MAXGRID+MAXGRID-1);

        *pot += v.SI;
		pot_1 = correct_IR(*pot);

        s.Cur[*idx] = RDCstep_reduced(s.NI, pot_1, m);
        s.Pot[(*idx)++] = *pot;
    } while (--xp);
}

void Step_CV_Ramp(int xp, int *pot, int *idx, MEC_EXP *m)
{
    int i;
    double pot_1;

    do {
        *pot += v.SI;
		pot_1 = correct_IR(*pot);
        s.Cur[*idx] = 0.0f;

        for (i=0; i<4; i++) {
            s.Cur[*idx] += RDCstep_reduced(s.NIDIV4, pot_1, m);
        }

        s.Cur[*idx] /= 4.0f;
        s.Pot[(*idx)++] = *pot;
    } while (--xp);
}


void Step_SWV(int xp, int *pot, int *idx, MEC_EXP *m)
{
    int tmp;
    double pot_1;

    do {
        *pot += v.SI;
		tmp = *pot - v.PH;
        pot_1 = (double)tmp + v.RU*icap0;
        if (v.IR == ir_dl) icap0 -= (double)(v.SI-v.PH)/v.RU;

        s.Cur[*idx] = RDCstep_reduced(s.NIDIV2, pot_1, m);
        s.Pot[(*idx)++] = tmp;

        tmp = *pot + v.PH;
		pot_1 = (double)tmp + v.RU*icap0;

		if (v.IR == ir_dl) icap0 -= (double)v.PH/v.RU;

        s.Cur[*idx] = RDCstep_reduced(s.NIDIV2, pot_1, m);
        s.Pot[(*idx)++] = tmp;
    } while (--xp);
}

double RDCstep(int num, int potE, MEC_EXP *m)
{
    double pot_1;

    if (v.IR == ir_dl) icap0 -= (double)(potE-oldpot)/v.RU;
    oldpot = potE;
	pot_1 = (double)potE + v.RU*icap0;

    return(RDCstep_reduced(num, pot_1, m));
}

double RDCstep_reduced(int num, int pot_1, MEC_EXP *m)
{
    int i;
	double cur1 = 0.0f;

	for (i=0; i<num; i++) {
        cur1 = Redox(pot_1, m);
        if (v.IR==ir_dl) cur1 += ((icap0*=Capacity)-intcp)/slope;
        Diffusione(m);
        if (m->nchim) Chimica(m);
    }

    return(cur1);
}

int correct_IR(int potE)
{
	int pot_1;
	pot_1 = (double)potE + v.RU*icap0;
	if (v.IR == ir_dl) icap0 -= (double)v.SI/v.RU;

	return(pot_1);
}

double Redox(double potE, MEC_EXP *m)
{
    const static double q1=2.05946309436;
    static int i;
    static double zz, pota, Kf, Kb, qRed, qOx;
    static double ReA[MAXREDOX], ReB[MAXREDOX], ReC[MAXREDOX], ReR[MAXREDOX], Flux[MAXREDOX];
    static double gam[MAXREDOX], bet;
	static uint64_t startTime;
	startTime = get_posix_clock_time();


    /* Redox MUST have the correct order, to obtain tridiagonalized matrix.
    NOTE: species cannot have more than 1 reduction or more than 1 oxidation.
    1)    In the first redox, Red(1) + n(1)e- <==> Ox(1), Red can never
          appear as Ox, Red(1) != Ox(i), 0<i<nredox.
    2)    If exist a redox having Red(i)=Ox(1), this redox must follow the
          previous redox.
          Otherwise any other is fine.
    3) goto 2
    */

    for (i=0; i < m->nredox; i++) {
        pota = (potE - m->Re[i].E) * s.FRT * m->Re[i].n; // normalized potential
        Kf = m->Re[i].Ke * exp(-m->Re[i].a * pota);
        Kb = m->Re[i].Ke * exp((1.0 - m->Re[i].a) * pota);
        qRed = q1 * m->D[m->Re[i].Red];
        qOx = q1 * m->D[m->Re[i].Ox];
        ReA[i] = (m->Re[i].Ox == m->Re[i-1].Red) ? -Kf/qOx : 0.0;
        ReB[i] = 1.0 + Kf/qOx + Kb/qRed;
        ReC[i] = (m->Re[i].Red == m->Re[i+1].Ox) ? -Kb/qRed : 0.0;
        ReR[i] = *(m->C+MAXGRID*m->Re[i].Ox)*Kf - *(m->C+MAXGRID*m->Re[i].Red)*Kb;
    }


    /* The tridiagonal matrix is stored into vectors ReA, ReB and ReC like:
       b1 c1  0  0 ...   The vector ReR are the known terms.
       a2 b2 c2  0 ...   The vector Flux store the solutions.
        0 a3 b3 c3 ...   nredox is the dimension. */

    bet = ReB[0];
    Flux[0] = ReR[0] / bet;

    for (i=1; i < m->nredox; i++) {
        gam[i] = ReC[i-1] / bet;
        bet = ReB[i] - ReA[i]*gam[i];
        Flux[i] = (ReR[i] - ReA[i]*Flux[i-1]) / bet;
    }

    for (i=m->nredox-1; i>0; i--) Flux[i-1] -= gam[i]*Flux[i];

    /* Compute the net current and mass Flux at electrode */
    for (i=0; i<m->nspec; i++) s.FluxJ[i]=0.0;
    for (i=zz=0; i < m->nredox; i++) {
        zz += Flux[i] * m->Re[i].n;
        s.FluxJ[m->Re[i].Ox] += Flux[i];
        s.FluxJ[m->Re[i].Red] -= Flux[i];
    }

    zz *= current_norm;

	elapsed_redox += get_posix_clock_time() - startTime;

    return ((double)zz);
}

void Chimica(MEC_EXP *m)
{
    static int i, iter, j, k, r2, p1, p2;
    static double c0[MAXSPEC], r[4][MAXCHIM], rb, *xc, *ptc, *rix;
    const static double r_con[4]={2.0,2.0,1.0,6.0};
	static uint64_t startTime;
	startTime = get_posix_clock_time();

    /* Runge-Kutta of order 4.
       rate(C) compute the homogeneous rate as function of C (concentrations)
       What we are doing here:

       c  = c0
       r1 = rate(c)
       c  = c0 + r1/2
       r2 = rate(c)
       c  = c0 + r2/2
       r3 = rate(c)
       c  = c0 + r3
       r4 = rate(c)
       c  = c0 + (r1 + 2 * (r2 + r3) + r4)/6
    */
    for (i=0; i <= s.NS; i++) {
        /* xc point to species 0, space grid i.
          (xc+k*MAXGRID) point to species k, space grid i. */
        xc = m->C + i;
        /* store original concentrations (grid i) of all species in c0[]: */
        for (k=0, ptc=xc; k < m->nspec; k++, ptc+=MAXGRID) c0[k] = *ptc;

        /* do the 4th order Runge-Kutta integration */
        for (iter=0; iter<4; iter++) {
            /* calculate reaction rates and store in r[iter][k] */
            for (k=0, rix=r[iter]; k < m->nchim; k++, rix++) {
                *rix = m->Ch[k].Kf * *(xc + m->Ch[k].a*MAXGRID);
                if ((r2=m->Ch[k].b) >= 0) *rix *= *(xc + r2*MAXGRID);
                if ((p1=m->Ch[k].c) >= 0) {
                    rb = m->Ch[k].Kb * *(xc + p1*MAXGRID);
                    if ((p2=m->Ch[k].d) >= 0) rb *= *(xc + p2*MAXGRID);
                    *rix -= rb;
                }
            }
            /* calculate next concentrations */
            for (j=0, ptc=xc; j < m->nspec; j++, ptc+=MAXGRID)
                for (k=0; k < m->nchim; k++)
                    *ptc=c0[j]+(double)s.c_idx[j][k]*r[iter][k]/r_con[iter];
        }

        for (j=0, ptc=xc; j < m->nspec; j++, ptc+=MAXGRID)
            for (k=0; k < m->nchim; k++)
                *ptc += (double)s.c_idx[j][k] * (r[0][k]+2.0*(r[1][k]+r[2][k]))/6.0;
    }

	elapsed_chim += get_posix_clock_time() - startTime;
}

void Diffusione(MEC_EXP *m)
{
    const static double p2=0.749153538438, q2=0.840896415254,
                        p3=0.594603557501, q3=0.667419927085,
                        p4=0.471937156342, q4=0.529731547180,
                        p1=0.943874312682;
    double ct0, ct1, ct2, ct3, *ct, d;
    int i,j,k=s.NS/3;
	static uint64_t startTime;
	startTime = get_posix_clock_time();

    for (i=0; i < m->nspec; i++) {
        ct = m->C + i*MAXGRID;
        ct0 = *(ct+1) - *ct;
        *ct++ += (d=m->D[i]) * p1 * ct0 - s.FluxJ[i];
        for (j=0; j < k; j++) {
            ct1 = *(ct+1) - *ct;
            ct2 = *(ct+2) - *(ct+1);
            ct3 = *(ct+3) - *(ct+2);
            *ct++ += d * (p2 * ct1 - q2 * ct0);
            *ct++ += d * (p3 * ct2 - q3 * ct1);
            *ct++ += d * (p4 * ct3 - q4 * ct2);
            ct0 = ct3;
            d /= 2.0;
        }
    }

	elapsed_diff += get_posix_clock_time() - startTime;
}

void Set_Simul(MEC_EXP *m)
{
    int i,j;
    double maxD, maxC, Time_exp, x, sg1;

    /*
       WARNING: Set_Simul modify many MEC_EXP variables. If you do more
       than one simulation (as in Best-Fitting procedure) you need to
       keep the old variables.
            Here are modified:

       m.Ch[].a m.Ch[].b    Because in simulation the species
       m.Ch[].c m.Ch[].d    start from 0 and not from 1.

       m.Ch[].Kf m.Ch[].Kb  To improve speed, some calculations are done here.
       m.Re[].Ke m.D[] m.C[]
       m.Re[].Ox m.Re[].Red                                    */

	// determine maximum concentration and diff coeff:
    maxD = (double)m->D[0]; maxC = C0[0];
    for (i=1; i < m->nspec; i++) {
        maxD = fmax(maxD, (double)m->D[i]);
        maxC = fmax(maxC, C0[i]);
    }

	// normalize concentration and diff coeff:
    for (i=0; i < m->nspec; i++) {
        m->D[i] *= DSTAR / maxD;             /* normalize Diff. coeff. */
        *(m->C + i*MAXGRID + MAXGRID - 1) = C0[i]/maxC;      /* normalize Conc. */
    }

	// calculate NI, which is the amounts of simulation steps per exp. step
    Calcola_NI(m);

	// calculate s.c_idx - this seems to be a matrix that holds all chemical reaction information
    if (m->nchim) {
        for (i=0; i < m->nchim; i++) {
            m->Ch[i].a--;
            m->Ch[i].b--;
            m->Ch[i].c--;
            m->Ch[i].d--;
        }
        for (i=0; i < m->nspec; i++)         /* compute the Chem index */
            for (j=0; j < m->nchim; j++)
                s.c_idx[i][j]= (m->Ch[j].c == i) + (m->Ch[j].d == i) - (m->Ch[j].a == i) - (m->Ch[j].b == i);
    }

	// s.Deltat is the step time in the simulation
    s.Deltat = v.ST / s.NI;

    /*
    		Space Grid
    */
    s.Deltax = sqrt(maxD * s.Deltat / DSTAR); // spacing of space grid

	// determine total experiment time for all the cycles:
    Time_exp=(double)v.NP*v.ST; // number of data points * step time
    switch (v.Mode) {
        case mode_SWV: Time_exp /= 2.0;
        case mode_CV: if (!v.AM) Time_exp /= 4.0;
        case mode_SDC: if (v.VD>0.0) Time_exp += v.VD * (v.ncycle - 1);
            Time_exp *= (double) v.NC;
        case mode_CA: break;
    }
	// add optional conditioning and equilibration time:
    oldpot=v.IP;
    if (v.CT>0.0) { Time_exp += v.CT; oldpot=v.CP; }
    if (v.ET>0.0) { Time_exp += v.ET; oldpot=v.IP; }

	// determine space grid size:
    x = 6.0 * sqrt(maxD * Time_exp);
    sg1 = log(2.0) / 6.0;
    s.NS=1;
    do (s.NS+=3); while ((s.Deltax * (exp(sg1 * ((double) s.NS - .5)) - 1.0) / (exp(sg1) - 1.0)) < x);
	// if space grid is too big, print warning and reduce in size:
	printf("Num points in space grid: %d\n", s.NS);
    if (s.NS >= MAXGRID) {
		printf("NS too large [%d]. Max space grid allowed: %d\r\n", s.NS, MAXGRID);
        while (s.NS>=MAXGRID) s.NS-=3;
    }

	// normalize Ke
    for (i=0; i < m->nredox; i++) {
        m->Re[i].Ke *= s.Deltat / s.Deltax;
        m->Re[i].Ox--; m->Re[i].Red--;
    }

	// do something to all _unused_ redox species...
    for (i = m->nredox; i < MAXREDOX; i++)
        m->Re[i].Red = m->Re[i].Ox = -1;

	// initialize (normalized!) concentrations space grid m->C
    for (i=0; i < m->nspec; i++) {
        s.FluxJ[i] = 0.0;
        for (j=0; j<=s.NS; j++) *(m->C+i*MAXGRID+j) = C0[i]/maxC;;
    }

	// normalize chemical rates:
    for (i=0; i < m->nchim; i++) {
        m->Ch[i].Kf *= s.Deltat;
        if (m->Ch[i].a >= 0 && m->Ch[i].b >= 0) m->Ch[i].Kf *= maxC; // bimolecular rxn
        m->Ch[i].Kb *= s.Deltat;
        if (m->Ch[i].c >= 0 && m->Ch[i].d >= 0) m->Ch[i].Kb *= maxC; // bimolecular rxn
    }

	// calculate some constants:
    s.NIDIV4 = s.NI / 4; s.NIDIV2 = s.NI / 2;
    s.FRT = 11.6045035066 / v.TE;

	// current normalization factor:
    current_norm = maxC * s.Deltax / s.Deltat * 0.0964846;
    /* the DME electrode is treated like flat electrode which is moving toward
    the bulk of the solution. The surface area is calculated as a growing
    sphere, which mass depend by the mercury flow MF and the time ST */
    if (v.WE==WE_DME) current_norm *= sqrt(7.0/3.0) * 4.0*M_PI*pow(3.0*v.AR*v.ST/(4.0*M_PI*13530.0), 2.0/3.0);
    else current_norm *= v.AR;

	// IR compensation:
    icap0=0.0f;
    if (v.IR==ir_none) v.RU=0.0;
    if (v.IR==ir_dl) Capacity=exp(-s.Deltat/(v.RU/1000.0*v.DL));
}

void Calcola_NI(MEC_EXP *m)
{
    double maxC=0.0, mK, Kmax=0.0;
    int i;
    /*      by using Runge-Kutta integration (of order 4), it is necessary,
      to get accurate results, that the product of step time (Deltat) and
      Kmax (maximum homogeneous rate constant) should be less than 0.33:

      Deltat * Kmax <= 0.33

      Deltat=v.ST/s.NI    ���>  v.ST * Kmax * 3 <= s.NI     */

    maxC = C0[0];
    for (i=1; i < m->nspec; i++) maxC = fmax(maxC, C0[i]);

	// determine maximum rate constant Kmax
    for (i=0; i < m->nchim; i++) {
        mK = m->Ch[i].Kf;
        if (m->Ch[i].a > 0 && m->Ch[i].b > 0) mK *= maxC;
        Kmax = fmax(Kmax, mK);
        mK = m->Ch[i].Kb;
        if (m->Ch[i].c > 0 && m->Ch[i].d > 0) mK *= maxC;
        Kmax = fmax(Kmax, mK);
    }

	// first, set s.NI to the minimum value necessary:
    if (s.NI < 4) {
        s.NI=4;
        switch (v.Mode) {
            case mode_CV:
                switch (v.AM) {
                    case 0: /* All */
                    case 5: /* Ramp */
                    case 1: s.NI*=2;
                    case 2: s.NI*=2;
                    case 3: s.NI*=2;
                    case 4: s.NI+=v.SI*4;
                }
                break;
            case mode_SWV: s.NI*=4;
            case mode_CA:
            case mode_SDC:
            default: s.NI+=v.SI*4; // huh? adding a potential to a time?
        }
    }

	// if s.NI cannot accomodate the highest rate, increase it:
    if ((int)(Kmax * v.ST * 3.0) > s.NI)
        s.NI = (int) (Kmax * v.ST * 3.0);

	// increase s.NI until it is a multiple of 4:
    while (s.NI % 4) (s.NI)++;
}

int CalcolaNP(void)
{
	int np, np2;

	if (v.Mode==mode_CA) {
		np=ROUND(v.VD/v.ST)+10;
		v.VD=v.ST*(np-10);
		if (v.T2>0.0) {np2=ROUND(v.T2/v.ST); v.T2=v.ST*np2; np+=np2;}
	} else {
		switch (v.ncycle) {
			case 1:
				np=abs(v.IP - v.FP) / v.SI + 1;
				break;
			case 2:
				np=abs(v.IP - v.V1) / v.SI + abs(v.V1 - v.FP) / v.SI + 1;
				break;
			case 3:
				np=abs(v.IP - v.V1) / v.SI + abs(v.V1 - v.V2) / v.SI + abs(v.V2 - v.FP) / v.SI + 1;
				break;
			default:
				np = 0;
				break;
		}
		if (v.Mode==mode_SWV) np = 2*np + v.ncycle;
		else if (v.Mode==mode_CV && !v.AM) np *= 4;
		if (v.NC>1 && v.SC!=v.NC) np*=2;
	}

	return(np);
}
