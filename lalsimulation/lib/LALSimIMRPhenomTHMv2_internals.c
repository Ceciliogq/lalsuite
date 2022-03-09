/*
 * Copyright (C) 2022 Hector Estelles
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


/**
 * \author Hector Estelles
 */

#include <math.h>

#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomTHMv2_fits.c"

/* LAL Header Files */
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>

/* GSL Header Files */
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>

/* Functions definition */
/* Function for computing orbital phase coefficients */

int IMRPhenomTSetPhase22v2Coefficients(IMRPhenomTPhase22Struct *pPhase, IMRPhenomTWaveformStruct *wf);
int IMRPhenomTSetHMAmplitudev2Coefficients(IMRPhenomTHMAmpStruct *pAmp, IMRPhenomTPhase22Struct *pPhase, IMRPhenomTWaveformStruct *wf);

/* Inspiral ansatz */

double IMRPhenomTInspiral_dtdx_TaylorT2(REAL8 x, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiral_dtdx_Calibrated_Early(REAL8 x, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiral_TofX_Calibrated_Early(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiral_PhiofX_Calibrated_Early(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);

double IMRPhenomTInspiral_dtdx_Calibrated_Late(REAL8 x, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiral_TofX_Calibrated_Late(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiral_PhiofX_Calibrated_Late(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);
double facNT4(REAL8 x, REAL8 eta);
double facT2(REAL8 x, REAL8 eta);

/* Merger ansatz */
double IMRPhenomTMerger_TofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTMerger_PhiofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase);

double IMRPhenomTv2_TofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTv2_PhiofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);
int IMRPhenomTv2_times_x_phase(
	REAL8Sequence **phaseorb,       /**< Values of the 22 phase for the waveform time array */
	REAL8Sequence **xorb,				/**< Values of the 22 frequency for the waveform time array */
	REAL8Sequence **times,
	IMRPhenomTWaveformStruct *wf,
	IMRPhenomTPhase22Struct *pPhase
);
COMPLEX16 IMRPhenomTInspiralAmpAnsatzHMv2(REAL8 x, IMRPhenomTHMAmpStruct *pAmp);
COMPLEX16 IMRPhenomTHMAmpv2(
  REAL8 t,
  UNUSED REAL8 x,
  IMRPhenomTHMAmpStruct *pAmp
);


/* ************************************* */
/* ******* STRUCTS INITIALIZATION ****** */
/* ************************************* */

int IMRPhenomTSetPhase22v2Coefficients(IMRPhenomTPhase22Struct *pPhase, IMRPhenomTWaveformStruct *wf){

	/* ***** Initialize parameters and variables ****** */

	REAL8 eta = wf->eta;     // Symmetric mass ratio
	REAL8 chi1L = wf->chi1L; // Dimensionless aligned spin of companion 1
	REAL8 chi2L = wf->chi2L; // Dimensionless aligned spin of companion 2
	REAL8 S = wf->Shat;      // Effective dimensionless spin parameter
	REAL8 dchi = wf->dchi;   // Dimensionless spin difference chi1L - chi2L
	REAL8 delta = wf->delta; // Asymmetry parameter
	pPhase->dtM = wf->dtM; // Dimensionless sampling rate

	REAL8 xCutInsp = 0.06 + 0.04*pow(1.0-4.0*eta,3);
	REAL8 xCutMerger = IMRPhenomTv2_xCutMergerFit(eta, S, dchi, delta);
	REAL8 xCutMerger2 = xCutMerger - 0.003;
  REAL8 xpeak = IMRPhenomTv2_xPeak(eta, S, dchi, delta);

	pPhase->xCutInsp = xCutInsp;
	pPhase->xCutMerger = xCutMerger2;
	pPhase->xpeak = xpeak;
	printf("---------------\n");
	printf("CASE: eta: %f chi1L: %f chi2L: %f\n", eta, chi1L, chi2L);
	printf("xcutI: %.8f xcutM: %.8f xpeak: %.8f\n", xCutInsp, xCutMerger, xpeak);

	/* Ringdown and damping frequency of final BH for the 22 mode (defined on LALSimIMRPhenomTHM_fits.c) */
	REAL8 fRING     = evaluate_QNMfit_fring22(wf->afinal) / (wf->Mfinal); // 22 mode ringdown frequency
	REAL8 fDAMP     = evaluate_QNMfit_fdamp22(wf->afinal) / (wf->Mfinal); // 22 mode damping frequency of ground state (n=1) QNM

	/* Angular ringdown and damping frequencies (omega =2*pi*f) */
	REAL8 alpha1RD  = 2*LAL_PI*fDAMP;
	REAL8 omegaRING = 2*LAL_PI*fRING;

	pPhase->omegaRING = omegaRING;
	pPhase->alpha1RD  = alpha1RD;

	pPhase->omegaRING_prec = 2*LAL_PI*evaluate_QNMfit_fring22(wf->afinal_prec) / (wf->Mfinal);
	pPhase->alpha1RD_prec  = 2*LAL_PI*evaluate_QNMfit_fdamp22(wf->afinal_prec) / (wf->Mfinal);

	pPhase->EulerRDslope = 2*LAL_PI*(evaluate_QNMfit_fring22(wf->afinal_prec) / (wf->Mfinal) - evaluate_QNMfit_fring21(wf->afinal_prec) / (wf->Mfinal)); // FIXME:
	if(wf->afinal<0)
	{
		pPhase->EulerRDslope = -pPhase->EulerRDslope;
	}

	/* Dimensionless minimum and reference frequencies */
	pPhase->MfRef = wf->MfRef;
	pPhase->Mfmin = wf->Mfmin;

  pPhase->xmin = pow(LAL_PI*pPhase->Mfmin,2./3);
  pPhase->xRef = pow(LAL_PI*pPhase->MfRef,2./3);

	/* GSL objects for solving system of equations via LU decomposition */
	gsl_vector *b, *sol;
	gsl_matrix *A;
	gsl_permutation *p;
	int s; /* Sign of permutation */

	/* ********************************** */
	/* *** RINGDOWN COEFFICIENTS ******** */
	/* ********************************** */

	/* For phase and frequency ringdown ansatz, coefficients are obtained from ringdown and damping frequencies, frequency at peak amplitude
	   and phenomenological fits of two free coefficients. No linear system solving is needed.
	   Coefficient c1 is defined in Eq. [9] of Damour&Nagar 2014 (10.1103/PhysRevD.90.024054, https://arxiv.org/abs/1406.0401).
	   Coefficients c2 and c3 are calibrated to NR/Teukolsky data.
	   Coefficient c4 is set to zero.
	   Explained also in Sec II.C, eq. 26 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */

	pPhase->omegaPeak = 2*pow(xpeak,3./2); // 22 Frequency at the peak 22 amplitude time (t=0 by definition)

	pPhase->c3 = IMRPhenomT_RD_Freq_D3_22(eta, S, dchi, delta);
	pPhase->c2 = IMRPhenomT_RD_Freq_D2_22(eta, S, dchi, delta);
	pPhase->c4 = 0.0;
	pPhase->c1 = (1 + pPhase->c3 + pPhase->c4)*(pPhase->omegaRING - pPhase->omegaPeak)/pPhase->c2/(pPhase->c3 + 2*pPhase->c4);
	pPhase->c1_prec = (1 + pPhase->c3 + pPhase->c4)*(pPhase->omegaRING_prec - pPhase->omegaPeak)/pPhase->c2/(pPhase->c3 + 2*pPhase->c4);

	/* ********************************** */
	/* *** INSPIRAL COEFFICIENTS ******** */
	/* ********************************** */

	/* PN coefficients of the TaylorT3 approximant employed in the inspiral ansatz are defined here.
	   Unknown coefficients for the higher order extension in the inspiral ansatz are obtained through fitted collocation points
	   solving a linear system of equations. */


	/*PN coefficients corresponding to the TaylorT3 implementation of the model,
	  as defined in Appendix A1 of Estelles et al 2020 (https://arxiv.org/pdf/2004.08302.pdf) */

	pPhase->dT2_1PN = 2.2113095238095237 + (11*eta)/4.;

	pPhase->dT2_1halfPN = (113*chi1L)/24. + (113*chi2L)/24. - (19*chi1L*eta)/6. - (19*chi2L*eta)/6. - 4*LAL_PI + (113*chi1L*pow(1 - 4*eta,0.5))/24. - (113*chi2L*pow(1 - 4*eta,0.5))/24.;

	pPhase->dT2_2PN = 3.010315295099521 + (5429*eta)/1008. - (79*chi1L*chi2L*eta)/8. - (81*pow(chi1L,2))/32. + (81*eta*pow(chi1L,2))/16. - (81*pow(chi2L,2))/32. + (81*eta*pow(chi2L,2))/16. - (81*pow(chi1L,2)*pow(1 - 4*eta,0.5))/32. + (81*pow(chi2L,2)*pow(1 - 4*eta,0.5))/32. + (617*pow(eta,2))/144.;

  pPhase->dT2_2halfPN = (146597*chi1L)/4032. + (146597*chi2L)/4032. - (1213*chi1L*eta)/36. - (1213*chi2L*eta)/36. - (7729*LAL_PI)/672. + (13*eta*LAL_PI)/8. + (146597*chi1L*pow(1 - 4*eta,0.5))/4032. - (146597*chi2L*pow(1 - 4*eta,0.5))/4032. + (7*chi1L*eta*pow(1 - 4*eta,0.5))/4. - (7*chi2L*eta*pow(1 - 4*eta,0.5))/4. - (17*chi1L*pow(eta,2))/4. - (17*chi2L*pow(eta,2))/4.;

  pPhase->dT2_3PN = 22.065018753916434 + (3147553127*eta)/1.2192768e7 - (6535*chi1L*chi2L*eta)/448. - (227*chi1L*LAL_PI)/12. - (227*chi2L*LAL_PI)/12. + 13*chi1L*eta*LAL_PI + 13*chi2L*eta*LAL_PI - (15103*pow(chi1L,2))/2304. + (21683*eta*pow(chi1L,2))/1152. - (15103*pow(chi2L,2))/2304. + (21683*eta*pow(chi2L,2))/1152. - (227*chi1L*LAL_PI*pow(1 - 4*eta,0.5))/12. + (227*chi2L*LAL_PI*pow(1 - 4*eta,0.5))/12. - (15103*pow(chi1L,2)*pow(1 - 4*eta,0.5))/2304. + (1645*eta*pow(chi1L,2)*pow(1 - 4*eta,0.5))/288. + (15103*pow(chi2L,2)*pow(1 - 4*eta,0.5))/2304. - (1645*eta*pow(chi2L,2)*pow(1 - 4*eta,0.5))/288. - (15211*pow(eta,2))/6912. - (1115*chi1L*chi2L*pow(eta,2))/72. + (613*pow(chi1L,2)*pow(eta,2))/144. + (613*pow(chi2L,2)*pow(eta,2))/144. + (25565*pow(eta,3))/5184. - (451*eta*pow(LAL_PI,2))/48.;

	pPhase->dT2_3PNlog = 856./105.;

  pPhase->dT2_3halfPN = (5030016755*chi1L)/2.4385536e7 + (5030016755*chi2L)/2.4385536e7 - (2113331119*chi1L*eta)/6.096384e6 - (2113331119*chi2L*eta)/6.096384e6 - (15419335*LAL_PI)/1.016064e6 - (75703*eta*LAL_PI)/6048. + 55*chi1L*chi2L*eta*LAL_PI - (2813*chi2L*eta*pow(chi1L,2))/64. + (57*LAL_PI*pow(chi1L,2))/4. - (57*eta*LAL_PI*pow(chi1L,2))/2. - (2917*pow(chi1L,3))/192. + (8819*eta*pow(chi1L,3))/192. - (2813*chi1L*eta*pow(chi2L,2))/64. + (57*LAL_PI*pow(chi2L,2))/4. - (57*eta*LAL_PI*pow(chi2L,2))/2. - (2917*pow(chi2L,3))/192. + (8819*eta*pow(chi2L,3))/192. + (5030016755*chi1L*pow(1 - 4*eta,0.5))/2.4385536e7 - (5030016755*chi2L*pow(1 - 4*eta,0.5))/2.4385536e7 - (5360987*chi1L*eta*pow(1 - 4*eta,0.5))/48384. + (5360987*chi2L*eta*pow(1 - 4*eta,0.5))/48384. - (2813*chi2L*eta*pow(chi1L,2)*pow(1 - 4*eta,0.5))/64. + (57*LAL_PI*pow(chi1L,2)*pow(1 - 4*eta,0.5))/4. - (2917*pow(chi1L,3)*pow(1 - 4*eta,0.5))/192. + (995*eta*pow(chi1L,3)*pow(1 - 4*eta,0.5))/64. + (2813*chi1L*eta*pow(chi2L,2)*pow(1 - 4*eta,0.5))/64. - (57*LAL_PI*pow(chi2L,2)*pow(1 - 4*eta,0.5))/4. + (2917*pow(chi2L,3)*pow(1 - 4*eta,0.5))/192. - (995*eta*pow(chi2L,3)*pow(1 - 4*eta,0.5))/64. + (208433*chi1L*pow(eta,2))/24192. + (208433*chi2L*pow(eta,2))/24192. + (14809*LAL_PI*pow(eta,2))/3024. - (chi2L*pow(chi1L,2)*pow(eta,2))/8. - (17*pow(chi1L,3)*pow(eta,2))/24. - (chi1L*pow(chi2L,2)*pow(eta,2))/8. - (17*pow(chi2L,3)*pow(eta,2))/24. + (397*chi1L*pow(1 - 4*eta,0.5)*pow(eta,2))/384. - (397*chi2L*pow(1 - 4*eta,0.5)*pow(eta,2))/384. - (1069*chi1L*pow(eta,3))/288. - (1069*chi2L*pow(eta,3))/288.;


  /* Solve inspiral coefficients */

  /* Inspiral reconstruction is based on the TaylorT3 approximant at 3.5PN with spin-orbit and spin-spin contributions, whose coefficients are defined above.
   TaylorT3 is parameterised in terms of a 'PN time parameter' theta=pow(eta*(tt0 - tEarly)/5,-1./8), which depends on an integration constant tt0 (coming from the TaylorT2 derivation) which has to be fixed
   in order to impose a determined frequency at a determined time. TaylorT3 will diverge when t=tt0, so this parameter can be understood as the merger time prediction of TaylorT3, which in general will not be
   the actual merger time. In order to overcome the divergence issue in cases where this occurs too early, we originally decided to set tt0 to the peak time of the 22. This produces a small frequency shift that
   over many cycles can reduce the accuracy of the model for long waveforms. To overcome this, an early collocation point at theta=0.33 is set to the actual value of TaylorT3 (through a tt0 fit) to impose
   the low frequency regime to follow the correct description. Along with this, the early inspiral region with theta<0.33 was selected to be modeled with pure TaylorT3 with the tt0 fit value.
   However, several tests suggest that actually imposing the theta=0.33 value is enough and two regions are not needed. We let as an option in the model to select the inspiral reconstruction with or without this split.
   Depending on the value of the LAL parameter 'PhenomTHMInspiralVersion' reconstruction of the inspiral phase can be splitted into two different regions. For default value (0) this will not be done.
   In default reconstruction, TaylorT3 with t0=0 will be employed, plus 6 additional higher order terms to be solved with the value of 5 fitted collocation points at fixed theta positions {0.45, 0.55, 0.65, 0.75, 0.82}
   + and early collocation point at theta=0.33 whose value will be the value of TaylorT3 with t0 computed from a fit.
   In non-default reconstruction, region with theta<0.33 will be modelled by TaylorT3 with the tt0 fit value and without extra coefficients.*/

	 /* initialize collocation point values for solving inspiral coefficient system */
	REAL8 dtinsppoints[3];
	dtinsppoints[0] = IMRPhenomTv2_Inspiral_CP1(eta, S, dchi, delta);
 	dtinsppoints[1] = IMRPhenomTv2_Inspiral_CP2(eta, S, dchi, delta);
 	dtinsppoints[2] = IMRPhenomTv2_Inspiral_CP3(eta, S, dchi, delta);

  /* Initialize higher order extension coefficients to 0. In this way, when calling the ansatz IMRPhenomTInspiralOmegaAnsatz22 it will return PN TaylorT3 value. */
	pPhase->dtInspC1 = 0.0;
	pPhase->dtInspC2 = 0.0;
	pPhase->dtInspC3 = 0.0;

 /* Set collocation points */
  REAL8 xpoints[4] = {0.045, 0.045 + (xCutInsp - 0.045)/2, xCutInsp};

	/*Set linear system, which is rank 6 */
	p = gsl_permutation_alloc(3);
	b = gsl_vector_alloc(3);
	sol = gsl_vector_alloc(3);
	A = gsl_matrix_alloc(3,3);

	/*Set A matrix and b vector*/

	REAL8 xx, xx9over2, xx5, xx11over2, dT2offset; // Initialize theta powers the diagonal A matrix
	/* theta is a weighted dimensionless time parameter defined in Eq. 315 of Blanchet 2014 (https://arxiv.org/abs/1310.1528).
	   Notice however that here is defined at the (-1/8) power, in order to represent the Post-Newtonian order (1/c). */

	/* Set up inspiral coefficient system.
	This system of equations is explained in eq. 10 of THM paper https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */
	for (UINT4 idx=0; idx<3; idx++)
	{
		/* Needed powers of x */
		xx = xpoints[idx];
		xx9over2 = pow(xx,4.5);
		xx5 = pow(xx,5);
		xx11over2 = pow(xx,5.5);

		//dT2offset = IMRPhenomTInspiral_dtdx_TaylorT2(xx, pPhase);

		gsl_vector_set(b,idx,(dtinsppoints[idx] /*- dT2offset*/));

		gsl_matrix_set(A,idx,0,xx9over2);
		gsl_matrix_set(A,idx,1,xx5);
		gsl_matrix_set(A,idx,2,xx11over2);

	}

	/* We now solve the system A x = b via an LU decomposition */
	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,sol);

	/* Set inspiral phenomenological coefficients from solution to A x = b */
	pPhase->dtInspC1 = gsl_vector_get(sol,0);
	pPhase->dtInspC2 = gsl_vector_get(sol,1);
	pPhase->dtInspC3 = gsl_vector_get(sol,2);

	printf("Insp coefs: %.8f %.8f %.8f\n", pPhase->dtInspC1, pPhase->dtInspC2, pPhase->dtInspC3);

	/* Boundary time between late inspiral and merger regions.
	This is selected to correspond to theta=0.81, earlier than the last collocation point at theta=0.82, because in this way
	the derivative at the boundary with the merger region, needed for the merger reconstruction, is less forced. */
	pPhase->tshiftT2 = 0.0;

	/* Tidy up in preparation for next GSL solve ... */
	gsl_vector_free(b);
	gsl_vector_free(sol);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

	/* ********************************** */
	/* *** Late Inspiral COEFFICIENTS ********** */
	/* ********************************** */

	/* initialize collocation point values for solving inspiral coefficient system */
 REAL8 tinsppoints[8];
 tinsppoints[0] = IMRPhenomTv2_LateInspiral_CP1(eta, S, dchi, delta);
 tinsppoints[1] = IMRPhenomTv2_LateInspiral_CP2(eta, S, dchi, delta);
 tinsppoints[2] = IMRPhenomTv2_LateInspiral_CP3(eta, S, dchi, delta);
 tinsppoints[3] = IMRPhenomTv2_LateInspiral_CP4(eta, S, dchi, delta);
 tinsppoints[4] = IMRPhenomTv2_LateInspiral_CP5(eta, S, dchi, delta);
 tinsppoints[5] = IMRPhenomTv2_LateInspiral_CP6(eta, S, dchi, delta);
 tinsppoints[6] = -100.0;
 tinsppoints[7] = IMRPhenomTInspiral_dtdx_Calibrated_Early(xCutInsp, pPhase);

 /* Initialize higher order extension coefficients to 0. In this way, when calling the ansatz IMRPhenomTInspiralOmegaAnsatz22 it will return PN TaylorT3 value. */
 pPhase->tLateC1 = 0.0;
 pPhase->tLateC2 = 0.0;
 pPhase->tLateC3 = 0.0;
 pPhase->tLateC4 = 0.0;
 pPhase->tLateC5 = 0.0;
 pPhase->tLateC6 = 0.0;
 pPhase->tLateC7 = 0.0;
 pPhase->tLateC8 = 0.0;

/* Set collocation points */
 REAL8 xpoints2[8] = {xCutInsp, xCutInsp + (xCutMerger - xCutInsp)/6, xCutInsp + 2.*(xCutMerger - xCutInsp)/6, xCutInsp + 3.*(xCutMerger - xCutInsp)/6, xCutInsp + 4.*(xCutMerger - xCutInsp)/6, xCutInsp + 5.*(xCutMerger - xCutInsp)/6, xCutMerger, xCutInsp};

 /*Set linear system, which is rank 6 */
 p = gsl_permutation_alloc(8);
 b = gsl_vector_alloc(8);
 sol = gsl_vector_alloc(8);
 A = gsl_matrix_alloc(8,8);

 /*Set A matrix and b vector*/

 REAL8 xx4, xx6, xx13over2, xx7, xx15over2, T2offset; // Initialize theta powers the diagonal A matrix
 /* theta is a weighted dimensionless time parameter defined in Eq. 315 of Blanchet 2014 (https://arxiv.org/abs/1310.1528).
		Notice however that here is defined at the (-1/8) power, in order to represent the Post-Newtonian order (1/c). */

 /* Set up inspiral coefficient system.
 This system of equations is explained in eq. 10 of THM paper https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */
 for (UINT4 idx=0; idx<7; idx++)
 {
	 /* Needed powers of x */
	 xx = xpoints2[idx];
	 xx4 = pow(xx,4);
	 xx9over2 = pow(xx,4.5);
	 xx5 = pow(xx,5);
	 xx11over2 = pow(xx,5.5);
	 xx6 = pow(xx,6);
	 xx13over2 = pow(xx,6.5);
	 xx7 = pow(xx,7);
	 xx15over2 = pow(xx,7.5);

	 T2offset = IMRPhenomTInspiral_TofX_Calibrated_Late(xx, wf, pPhase);

	 gsl_vector_set(b,idx,(tinsppoints[idx] - T2offset)/facT2(xx,eta));

	 gsl_matrix_set(A,idx,0,xx4);
	 gsl_matrix_set(A,idx,1,xx9over2);
	 gsl_matrix_set(A,idx,2,xx5);
	 gsl_matrix_set(A,idx,3,xx11over2);
	 gsl_matrix_set(A,idx,4,xx6);
	 gsl_matrix_set(A,idx,5,xx13over2);
	 gsl_matrix_set(A,idx,6,xx7);
	 gsl_matrix_set(A,idx,7,xx15over2);

 }
 xx = xCutInsp;
 xx4 = pow(xx,4);
 xx9over2 = pow(xx,4.5);
 xx5 = pow(xx,5);
 xx11over2 = pow(xx,5.5);
 xx6 = pow(xx,6);
 xx13over2 = pow(xx,6.5);
 xx7 = pow(xx,7);
 xx15over2 = pow(xx,7.5);

 dT2offset = IMRPhenomTInspiral_dtdx_TaylorT2(xCutInsp, pPhase);
 gsl_vector_set(b,7,(tinsppoints[7] - dT2offset));
 gsl_matrix_set(A,7,0,0);
 gsl_matrix_set(A,7,1,-0.125*xx9over2);
 gsl_matrix_set(A,7,2,-0.25*xx5);
 gsl_matrix_set(A,7,3,-0.375*xx11over2);
 gsl_matrix_set(A,7,4,-0.5*xx6);
 gsl_matrix_set(A,7,5,-0.625*xx13over2);
 gsl_matrix_set(A,7,6,-0.75*xx7);
 gsl_matrix_set(A,7,7,-0.875*xx15over2);

 /* We now solve the system A x = b via an LU decomposition */
 gsl_linalg_LU_decomp(A,p,&s);
 gsl_linalg_LU_solve(A,p,b,sol);

 /* Set inspiral phenomenological coefficients from solution to A x = b */
 pPhase->tLateC1 = gsl_vector_get(sol,0);
 pPhase->tLateC2 = gsl_vector_get(sol,1);
 pPhase->tLateC3 = gsl_vector_get(sol,2);
 pPhase->tLateC4 = gsl_vector_get(sol,3);
 pPhase->tLateC5 = gsl_vector_get(sol,4);
 pPhase->tLateC6 = gsl_vector_get(sol,5);
 pPhase->tLateC7 = gsl_vector_get(sol,6);
 pPhase->tLateC8 = gsl_vector_get(sol,7);

 printf("Int coefs: %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", pPhase->tLateC1, pPhase->tLateC2, pPhase->tLateC3, pPhase->tLateC4, pPhase->tLateC5, pPhase->tLateC6, pPhase->tLateC7, pPhase->tLateC8);

 /* Tidy up in preparation for next GSL solve ... */
 gsl_vector_free(b);
 gsl_vector_free(sol);
 gsl_matrix_free(A);
 gsl_permutation_free(p);

  /* ********************************** */
	/* *** MERGER COEFFICIENTS ********** */
	/* ********************************** */
	REAL8 invx;
  REAL8 xpoints3[5] = {xCutMerger2, xCutMerger + (xpeak - xCutMerger)/4, xCutMerger + 2*(xpeak - xCutMerger)/4, xCutMerger + 3*(xpeak - xCutMerger)/4, xpeak};

  REAL8 tmergerpoints[5];
	tmergerpoints[0] = IMRPhenomTInspiral_TofX_Calibrated_Late(xCutMerger2, wf, pPhase);
 	tmergerpoints[1] = IMRPhenomTv2_Merger_CP1(eta, S, dchi, delta);
 	tmergerpoints[2] = IMRPhenomTv2_Merger_CP2(eta, S, dchi, delta);
  tmergerpoints[3] = IMRPhenomTv2_Merger_CP3(eta, S, dchi, delta);
  tmergerpoints[4] = 0.0;

	/*Set linear system for solving merger coefficients*/
	p = gsl_permutation_alloc(6);
	b = gsl_vector_alloc(6);
	sol = gsl_vector_alloc(6);
	A = gsl_matrix_alloc(6,6);

  for (UINT4 idx=0; idx<5; idx++)
	{
    xx = xpoints3[idx];
		invx = 1./xx;

		gsl_vector_set(b,idx,tmergerpoints[idx]);

		gsl_matrix_set(A,idx,0,1.0);
		gsl_matrix_set(A,idx,1,invx);
		gsl_matrix_set(A,idx,2,invx*invx);
    gsl_matrix_set(A,idx,3,invx*invx*invx);
    gsl_matrix_set(A,idx,4,invx*invx*invx*invx);
		gsl_matrix_set(A,idx,5,invx*invx*invx*invx*invx);

	}
	invx = 1./xCutMerger2;
	gsl_vector_set(b,5,IMRPhenomTInspiral_dtdx_Calibrated_Late(xCutMerger2, pPhase)/facNT4(xCutMerger2,eta));

	gsl_matrix_set(A,5,0,0);
	gsl_matrix_set(A,5,1,-invx*invx);
	gsl_matrix_set(A,5,2,-2*invx*invx*invx);
	gsl_matrix_set(A,5,3,-3*invx*invx*invx*invx);
	gsl_matrix_set(A,5,4,-4*invx*invx*invx*invx*invx);
	gsl_matrix_set(A,5,5,-5*invx*invx*invx*invx*invx*invx);

  gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,sol);

  pPhase->cc0M = gsl_vector_get(sol,0);
	pPhase->cc1M = gsl_vector_get(sol,1);
	pPhase->cc2M = gsl_vector_get(sol,2);
  pPhase->cc3M = gsl_vector_get(sol,3);
  pPhase->cc4M = gsl_vector_get(sol,4);
	pPhase->cc5M = gsl_vector_get(sol,5);

	printf("Merger coefs: %.8f %.8f %.8f %.8f %.8f %.8f\n", pPhase->cc0M, pPhase->cc1M, pPhase->cc2M, pPhase->cc3M, pPhase->cc4M, pPhase->cc5M);

	// Free the gsl objects employed in solving the coefficient system
	gsl_vector_free(b);
	gsl_vector_free(sol);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

	/* *********************************************** */
	/* *** PHASE CONTINUITY BETWEEN REGIONS ********** */
	/* *********************************************** */

	/* Phase of the different regions correspond to the analytical integration of the frequency ansatz.
	Here we adapt the corresponding offsets for each phase ansatz in order to obtain continuity between regions */

	/* First we set the offsets to zero in order to call the phase ansatzs without the offsets */
	pPhase->phOffInsp = 0.;
	pPhase->phOffInt = 0.;
	pPhase->phOffMerger = 0.;
	pPhase->phOffRD = 0.;

	REAL8 phCutInsp = IMRPhenomTInspiral_PhiofX_Calibrated_Early(xCutInsp, wf, pPhase); //Value of inspiral phase at inspiral-merger boundary
	REAL8 phCutInt = IMRPhenomTInspiral_PhiofX_Calibrated_Late(xCutInsp, wf, pPhase); //Value of merger phase at merger-ringdown boundary

	pPhase->phOffInt = phCutInsp - phCutInt;

  REAL8 phCut2Int = IMRPhenomTInspiral_PhiofX_Calibrated_Late(xCutMerger2, wf, pPhase); //Value of inspiral phase at inspiral-merger boundary
	REAL8 phCut2Merger = IMRPhenomTMerger_PhiofX_Calibrated(xCutMerger2, pPhase); //Value of merger phase at merger-ringdown boundary

	// Needed offset for merger ansatz in order to match inspiral phase value at the boundary
  pPhase->phOffMerger = phCut2Int - phCut2Merger;

	/*REAL8 phCut3Merger = IMRPhenomTMerger_PhiofX_Calibrated(xpeak, pPhase);
	REAL8 phCutRD = 0.5*IMRPhenomTRDPhaseAnsatz22*/
	pPhase->phOffRD = 2*IMRPhenomTMerger_PhiofX_Calibrated(xpeak, pPhase); //Neded offset for ringdown ansatz in order to match merger phase value at the boundary


  /* *********************************************** */
	/* *** EPOCH AND WAVEFORM LENGHT ***************** */
	/* *********************************************** */
	REAL8 tCutT2 = IMRPhenomTInspiral_TofX_Calibrated_Early(xCutInsp, wf, pPhase);
  pPhase->tshiftT2 = -tCutT2 + IMRPhenomTInspiral_TofX_Calibrated_Late(xCutInsp, wf, pPhase);

  if(pPhase->xmin <= xCutInsp){
    pPhase->tmin = IMRPhenomTInspiral_TofX_Calibrated_Early(pPhase->xmin, wf, pPhase);
  }
  else if(pPhase->xmin > xCutMerger2){
    pPhase->tmin = IMRPhenomTMerger_TofX_Calibrated(pPhase->xmin, pPhase);
  }
  else{
    pPhase->tmin = IMRPhenomTInspiral_TofX_Calibrated_Late(pPhase->xmin, wf, pPhase);
  }


		//printf("pPhase->tmin: %.8f\n", pPhase->tmin);

    /* Now we repeat the same procedure for determining the time at which the specified reference frequency occurs. This is needed to set the reference phase */

    /* First we check if fmin and fref coincide, to save the computation */
    if(wf->fRef == wf->fmin)
    {
    	pPhase->tRef = pPhase->tmin;
    }
    /* If not, we repeat the same operation to obtain the reference time */
    else
    {
      if(pPhase->xRef <= xCutInsp){
        pPhase->tRef = IMRPhenomTInspiral_TofX_Calibrated_Early(pPhase->xRef, wf, pPhase);
      }
      else if(pPhase->tRef > xCutMerger2){
        pPhase->tRef = IMRPhenomTMerger_TofX_Calibrated(pPhase->xRef, pPhase);
      }
      else{
        pPhase->tRef = IMRPhenomTInspiral_TofX_Calibrated_Late(pPhase->xRef, wf, pPhase);
      }
    }

    pPhase->tminSec = wf->M_sec*pPhase->tmin;

    /* Required waveform length. We select maximum time of 500M since it is enough to contain the non-negligible ringdown content of the modes up to m=5.
    Length of early and late inspiral regions are computed since in the case of non-default reconstruction (4 regions) this is needed to compute frequencies
    in both regimes with the two different TaylorT3 implementations. If default reconstruction, it is harmless. */

    pPhase->wflength = floor((tEnd - pPhase->tmin)/wf->dtM);
		//printf("pPhase->wflength: %zu\n", pPhase->wflength);

	return XLAL_SUCCESS;
}

int IMRPhenomTSetHMAmplitudev2Coefficients(IMRPhenomTHMAmpStruct *pAmp, IMRPhenomTPhase22Struct *pPhase, IMRPhenomTWaveformStruct *wf)
{
	REAL8 eta   = wf->eta;   // Symmetric mass ratio
	REAL8 chi1 = wf->chi1L; // Dimensionless aligned spin of companion 1
	REAL8 chi2 = wf->chi2L; // Dimensionless aligned spin of companion 2
	REAL8 S     = wf->Shat;  // Dimensionless effective spin parameters employed for fit evaluation
	REAL8 dchi  = wf->dchi;  // Dimensionless spin difference chi1 - chi2
	REAL8 delta = wf->delta; // Mass asymmetry parameter

	pAmp->fac0 = 2*eta*sqrt(16*LAL_PI/5); // Precomputed amplitude global factor (see eq. 12-14 of PhenomTHM paper https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf)

	/* GSL objects for solving system of equations via LU decomposition */
	gsl_vector *b, *x;
	gsl_matrix *A;
	gsl_permutation *p;
	int s;

	/* For the selected mode, all needed fit values, PN coefficients and ringdown coefficients to fully determine the amplitude ansatz are computed. */

	REAL8 ampInspCP[6]; // Inspiral collocation point values.
	REAL8 ampRDC3, fDAMP, fDAMPn2, fDAMP_prec, fDAMPn2_prec, coshc3, tanhc3; // Ringdown ansatz c3 free coefficient, damping frequencies and needed quantities to compute ringdown ansatz coefficients.
	/* Ringdown ansatz coefficients c_{1,2,3,4} as defined in eq. [6-8] of Damour&Nagar 2014 (https://arxiv.org/pdf/1406.0401.pdf), also explained in eq.26 of THM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf)
	 and alpha_1, alpha_2, alpha_21 as defined in Sec IV of the same reference are stored in the amplitude struct since they are directly passed to the ringdown ansatz.
	 tshift, also passed to the amplitude struct, corresponds to the peak amplitude time of the (l,m) mode (0 for the 22 by construction). */
	gsl_complex phi; // Argument for the gsl hyperbolic functions. It is a real quantity but gsl functions expect a complex number, so it is constructed as (phi,0).

	/* The PN amplitude coefficients of the inspiral ansatz are defined in Appendix A of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf).
	   They are constructed from:
	   - 3PN expressions from Blanchet et al. 2008 (Class.Quant.Grav.25:165003,2008, https://arxiv.org/abs/0802.1249)
	   - 2PN spin corrections from Buonanno et al. 2013 (Phys. Rev D87, 044009, 2013, https://arxiv.org/abs/1209.6349)
	   - 1.5PN contributions from Arun et al. 2008 (Phys.Rev.D79:104023,2009; Erratum-ibid.D84:049901,2011, https://arxiv.org/abs/0810.5336) */

		/* Needed fits for collocation points and ringdown coefficient */
		ampInspCP[0] = IMRPhenomT_Inspiral_Amp_CP1_22v2(eta, S, dchi, delta);
		ampInspCP[1] = IMRPhenomT_Inspiral_Amp_CP2_22v2(eta, S, dchi, delta);
		ampInspCP[2] = IMRPhenomT_Inspiral_Amp_CP3_22v2(eta, S, dchi, delta);
		ampInspCP[3] = IMRPhenomT_Inspiral_Amp_CP4_22v2(eta, S, dchi, delta);
		ampInspCP[4] = IMRPhenomT_Inspiral_Amp_CP5_22v2(eta, S, dchi, delta);

		ampRDC3 = IMRPhenomT_RD_Amp_C3_22(eta, S);

		fDAMP     = evaluate_QNMfit_fdamp22(wf->afinal) / (wf->Mfinal); //damping frequency of 122 QNM
		fDAMPn2   = evaluate_QNMfit_fdamp22n2(wf->afinal) / (wf->Mfinal); //damping frequency of 222 QNM

		fDAMP_prec = evaluate_QNMfit_fdamp22(wf->afinal_prec) / (wf->Mfinal);
		fDAMPn2_prec   = evaluate_QNMfit_fdamp22n2(wf->afinal_prec) / (wf->Mfinal);

		pAmp->tshift = 0.0; //Peak time, by construction 0 for l=2, m=2

		/* PN Amplitude coefficients, defined in eq.A1 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf (before global rotation applied, see eq.13 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf) */

		pAmp->ampN = 1.0;
		pAmp->amp0halfPNreal = 0.0;
		pAmp->amp0halfPNimag = 0.0;

		pAmp->amp1PNreal = -2.5476190476190474 + (55*eta)/42.;
		pAmp->amp1PNimag = 0.0;

		pAmp->amp1halfPNreal = (-2*chi1)/3. - (2*chi2)/3. - (2*chi1*delta)/(3.*((1 - delta)/2. + (1 + delta)/2.)) + (2*chi2*delta)/(3.*((1 - delta)/2. + (1 + delta)/2.)) + (2*chi1*eta)/3. + (2*chi2*eta)/3. + 2*LAL_PI;
    	pAmp->amp1halfPNimag = 0.0;

    	pAmp->amp2PNreal = -1.437169312169312 + pow(chi1,2)/2. + pow(chi2,2)/2. + (pow(chi1,2)*delta)/2. - (pow(chi2,2)*delta)/2. - (1069*eta)/216. - pow(chi1,2)*eta + 2*chi1*chi2*eta - pow(chi2,2)*eta + (2047*pow(eta,2))/1512.;
    	pAmp->amp2PNimag = 0.0;

    	pAmp->amp2halfPNreal = - (107*LAL_PI)/21. + (34*eta*LAL_PI)/21.;
    	pAmp->amp2halfPNimag =  -24.*eta;

    	pAmp->amp3PNreal = 41.78634662956092 - (278185*eta)/33264. - (20261*pow(eta,2))/2772. + (114635*pow(eta,3))/99792. - (856*LAL_GAMMA)/105. + (2*pow(LAL_PI,2))/3. + (41*eta*pow(LAL_PI,2))/96.;
    	pAmp->amp3PNimag = (428./105)*LAL_PI;

    	pAmp->amp3halfPNreal = (-2173*LAL_PI)/756. - (2495*eta*LAL_PI)/378. + (40*pow(eta,2)*LAL_PI)/27.;
    	pAmp->amp3halfPNimag = (14333*eta)/162. - (4066*pow(eta,2))/945.;

    	pAmp->amplog = -428/105.;

	REAL8 ampPeak = ampInspCP[4]*pAmp->fac0*pPhase->xpeak;

	/***********************************************************/
	/************** RINGDOWN ANSATZ COEFFICIENTS ***************/
	/***********************************************************/

	/* Ringdown ansatz coefficients as defined in in eq. [6-8] of Damour&Nagar 2014 (https://arxiv.org/pdf/1406.0401.pdf).
	See also eq.26c-e of THM paper: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf
	Essentially, c3 is the only calibrated coefficient. c2 accounts for the effect of the first overtone in the early ringdown, and c1 and c4 fix the peak amplitude and null derivative at peak */

	pAmp->alpha1RD = 2*LAL_PI*fDAMP;
	pAmp->alpha2RD = 2*LAL_PI*fDAMPn2;
	pAmp->alpha21RD = 0.5*(pAmp->alpha2RD - pAmp->alpha1RD); //Coefficient c_2 of ringdown amplitude ansatz as defined in equation 7 of Damour&Nagar 2014 (https://arxiv.org/pdf/1406.0401.pdf) and eq.26d of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf

	pAmp->alpha1RD_prec = 2*LAL_PI*fDAMP_prec;
	pAmp->alpha2RD_prec = 2*LAL_PI*fDAMPn2_prec;
	pAmp->alpha21RD_prec = 0.5*(pAmp->alpha2RD_prec - pAmp->alpha1RD_prec);

	pAmp->c3 = ampRDC3;
	pAmp->c2 = 0.5*(pAmp->alpha2RD - pAmp->alpha1RD);
	pAmp->c2_prec = 0.5*(pAmp->alpha2RD_prec - pAmp->alpha1RD_prec);

	phi = gsl_complex_rect(pAmp->c3,0); // Needed complex parameter for gsl hyperbolic functions
	gsl_complex coshphi = gsl_complex_cosh(phi);
	coshc3 =  GSL_REAL(coshphi);
	gsl_complex tanhphi = gsl_complex_tanh(phi);
	tanhc3 =  GSL_REAL(tanhphi);

	/* This condition ensures that the second derivative of the amplitude at the mode peak is always zero or negative, not producing then a second peak in the ringdown */
	if(fabs(pAmp->c2) > fabs(0.5*pAmp->alpha1RD/tanhc3))
	{
		pAmp->c2 = -0.5*pAmp->alpha1RD/tanhc3;
	}
	if(fabs(pAmp->c2_prec) > fabs(0.5*pAmp->alpha1RD_prec/tanhc3))
	{
		pAmp->c2_prec = -0.5*pAmp->alpha1RD_prec/tanhc3;
	}


	pAmp->c1 = ampPeak*pAmp->alpha1RD*coshc3*coshc3/pAmp->c2;
	pAmp->c1_prec = ampPeak*pAmp->alpha1RD_prec*coshc3*coshc3/pAmp->c2_prec;
	pAmp->c4 = ampPeak - pAmp->c1*tanhc3;
	pAmp->c4_prec = ampPeak - pAmp->c1_prec*tanhc3;

	/***********************************************************/
	/************** INSPIRAL COEFFICIENTS SOLUTION *************/
	/***********************************************************/

	/* In order to obtain the value of the 3 unknown extra coefficients of the inspiral amplitude ansatz, we need to solve a linear system
      where the ansatz is imposed to match collocation point values at the corresponding collocation point times.
      See equation 15 of THM paper: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */

	/* Initialise the extra coefficients to zero, so calls to the amplitude ansatz function returns pure PN */
	pAmp->inspC1 = 0.0;
	pAmp->inspC2 = 0.0;
	pAmp->inspC3 = 0.0;
	pAmp->inspC4 = 0.0;
	pAmp->inspC5 = 0.0;
	pAmp->inspC6 = 0.0;

	REAL8 xcut = pPhase->xCutInsp;
	REAL8 xpeak = pPhase->xpeak;
	REAL8 xinsppoints[5]     = {xcut, xcut + (xpeak-xcut)/4, xcut + 2*(xpeak-xcut)/4, xcut + 3*(xpeak-xcut)/4, xpeak};

	/* We allocate a rank three linear system to solve the coefficients */
	p = gsl_permutation_alloc(6);
	b = gsl_vector_alloc(6);
	x = gsl_vector_alloc(6);
	A = gsl_matrix_alloc(6,6);

	REAL8 xx, x4, x4half, x5, x5half, x6, x6half; // Needed powers of PN parameter x=v^2=(\omega_orb)^(2/3)=(0.5\omega_22)^(2/3)
	REAL8 ampoffset; // Known PN part of the amplitude
	REAL8 bi; // CP value - known PN amplitude vector

	/* In this loop over collocation points, the components of the solution vector b and the basis matrix A are established.
	See equation 15 of THM paper: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */
	for (UINT4 idx=0; idx<5; idx++)
	{
		xx = xinsppoints[idx]; // PN expansion parameter of the amplitude
		x4 = xx*xx*xx*xx; // Needed powers
		x4half = x4*sqrt(xx);
		x5 = x4*xx;
		x5half = x5*sqrt(xx); // Needed powers
		x6 = x5*xx;
		x6half = x6*sqrt(xx);

		ampoffset = creal(IMRPhenomTInspiralAmpAnsatzHMv2(xx, pAmp)); // Real part of the known PN contribution
		bi = ampInspCP[idx] - (1./pAmp->fac0/xx)*(ampoffset);// Solution vector: collocation point value minus the know PN part of the ansatz, factored by the amplitude factor to not include it in each basis function (the powers of x)

		gsl_vector_set(b,idx,bi); // Set b vector

		/*Set basis matrix elements, Basis functions are the higher order powers of x that we add to the PN ansatz */
		gsl_matrix_set(A,idx,0,x4);
		gsl_matrix_set(A,idx,1,x4half);
		gsl_matrix_set(A,idx,2,x5);
		gsl_matrix_set(A,idx,3,x5half);
		gsl_matrix_set(A,idx,4,x6);
		gsl_matrix_set(A,idx,5,x6half);
	}
	bi =  1.0 + 2*pAmp->amp1PNreal*xx + 2.5*pAmp->amp1halfPNreal*xx*sqrt(xx) + 3*pAmp->amp2PNreal*xx*xx  + 3.5*pAmp->amp2halfPNreal*xx*xx*sqrt(xx) + 4*pAmp->amp3PNreal*xx*xx*xx + 4.5*pAmp->amp3halfPNreal*xx*xx*xx*sqrt(xx) + 4.0*pAmp->amplog*log(16*xx)*xx*xx*xx + pAmp->amplog*xx*xx*xx;

	gsl_vector_set(b,5,-bi);

	gsl_matrix_set(A,5,0,5*x4);
	gsl_matrix_set(A,5,1,5.5*x4half);
	gsl_matrix_set(A,5,2,6.0*x5);
	gsl_matrix_set(A,5,3,6.5*x5half);
	gsl_matrix_set(A,5,4,7.0*x6);
	gsl_matrix_set(A,5,5,7.5*x6half);

	/* We now solve the system A x = b via an LU decomposition */
	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,x);

	/* Set the extra pseudo-PN coefficients with the solutions of the system */
	pAmp->inspC1 = gsl_vector_get(x,0);
	pAmp->inspC2 = gsl_vector_get(x,1);
	pAmp->inspC3 = gsl_vector_get(x,2);
	pAmp->inspC4 = gsl_vector_get(x,3);
	pAmp->inspC5 = gsl_vector_get(x,4);
	pAmp->inspC6 = gsl_vector_get(x,5);

	/* Free the gsl solver */
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

	return XLAL_SUCCESS;
}

/* ************************************ */
        /* ANSATZ for the orbital phase */
/* ************************************ */

/* Inspiral */

/* Orbital TaylorT2 quantities */
double IMRPhenomTInspiral_dtdx_TaylorT2(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 xsqrt = sqrt(x);
	REAL8 out = 1. + x*pPhase->dT2_1PN + x*xsqrt*pPhase->dT2_1halfPN + x*x*pPhase->dT2_2PN + x*x*xsqrt*pPhase->dT2_2halfPN + x*x*x*(pPhase->dT2_3PN + log(x)*pPhase->dT2_3PNlog) + x*x*x*xsqrt*pPhase->dT2_3halfPN;

  return out;
}

/* Early inspiral */
double IMRPhenomTInspiral_dtdx_Calibrated_Early(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 xsqrt = sqrt(x);
	REAL8 out = 1. + x*pPhase->dT2_1PN + x*xsqrt*pPhase->dT2_1halfPN + x*x*pPhase->dT2_2PN + x*x*xsqrt*pPhase->dT2_2halfPN + x*x*x*(pPhase->dT2_3PN + log(x)*pPhase->dT2_3PNlog) + x*x*x*xsqrt*pPhase->dT2_3halfPN + xsqrt*x*x*x*x*pPhase->dtInspC1 + x*x*x*x*x*pPhase->dtInspC2 + xsqrt*x*x*x*x*x*pPhase->dtInspC3;

  return out;
}

double IMRPhenomTInspiral_TofX_Calibrated_Early(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase){

	REAL8 eta = wf->eta;
	REAL8 out = -0.0013020833333333333*(pow(eta,-1)*pow(x,-4)*(15 + 20*pPhase->dT2_1PN*x + 24*pPhase->dT2_1halfPN*pow(x,1.5) + 30*pPhase->dT2_2PN*pow(x,2) + 40*pPhase->dT2_2halfPN*pow(x,2.5) + 60*pPhase->dT2_3PNlog*pow(x,3) + 60*(pPhase->dT2_3PN + pPhase->dT2_3PNlog*log(x))*pow(x,3) + 120*pPhase->dT2_3halfPN*pow(x,3.5) - 120*pPhase->dtInspC1*pow(x,4.5) - 60*pPhase->dtInspC2*pow(x,5) - 40*pPhase->dtInspC3*pow(x,5.5)));

  return out + pPhase->tshiftT2;
}

double IMRPhenomTInspiral_PhiofX_Calibrated_Early(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase){

	REAL8 eta = wf->eta;
	REAL8 out = (5*pow(eta,-1)*(pPhase->dT2_3halfPN*x + pPhase->dT2_2halfPN*log(x) - (2*pow(x,-2.5))/5. - (2*pPhase->dT2_1PN*pow(x,-1.5))/3. - pPhase->dT2_1halfPN*pow(x,-1) - 2*pPhase->dT2_2PN*pow(x,-0.5) - 4*pPhase->dT2_3PNlog*pow(x,0.5) + 2*(pPhase->dT2_3PN + pPhase->dT2_3PNlog*log(x))*pow(x,0.5) + (pPhase->dtInspC1*pow(x,2))/2. + (2*pPhase->dtInspC2*pow(x,2.5))/5. + (pPhase->dtInspC3*pow(x,3))/3.))/64.;
  return out;
}

/* Late inspiral */
double IMRPhenomTInspiral_dtdx_Calibrated_Late(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 out = 1. + pPhase->dT2_1PN*x + pPhase->dT2_1halfPN*pow(x,1.5) + pPhase->dT2_2PN*pow(x,2) + pPhase->dT2_2halfPN*pow(x,2.5) + pPhase->dT2_3PN*pow(x,3) + pPhase->dT2_3PNlog*log(x)*pow(x,3) + pPhase->dT2_3halfPN*pow(x,3.5) - 0.125*pPhase->tLateC2*pow(x,4.5) - 0.25*pPhase->tLateC3*pow(x,5) - 0.375*pPhase->tLateC4*pow(x,5.5) - 0.5*pPhase->tLateC5*pow(x,6) - 0.625*pPhase->tLateC6*pow(x,6.5) - 0.75*pPhase->tLateC7*pow(x,7) - 0.875*pPhase->tLateC8*pow(x,7.5);

  return out;
}

double IMRPhenomTInspiral_TofX_Calibrated_Late(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase){

	REAL8 eta = wf->eta;
	REAL8 out = 1 + (4*pPhase->dT2_1PN*x)/3. + (8*pPhase->dT2_1halfPN*pow(x,1.5))/5. + 2*pPhase->dT2_2PN*pow(x,2) + (8*pPhase->dT2_2halfPN*pow(x,2.5))/3. + 4*pPhase->dT2_3PN*pow(x,3) + 4*pPhase->dT2_3PNlog*pow(x,3) + 4*pPhase->dT2_3PNlog*log(x)*pow(x,3) + 8*pPhase->dT2_3halfPN*pow(x,3.5) + pPhase->tLateC1*pow(x,4) + pPhase->tLateC2*pow(x,4.5) + pPhase->tLateC3*pow(x,5) +
   pPhase->tLateC4*pow(x,5.5) + pPhase->tLateC5*pow(x,6) + pPhase->tLateC6*pow(x,6.5) + pPhase->tLateC7*pow(x,7) + pPhase->tLateC8*pow(x,7.5);

  return (-5*pow(eta,-1)*pow(x,-4))/256.*out;
}

double IMRPhenomTInspiral_PhiofX_Calibrated_Late(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase){

	REAL8 eta = wf->eta;
	REAL8 out = pow(eta,-1)*pow(x,-2.5)*(-0.031249999999999997 - 0.05208333333333333*pPhase->dT2_1PN*x - 0.078125*pPhase->dT2_1halfPN*pow(x,1.5) - 0.15625*pPhase->dT2_2PN*pow(x,2) + 0.15625*pPhase->dT2_3PN*pow(x,3) - 0.3125*pPhase->dT2_3PNlog*pow(x,3) + log(x)*(0.078125*pPhase->dT2_2halfPN*pow(x,2.5) + 0.15625*pPhase->dT2_3PNlog*pow(x,3)) + 0.078125*pPhase->dT2_3halfPN*pow(x,3.5) - 0.004882812500000002*pPhase->tLateC2*pow(x,4.5) - 0.007812499999999999*pPhase->tLateC3*pow(x,5) - 0.009765625*pPhase->tLateC4*pow(x,5.5) - 0.011160714285714284*pPhase->tLateC5*pow(x,6) - 0.012207031249999998*pPhase->tLateC6*pow(x,6.5) - 0.013020833333333332*pPhase->tLateC7*pow(x,7) - 0.013671875000000002*pPhase->tLateC8*pow(x,7.5));

  return out + pPhase->phOffInt;
}

double facNT4(REAL8 x, REAL8 eta){return 64./5*x*x*x*x*x*eta;}
double facT2(REAL8 x, REAL8 eta){return (-5*pow(eta,-1)*pow(x,-4))/256.;}


/* Merger */
double IMRPhenomTMerger_TofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 out = pPhase->cc0M + pPhase->cc5M*pow(x,-5) + pPhase->cc4M*pow(x,-4) + pPhase->cc3M*pow(x,-3) + pPhase->cc2M*pow(x,-2) + pPhase->cc1M*pow(x,-1);
  return out;
}

double IMRPhenomTMerger_PhiofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 out = (pow(x,-3.5)*(50*pPhase->cc5M + 14*x*(4*pPhase->cc4M + 5*x*(pPhase->cc3M + 2*pPhase->cc2M*x - pPhase->cc1M*pow(x,2)))))/35.;
  return out + pPhase->phOffMerger;
}

/* Full */
double IMRPhenomTv2_TofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase)
{
  REAL8 tt;

  if(x <= pPhase->xCutInsp)
  {
        tt = IMRPhenomTInspiral_TofX_Calibrated_Early(x, wf, pPhase);
  }
  else if(x > pPhase->xCutMerger)
  {
        tt = IMRPhenomTMerger_TofX_Calibrated(x, pPhase);
  }
  else
  {
        tt = IMRPhenomTInspiral_TofX_Calibrated_Late(x, wf, pPhase);
  }

  return tt;
}

double IMRPhenomTv2_PhiofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase)
{
  REAL8 ph;

  if(x <= pPhase->xCutInsp)
  {
        ph = IMRPhenomTInspiral_PhiofX_Calibrated_Early(x, wf, pPhase);
  }
  else if(x > pPhase->xCutMerger)
  {
        ph = IMRPhenomTMerger_PhiofX_Calibrated(x, pPhase);
  }
  else
  {
        ph = IMRPhenomTInspiral_PhiofX_Calibrated_Late(x, wf, pPhase);
  }

  return ph;
}

int IMRPhenomTv2_times_x_phase(
	REAL8Sequence **phaseorb,       /**< Values of the 22 phase for the waveform time array */
	REAL8Sequence **xorb,				/**< Values of the 22 frequency for the waveform time array */
	REAL8Sequence **times,
	IMRPhenomTWaveformStruct *wf,
	IMRPhenomTPhase22Struct *pPhase
){
	REAL8 dx = 0.00005;
	REAL8 dt = 1.0;
	REAL8 xx, tt, ph22;

	size_t insp_length = floor((pPhase->xpeak - pPhase->xmin)/dx)+1;
  size_t rd_length = floor((tEnd)/dt)+3;
	size_t length_coarse = insp_length + rd_length;

	*times = XLALCreateREAL8Sequence(length_coarse);
	*xorb = XLALCreateREAL8Sequence(length_coarse);
	*phaseorb = XLALCreateREAL8Sequence(length_coarse);

	// INSPIRAL
	for(UINT4 jdx = 0; jdx < insp_length; jdx++){

		xx = pPhase->xmin + jdx*dx;
		(*xorb)->data[jdx] = xx;
		(*times)->data[jdx] = IMRPhenomTv2_TofX_Calibrated(xx, wf, pPhase);
		(*phaseorb)->data[jdx] = IMRPhenomTv2_PhiofX_Calibrated(xx, wf, pPhase);
	}

  for(UINT4 jdx = 0; jdx < rd_length; jdx++){

		tt = jdx*dt;
		(*times)->data[insp_length +jdx] = tt;

    (*xorb)->data[insp_length + jdx] = (*xorb)->data[insp_length - 1]+0.00001*jdx;
    ph22 = IMRPhenomTRDPhaseAnsatz22(tt, pPhase);
    (*phaseorb)->data[insp_length + jdx] = 0.5*ph22;
	}

	return XLAL_SUCCESS;
}

/* Amplitudes */
COMPLEX16 IMRPhenomTInspiralAmpAnsatzHMv2(REAL8 x, IMRPhenomTHMAmpStruct *pAmp)
{
	REAL8 fac = pAmp->fac0*x;

	REAL8 xhalf = sqrt(x);
	REAL8 x1half = x*xhalf;
	REAL8 x2 = x*x;
	REAL8 x2half = x2*xhalf;
	REAL8 x3 = x2*x;
	REAL8 x3half = x3*xhalf;
	REAL8 x4 = x2*x2;
	REAL8 x4half = x4*xhalf;
	REAL8 x5 = x3*x2;
	REAL8 x5half = x5*xhalf;
	REAL8 x6 = x4*x2;
	REAL8 x6half = x6*xhalf;

	REAL8 ampreal = pAmp->ampN + pAmp->amp0halfPNreal*xhalf + pAmp->amp1PNreal*x + pAmp->amp1halfPNreal*x1half + pAmp->amp2PNreal*x2  + pAmp->amp2halfPNreal*x2half + pAmp->amp3PNreal*x3 + pAmp->amp3halfPNreal*x3half + pAmp->amplog*log(16*x)*x3;
	REAL8 ampimag = pAmp->amp0halfPNimag*xhalf + pAmp->amp1PNimag*x + pAmp->amp1halfPNimag*x1half + pAmp->amp2PNimag*x2  + pAmp->amp2halfPNimag*x2half + pAmp->amp3PNimag*x3 + pAmp->amp3halfPNimag*x3half;
	COMPLEX16 amp = crect(ampreal + pAmp->inspC1*x4 + pAmp->inspC2*x4half + pAmp->inspC3*x5 + pAmp->inspC4*x5half + pAmp->inspC5*x6 + pAmp->inspC6*x6half, ampimag);

	return fac*amp;
}

COMPLEX16 IMRPhenomTHMAmpv2(
  REAL8 t,
  UNUSED REAL8 x,
  IMRPhenomTHMAmpStruct *pAmp
)
{
    COMPLEX16 amp;

    if(t <= 0.0)
        {
          amp = IMRPhenomTInspiralAmpAnsatzHMv2(x, pAmp);
        }
        else
        {
          amp = IMRPhenomTRDAmpAnsatzHM(t, pAmp);
        }

    return amp;
}
