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

/* Inspiral ansatz */

double IMRPhenomTInspiral_dtdx_TaylorT2(REAL8 x, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiral_dtdx_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiral_TofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTInspiral_PhiofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase);
double facNT4(REAL8 x, REAL8 eta);

/* Intermediate ansatz */
double IMRPhenomTIntermediate_TofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase);
double IMRPhenomTIntermediate_PhiofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase);
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

	REAL8 xmeco = pow(LAL_PI*XLALSimIMRPhenomXfMECO(eta,chi1L,chi2L),2./3);
  pPhase->xmeco = xmeco;
  REAL8 xpeak = IMRPhenomTv2_xPeak(eta, S, dchi, delta);
  pPhase->xpeak = xpeak;
	printf("---------------\n");
	printf("CASE: eta: %f chi1L: %f chi2L: %f\n", eta, chi1L, chi2L);
	printf("xmeco: %.8f xpeak: %.8f\n", xmeco, xpeak);

  REAL8 xCutInsp = 0.9*xmeco;
  REAL8 xCutMerger = 1.2*xmeco;

  pPhase->xCutInsp = xCutInsp;
  pPhase->xCutMerger = xCutMerger;

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
	REAL8 dtinsppoints[4];
	dtinsppoints[0] = IMRPhenomTv2_Inspiral_CP1(eta, S, dchi, delta);
 	dtinsppoints[1] = IMRPhenomTv2_Inspiral_CP2(eta, S, dchi, delta);
 	dtinsppoints[2] = IMRPhenomTv2_Inspiral_CP3(eta, S, dchi, delta);
  dtinsppoints[3] = IMRPhenomTv2_Inspiral_CP4(eta, S, dchi, delta);

  /* Initialize higher order extension coefficients to 0. In this way, when calling the ansatz IMRPhenomTInspiralOmegaAnsatz22 it will return PN TaylorT3 value. */
	pPhase->dtInspC1 = 0.0;
	pPhase->dtInspC2 = 0.0;
	pPhase->dtInspC3 = 0.0;
  pPhase->dtInspC4 = 0.0;

 /* Set collocation points */

  REAL8 xpoints[4] = {0.045, 0.045 + (xCutInsp - 0.045)/3, 0.045 + 2*(xCutInsp - 0.045)/3, xCutInsp};

	/*Set linear system, which is rank 6 */
	p = gsl_permutation_alloc(4);
	b = gsl_vector_alloc(4);
	sol = gsl_vector_alloc(4);
	A = gsl_matrix_alloc(4,4);

	/*Set A matrix and b vector*/

	REAL8 xx, xx9over2, xx5, xx11over2, xx6, dT2offset; // Initialize theta powers the diagonal A matrix
	/* theta is a weighted dimensionless time parameter defined in Eq. 315 of Blanchet 2014 (https://arxiv.org/abs/1310.1528).
	   Notice however that here is defined at the (-1/8) power, in order to represent the Post-Newtonian order (1/c). */

	/* Set up inspiral coefficient system.
	This system of equations is explained in eq. 10 of THM paper https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */
	for (UINT4 idx=0; idx<4; idx++)
	{
		/* Needed powers of x */
		xx = xpoints[idx];
		xx9over2 = pow(xx,4.5);
		xx5 = pow(xx,5);
		xx11over2 = pow(xx,5.5);
    xx6 = pow(xx,6);

		dT2offset = IMRPhenomTInspiral_dtdx_TaylorT2(xx, pPhase);

		gsl_vector_set(b,idx,(dtinsppoints[idx] - dT2offset));

		gsl_matrix_set(A,idx,0,xx9over2);
		gsl_matrix_set(A,idx,1,xx5);
		gsl_matrix_set(A,idx,2,xx11over2);
    gsl_matrix_set(A,idx,3,xx6);

	}

	/* We now solve the system A x = b via an LU decomposition */
	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,sol);

	/* Set inspiral phenomenological coefficients from solution to A x = b */
	pPhase->dtInspC1 = gsl_vector_get(sol,0);
	pPhase->dtInspC2 = gsl_vector_get(sol,1);
	pPhase->dtInspC3 = gsl_vector_get(sol,2);
  pPhase->dtInspC4 = gsl_vector_get(sol,3);

	printf("Insp coefs: %.8f %.8f %.8f %.8f\n", pPhase->dtInspC1, pPhase->dtInspC2, pPhase->dtInspC3, pPhase->dtInspC4);

	/* Boundary time between late inspiral and merger regions.
	This is selected to correspond to theta=0.81, earlier than the last collocation point at theta=0.82, because in this way
	the derivative at the boundary with the merger region, needed for the merger reconstruction, is less forced. */
	pPhase->tshiftT2 = 0.0;
	REAL8 tCutT2 = IMRPhenomTInspiral_TofX_Calibrated(xCutInsp, wf, pPhase);

	/* Tidy up in preparation for next GSL solve ... */
	gsl_vector_free(b);
	gsl_vector_free(sol);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

  /* ********************************** */
	/* *** INTERMEDIATE COEFFICIENTS ********** */
	/* ********************************** */
  REAL8 invx;
  REAL8 xpoints2[6] = {xCutInsp, xCutInsp + (xCutMerger - xCutInsp)/5, xCutInsp + 2*(xCutMerger - xCutInsp)/5, xCutInsp + 3*(xCutMerger - xCutInsp)/5, xCutInsp + 4*(xCutMerger - xCutInsp)/5, xCutMerger};

  REAL8 tintpoints[6];
	tintpoints[0] = IMRPhenomTv2_Intermediate_CP1(eta, S, dchi, delta);
 	tintpoints[1] = IMRPhenomTv2_Intermediate_CP2(eta, S, dchi, delta);
 	tintpoints[2] = IMRPhenomTv2_Intermediate_CP3(eta, S, dchi, delta);
  tintpoints[3] = IMRPhenomTv2_Intermediate_CP4(eta, S, dchi, delta);
  tintpoints[4] = IMRPhenomTv2_Intermediate_CP5(eta, S, dchi, delta);
  tintpoints[5] = IMRPhenomTv2_Intermediate_CP6(eta, S, dchi, delta);

	/*Set linear system for solving merger coefficients*/
	p = gsl_permutation_alloc(6);
	b = gsl_vector_alloc(6);
	sol = gsl_vector_alloc(6);
	A = gsl_matrix_alloc(6,6);

  for (UINT4 idx=0; idx<6; idx++)
	{
    xx = xpoints2[idx];
		invx = 1./xx;

		gsl_vector_set(b,idx,tintpoints[idx]);

		gsl_matrix_set(A,idx,0,1.0);
		gsl_matrix_set(A,idx,1,invx);
		gsl_matrix_set(A,idx,2,invx*invx);
    gsl_matrix_set(A,idx,3,invx*invx*invx);
    gsl_matrix_set(A,idx,4,invx*invx*invx*invx);
    gsl_matrix_set(A,idx,5,invx*invx*invx*invx*invx);

	}
  gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,sol);

  pPhase->bb0Int = gsl_vector_get(sol,0);
	pPhase->bb1Int = gsl_vector_get(sol,1);
	pPhase->bb2Int = gsl_vector_get(sol,2);
  pPhase->bb3Int = gsl_vector_get(sol,3);
  pPhase->bb4Int = gsl_vector_get(sol,4);
  pPhase->bb5Int = gsl_vector_get(sol,5);

	printf("Int coefs: %.8f %.8f %.8f %.8f %.8f %.8f\n", pPhase->bb0Int, pPhase->bb1Int, pPhase->bb2Int, pPhase->bb3Int, pPhase->bb4Int, pPhase->bb5Int);

	// Free the gsl objects employed in solving the coefficient system
	gsl_vector_free(b);
	gsl_vector_free(sol);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

  /* ********************************** */
	/* *** INTERMEDIATE COEFFICIENTS ********** */
	/* ********************************** */

  REAL8 xpoints3[5] = {xCutMerger, xCutMerger + (xpeak - xCutMerger)/4, xCutMerger + 2*(xpeak - xCutMerger)/4, xCutMerger + 3*(xpeak - xCutMerger)/4, xpeak};

  REAL8 tmergerpoints[5];
	tmergerpoints[0] = IMRPhenomTv2_Intermediate_CP6(eta, S, dchi, delta);
 	tmergerpoints[1] = IMRPhenomTv2_Merger_CP7(eta, S, dchi, delta);
 	tmergerpoints[2] = IMRPhenomTv2_Merger_CP8(eta, S, dchi, delta);
  tmergerpoints[3] = IMRPhenomTv2_Merger_CP9(eta, S, dchi, delta);
  tmergerpoints[4] = 0.0;

	/*Set linear system for solving merger coefficients*/
	p = gsl_permutation_alloc(5);
	b = gsl_vector_alloc(5);
	sol = gsl_vector_alloc(5);
	A = gsl_matrix_alloc(5,5);

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

	}
  gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,sol);

  pPhase->cc0M = gsl_vector_get(sol,0);
	pPhase->cc1M = gsl_vector_get(sol,1);
	pPhase->cc2M = gsl_vector_get(sol,2);
  pPhase->cc3M = gsl_vector_get(sol,3);
  pPhase->cc4M = gsl_vector_get(sol,4);

	printf("Int coefs: %.8f %.8f %.8f %.8f %.8f\n", pPhase->cc0M, pPhase->cc1M, pPhase->cc2M, pPhase->cc3M, pPhase->cc4M);

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

	REAL8 phCutInsp = IMRPhenomTInspiral_PhiofX_Calibrated(xCutInsp, wf, pPhase); //Value of inspiral phase at inspiral-merger boundary
	REAL8 phCutInt = IMRPhenomTIntermediate_PhiofX_Calibrated(xCutInsp, pPhase); //Value of merger phase at merger-ringdown boundary

	pPhase->phOffInt = phCutInsp - phCutInt;

  REAL8 phCut2Int = IMRPhenomTIntermediate_PhiofX_Calibrated(xCutMerger, pPhase); //Value of inspiral phase at inspiral-merger boundary
	REAL8 phCut2Merger = IMRPhenomTMerger_PhiofX_Calibrated(xCutMerger, pPhase); //Value of merger phase at merger-ringdown boundary

	// Needed offset for merger ansatz in order to match inspiral phase value at the boundary
  pPhase->phOffMerger = phCut2Int - phCut2Merger;

	/*REAL8 phCut3Merger = IMRPhenomTMerger_PhiofX_Calibrated(xpeak, pPhase);
	REAL8 phCutRD = 0.5*IMRPhenomTRDPhaseAnsatz22*/
	pPhase->phOffRD = 2*IMRPhenomTMerger_PhiofX_Calibrated(xpeak, pPhase); //Neded offset for ringdown ansatz in order to match merger phase value at the boundary


  /* *********************************************** */
	/* *** EPOCH AND WAVEFORM LENGHT ***************** */
	/* *********************************************** */

  pPhase->tshiftT2 = -tCutT2 + IMRPhenomTIntermediate_TofX_Calibrated(xCutInsp, pPhase);

  if(pPhase->xmin <= xCutInsp){
    pPhase->tmin = IMRPhenomTInspiral_TofX_Calibrated(pPhase->xmin, wf, pPhase);
  }
  else if(pPhase->xmin > xCutMerger){
    pPhase->tmin = IMRPhenomTMerger_TofX_Calibrated(pPhase->xmin, pPhase);
  }
  else{
    pPhase->tmin = IMRPhenomTIntermediate_TofX_Calibrated(pPhase->xmin, pPhase);
  }


		printf("pPhase->tmin: %.8f\n", pPhase->tmin);

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
        pPhase->tRef = IMRPhenomTInspiral_TofX_Calibrated(pPhase->xRef, wf, pPhase);
      }
      else if(pPhase->tRef > xCutMerger){
        pPhase->tRef = IMRPhenomTMerger_TofX_Calibrated(pPhase->xRef, pPhase);
      }
      else{
        pPhase->tRef = IMRPhenomTIntermediate_TofX_Calibrated(pPhase->xRef, pPhase);
      }
    }

    pPhase->tminSec = wf->M_sec*pPhase->tmin;

    /* Required waveform length. We select maximum time of 500M since it is enough to contain the non-negligible ringdown content of the modes up to m=5.
    Length of early and late inspiral regions are computed since in the case of non-default reconstruction (4 regions) this is needed to compute frequencies
    in both regimes with the two different TaylorT3 implementations. If default reconstruction, it is harmless. */

    pPhase->wflength = floor((tEnd - pPhase->tmin)/wf->dtM);
		printf("pPhase->wflength: %zu\n", pPhase->wflength);

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

double IMRPhenomTInspiral_dtdx_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 xsqrt = sqrt(x);
	REAL8 out = 1. + x*pPhase->dT2_1PN + x*xsqrt*pPhase->dT2_1halfPN + x*x*pPhase->dT2_2PN + x*x*xsqrt*pPhase->dT2_2halfPN + x*x*x*(pPhase->dT2_3PN + log(x)*pPhase->dT2_3PNlog) + x*x*x*xsqrt*pPhase->dT2_3halfPN + xsqrt*x*x*x*x*pPhase->dtInspC1 + x*x*x*x*x*pPhase->dtInspC2 + xsqrt*x*x*x*x*x*pPhase->dtInspC3 + x*x*x*x*x*x*pPhase->dtInspC4;

  return out;
}

double IMRPhenomTInspiral_TofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase){

	REAL8 eta = wf->eta;
	REAL8 out = -0.0013020833333333333*(pow(eta,-1)*pow(x,-4)*(15 + 20*pPhase->dT2_1PN*x + 24*pPhase->dT2_1halfPN*pow(x,1.5) + 30*pPhase->dT2_2PN*pow(x,2) + 40*pPhase->dT2_2halfPN*pow(x,2.5) + 60*pPhase->dT2_3PNlog*pow(x,3) + 60*(pPhase->dT2_3PN + pPhase->dT2_3PNlog*log(x))*pow(x,3) + 120*pPhase->dT2_3halfPN*pow(x,3.5) - 120*pPhase->dtInspC1*pow(x,4.5) - 60*pPhase->dtInspC2*pow(x,5) - 40*pPhase->dtInspC3*pow(x,5.5) - 30*pPhase->dtInspC4*pow(x,6)));

  return out + pPhase->tshiftT2;
}

double IMRPhenomTInspiral_PhiofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase){

	REAL8 eta = wf->eta;
	REAL8 out = (5*pow(eta,-1)*(pPhase->dT2_3halfPN*x + pPhase->dT2_2halfPN*log(x) - (2*pow(x,-2.5))/5. - (2*pPhase->dT2_1PN*pow(x,-1.5))/3. - pPhase->dT2_1halfPN*pow(x,-1) - 2*pPhase->dT2_2PN*pow(x,-0.5) - 4*pPhase->dT2_3PNlog*pow(x,0.5) + 2*(pPhase->dT2_3PN + pPhase->dT2_3PNlog*log(x))*pow(x,0.5) + (pPhase->dtInspC1*pow(x,2))/2. + (2*pPhase->dtInspC2*pow(x,2.5))/5. + (pPhase->dtInspC3*pow(x,3))/3. + (2*pPhase->dtInspC4*pow(x,3.5))/7.))/64.;
  return out;
}

double facNT4(REAL8 x, REAL8 eta){return 64./5*x*x*x*x*x*eta;}

/* Intermediate */

double IMRPhenomTIntermediate_TofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 out = pPhase->bb0Int + pPhase->bb5Int*pow(x,-5) + pPhase->bb4Int*pow(x,-4) + pPhase->bb3Int*pow(x,-3) + pPhase->bb2Int*pow(x,-2) + pPhase->bb1Int*pow(x,-1);
  return out;
}

double IMRPhenomTIntermediate_PhiofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 out = (pow(x,-3.5)*(50*pPhase->bb5Int + 14*x*(4*pPhase->bb4Int + 5*x*(pPhase->bb3Int + 2*pPhase->bb2Int*x - pPhase->bb1Int*pow(x,2)))))/35.;
  return out + pPhase->phOffInt;
}

/* Merger */

double IMRPhenomTMerger_TofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 out = pPhase->cc0M + pPhase->cc4M*pow(x,-4) + pPhase->cc3M*pow(x,-3) + pPhase->cc2M*pow(x,-2) + pPhase->cc1M*pow(x,-1);
  return out;
}

double IMRPhenomTMerger_PhiofX_Calibrated(REAL8 x, IMRPhenomTPhase22Struct *pPhase){

	REAL8 out = ((8*pPhase->cc4M + 10*x*(pPhase->cc3M + x*(2*pPhase->cc2M - pPhase->cc1M*x)))*pow(x,-2.5))/5.;
  return out + pPhase->phOffMerger;
}

double IMRPhenomTv2_TofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase)
{
  REAL8 tt;

  if(x <= pPhase->xCutInsp)
  {
        tt = IMRPhenomTInspiral_TofX_Calibrated(x, wf, pPhase);
  }
  else if(x > pPhase->xCutMerger)
  {
        tt = IMRPhenomTMerger_TofX_Calibrated(x, pPhase);
  }
  else
  {
        tt = IMRPhenomTIntermediate_TofX_Calibrated(x, pPhase);
  }

  return tt;
}

double IMRPhenomTv2_PhiofX_Calibrated(REAL8 x, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase)
{
  REAL8 ph;

  if(x <= pPhase->xCutInsp)
  {
        ph = IMRPhenomTInspiral_PhiofX_Calibrated(x, wf, pPhase);
  }
  else if(x > pPhase->xCutMerger)
  {
        ph = IMRPhenomTMerger_PhiofX_Calibrated(x, pPhase);
  }
  else
  {
        ph = IMRPhenomTIntermediate_PhiofX_Calibrated(x, pPhase);
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
	REAL8 dx = 0.0001;
	REAL8 dt = 1.0;
	REAL8 xx, tt, ph22;

	size_t insp_length = floor((pPhase->xpeak - pPhase->xmin)/dx);
  size_t rd_length = floor((tEnd)/dt) + 3;
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

    (*xorb)->data[insp_length + jdx] = (*xorb)->data[insp_length - 1]+0.0001*jdx;
    ph22 = IMRPhenomTRDPhaseAnsatz22(tt, pPhase);
    (*phaseorb)->data[insp_length + jdx] = 0.5*ph22;
	}

	return XLAL_SUCCESS;
}
