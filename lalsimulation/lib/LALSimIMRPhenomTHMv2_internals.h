#ifndef _LALSIM_IMR_PHENOMTHMV2_INTERNALS_H
#define _LALSIM_IMR_PHENOMTHMV2_INTERNALS_H

/*
 * Copyright (C) 2020 Hector Estelles
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


/*
 * \author Hector Estelles
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

/* Standard libraries */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* LAL */
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>

/* IMRPhenomT */
#include <lal/LALSimInspiral.h>

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

COMPLEX16 IMRPhenomTHMAmpv2(REAL8 t, REAL8 x,IMRPhenomTHMAmpStruct *pAmp);

COMPLEX16 IMRPhenomTInspiralAmpAnsatzHMv2(REAL8 x, IMRPhenomTHMAmpStruct *pAmp);

#ifdef __cplusplus
}
#endif

#endif
