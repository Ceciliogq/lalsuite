#ifndef _LALSIM_IMR_PHENOMX_SHARED_H
#define _LALSIM_IMR_PHENOMX_SHARED_H

#ifdef __cplusplus
extern "C" {
#endif
    
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"

void IMRPhenomX_FillArray(int n, float *x, float *y);
    
int IMRPhenomX_Frequency_Loop(COMPLEX16FrequencySeries **htilde22, 
                     REAL8Sequence *freqs, 
                     IMRPhenomXWaveformStruct *pWF, 
                     IMRPhenomXAmpCoefficients *pAmp22, 
                     IMRPhenomXPhaseCoefficients *pPhase22,
                     UINT4 offset);


#ifdef __cplusplus
}
#endif

#endif
