#define PHENOMXDEBUG 0


#include "LALSimIMRPhenomX_shared.h"


void IMRPhenomX_FillArray(int n, float *x, float *y)
{
  for (int i = 0; i < n; i++)
      y[i] = x[i] + y[i];
}



int IMRPhenomX_Frequency_Loop(COMPLEX16FrequencySeries **htilde22,
                     UNUSED  REAL8Sequence *freqs,
                         IMRPhenomXWaveformStruct *pWF,
                         IMRPhenomXAmpCoefficients *pAmp22,
                         IMRPhenomXPhaseCoefficients *pPhase22,
                         UINT4 offset,
                         UINT4 N)
{

  printf("\nInside PhXFreqLoop CUDA file");

  #ifdef __cplusplus
  printf("\nC++\n");
  #else
  printf("\nNo C++\n");
  #endif

  #ifdef __NVCC__
  printf("NVCC defined\n");
  #else
  printf("NVCC NOT defined\n");
  #endif

  #if defined(CUDA)
  printf("CUDA defined\n");
  #else
  printf("CUDA NOT defined\n");
  #endif

  #if defined(LALSIMULATION_CUDA_ENABLED)
  printf("LALSIMULATION_CUDA_ENABLED defined\n");
  #else
  printf("LALSIMULATION_CUDA_ENABLED NOT defined\n");
  #endif

  N = 100; //freqs->length;
  float *x = NULL;// *y = NULL;

  // Allocate Unified Memory – accessible from CPU or GPU
  cudaMallocManaged((void**)&x, N*sizeof(float), cudaMemAttachGlobal);
  //cudaMallocManaged(&y, N*sizeof(float));

  // initialize x and y arrays on the host
 // for (int i = 0; i < N; i++) {
 //   x[i] = 1.0f;
 //   y[i] = 2.0f;
 // }

  // Run kernel on 1M elements on the GPU
  //IMRPhenomX_FillArray<<<1, 1>>>(N, x, y);

  // Wait for GPU to finish before accessing on host
  //cudaDeviceSynchronize();

  for (UINT4 idx = 0; idx < N; idx++)
  {
	  /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
      ((*htilde22)->data->data)[idx] = 0.*(pWF->amp0 + pAmp22->fAmpMatchIN + pPhase22->C1Int + offset);
  }


  // Free memory
  //cudaFree(x);
  //cudaFree(y);

  return XLAL_SUCCESS;
}

// int IMRPhenomX_Frequency_Loop(COMPLEX16FrequencySeries **htilde22,
//                          REAL8Sequence *freqs,
//                          IMRPhenomXWaveformStruct *pWF,
//                          IMRPhenomXAmpCoefficients *pAmp22,
//                          IMRPhenomXPhaseCoefficients *pPhase22,
//                          UINT4 offset)
// {
//
//   /* initial_status used to track  */
//   UINT4 initial_status = XLAL_SUCCESS;
//   int status = initial_status;
//
//   /* Linear time and phase shifts so that model peaks near t ~ 0 */
//   REAL8 lina = 0;
//    double linb=IMRPhenomX_TimeShift_22(pPhase22, pWF);
//
//   /* 1/eta is used to re-scale phase */
//   REAL8 inveta    = (1.0 / pWF->eta);
//
//   /*
//       Here we declare explicit REAL8 variables for main loop in order to avoid numerous
//       pointer calls.
//   */
//   //REAL8 MfRef     = pWF->MfRef;
//   REAL8 Msec      = pWF->M_sec;
//
//   REAL8 C1IM      = pPhase22->C1Int;
//   REAL8 C2IM      = pPhase22->C2Int;
//   REAL8 C1RD      = pPhase22->C1MRD;
//   REAL8 C2RD      = pPhase22->C2MRD;
//
//   REAL8 fPhaseIN  = pPhase22->fPhaseMatchIN;
//   REAL8 fPhaseIM  = pPhase22->fPhaseMatchIM;
//   REAL8 fAmpIN    = pAmp22->fAmpMatchIN;
//   REAL8 fAmpIM    = pAmp22->fAmpRDMin;
//
//   REAL8 Amp0      = pWF->amp0 * pWF->ampNorm;
//
//
//
//   for (UINT4 idx = 0; idx < freqs->length; idx++)
//   {
//     double Mf    = Msec * freqs->data[idx];   // Mf is declared locally inside the loop
//     UINT4 jdx    = idx  + offset;             // jdx is declared locally inside the loop
//
//     /* Initialize a struct containing useful powers of Mf */
//     IMRPhenomX_UsefulPowers powers_of_Mf;
//     initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
//     if(initial_status != XLAL_SUCCESS)
//     {
//       status = initial_status;
//       XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
//     }
//     else
//     {
//       /* Generate amplitude and phase at MfRef */
//       REAL8 amp = 0.0;
//       REAL8 phi = 0.0;
//
//       /* The functions in this routine are inlined to help performance. */
//       /* Construct phase */
//       if(Mf < fPhaseIN)
//       {
//         phi = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pPhase22);
//       }
//       else if(Mf > fPhaseIM)
//       {
//         phi = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1RD + (C2RD * Mf);
//       }
//       else
//       {
//         phi = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1IM + (C2IM * Mf);
//       }
//
// 	  /* Scale phase by 1/eta */
// 	  phi  *= inveta;
//       phi  += linb*Mf + lina + pWF->phifRef;
//
// 	  /* Construct amplitude */
// 	  if(Mf < fAmpIN)
// 	  {
// 		  amp = IMRPhenomX_Inspiral_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
// 	  }
// 	  else if(Mf > fAmpIM)
// 	  {
// 		  amp = IMRPhenomX_Ringdown_Amp_22_Ansatz(Mf, pWF, pAmp22);
// 	  }
// 	  else
// 	  {
//         amp = IMRPhenomX_Intermediate_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
//       }
//
// 	  /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
//       ((*htilde22)->data->data)[jdx] = Amp0 * 0. * powers_of_Mf.m_seven_sixths * amp * cexp(I * phi);
//     }
//   }
//
//     return status;
// }
