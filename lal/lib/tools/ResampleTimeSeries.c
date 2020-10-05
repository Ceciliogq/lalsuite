/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Patrick Brady
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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/ResampleTimeSeries.h>

#if __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * \defgroup ResampleTimeSeries_c Module ResampleTimeSeries.c
 * \ingroup ResampleTimeSeries_h
 *
 * \author Brown, D. A., Brady, P. R., Charlton, P.
 *
 * \brief Downsamples a time series in place by an integer power of two.
 *
 * The routine LALResampleREAL4TimeSeries() provided functionality to
 * downsample a time series in place by an integer factor which is a power of
 * two. Upsampling, non-integer resampling and resampling by a factor which is
 * not a power of two are currently unsupported. Attempts to use this methods
 * will cause the function to abort with an error message.
 *
 * On entry the input time series \c ts should contain the data to be
 * resampled, with the data, length and sample interbal of the time series
 * populated. The parameter structure \c params should have the value of
 * \c deltaT set to the desired value sample interval of the output data.
 * The type of filter used to perform the low pass filter should be set in
 * the parameter structure field \c filterType. The parameter structure
 * field \c filterParams may be ignored at present. It is designed to
 * allow user defined low pass filters and is not currently implemented. The
 * resampling function will behave correctly even if this contains garbage on
 * entry.
 *
 * On exit, the time series \c ts will contain the resampled time series.
 * The length of the time series will be reduced by the resampling ratio and
 * the sample interval will be set correctly.
 *
 * <em>There is no time shift in the output time series</em> when either the
 * #defaultButterworth or #LDASfirLP low pass filter types are
 * used. The timestamp of each point in the output time series is the same as
 * the timestamp of the corresponding point in input time series.
 *
 * <em>There will be corrupted data at the start and end of the time series</em>
 * when either the #defaultButterworth or #LDASfirLP low pass filter
 * types are used. This is caused by corruption due to the low pass filter.
 * Users should take care to truncate these points from the time series before
 * using the output in subsequent filtering applications.
 *
 * ### Algorithm ###
 *
 * The input time series is first low passed to remove any power above the new
 * Nyquist frequency. There are two available low pass filters:
 * <ol>
 * <li> #defaultButterworth: The input time series has a time domain low
 * pass filter applied by the LALDButterworthREAL4TimeSeries() function
 * from the tdfilters package. The filter order is \f$20\f$ and the attenuation is
 * \f$0.1\f$ at the new Nyquist frequency. Since the butterworth filter is applied
 * forwards and backwards there is no time delay or phase shift of the output
 * data. The disadvantage of the butterworth filter is that since it is an IIR
 * filter, it is not possible to determine the exact length of corrupted data at
 * the start and end of the output time series. Care should be taked to avoid
 * using these regions of data in subsequent filtering.</li>
 *
 * <li> #LDASfirLP: The input time series has a time domain low
 * pass filter applied by the LALDIIRFilterREAL4Vector() function
 * from the tdfilters package. This applies an FIR filter with coefficents
 * generated by the LDAS <tt>firlp()</tt> dataconditioning API action. FIR
 * coefficents are available for downsampling by a factor of 2, 4 or 8. An
 * attempt to downsample by any other ratio will cause an error. The FIR
 * coeffients are equivalent to the default <tt>resample()</tt> action in the
 * dataconditioning API.
 *
 * The FIR filter causes the a time delay of 10 points in the output time series
 * and corruption of the first \f$n\f$ points. \f$n\f$ is given by the order of the
 * filter by \f$n = 2 \times 10 \times q\f$, where \f$q\f$ is the resampling ratio. To
 * account for this, we do the following:
 * <ol>
 * <li> The first \f$n/2\f$ points of the time series are deleted.</li>
 * <li> The whole time series is shifted \f$n/2\f$ points "left" to remove the
 * time delay.</li>
 * <li> The first and last \f$n/2\f$ points of the output time series are set to
 * zero to make sure that the corrupted data contains sensible numbers and not
 * garbage or \c Inf.</li>
 * </ol>
 * This means that there is no time delay in the output time series, but the
 * first and last \f$n/2q\f$ points are set to zero. Care should be taked to avoid
 * using these regions of data in subsequent filtering. If the debug level is
 * set to LALWARNING, a message is printed reporting how many points at the
 * start and end of time series are corrupted.
 *
 * The filter coefficents used were produced by LDAS-CIT running version 0.7.0 of
 * LDAS. See the LDAS dataconditioning API documentation for more information.
 * </ol>
 *
 */
/** @{ */

/** \see See \ref ResampleTimeSeries_c for documentation */
int XLALResampleREAL4TimeSeries( REAL4TimeSeries *series, REAL8 dt )
{
  const INT4 filterOrder = 20;
  const REAL8 newNyquistAmplitude = 0.1;
  REAL8 newNyquistFrequency;
  UINT4 resampleFactor;
  REAL4 *dataPtr = NULL;
  UINT4 j;

  resampleFactor = floor( dt / series->deltaT + 0.5 );
  newNyquistFrequency = 0.5 / dt;

  /* check that the resampling factor is valid */
  if ( resampleFactor < 1 ||
      fabs( dt - resampleFactor * series->deltaT ) > 1e-3 * series->deltaT )
    XLAL_ERROR( XLAL_EINVAL );

  /* just return if no resampling is required */
  if ( resampleFactor == 1 )
  {
    XLALPrintInfo( "XLAL Info - %s: No resampling required", __func__ );
    return 0;
  }

  /* check that we are resampling by a power of two */
  if ( resampleFactor & (resampleFactor - 1) )
    XLAL_ERROR( XLAL_EINVAL );

  if ( XLALLowPassREAL4TimeSeries( series, newNyquistFrequency,
        newNyquistAmplitude, filterOrder ) < 0 )
    XLAL_ERROR( XLAL_EFUNC );

  /* decimate the time series */
  series->deltaT = dt;
  series->data->length /= resampleFactor;
  dataPtr = series->data->data;
  for ( j = 0; j < series->data->length; ++j )
  {
    series->data->data[j] = *dataPtr;
    dataPtr += resampleFactor;
  }
  series->data->data = LALRealloc( series->data->data,
      series->data->length * sizeof( *series->data->data ) );

  return 0;
}

/** \see See \ref ResampleTimeSeries_c for documentation */
int XLALResampleREAL8TimeSeries( REAL8TimeSeries *series, REAL8 dt )
{
  const INT4 filterOrder = 20;
  const REAL8 newNyquistAmplitude = 0.1;
  REAL8 newNyquistFrequency;
  UINT4 resampleFactor;
  REAL8 *dataPtr = NULL;
  UINT4 j;

  resampleFactor = floor( dt / series->deltaT + 0.5 );
  newNyquistFrequency = 0.5 / dt;

  /* check that the resampling factor is valid */
  if ( resampleFactor < 1 ||
      fabs( dt - resampleFactor * series->deltaT ) > 1e-3 * series->deltaT )
    XLAL_ERROR( XLAL_EINVAL );

  /* just return if no resampling is required */
  if ( resampleFactor == 1 )
  {
    XLALPrintInfo( "XLAL Info - %s: No resampling required", __func__ );
    return 0;
  }

  /* check that we are resampling by a power of two */
  if ( resampleFactor & (resampleFactor - 1) )
    XLAL_ERROR( XLAL_EINVAL );

  if ( XLALLowPassREAL8TimeSeries( series, newNyquistFrequency,
        newNyquistAmplitude, filterOrder ) < 0 )
    XLAL_ERROR( XLAL_EFUNC );

  /* decimate the time series */
  series->deltaT = dt;
  series->data->length /= resampleFactor;
  dataPtr = series->data->data;
  for ( j = 0; j < series->data->length; ++j )
  {
    series->data->data[j] = *dataPtr;
    dataPtr += resampleFactor;
  }
  series->data->data = LALRealloc( series->data->data,
      series->data->length * sizeof( *series->data->data ) );

  return 0;
}


/**
 * \deprecated Use XLALResampleREAL4TimeSeries() instead.
 */
void
LALResampleREAL4TimeSeries(
    LALStatus          *status,
    REAL4TimeSeries    *ts,
    ResampleTSParams   *params
    )
{
  UINT4       resampleFactor;
  UINT4       j;
  REAL4      *dataPtr = NULL;
  /* INT8        startTimeNS = 0; */
  CHAR        warnMsg[256];
  UINT4       corrupted = 0;
  UINT4       filterOrder = 0;
  REAL8Vector recursCoef;
  REAL8Vector directCoef;
  REAL8       recursd0 = -1;

  /* LDAS low pass FIR filter coefficents for resampling by 2, 4 and 8      */
  /* thse were generated using the firlp action to get the FIR coeffs       */
  /* used by the resample action. FIR coeffs are provided for the default   */
  /* resample action that used a Kaiser window with beta = 5 and a filter   */
  /* order parameter n = 10. The order of the filter is 2 * n * resampRatio */
  /* UINT4 ldasFilterOrderParameter = 10; */
  UINT4 ldasByTwoOrder = 40;
  UNUSED CHAR  ldasByTwoMsg[] =
    "LDAS FIR filter generated by firlp(0.5,40) using LDAS-CIT 0.7.0";
  REAL8 ldasByTwo[] =
  {
    -7.1587345178999983e-19, -1.0514587726686521e-03, 1.8542385840153159e-18,
    2.5089668121365313e-03, -3.4940298359108189e-18, -4.8948343392593557e-03,
    5.6062660030794468e-18, 8.5565590024681924e-03, -8.0907458797357440e-18,
    -1.3989973238439804e-02, 1.0780775955476128e-17, 2.2023120574050724e-02,
    -1.3459232902020823e-17, -3.4340179881765416e-02, 1.5884076076327108e-17,
    5.5288283911880384e-02, -1.7819618555760906e-17, -1.0093019739775573e-01,
    1.9068667053692095e-17, 3.1670034568032440e-01, 5.0025873529805742e-01,
    3.1670034568032440e-01, 1.9068667053692095e-17, -1.0093019739775573e-01,
    -1.7819618555760906e-17, 5.5288283911880384e-02, 1.5884076076327108e-17,
    -3.4340179881765416e-02, -1.3459232902020823e-17, 2.2023120574050724e-02,
    1.0780775955476128e-17, -1.3989973238439804e-02, -8.0907458797357440e-18,
    8.5565590024681924e-03, 5.6062660030794468e-18, -4.8948343392593557e-03,
    -3.4940298359108189e-18, 2.5089668121365313e-03, 1.8542385840153159e-18,
    -1.0514587726686521e-03, -7.1587345178999983e-19
  };
  UINT4 ldasByFourOrder = 80;
  UNUSED CHAR  ldasByFourMsg[] =
    "LDAS FIR filter generated by firlp(0.25,80) using LDAS-CIT 0.7.0";
  REAL8 ldasByFour[] =
  {
    -3.5797933214045194e-19, -2.8264939479410322e-04, -5.2579196542625766e-04,
    -4.7541698372288916e-04, 9.2722964970291765e-19, 7.3162852724022317e-04,
    1.2546327308624197e-03, 1.0630797667142953e-03, -1.7472228702023229e-18,
    -1.4830129609431080e-03, -2.4477084927857239e-03, -2.0065313653173720e-03,
    2.8034666665817767e-18, 2.6512725192001795e-03, 4.2787887572376141e-03,
    3.4384897676469164e-03, -4.0458544723286438e-18, -4.3947913127916887e-03,
    -6.9958192527421722e-03, -5.5549713919258352e-03, 5.3910296112353999e-18,
    6.9667290898765182e-03, 1.1012871024947616e-02, 8.6988136875756853e-03,
    -6.7304174967527627e-18, -1.0855771610536039e-02, -1.7172133746431371e-02,
    -1.3606372767203865e-02, 7.9429834019598979e-18, 1.7243084699552668e-02,
    2.7647432518244298e-02, 2.2320837133020830e-02, -8.9108698382911643e-18,
    -3.0033587538646083e-02, -5.0471105705777050e-02, -4.3494435742894542e-02,
    9.5354684261874058e-18, 7.4135704901104452e-02, 1.5836902171998732e-01,
    2.2490814275257559e-01, 2.5015914127230293e-01, 2.2490814275257559e-01,
    1.5836902171998732e-01, 7.4135704901104452e-02, 9.5354684261874058e-18,
    -4.3494435742894542e-02, -5.0471105705777050e-02, -3.0033587538646083e-02,
    -8.9108698382911643e-18, 2.2320837133020830e-02, 2.7647432518244298e-02,
    1.7243084699552668e-02, 7.9429834019598979e-18, -1.3606372767203865e-02,
    -1.7172133746431371e-02, -1.0855771610536039e-02, -6.7304174967527627e-18,
    8.6988136875756853e-03, 1.1012871024947616e-02, 6.9667290898765182e-03,
    5.3910296112353999e-18, -5.5549713919258352e-03, -6.9958192527421722e-03,
    -4.3947913127916887e-03, -4.0458544723286438e-18, 3.4384897676469164e-03,
    4.2787887572376141e-03, 2.6512725192001795e-03, 2.8034666665817767e-18,
    -2.0065313653173720e-03, -2.4477084927857239e-03, -1.4830129609431080e-03,
    -1.7472228702023229e-18, 1.0630797667142953e-03, 1.2546327308624197e-03,
    7.3162852724022317e-04, 9.2722964970291765e-19, -4.7541698372288916e-04,
    -5.2579196542625766e-04, -2.8264939479410322e-04, -3.5797933214045194e-19
  };
  UINT4 ldasByEightOrder = 160;
  UNUSED CHAR  ldasByEightMsg[] =
    "LDAS FIR filter generated by firlp(0.125,160) using LDAS-CIT 0.7.0";
  REAL8 ldasByEight[] =
  {
    -1.7899485045886187e-19, -6.5785565693739621e-05, -1.4132879082976897e-04,
    -2.1264395052204678e-04, -2.6290359742617416e-04, -2.7550349276906565e-04,
    -2.3771537702543715e-04, -1.4425544130102367e-04, 4.6362825333301403e-19,
    1.7887185358663154e-04, 3.6582485933411455e-04, 5.2717915649190441e-04,
    6.2733453548486414e-04, 6.3531304436406394e-04, 5.3155527927017054e-04,
    3.1369007704874736e-04, -8.7363673902677791e-19, -3.7044471629615059e-04,
    -7.4152795801187910e-04, -1.0475948732225806e-03, -1.2238896950094572e-03,
    -1.2184233937642017e-03, -1.0032947419855074e-03, -5.8331646493291016e-04,
    1.4017739341284856e-18, 6.7046275382674411e-04, 1.3256746563059471e-03,
    1.8513153796509275e-03, 2.1394563456147114e-03, 2.1081914808740686e-03,
    1.7192946812996687e-03, 9.9055884251036106e-04, -2.0229858297200544e-18,
    -1.1198042255458091e-03, -2.1974593033835155e-03, -3.0470761300287886e-03,
    -3.4980109424040981e-03, -3.4255709083983424e-03, -2.7775661451061263e-03,
    -1.5917261396652274e-03, 2.6955928804956133e-18, 1.7824568219118389e-03,
    3.4834654396768100e-03, 4.8124727581638303e-03, 5.5065950049312329e-03,
    5.3772206622082304e-03, 4.3495328232139681e-03, 2.4876883619736412e-03,
    -3.3653062207633390e-18, -2.7788423673236681e-03, -5.4280430225538204e-03,
    -7.4993135801407406e-03, -8.5863155663860897e-03, -8.3949478976005458e-03,
    -6.8033834361075291e-03, -3.9013002754523102e-03, 3.9716067341932858e-18,
    4.3911559143580796e-03, 8.6217920704845935e-03, 1.1985707616989357e-02,
    1.3824116659430464e-02, 1.3633243414968320e-02, 1.1160741825101029e-02,
    6.4757499732487761e-03, -4.4555639696470439e-18, -7.5079556703908984e-03,
    -1.5017228726807877e-02, -2.1332840162384223e-02, -2.5236283794044537e-02,
    -2.5641711624359045e-02, -2.1747847773897343e-02, -1.3165154002435106e-02,
    4.7678723092621354e-18, 1.7117561087684679e-02, 3.7068926111156336e-02,
    5.8343464299929801e-02, 7.9186804418539564e-02, 9.7784206656366959e-02,
    1.1245732857891018e-01, 1.2185045954177413e-01, 1.2508319353304170e-01,
    1.2185045954177413e-01, 1.1245732857891018e-01, 9.7784206656366959e-02,
    7.9186804418539564e-02, 5.8343464299929801e-02, 3.7068926111156336e-02,
    1.7117561087684679e-02, 4.7678723092621354e-18, -1.3165154002435106e-02,
    -2.1747847773897343e-02, -2.5641711624359045e-02, -2.5236283794044537e-02,
    -2.1332840162384223e-02, -1.5017228726807877e-02, -7.5079556703908984e-03,
    -4.4555639696470439e-18, 6.4757499732487761e-03, 1.1160741825101029e-02,
    1.3633243414968320e-02, 1.3824116659430464e-02, 1.1985707616989357e-02,
    8.6217920704845935e-03, 4.3911559143580796e-03, 3.9716067341932858e-18,
    -3.9013002754523102e-03, -6.8033834361075291e-03, -8.3949478976005458e-03,
    -8.5863155663860897e-03, -7.4993135801407406e-03, -5.4280430225538204e-03,
    -2.7788423673236681e-03, -3.3653062207633390e-18, 2.4876883619736412e-03,
    4.3495328232139681e-03, 5.3772206622082304e-03, 5.5065950049312329e-03,
    4.8124727581638303e-03, 3.4834654396768100e-03, 1.7824568219118389e-03,
    2.6955928804956133e-18, -1.5917261396652274e-03, -2.7775661451061263e-03,
    -3.4255709083983424e-03, -3.4980109424040981e-03, -3.0470761300287886e-03,
    -2.1974593033835155e-03, -1.1198042255458091e-03, -2.0229858297200544e-18,
    9.9055884251036106e-04, 1.7192946812996687e-03, 2.1081914808740686e-03,
    2.1394563456147114e-03, 1.8513153796509275e-03, 1.3256746563059471e-03,
    6.7046275382674411e-04, 1.4017739341284856e-18, -5.8331646493291016e-04,
    -1.0032947419855074e-03, -1.2184233937642017e-03, -1.2238896950094572e-03,
    -1.0475948732225806e-03, -7.4152795801187910e-04, -3.7044471629615059e-04,
    -8.7363673902677791e-19, 3.1369007704874736e-04, 5.3155527927017054e-04,
    6.3531304436406394e-04, 6.2733453548486414e-04, 5.2717915649190441e-04,
    3.6582485933411455e-04, 1.7887185358663154e-04, 4.6362825333301403e-19,
    -1.4425544130102367e-04, -2.3771537702543715e-04, -2.7550349276906565e-04,
    -2.6290359742617416e-04, -2.1264395052204678e-04, -1.4132879082976897e-04,
    -6.5785565693739621e-05, -1.7899485045886187e-19
  };

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( ts, status,
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( ts->deltaT > 0, status,
      RESAMPLETIMESERIESH_ERATE, RESAMPLETIMESERIESH_MSGERATE );
  ASSERT( ts->data, status,
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( ts->data->length, status,
      RESAMPLETIMESERIESH_EZERO, RESAMPLETIMESERIESH_MSGEZERO );
  ASSERT( ts->data->data, status,
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( params, status,
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( params->deltaT > 0, status,
      RESAMPLETIMESERIESH_ERATE, RESAMPLETIMESERIESH_MSGERATE );

  /* determine the factor by which to resample */
  resampleFactor = floor( params->deltaT / ts->deltaT + 0.5 );

  /* check that the resampling factor is valid */
  if ( resampleFactor < 1 ||
      fabs( params->deltaT - resampleFactor * ts->deltaT )
      > 1e-3 * ts->deltaT )
  {
    ABORT( status, RESAMPLETIMESERIESH_EINVD, RESAMPLETIMESERIESH_MSGEINVD );
  }

  /* just return if no resampling is required */
  if ( resampleFactor == 1 )
  {
    LALInfo( status, "No resampling required" );
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

  /* check that we are resampling by a power of two */
  if ( resampleFactor & (resampleFactor - 1) )
  {
    ABORT( status, RESAMPLETIMESERIESH_ELOG2, RESAMPLETIMESERIESH_MSGELOG2 );
  }

  if ( params->filterType == defaultButterworth )
  {
    params->filterParams.butterworth.nMax = 20;
    params->filterParams.butterworth.f1   = 0.5 / params->deltaT;
    params->filterParams.butterworth.a1   = 0.1;
    params->filterParams.butterworth.f2   = LAL_REAL4_MAX;
    params->filterParams.butterworth.a1   = 0.0;

    LALDButterworthREAL4TimeSeries( status->statusPtr, ts,
        &(params->filterParams.butterworth) );
    CHECKSTATUSPTR( status );
  }
  else if ( params->filterType == LDASfirLP )
  {
    recursCoef.length = 1;
    recursCoef.data = &recursd0;
    params->filterParams.iirfilter.name = "FIR filter";
    params->filterParams.iirfilter.deltaT = ts->deltaT;
    params->filterParams.iirfilter.recursCoef = &recursCoef;
    params->filterParams.iirfilter.directCoef = &directCoef;
    params->filterParams.iirfilter.history = NULL;

    if ( resampleFactor == 2 )
    {
      directCoef.data = ldasByTwo;
      filterOrder = ldasByTwoOrder;
      LALInfo( status, ldasByTwoMsg );
    }
    else if ( resampleFactor == 4 )
    {
      directCoef.data = ldasByFour;
      filterOrder = ldasByFourOrder;
      LALInfo( status, ldasByFourMsg );
    }
    else if ( resampleFactor == 8 )
    {
      directCoef.data = ldasByEight;
      filterOrder = ldasByEightOrder;
      LALInfo( status, ldasByEightMsg );
    }
    else
    {
      ABORT( status, RESAMPLETIMESERIESH_ELDAS, RESAMPLETIMESERIESH_MSGELDAS );
    }

    directCoef.length = filterOrder + 1;

    LALDCreateVector( status->statusPtr,
        &(params->filterParams.iirfilter.history), filterOrder );
    CHECKSTATUSPTR( status );

    LALDIIRFilterREAL4Vector( status->statusPtr, ts->data,
        &(params->filterParams.iirfilter) );
    CHECKSTATUSPTR( status );

    /* account for the corruption of the data by the fir filter */
    corrupted = filterOrder;
    snprintf( warnMsg, 256 * sizeof(CHAR),
        "Corrupted %d points at start and end of time series", corrupted / 2 );
    LALWarning( status, warnMsg );
    for ( j = 0; j < corrupted; ++j )
    {
      ts->data->data[j] = 0.0;
    }
    memmove( ts->data->data, ts->data->data + corrupted / 2,
        (ts->data->length - corrupted / 2) * sizeof(REAL4) );
    for ( j = ts->data->length - corrupted / 2; j < ts->data->length; ++j )
    {
      ts->data->data[j] = 0.0;
    }

    LALDDestroyVector( status->statusPtr,
        &(params->filterParams.iirfilter.history) );
    CHECKSTATUSPTR( status );
  }

  /* decimate the time series */
  ts->deltaT = params->deltaT;
  ts->data->length /= resampleFactor;
  dataPtr = ts->data->data;
  for ( j = 0; j < ts->data->length; ++j )
  {
    ts->data->data[j] = *dataPtr;
    dataPtr += resampleFactor;
  }
  if ( LALRealloc( ts->data->data, ts->data->length * sizeof(REAL4) ) == NULL ) {
    XLALPrintError ("LALRealloc() failed!\n");
    ABORT( status, RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
/** @} */
