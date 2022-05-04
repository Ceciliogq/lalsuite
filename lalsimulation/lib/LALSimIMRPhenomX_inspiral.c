/*
 * Copyright (C) 2018 Geraint Pratten
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
 * \author Geraint Pratten
 */

#include <lal/XLALError.h>

#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomX_internals.h"

/******************************* IMRPhenomX Amplitude Functions: Inspiral *******************************/
/*
	Phenomenological coefficients for pseudo-PN corrections to TaylorF2 insiral:
	Amp(f,eta,chi1,chi2) = TF2(f,eta,chi1,chi2) + A0 * ( delta1*f**(7/3) + delta2*f**(8/3) + delta3*f**(9/3) )
*/


/*
    Value of amplitude collocation point at 2/4 f^A_T, defined in Eq. 5.7 of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Amp_22_v2(double eta, double S, double dchi, double delta, int InsAmpFlag){

  double eta2 = (eta*eta);
  double eta3 = (eta2*eta);
  double eta4 = (eta3*eta);

  double S2 = (S*S);
  double S3 = (S2*S);

  double noSpin, eqSpin, uneqSpin;

  switch ( InsAmpFlag )
	{
		case 103:
		{
      noSpin = (-0.015178276424448592 - 0.06098548699809163*eta + 0.4845148547154606*eta2)/(1. + 0.09799277215675059*eta);

      eqSpin = ((0.02300153747158323 + 0.10495263104245876*eta2)*S + (0.04834642258922544 - 0.14189350657140673*eta)*eta*S3
                    + (0.01761591799745109 - 0.14404522791467844*eta2)*S2)/(1. - 0.7340448493183307*S);

      uneqSpin = dchi*delta*eta4*(0.0018724905795891192 + 34.90874132485147*eta);

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_v2: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
	}

  return (noSpin + eqSpin + uneqSpin);

}

/*
    Value of amplitude collocation point at 3/4 f^A_T, defined in Eq. 5.7 of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Amp_22_v3(double eta, double S, double dchi, double delta, int InsAmpFlag){

  double eta2 = (eta*eta);
  double eta3 = (eta2*eta);
  double eta4 = (eta3*eta);

  double noSpin, eqSpin, uneqSpin;

  switch ( InsAmpFlag )
	{
		case 103:
		{
      noSpin = (-0.058572000924124644 - 1.1970535595488723*eta + 8.4630293045015*eta2)/(1. + 15.430818840453686*eta);

      eqSpin = ((-0.08746408292050666 + eta*(-0.20646621646484237 - 0.21291764491897636*S)
                  + eta2*(0.788717372588848 + 0.8282888482429105*S) - 0.018924013869130434*S)*S)/(-1.332123330797879 + 1.*S);

      uneqSpin = dchi*delta*eta4*(0.004389995099201855 + 105.84553997647659*eta);

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_v3: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
	}

  return (noSpin + eqSpin + uneqSpin);

}

/*
    Value of amplitude collocation point at 4/4 f^A_T, defined in Eq. 5.7 of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Amp_22_v4(double eta, double S, double dchi, double delta, int InsAmpFlag){

  double eta2 = (eta*eta);
  double eta3 = (eta2*eta);
  double eta4 = (eta3*eta);

  double S2 = (S*S);

  double noSpin, eqSpin, uneqSpin;

  switch ( InsAmpFlag )
	{
		case 103:
		{
      noSpin = (-0.16212854591357853 + 1.617404703616985*eta - 3.186012733446088*eta2 + 5.629598195000046*eta3)/(1. + 0.04507019231274476*eta);

      eqSpin = (S*(1.0055835408962206 + eta2*(18.353433894421833 - 18.80590889704093*S) - 0.31443470118113853*S
                  + eta*(-4.127597118865669 + 5.215501942120774*S) + eta3*(-41.0378120175805
                    + 19.099315016873643*S)))/(5.852706459485663 - 5.717874483424523*S + 1.*S2);

      uneqSpin = dchi*delta*eta4*(0.05575955418803233 + 208.92352600701068*eta);

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_v4: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
	}

  return (noSpin + eqSpin + uneqSpin);

}

/*
    f1 = freq. at 2/4 f^A_T
    f2 = freq. at 3/4 f^A_T
    f3 = freq. at 4/4 f^A_T

    v1 = Value of collocation point at f1
    v2 = Value of collocation point at f2
    v3 = Value of collocation point at f3

    These expressions do *not* depend on calibration if model assumes:
    TF2 + rho1 f^(7/3) + rho2 f^(8/3) + rho3 f^(9/3)
*/

/* Get Pseudo PN coefficient at f^(7/3) */
static double IMRPhenomX_Inspiral_Amp_22_rho1(double v1, double v2, double v3, double F1, double F2, double F3, int InsAmpFlag)
{
  double retVal;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      double f1p1o3, f2p1o3, f3p1o3;
      double f1p7o3, f1p8o3, f2p7o3, f2p8o3, f3p7o3, f3p8o3;

      f1p1o3 = cbrt(F1); // f1^(1/3)
      f2p1o3 = cbrt(F2); // f2^(1/3)
      f3p1o3 = cbrt(F3); // f3^(1/3)

      f1p7o3 = pow_7_of(f1p1o3); // f1^(7/3)
      f1p8o3 = f1p7o3 * f1p1o3;  // f1^(8/3)

      f2p7o3 = pow_7_of(f2p1o3); // f2^(7/3)
      f2p8o3 = f2p7o3 * f2p1o3;  // f2^(8/3)

      f3p7o3 = pow_7_of(f3p1o3); // f3^(7/3)
      f3p8o3 = f3p7o3 * f3p1o3;  // f3^(8/3)

      retVal = (-(f2p8o3*(F3*F3*F3)*v1) + F2*F2*F2*f3p8o3*v1 + f1p8o3*(F3*F3*F3)*v2 - F1*F1*F1*f3p8o3*v2 - f1p8o3*(F2*F2*F2)*v3 + F1*F1*F1*f2p8o3*v3)
                      / (f1p7o3*(f1p1o3 - f2p1o3)*f2p7o3*(f1p1o3 - f3p1o3)*(f2p1o3 - f3p1o3)*f3p7o3);

      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_rho1: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
  }

  return retVal;
}

/* Get Pseudo PN coefficient at f^(8/3) */
static double IMRPhenomX_Inspiral_Amp_22_rho2(double v1, double v2, double v3, double F1, double F2, double F3, int InsAmpFlag)
{
  double retVal;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      double f1p1o3, f2p1o3, f3p1o3;
      double f1p7o3, f2p7o3, f3p7o3;

      f1p1o3 = cbrt(F1); // f1^(1/3)
      f2p1o3 = cbrt(F2); // f2^(1/3)
      f3p1o3 = cbrt(F3); // f3^(1/3)

      f1p7o3 = pow_7_of(f1p1o3); // f1^(7/3)

      f2p7o3 = pow_7_of(f2p1o3); // f2^(7/3)

      f3p7o3 = pow_7_of(f3p1o3); // f3^(7/3)

      retVal = ( f2p7o3*(F3*F3*F3)*v1 - F2*F2*F2*f3p7o3*v1 - f1p7o3*(F3*F3*F3)*v2 + F1*F1*F1*f3p7o3*v2 + f1p7o3*(F2*F2*F2)*v3 - F1*F1*F1*f2p7o3*v3 )
                        / (f1p7o3*(f1p1o3 - f2p1o3)*f2p7o3*(f1p1o3 - f3p1o3)*(f2p1o3 - f3p1o3)*f3p7o3);

      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_rho2: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
  }

  return retVal;
}

/* Get Pseudo PN coefficient at f^(9/3) */
static double IMRPhenomX_Inspiral_Amp_22_rho3(double v1, double v2, double v3, double F1, double F2, double F3, int InsAmpFlag)
{
  double retVal;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      double f1p1o3, f2p1o3, f3p1o3;
      double f1p7o3, f1p8o3, f2p7o3, f2p8o3, f3p7o3, f3p8o3;

      f1p1o3 = cbrt(F1); // f1^(1/3)
      f2p1o3 = cbrt(F2); // f2^(1/3)
      f3p1o3 = cbrt(F3); // f3^(1/3)

      f1p7o3 = pow_7_of(f1p1o3); // f1^(7/3)
      f1p8o3 = f1p7o3 * f1p1o3;  // f1^(8/3)

      f2p7o3 = pow_7_of(f2p1o3); // f2^(7/3)
      f2p8o3 = f2p7o3 * f2p1o3;  // f2^(8/3)

      f3p7o3 = pow_7_of(f3p1o3); // f3^(7/3)
      f3p8o3 = f3p7o3 * f3p1o3;  // f3^(8/3)

      retVal = ( f2p8o3*f3p7o3*v1 - f2p7o3*f3p8o3*v1 - f1p8o3*f3p7o3*v2 + f1p7o3*f3p8o3*v2 + f1p8o3*f2p7o3*v3 - f1p7o3*f2p8o3*v3 )
                        / ( f1p7o3*(f1p1o3 - f2p1o3)*f2p7o3*(f1p1o3 - f3p1o3)*(f2p1o3 - f3p1o3)*f3p7o3 );

      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_rho3: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
  }

  return retVal;
}


/*
 *  TaylorF2 PN Amplitude + pseudo-PN coefficients
 */
static double IMRPhenomX_Inspiral_Amp_22_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp)
{
  double pnAmp;

  int InsAmpFlag = pWF->IMRPhenomXInspiralAmpVersion;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      // Re-factor expression
      pnAmp = (
          pAmp->pnInitial // 1.0
        + powers_of_Mf->one_third      * pAmp->pnOneThird
        + powers_of_Mf->two_thirds     * pAmp->pnTwoThirds
        + Mf                           * pAmp->pnThreeThirds
        + Mf*(
          + powers_of_Mf->one_third    * pAmp->pnFourThirds
          + powers_of_Mf->two_thirds   * pAmp->pnFiveThirds
          + Mf                         * pAmp->pnSixThirds
          + Mf*(
            + powers_of_Mf->one_third  * pAmp->rho1
            + powers_of_Mf->two_thirds * pAmp->rho2
            + Mf                       * pAmp->rho3
            )
          )
      );
      break;
    }
    default :
    {
        pnAmp = 0.0;
    }
  }

  return pnAmp;
}

/*
 *  Derivative of TaylorF2 PN Amplitude + pseudo-PN coefficients
 */
static double IMRPhenomX_Inspiral_Amp_22_DAnsatz(double Mf, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp) {

  double DAmpIns;

  int InsAmpFlag = pWF->IMRPhenomXInspiralAmpVersion;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      double eta   = pWF->eta;
      double chi1L = pWF->chi1L;
      double chi2L = pWF->chi2L;
      double rho1  = pAmp->rho1;
      double rho2  = pAmp->rho2;
      double rho3  = pAmp->rho3;

      double chi1L2 = chi1L  * chi1L;
      double chi1L3 = chi1L2 * chi1L;
      double chi2L2 = chi2L  * chi2L;

      double eta2   = eta*eta;
      double eta3   = eta*eta2;
      double Mf2    = Mf*Mf;
      double LALPi  = LAL_PI;

      double delta  = pWF->delta;

      DAmpIns =
      ( ((chi2L*(81 - 81*delta - 44*eta) + chi1L*(81*(1 + delta) - 44*eta))*LALPi)/48.
      + ((-969 + 1804*eta)*pow(LALPi,2./3.))/(1008.*pow(Mf,1/3.))
      + ((-27312085 - 10287648*chi2L2 + 10287648*chi2L2*delta - 10287648*chi1L2*(1 + delta)
         + 24*(-1975055 + 857304*chi1L2 - 994896*chi1L*chi2L + 857304*chi2L2)*eta + 35371056*eta2)*pow(LALPi,4./3.)*pow(Mf,1./3.))/6.096384e6
      + (5*pow(LALPi,5./3.)*(-6048*chi1L3*(-1 - delta + (3 + delta)*eta) + chi1L*(287213*(1 + delta) - 4*(93414 + 2083*delta)*eta - 35632*eta2)
               + chi2L*(-((287213 + 6048*chi2L2)*(-1 + delta)) + 4*(-93414 + 1512*chi2L2*(-3 + delta) + 2083*delta)*eta - 35632*eta2)
               + 42840*(-1 + 4*eta)*LALPi)*pow(Mf,2./3.))/96768.
       - (pow(LALPi,2.0)*(-336*(-3248849057 + 1809550512*chi1L2 - 2954929824*chi1L*chi2L + 1809550512*chi2L2)*eta2 - 324322727232*eta3
             + 7*(177520268561 + 29362199328*chi2L2 - 29362199328*chi2L2*delta + 29362199328*chi1L2*(1 + delta)
             + 12160253952*(chi1L + chi2L + chi1L*delta - chi2L*delta)*LALPi)
             + 12*eta*(-545384828789 + 49568837472*chi1L*chi2L - 12312458928*chi2L2 - 21943440288*chi2L2*delta
             + 77616*chi1L2*(-158633 + 282718*delta) - 8345272320*(chi1L + chi2L)*LALPi
             + 21384760320*pow(LALPi,2.0)))*Mf)/3.0042980352e10
       + (7.0/3.0)*pow(Mf,4.0/3.0)*rho1 + (8.0/3.0)*pow(Mf,5.0/3.0)*rho2 + 3*Mf2*rho3 );

      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_DAnsatz: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
  }

  return DAmpIns;
}

/******************************* Phase Functions: Inspiral *******************************/
/*
	Phenomenological coefficients for pseudo-PN corrections to TaylorF2 insiral:
	Psi(f,eta,chi1,chi2) = TF2(f,eta,chi1,chi2) + (1/eta) * ( alpha0 + alpha1*f**(3/3) + (3/4)*alpha2*f**(4/3) + (3/5)*alpha3*f**(5/3) + (1/2)*alpha4*f**(6/3) )
*/

/*
	This code is designed and intended to be modular. If a new PN TaylorF2 approximant is produced or a new fit against NR
	hybrids generated, then you can add the fit to the collocation points by simply adding in the extra case.

	This is intended to avoid a massive over-duplication of near-identical functions.
*/

/*
    Value of phase collocation point at v_3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_v3(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2  = eta*eta;
	double eta3  = eta2*eta;
	double eta4  = eta3*eta;
	double eta5  = eta4*eta;

	double S2    = S*S;
	double S3    = S2*S;

  double dchi2 = dchi*dchi;

	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* Canonical, 3 pseudo PN terms */
		{
      noSpin = (15415.000000000007 + 873401.6255736464*eta + 376665.64637025696*eta2 - 3.9719980569125614e6*eta3 + 8.913612508054944e6*eta4)/(1. + 46.83697749859996*eta);

      eqSpin = (S*(397951.95299014193 - 207180.42746987*S + eta3*(4.662143741417853e6 - 584728.050612325*S - 1.6894189124921719e6*S2) + eta*(-1.0053073129700898e6 + 1.235279439281927e6*S - 174952.69161683554*S2) - 130668.37221912303*S2 + eta2*(-1.9826323844247842e6 + 208349.45742548333*S + 895372.155565861*S2)))/(-9.675704197652225 + 3.5804521763363075*S + 2.5298346636273306*S2 + 1.*S3);

      uneqSpin = -1296.9289110696955*dchi2*eta + dchi*delta*eta*(-24708.109411857182 + 24703.28267342699*eta + 47752.17032707405*S);

  		break;
		}
		case 105:				/* Canonical, 4 pseudo PN terms */
		{
      noSpin = (11717.332402222377 + 4.972361134612872e6*eta + 2.137585030930089e7*eta2 - 8.882868155876668e7*eta3 + 2.4104945956043008e8*eta4 - 2.3143426271719798e8*eta5)/(1. + 363.5524719849582*eta);

      eqSpin = (S*(52.42001436116159 - 50.547943589389966*S + eta3*S*(-15355.56020802297 + 20159.588079899433*S) + eta2*(-286.9576245212502 + 2795.982637986682*S - 2633.1870842242447*S2) - 1.0824224105690476*S2 + eta*(-123.78531181532225 + 136.1961976556154*S - 7.534890781638552*S3) + 5.973206330664007*S3 + eta4*(1777.2176433016125 + 24069.288079063674*S - 44794.9522164669*S2 + 1584.1515277998406*S3)))/(-0.0015307616935628491 + (0.0010676159178395538 - 0.25*eta3 + 1.*eta4)*S);

      uneqSpin = -1357.9794908614106*dchi2*eta + dchi*delta*eta*(-23093.829989687543 + 21908.057881789653*eta + 49493.91485992256*S);

  		break;
		}
		case 114:				/* Extended, 3 pseudo PN terms */
		{
      noSpin = (68014.00000000003 + 1.1513072539654972e6*eta - 2.725589921577228e6*eta2 + 312571.92531733884*eta3)/(1. + 17.48539665509149*eta);

      eqSpin = (S*(-34467.00643820664 + 99693.81839115614*eta + 144345.24343461913*eta4 + (23618.044919850676 - 89494.69555164348*eta + 725554.5749749158*eta4 - 103449.15865381068*eta2)*S + (10350.863429774612 - 73238.45609787296*eta + 3.559251543095961e6*eta4 + 888228.5439003729*eta2 - 3.4602940487291473e6*eta3)*S2))/(1. - 0.056846656084188936*S - 0.32681474740130184*S2 - 0.30562055811022015*S3);

      uneqSpin = -1182.4036752941936*dchi2*eta + dchi*delta*eta*(-0.39185419821851025 - 99764.21095663306*eta + 41826.177356107364*S);

  		break;
		}
		case 115:				/* Extended, 4 pseudo PN terms */
		{
      noSpin = (60484.00000000003 + 4.370611564781374e6*eta - 5.785128542827255e6*eta2 - 8.82141241633613e6*eta3 + 1.3406996319926713e7*eta4)/(1. + 70.89393713617065*eta);

      eqSpin = (S*(21.91241092620993 - 32.57779678272056*S + eta2*(-102.4982890239095 + 2570.7566494633033*S - 2915.1250015652076*S2) + 8.130585173873232*S2 + eta*(-28.158436727309365 + 47.42665852824481*S2) + eta3*(-1635.6270690726785 - 13745.336370568011*S + 19539.310912464192*S2) + 1.2792283911312285*S3 + eta4*(5558.382039622131 + 21398.7730201213*S - 37773.40511355719*S2 + 768.6183344184254*S3)))/(-0.0007758753818017038 + (0.0005304023864415552 - 0.25000000000000006*eta3 + 1.*eta4)*S);

      uneqSpin = -1223.769262298142*dchi2*eta + dchi*delta*eta*(-16.705471562129436 - 93771.93750060834*eta + 43675.70151058481*S);

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_v3: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/*
    Value of phase collocation point for d13 = v1 - v3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_d13(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double eta4 = eta3*eta;
	double eta5 = eta4*eta;

	double S2   = S*S;
	double S3   = S2*S;
	double S4   = S3*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* Canonical, 3 pseudo PN terms */
		{
      noSpin = (-17294.000000000007 - 19943.076428555978*eta + 483033.0998073767*eta2)/(1. + 4.460294035404433*eta);

      eqSpin = (S*(68384.62786426462 + 67663.42759836042*S - 2179.3505885609297*S2 + eta*(-58475.33302037833 + 62190.404951852535*S + 18298.307770807573*S2 - 303141.1945565486*S3) + 19703.894135534803*S3 + eta2*(-148368.4954044637 - 758386.5685734496*S - 137991.37032619823*S2 + 1.0765877367729193e6*S3) + 32614.091002011017*S4))/(2.0412979553629143 + 1.*S);

      uneqSpin = 12017.062595934838*dchi*delta*eta;
  		break;
		}
		case 105:				/* Canonical, 4 pseudo PN terms */
		{
      noSpin = (-14234.000000000007 + 16956.107542097994*eta + 176345.7518697656*eta2)/(1. + 1.294432443903631*eta);

      eqSpin = (S*(814.4249470391651 + 539.3944162216742*S + 1985.8695471257474*S2 + eta*(-918.541450687484 + 2531.457116826593*S - 14325.55532310136*S2 - 19213.48002675173*S3) + 1517.4561843692861*S3 + eta2*(-517.7142591907573 - 14328.567448748548*S + 21305.033147575057*S2 + 50143.99945676916*S3)))/(0.03332712934306297 + 0.0025905919215826172*S + (0.07388087063636305 - 0.6397891808905215*eta + 1.*eta2)*S2);

      uneqSpin = dchi*delta*eta*(0.09704682517844336 + 69335.84692284222*eta);

  		break;
		}
		case 114:				/* Extended, 3 pseudo PN terms */
		{
      noSpin = (-36664.000000000015 + 277640.10051158903*eta - 581120.4916255298*eta2 + 1.415628418251648e6*eta3 - 7.640937162029471e6*eta4 + 1.1572710625359124e7*eta5)/(1. - 4.011038704323779*eta);

      eqSpin = (S*(-38790.01253014577 - 50295.77273512981*S + 15182.324439704937*S2 + eta2*(57814.07222969789 + 344650.11918139807*S + 17020.46497164955*S2 - 574047.1384792664*S3) + 24626.598127509922*S3 + eta*(23058.264859112394 - 16563.935447608965*S - 36698.430436426395*S2 + 105713.91549712936*S3)))/(-1.5445637219268247 - 0.24997068896075847*S + 1.*S2);

      uneqSpin = 74115.77361380383*dchi*delta*eta2;

  		break;
		}
		case 115:				/* Extended, 4 pseudo PN terms */
		{
      noSpin = (-29240.00000000001 - 12488.41035199958*eta + 1.3911845288427814e6*eta2 - 3.492477584609041e6*eta3)/(1. + 2.6711462529779824*eta - 26.80876660227278*eta2);

      eqSpin = (S*(-29536.155624432842 - 40024.5988680615*S + 11596.401177843705*S2 + eta2*(122185.06346551726 + 351091.59147835104*S - 37366.6143666202*S2 - 505834.54206320125*S3) + 20386.815769841945*S3 + eta*(-9638.828456576934 - 30683.453790630676*S - 15694.962417099561*S2 + 91690.51338194775*S3)))/(-1.5343852108869265 - 0.2215651087041266*S + 1.*S2);

      uneqSpin = 68951.75344813892*dchi*delta*eta2;

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_d13: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/*
    Value of phase collocation point for d23 = v2 - v3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_d23(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2  = eta*eta;
	double eta3  = eta2*eta;
	
	double dchi2 = dchi * dchi;

	double S2    = S*S;
	double S3    = S2*S;
	double S4    = S3*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* Canonical, 4 pseudo PN coefficients */
		{
      noSpin = (-7579.300000000004 - 120297.86185566607*eta + 1.1694356931282217e6*eta2 - 557253.0066989232*eta3)/(1. + 18.53018618227582*eta);

      eqSpin = (S*(-27089.36915061857 - 66228.9369155027*S + eta2*(150022.21343386435 - 50166.382087278434*S - 399712.22891153296*S2) - 44331.41741405198*S2 + eta*(50644.13475990821 + 157036.45676788126*S + 126736.43159783827*S2) + eta3*(-593633.5370110178 - 325423.99477314285*S + 847483.2999508682*S2)))/(-1.5232497464826662 - 3.062957826830017*S - 1.130185486082531*S2 + 1.*S3);

      uneqSpin = 3843.083992827935*dchi*delta*eta;

  		break;
		}
		case 105:				/* Canonical, 5 pseudo PN coefficients */
		{
      noSpin = (-7520.900000000003 - 49463.18828584058*eta + 437634.8057596484*eta2)/(1. + 9.10538019868398*eta);

      eqSpin = (S*(25380.485895523005 + 30617.553968012628*S + 5296.659585425608*S2 + eta*(-49447.74841021405 - 94312.78229903466*S - 5731.614612941746*S3) + 2609.50444822972*S3 + 5206.717656940992*S4 + eta2*(54627.758819129864 + 157460.98527210607*S - 69726.85196686552*S2 + 4674.992397927943*S3 + 20704.368650323784*S4)))/(1.5668927528319367 + 1.*S);

      uneqSpin = -95.38600275845481*dchi2*eta + dchi*delta*eta*(3271.8128884730654 + 12399.826237672185*eta + 9343.380589951552*S);

  		break;
		}
		case 114:				/* Extended, 4 pseudo PN coefficients */
		{
      noSpin = (-17762.000000000007 - 1.6929191194109183e6*eta + 8.420903644926643e6*eta2)/(1. + 98.061533474615*eta);

      eqSpin = (S*(-46901.6486082098 - 83648.57463631754*S + eta2*(1.2502334322912344e6 + 1.4500798116821344e6*S - 1.4822181506831646e6*S2) - 41234.966418619966*S2 + eta*(-24017.33452114588 - 15241.079745314566*S + 136554.48806839858*S2) + eta3*(-3.584298922116994e6 - 3.9566921791790277e6*S + 4.357599992831832e6*S2)))/(-3.190220646817508 - 3.4308485421201387*S - 0.6347932583034377*S2 + 1.*S3);

      uneqSpin = 24906.33337911219*dchi*delta*eta2;

  		break;
		}
		case 115:				/* Extended, 5 pseudo PN coefficients */
		{
      noSpin = (-18482.000000000007 - 1.2846128476247871e6*eta + 4.894853535651343e6*eta2 + 3.1555931338015324e6*eta3)/(1. + 82.79386070797756*eta);

      eqSpin = (S*(-19977.10130179636 - 24729.114403562427*S + 10325.295899053815*S2 + eta*(30934.123894659646 + 58636.235226102894*S - 32465.70372990005*S2 - 38032.16219587224*S3) + 15485.725898689267*S3 + eta2*(-38410.1127729419 - 87088.84559983511*S + 61286.73536122325*S2 + 42503.913487705235*S3)))/(-1.5148031011828222 - 0.24267195338525768*S + 1.*S2);

      uneqSpin = 5661.027157084334*dchi*delta*eta;

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_d23: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/*
    Value of phase collocation point for d43 = v4 - v3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_d43(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double eta4 = eta3*eta;

	double S2   = S*S;
	double S3   = S2*S;
	
	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* Canonical, 3 pseudo PN coefficients */
		{
      noSpin = (2439.000000000001 - 31133.52170083207*eta + 28867.73328134167*eta2)/(1. + 0.41143032589262585*eta);

      eqSpin = (S*(16116.057657391262 + eta3*(-375818.0132734753 - 386247.80765802023*S) + eta*(-82355.86732027541 - 25843.06175439942*S) + 9861.635308837876*S + eta2*(229284.04542668918 + 117410.37432997991*S)))/(-3.7385208695213668 + 0.25294420589064653*S + 1.*S2);

      uneqSpin = 194.5554531509207*dchi*delta*eta;

  		break;
		}
		case 105:				/* Canonical, 4 pseudo PN coefficients */
		{
      noSpin = (4085.300000000002 + 62935.7755506329*eta - 1.3712743918777364e6*eta2 + 5.024685134555112e6*eta3 - 3.242443755025284e6*eta4)/(1. + 20.889132970603523*eta - 99.39498823723363*eta2);

      eqSpin = (S*(-299.6987332025542 - 106.2596940493108*S + eta3*(2383.4907865977148 - 13637.11364447208*S - 14808.138346145908*S2) + eta*(1205.2993091547498 - 508.05838536573464*S - 1453.1997617403304*S2) + 132.22338129554674*S2 + eta2*(-2438.4917103042208 + 5032.341879949591*S + 7026.9206794027405*S2)))/(0.03089183275944264 + 1.*eta3*S - 0.010670764224621489*S2);

      uneqSpin = -1392.6847421907178*dchi*delta*eta;

  		break;
		}
		case 114:				/* Extended, 3 pseudo PN coefficients */
		{
      noSpin = (5749.000000000003 - 37877.95816426952*eta)/(1. + 1.1883386102990128*eta);

      eqSpin = ((-4285.982163759047 + 24558.689969419473*eta - 49270.2296311733*eta2)*S + eta*(-24205.71407420114 + 70777.38402634041*eta)*S2 + (2250.661418551257 + 187.95136178643946*eta - 11976.624134935797*eta2)*S3)/(1. - 0.7220334077284601*S);

      uneqSpin = dchi*delta*eta*(339.69292150803585 - 3459.894150148715*S);

  		break;
		}
		case 115:				/* Extended, 4 pseudo PN coefficients */
		{
      noSpin = (9760.400000000005 + 9839.852773121198*eta - 398521.0434645335*eta2 + 267575.4709475981*eta3)/(1. + 6.1355249449135005*eta);

      eqSpin = (S*(-1271.406488219572 + eta2*(-9641.611385554736 - 9620.333878140807*S) - 1144.516395171019*S + eta*(5155.337817255137 + 4450.755534615418*S)))/(0.1491519640750958 + (-0.0008208549820159909 - 0.15468508831447628*eta + 0.7266887643762937*eta2)*S + (0.02282542856845755 - 0.445924460572114*eta + 1.*eta2)*S2);

      uneqSpin = -1366.7949288045616*dchi*delta*eta;

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_d43: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/*
    Value of phase collocation point for d53 = v5 - v3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_d53(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2 = eta*eta;
	double eta3 = eta2*eta;

	double S2   = S*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* This should not be called for 4 pseudo-PN coefficients. Return 0 and print warning just in case... */
		{
			XLAL_ERROR_REAL8(XLAL_EINVAL, "Calling IMRPhenomX_Inspiral_Phase_22_d53 but trying to pass InspPhaseFlag for 4 pseudo-PN coefficients. Check this.\n");

  		break;
		}
		case 105:				/* Canonical, 5 pseudo PN coefficients */
		{
      noSpin = (5474.400000000003 + 131008.0112992443*eta - 1.9692364337640922e6*eta2 + 1.8732325307375633e6*eta3)/(1. + 32.90929274981482*eta);

      eqSpin = (S*(18609.016486281424 - 1337.4947536109685*S + eta2*(219014.98908698096 - 307162.33823247004*S - 124111.02067626518*S2) - 7394.9595046977365*S2 + eta*(-87820.37490863055 + 53564.4178831741*S + 34070.909093771494*S2) + eta3*(-248096.84257893753 + 536024.5354098587*S + 243877.65824670633*S2)))/(-1.5282904337787517 + 1.*S);

      uneqSpin = -429.1148607925461*dchi*delta*eta;

  		break;
		}
		case 114:				/* Extended, 4 pseudo PN coefficients */
		{
			XLAL_ERROR_REAL8(XLAL_EINVAL, "Calling IMRPhenomX_Inspiral_Phase_22_d53 but trying to pass InspPhaseFlag for 4 pseudo-PN coefficients. Check this.\n");
		}
		case 115:				/* Extended, 5 pseudo PN terms */
		{
      noSpin = (12971.000000000005 - 93606.05144508784*eta + 102472.4473167639*eta2)/(1. - 0.8909484992212859*eta);

      eqSpin = (S*(16182.268123259992 + 3513.8535400032874*S + eta2*(343388.99445324624 - 240407.0282222587*S - 312202.59917289804*S2) - 10814.056847109632*S2 + eta*(-94090.9232151429 + 35305.66247590705*S + 65450.36389642103*S2) + eta3*(-484443.15601144277 + 449511.3965208116*S + 552355.592066788*S2)))/(-1.4720837917195788 + 1.*S);

      uneqSpin = -494.2754225110706*dchi*delta*eta;

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_d53: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

static double IMRPhenomX_Inspiral_Phase_22_p1(IMRPhenomXWaveformStruct *pWF)
{
    double total;
    double S = pWF->STotR;
    double chidiff = pWF->dchi/2.;
    double eta = pWF->eta;
    double eta1 = eta;
    double eta2 = eta * eta1;
    double eta3 = eta * eta2;
    double eta4 = eta * eta3;
    double eta5 = eta * eta4;
    double eta6 = eta * eta5;
    double eta7 = eta * eta6;
    double eta8 = eta * eta7;
    double S1 = S;
    double S2 = S * S1;
    double S3 = S * S2;
    double S4 = S * S3;
    double S5 = S * S4;
    double chidiff1 = chidiff;
    double chidiff2 = chidiff * chidiff1;
    double chidiff3 = chidiff * chidiff2;
    total = 4519.7565316614055 - 2611.2832356628037*eta1 + chidiff2*(-9.64898875603356e6 + 4.774136039413815e8*eta1 - 1.0126370619124613e10*eta2 + 1.204402186677466e11*eta3 - 8.79992926696473e11*eta4 + 4.0508259399700923e12*eta5 - 1.14882207823019e13*eta6 + 1.8373732195462902e13*eta7 - 1.2700812276830186e13*eta8) + chidiff1*(-9.955334580156608e6 + 4.575863833419086e8*eta1 - 9.030596017354996e9*eta2 + 1.0007566651361551e11*eta3 - 6.818969201923483e11*eta4 + 2.9282860078355347e12*eta5 - 7.746100567095989e12*eta6 + 1.1548617002143367e13*eta7 - 7.434375182981527e12*eta8) + chidiff1*(-1.504115426754405e8 + 7.135562478407768e9*eta1 - 1.4549737691426587e11*eta2 + 1.6674307048636768e12*eta3 - 1.1760830502543629e13*eta4 + 5.233877728223997e13*eta5 - 1.436720050382833e14*eta6 + 2.2263685197595975e14*eta7 - 1.492479381758081e14*eta8)*S1 + 457.6475788624141*(-9736.866761195219 + 512524.1159527907*eta1 - 1.1594981410879988e7*eta2 + 1.4713376353356564e8*eta3 - 1.1448448612454278e9*eta4 + 5.593484167829484e9*eta5 - 1.6765211583910353e10*eta6 + 2.8204706254398148e10*eta7 - 2.040908621640001e10*eta8)*S1 + 362.750677372211*(-15831.954639401029 + 1.1893534887694176e6*eta1 - 3.361363396032276e7*eta2 + 4.9650675310588276e8*eta3 - 4.308304139547748e9*eta4 + 2.282668302527966e10*eta5 - 7.279124941389766e10*eta6 + 1.2854426772897418e11*eta7 - 9.668914371090125e10*eta8)*S2 + chidiff2*(1.1084015105853367e9 - 5.361906297327003e10*eta1 + 1.1154944032194292e12*eta2 - 1.3047251923934688e13*eta3 + 9.393059743485292e13*eta4 - 4.266159799375811e14*eta5 + 1.1947920183727885e15*eta6 - 1.8880381939048675e15*eta7 + 1.289845616442401e15*eta8)*S2 + chidiff3*(-1.5740313274156768e9 + 7.628704712208502e10*eta1 - 1.589917830400508e12*eta2 + 1.8626212624240574e13*eta3 - 1.3427459154458231e14*eta4 + 6.104503178016684e14*eta5 - 1.7105787737022912e15*eta6 + 2.7032188869370615e15*eta7 - 1.8458095196563602e15*eta8)*S3 - 215.38668666551237*(-519726.4172960017 + 2.5301901885822758e7*eta1 - 5.309393513255078e8*eta2 + 6.275179186625831e9*eta3 - 4.5706965199869415e10*eta4 + 2.101786959153881e11*eta5 - 5.961254791697051e11*eta6 + 9.539411658329324e11*eta7 - 6.597518525876953e11*eta8)*S3 - 767.3626885121919*(-60305.70554946897 + 3.2145141834742404e6*eta1 - 7.313432751189949e7*eta2 + 9.28806632527406e8*eta3 - 7.211289192911002e9*eta4 + 3.509681545175549e10*eta5 - 1.0470638252660025e11*eta6 + 1.7529838931096756e11*eta7 - 1.2625171869187671e11*eta8)*S4 + 569.1806585868325*(-288360.497063665 + 1.3983058486748807e7*eta1 - 2.922538375745895e8*eta2 + 3.4402968857475986e9*eta3 - 2.4958463209286537e10*eta4 + 1.1432130353455322e11*eta5 - 3.230277782365166e11*eta6 + 5.1506978208867584e11*eta7 - 3.550241846798408e11*eta8)*S5;

    return total;
}

static double IMRPhenomX_Inspiral_Phase_22_p2(IMRPhenomXWaveformStruct *pWF)
{
    double total;
    double S = pWF->STotR;
    double chidiff = pWF->dchi/2.;
    double eta = pWF->eta;
    double eta1 = eta;
    double eta2 = eta * eta1;
    double eta3 = eta * eta2;
    double eta4 = eta * eta3;
    double eta5 = eta * eta4;
    double eta6 = eta * eta5;
    double eta7 = eta * eta6;
    double eta8 = eta * eta7;
    double S1 = S;
    double S2 = S * S1;
    double S3 = S * S2;
    double S4 = S * S3;
    double S5 = S * S4;
    double chidiff1 = chidiff;
    double chidiff2 = chidiff * chidiff1;
    total = 6802.504666171843 - 31242.5773783971*eta1 + 130171.60395748973*eta2 - 198648.2553218823*eta3 + chidiff1*(-302084.6017781467 + 1.1196113279842721e7*eta1 - 1.6233637743294013e8*eta2 + 1.0802891466572094e9*eta3 - 1.940000965215533e9*eta4 - 1.746635758618033e10*eta5 + 1.2192111490158032e11*eta6 - 3.072088948587896e11*eta7 + 2.8765887483424e11*eta8) + chidiff2*(3.720509793121843e6 - 1.7545591212178227e8*eta1 + 3.564380716219673e9*eta2 - 4.077818520452774e10*eta3 + 2.8757278227808997e11*eta4 - 1.2809928969555525e12*eta5 + 3.5221630149435293e12*eta6 - 5.46876580943075e12*eta7 + 3.6733680554505645e12*eta8) + chidiff1*(-1.5410638688588323e7 + 7.284017193348346e8*eta1 - 1.480480907816028e10*eta2 + 1.692151005787323e11*eta3 - 1.191016579151948e12*eta4 + 5.292126247454327e12*eta5 - 1.4512131115663654e13*eta6 + 2.247601632632997e13*eta7 - 1.5065749814901932e13*eta8)*S1 + 172.08667078017876*(2684.349729739686 - 137924.910618581*eta1 + 3.0174092968270625e6*eta2 - 3.680978709529907e7*eta3 + 2.743997505577885e8*eta4 - 1.2820158575424685e9*eta5 + 3.67133005678073e9*eta6 - 5.899687945986301e9*eta7 + 4.077968848846874e9*eta8)*S1 - 109.84488523248508*(-42121.62182276369 + 1.966827779663158e6*eta1 - 3.957383217356902e7*eta2 + 4.4873516764261794e8*eta3 - 3.140492143886465e9*eta4 + 1.3906992694297485e10*eta5 - 3.809092579665219e10*eta6 + 5.904744565279731e10*eta7 - 3.969086674618475e10*eta8)*S2 + chidiff2*(2.785044661690826e7 - 1.296457743550634e9*eta1 + 2.5879054969266613e10*eta2 - 2.8966554097735645e11*eta3 + 1.990924091243516e12*eta4 - 8.614620450719862e12*eta5 + 2.2941335072408844e13*eta6 - 3.441259823935184e13*eta7 + 2.2281016234053117e13*eta8)*S2 + 117.01658246295077*(12384.359275654264 - 464724.2776273711*eta1 + 6.989807299428555e6*eta2 - 5.188223801938274e7*eta3 + 1.6912599218518743e8*eta4 + 1.1176940840874045e8*eta5 - 2.6439222229512177e9*eta6 + 7.65583065616598e9*eta7 - 7.470415430465679e9*eta8)*S3 - 221.21888692433802*(23918.107710077144 - 1.109553466088309e6*eta1 + 2.2208638666154485e7*eta2 - 2.5086480785159123e8*eta3 + 1.7514759872484207e9*eta4 - 7.748334405805534e9*eta5 + 2.1229825304676384e10*eta6 - 3.2961688565789448e10*eta7 + 2.2215570305257996e10*eta8)*S4 + 109.14396124225325*(-40092.63009167946 + 1.7630640087484568e6*eta1 - 3.324141229396693e7*eta2 + 3.511553602436209e8*eta3 - 2.2737451928886595e9*eta4 + 9.240581379271824e9*eta5 - 2.301096228969979e10*eta6 + 3.2078876553247677e10*eta7 - 1.9143277182719604e10*eta8)*S5;

    return total;
}

static double IMRPhenomX_Inspiral_Phase_22_p3(IMRPhenomXWaveformStruct *pWF)
{
    double total;
    double S = pWF->STotR;
    double chidiff = pWF->dchi/2.;
    double eta = pWF->eta;
    double eta1 = eta;
    double eta2 = eta * eta1;
    double eta3 = eta * eta2;
    double eta4 = eta * eta3;
    double eta5 = eta * eta4;
    double eta6 = eta * eta5;
    double eta7 = eta * eta6;
    double eta8 = eta * eta7;
    double S1 = S;
    double S2 = S * S1;
    double S3 = S * S2;
    double S4 = S * S3;
    double S5 = S * S4;
    double chidiff1 = chidiff;
    double chidiff2 = chidiff * chidiff1;
    double chidiff3 = chidiff * chidiff2;
    total = 9899.499417739311 - 107719.5455984343*eta1 + 934013.2337711933*eta2 - 4.443779096991549e6*eta3 + 1.0957463175309092e7*eta4 - 1.095294229295704e7*eta5 + chidiff1*(-988648.3586590848 + 4.60475111288568e7*eta1 - 9.210729484002321e8*eta2 + 1.0341197701078651e10*eta3 - 7.132957709954004e10*eta4 + 3.0972093544385284e11*eta5 - 8.271946086040485e11*eta6 + 1.2429238856536748e12*eta7 - 8.046719854747562e11*eta8) + chidiff2*(185985.05896009022 - 9.639275431068508e6*eta1 + 2.192333669229349e8*eta2 - 2.8339468128717933e9*eta3 + 2.2620744877723763e10*eta4 - 1.1370090500575662e11*eta5 + 3.507197508431647e11*eta6 - 6.065623941847914e11*eta7 + 4.5040131423325635e11*eta8) + chidiff3*(2.8071215648028813e6 - 1.3401882714464022e8*eta1 + 2.7508051702240953e9*eta2 - 3.1735495062903145e10*eta3 + 2.2531933274643768e11*eta4 - 1.0091720449027404e12*eta5 + 2.7871305773713247e12*eta6 - 4.3433893711411846e12*eta7 + 2.9263758559401494e12*eta8) + chidiff1*(-7.013415292991192e6 + 3.395517460047173e8*eta1 - 7.079818297069201e9*eta2 + 8.311115317869359e10*eta3 - 6.01319954990838e11*eta4 + 2.7479823205279746e12*eta5 - 7.752032891019225e12*eta6 + 1.2350985825268125e13*eta7 - 8.514846723152369e12*eta8)*S1 - 28.324506994041492*(13459.105129388166 - 704384.7075543567*eta1 + 1.584510336886581e7*eta2 - 1.9940674548011428e8*eta3 + 1.5356788974353833e9*eta4 - 7.418737718267268e9*eta5 + 2.1982987757507126e10*eta6 - 3.65780745764341e10*eta7 + 2.619952305484308e10*eta8)*S1 - 268.53601259094523*(-115.77728690085621 - 11375.214367283645*eta1 + 601374.507548577*eta2 - 1.1383577066962186e7*eta3 + 1.1349654961811613e8*eta4 - 6.581216223324171e8*eta5 + 2.2367565707548833e9*eta6 - 4.1417783868551927e9*eta7 + 3.23125310449782e9*eta8)*S2 + chidiff2*(1.9212396436131816e7 - 9.174680503041137e8*eta1 + 1.8837956899313965e10*eta2 - 2.174528076950341e11*eta3 + 1.5451948505045312e12*eta4 - 6.928616480006358e12*eta5 + 1.9163666163611164e13*eta6 - 2.9919055029144855e13*eta7 + 2.0203223404606586e13*eta8)*S2 + 193.54188880524111*(20555.639016138994 - 1.0362721045181333e6*eta1 + 2.2423304792193398e7*eta2 - 2.720846997099407e8*eta3 + 2.0262362610062382e9*eta4 - 9.491467630201742e9*eta5 + 2.733696207030013e10*eta6 - 4.4304520805950554e10*eta7 + 3.096422502633316e10*eta8)*S3 - 89.10224127671712*(-5056.06601393173 + 313297.9246351895*eta1 - 8.000633986385037e6*eta2 + 1.1132117091038497e8*eta3 - 9.303512952125123e8*eta4 + 4.809245698562199e9*eta5 - 1.5081090117351822e10*eta6 + 2.6319697662375946e10*eta7 - 1.9627077428387222e10*eta8)*S4 - 50.879452255451945*(61956.50145921655 - 3.2289520907671605e6*eta1 + 7.190739223198949e7*eta2 - 8.943120198922313e8*eta3 + 6.801440541734575e9*eta4 - 3.243274966490868e10*eta5 + 9.482880245917538e10*eta6 - 1.5564762671448236e11*eta7 + 1.0994253248014148e11*eta8)*S5;
    
    return total;
}

static double IMRPhenomX_Inspiral_Phase_22_p4(IMRPhenomXWaveformStruct *pWF)
{
    double total;
    double S = pWF->STotR;
    double chidiff = pWF->dchi/2.;
    double eta = pWF->eta;
    double eta1 = eta;
    double eta2 = eta * eta1;
    double eta3 = eta * eta2;
    double eta4 = eta * eta3;
    double eta5 = eta * eta4;
    double eta6 = eta * eta5;
    double eta7 = eta * eta6;
    double eta8 = eta * eta7;
    double S1 = S;
    double S2 = S * S1;
    double S3 = S * S2;
    double S4 = S * S3;
    double S5 = S * S4;
    double chidiff1 = chidiff;
    double chidiff2 = chidiff * chidiff1;
    double chidiff3 = chidiff * chidiff2;
    total = 12223.303518621497 - 184338.04580837383*eta1 + 2.0795605573727477e6*eta2 - 1.3811902067228772e7*eta3 + 5.37847594692131e7*eta4 - 1.1369469743879686e8*eta5 + 1.0074709551638646e8*eta6 + chidiff1*(-1.1062767929994464e6 + 5.153451139512685e7*eta1 - 1.0321863979806923e9*eta2 + 1.1618163335156235e10*eta3 - 8.044223746720557e10*eta4 + 3.5107587299645935e11*eta5 - 9.437487225779009e11*eta6 + 1.4294304928295515e12*eta7 - 9.343897515278184e11*eta8) + chidiff2*(338489.6706022764 - 1.7403953679292176e7*eta1 + 3.888546032971283e8*eta2 - 4.909462300105638e9*eta3 + 3.818682677313182e10*eta4 - 1.8703282431129843e11*eta5 + 5.628371651307019e11*eta6 - 9.5133852541667e11*eta7 + 6.917343025741147e11*eta8) + chidiff3*(3.257645582852421e6 - 1.5353490184412998e8*eta1 + 3.1150232110411897e9*eta2 - 3.556689530407699e10*eta3 + 2.5020255702663123e11*eta4 - 1.1114567674973105e12*eta5 + 3.047259450265852e12*eta6 - 4.71790789573228e12*eta7 + 3.1602473697743843e12*eta8) + chidiff1*(-6.370941789590047e6 + 3.09062995881748e8*eta1 - 6.460544109035304e9*eta2 + 7.606853853795491e10*eta3 - 5.521986725623206e11*eta4 + 2.532506948690769e12*eta5 - 7.17078061832084e12*eta6 + 1.1468438592054572e13*eta7 - 7.936805820999974e12*eta8)*S1 - 153.05192275164762*(4588.337424728516 - 232093.4389642095*eta1 + 5.062386361215799e6*eta2 - 6.2011204602591544e7*eta3 + 4.6650652585302854e8*eta4 - 2.2083886908995223e9*eta5 + 6.429627028092716e9*eta6 - 1.0535736665795454e10*eta7 + 7.446041948928902e9*eta8)*S1 - 325.8623698635016*(121.042325889043 - 14854.404700615993*eta1 + 504772.76166784606*eta2 - 8.2497217351480415e6*eta3 + 7.64555553140091e7*eta4 - 4.246522171454323e8*eta5 + 1.403745720841024e9*eta6 - 2.550837001463362e9*eta7 + 1.9642095070639286e9*eta8)*S2 + chidiff2*(1.692043615030946e7 - 8.068812213604535e8*eta1 + 1.6554380788949316e10*eta2 - 1.9106073046419745e11*eta3 + 1.3582315075121965e12*eta4 - 6.096305224191475e12*eta5 + 1.688737668295182e13*eta6 - 2.641902990192917e13*eta7 + 1.7884844781400465e13*eta8)*S2 + 208.79096253610817*(27955.9125110607 - 1.3983738865010922e6*eta1 + 3.0030726686248127e7*eta2 - 3.6181435942691624e8*eta3 + 2.676894586837965e9*eta4 - 1.24651815170199e10*eta5 + 3.5711646909074165e10*eta6 - 5.760624204627082e10*eta7 + 4.009617265505005e10*eta8)*S3 - 82.58614064808282*(-1934.3130065800804 + 139116.50902279763*eta1 - 3.8903565191926057e6*eta2 + 5.7498116206102274e7*eta3 - 5.01368996417243e8*eta4 + 2.6738645677699804e9*eta5 - 8.586972651214526e9*eta6 + 1.527056655297533e10*eta7 - 1.1563195949999542e10*eta8)*S4 - 81.78172545654509*(69265.16672364155 - 3.5188428107304987e6*eta1 + 7.658062025122945e7*eta2 - 9.331445106421219e8*eta3 + 6.96997553516923e9*eta4 - 3.271583963513568e10*eta5 + 9.434986620707211e10*eta6 - 1.5302603147043164e11*eta7 + 1.0698544874303833e11*eta8)*S5;

    return total;
}

static double IMRPhenomX_Inspiral_Phase_22_p5(IMRPhenomXWaveformStruct *pWF)
{
    double total;
    double S = pWF->STotR;
    double chidiff = pWF->dchi/2.;
    double eta = pWF->eta;
    double eta1 = eta;
    double eta2 = eta * eta1;
    double eta3 = eta * eta2;
    double eta4 = eta * eta3;
    double eta5 = eta * eta4;
    double eta6 = eta * eta5;
    double eta7 = eta * eta6;
    double eta8 = eta * eta7;
    double S1 = S;
    double S2 = S * S1;
    double S3 = S * S2;
    double S4 = S * S3;
    double S5 = S * S4;
    double chidiff1 = chidiff;
    double chidiff2 = chidiff * chidiff1;
    double chidiff3 = chidiff * chidiff2;
    total = 13068.869405273277 - 208719.014368122*eta1 + 2.393749544104357e6*eta2 - 1.6118106502561683e7*eta3 + 6.35149187163301e7*eta4 - 1.3563422722574383e8*eta5 + 1.2122625655708757e8*eta6 + chidiff1*(-1.1346912335370549e6 + 5.2901824285784006e7*eta1 - 1.060843330338188e9*eta2 + 1.1960143771000502e10*eta3 - 8.29828482129474e10*eta4 + 3.6310012727119434e11*eta5 - 9.791268332829137e11*eta6 + 1.4885483708775376e12*eta7 - 9.77316464156896e11*eta8) + chidiff2*(110961.07593861417 - 6.057307212497529e6*eta1 + 1.4516431101961306e8*eta2 - 1.9644772970515106e9*eta3 + 1.6275357488759272e10*eta4 - 8.421545338799582e10*eta5 + 2.6556883978961383e11*eta6 - 4.669434283641991e11*eta7 + 3.509703243263938e11*eta8) + chidiff3*(3.267444776833002e6 - 1.536363335586708e8*eta1 + 3.1105567336772456e9*eta2 - 3.5450232769343575e10*eta3 + 2.4897906193522418e11*eta4 - 1.1044788269101484e12*eta5 + 3.0245277040645454e12*eta6 - 4.67804383096105e12*eta7 + 3.130977498578685e12*eta8) + chidiff1*(-6.144518388941518e6 + 2.978794405469369e8*eta1 - 6.223125084375237e9*eta2 + 7.323648049789854e10*eta3 - 5.3141814142908527e11*eta4 + 2.4363811842103896e12*eta5 - 6.896826301669711e12*eta6 + 1.1028336151896166e13*eta7 - 7.631480580680725e12*eta8)*S1 - 201.29625853950463*(3426.151369152376 - 173933.4333609392*eta1 + 3.810139759079006e6*eta2 - 4.687278222812959e7*eta3 + 3.5406957757174265e8*eta4 - 1.6825618379216716e9*eta5 + 4.916074094298173e9*eta6 - 8.081732298988556e9*eta7 + 5.728540656278275e9*eta8)*S1 - 356.6790175246098*(-185.03847226166693 + 2180.694287932184*eta1 + 102754.74212872074*eta2 - 2.975693343703083e6*eta3 + 3.429586226575249e7*eta4 - 2.1398189787438172e8*eta5 + 7.600228570436697e8*eta6 - 1.4493949576746433e9*eta7 + 1.155005124225028e9*eta8)*S2 + chidiff2*(1.6494383281695059e7 - 7.863349867866403e8*eta1 + 1.6129603154777107e10*eta2 - 1.8613300100383533e11*eta3 + 1.3230749891095264e12*eta4 - 5.938042294598124e12*eta5 + 1.6447603202969436e13*eta6 - 2.572829257982972e13*eta7 + 1.741469831553134e13*eta8)*S2 + 202.55782307766498*(30514.903074204405 - 1.5288846331143344e6*eta1 + 3.2889330918312784e7*eta2 - 3.969483667156421e8*eta3 + 2.94205122660503e9*eta4 - 1.3724275191813652e10*eta5 + 3.938799187216622e10*eta6 - 6.3645791021214905e10*eta7 + 4.437370714786865e10*eta8)*S3 - 90.29446915177176*(414.6844618996613 + 13345.895499183154*eta1 - 1.008696578842175e6*eta2 + 2.0578080298582397e7*eta3 - 2.1187452804073665e8*eta4 + 1.25007431914817e9*eta5 - 4.2935322044232655e9*eta6 + 8.005217704685123e9*eta7 - 6.275336210587044e9*eta8)*S4 - 83.91524653031406*(74302.85593826405 - 3.7746379587366753e6*eta1 + 8.216700673765357e7*eta2 - 1.0017035574293776e9*eta3 + 7.487299149136281e9*eta4 - 3.517456578435511e10*eta5 + 1.0154159570054008e11*eta6 - 1.648692678154233e11*eta7 + 1.1539777316831628e11*eta8)*S5;

    return total;
}

static double IMRPhenomX_Inspiral_Phase_22_p6(IMRPhenomXWaveformStruct *pWF)
{
    double total;
    double S = pWF->STotR;
    double chidiff = pWF->dchi/2.;
    double eta = pWF->eta;
    double eta1 = eta;
    double eta2 = eta * eta1;
    double eta3 = eta * eta2;
    double eta4 = eta * eta3;
    double eta5 = eta * eta4;
    double eta6 = eta * eta5;
    double eta7 = eta * eta6;
    double eta8 = eta * eta7;
    double S1 = S;
    double S2 = S * S1;
    double S3 = S * S2;
    double S4 = S * S3;
    double S5 = S * S4;
    double chidiff1 = chidiff;
    double chidiff2 = chidiff * chidiff1;
    double chidiff3 = chidiff * chidiff2;
    total = 13439.50771961634 - 220872.47653260114*eta1 + 2.5548277051871163e6*eta2 - 1.7300242657260023e7*eta3 + 6.845441500889009e7*eta4 - 1.4662098776995435e8*eta5 + 1.3132341927541855e8*eta6 + chidiff1*(-1.2013628595408874e6 + 5.609817047716815e7*eta1 - 1.1266635819209006e9*eta2 + 1.2721866044583317e10*eta3 - 8.840980616915173e10*eta4 + 3.8750997758681213e11*eta5 - 1.046905346309988e12*eta6 + 1.5948956776532112e12*eta7 - 1.0495820645853236e12*eta8) + chidiff2*(6396.90096730777 - 923038.2951488711*eta1 + 3.629072690093656e7*eta2 - 6.624095800960374e8*eta3 + 6.670112754228356e9*eta4 - 3.945761619533916e10*eta5 + 1.3689827209669908e11*eta6 - 2.5824330570565085e11*eta7 + 2.0470180601938974e11*eta8) + chidiff3*(3.467599999648828e6 - 1.63320987363281e8*eta1 + 3.3119284747131176e9*eta2 - 3.780317167321989e10*eta3 + 2.6589809833659833e11*eta4 - 1.1812335178012505e12*eta5 + 3.2392943216168145e12*eta6 - 5.01724988063893e12*eta7 + 3.3627094934581533e12*eta8) + chidiff1*(-6.200497967443122e6 + 3.005745873773965e8*eta1 - 6.278962303877359e9*eta2 + 7.3887758845549e10*eta3 - 5.361021642246001e11*eta4 + 2.4576788028873384e12*eta5 - 6.956683552407369e12*eta6 + 1.1123512719960342e13*eta7 - 7.697096989874431e12*eta8)*S1 - 210.73432339473712*(3636.8440955438455 - 183826.26615613265*eta1 + 4.0098155561561515e6*eta2 - 4.913945623175002e7*eta3 + 3.6991874096017885e8*eta4 - 1.7525357732980201e9*eta5 + 5.106725546734686e9*eta6 - 8.375012375390846e9*eta7 + 5.923662088428856e9*eta8)*S1 - 376.94413516743265*(-412.45143134115676 + 13385.832604225112*eta1 - 135199.3146009821*eta2 - 128778.57632353477*eta3 + 1.3304413703968465e7*eta4 - 1.1629336502929439e8*eta5 + 4.7976219413443136e8*eta6 - 9.96073120475827e8*eta7 + 8.383709557889917e8*eta8)*S2 + chidiff2*(1.6642663727385882e7 - 7.93873343126672e8*eta1 + 1.6294381624893711e10*eta2 - 1.8815470082886902e11*eta3 + 1.3383019198791726e12*eta4 - 6.010154595433331e12*eta5 + 1.6657416387928002e13*eta6 - 2.6071416735117438e13*eta7 + 1.765635513993481e13*eta8)*S2 + 199.3214453830331*(32516.167842364277 - 1.6248112715368238e6*eta1 + 3.487002801907473e7*eta2 - 4.1998179155686307e8*eta3 + 3.1071969268181114e9*eta4 - 1.4472413158551647e10*eta5 + 4.14806331621096e10*eta6 - 6.695213262700122e10*eta7 + 4.663414241736005e10*eta8)*S3 - 93.92008564699381*(975.0840496947777 - 13617.43574625765*eta1 - 443283.332469991*eta2 + 1.3850990278276674e7*eta3 - 1.623307981767618e8*eta4 + 1.0192013487414455e9*eta5 - 3.629540777539328e9*eta6 + 6.928424275723381e9*eta7 - 5.521701668638111e9*eta8)*S4 - 85.24641836595326*(75110.99078020648 - 3.8004290944280396e6*eta1 + 8.244389795336665e7*eta2 - 1.0021240755327436e9*eta3 + 7.471664079038361e9*eta4 - 3.502605963032926e10*eta5 + 1.0092806085411728e11*eta6 - 1.636160830091252e11*eta7 + 1.1436583973375378e11*eta8)*S5;

    return total;
}

/*
		See section VII.A of arXiv:2001.11412

		This function solves for the pseudo-PN coefficients. The structure of the equations is:

			c + alpha * x + beta * x^2 + gamma * x^3 + xi * x^4

		where x = f^(1/3).

		For 3 pseudo-PN parameters, we solve for: (c,alpha,beta,gamma).
		For 4 pseudo-PN parameters, we solve for: (c,alpha,beta,gamma,xi).

		Phase Derivative: TaylorF2 + pseudo-PN coefficients
*/
static double IMRPhenomX_Inspiral_Phase_22_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXPhaseCoefficients *pPhase)
{
    double phaseIN;

    // Assemble PN phase derivative series
  	phaseIN  = pPhase->dphi0; 																						    // f^{0/3}
  	phaseIN += pPhase->dphi1 	* powers_of_Mf->one_third; 								      // f^{1/3}
  	phaseIN += pPhase->dphi2 	* powers_of_Mf->two_thirds; 									    // f^{2/3}
  	phaseIN += pPhase->dphi3 	* Mf; 									        // f^{3/3}
  	phaseIN += pPhase->dphi4 	* powers_of_Mf->four_thirds; 							      // f^{4/3}
  	phaseIN += pPhase->dphi5 	* powers_of_Mf->five_thirds; 							      // f^{5/3}
  	phaseIN += pPhase->dphi6  * powers_of_Mf->two;													    // f^{6/3}
  	phaseIN += pPhase->dphi6L * powers_of_Mf->two * powers_of_Mf->log;			    // f^{6/3}, Log[f]
  	phaseIN += pPhase->dphi7  * powers_of_Mf->seven_thirds;									  // f^{7/3}
  	phaseIN += pPhase->dphi8  * powers_of_Mf->eight_thirds;									  // f^{8/3}
  	phaseIN += pPhase->dphi8L * powers_of_Mf->eight_thirds * powers_of_Mf->log;	// f^{8/3}
  	phaseIN += pPhase->dphi9  * powers_of_Mf->three;														// f^{9/3}
  	phaseIN += pPhase->dphi9L * powers_of_Mf->three * powers_of_Mf->log;				// f^{9/3}

  	// Add pseudo-PN Coefficient
  	phaseIN += ( 		pPhase->a0 * powers_of_Mf->eight_thirds
  								+ pPhase->a1 * powers_of_Mf->three
  								+ pPhase->a2 * powers_of_Mf->eight_thirds * powers_of_Mf->two_thirds
  								+ pPhase->a3 * powers_of_Mf->eight_thirds * powers_of_Mf->itself
  								+ pPhase->a4 * powers_of_Mf->eight_thirds * powers_of_Mf->four_thirds
  							);

  	phaseIN  = phaseIN * powers_of_Mf->m_eight_thirds * (5.0 / (128.0 * powers_of_lalpi.five_thirds));

    return phaseIN;
}

/**
 * Ansatz for the inspiral phase.
 * The TaylorF2 coefficients are defined elsewhere.
 */
static double IMRPhenomX_Inspiral_Phase_22_AnsatzInt(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXPhaseCoefficients *pPhase)
{

  // Assemble PN phasing series
  //const double v    = powers_of_Mf->one_third * powers_of_lalpi.one_third; // v = (\pi M f)^{1/3}
  //const double logv = log(v);

  // Sum up phasing contributions
  double phasing    = 0.0;

  /* The PN Phasing series is normalised by: 3 / (128 * eta * pi^{5/3} ) */
  /* 0  */ phasing += pPhase->phi0;                                                   // f^{-5/3}, v = 0;  Newt.
  /* 1  */ phasing += pPhase->phi1  * powers_of_Mf->one_third;                        // f^{-4/3}, v = 1;  0.5PN
  /* 2  */ phasing += pPhase->phi2  * powers_of_Mf->two_thirds;                       // f^{-3/3}, v = 2;  1.0PN
  /* 3  */ phasing += pPhase->phi3  * Mf;                                             // f^{-2/3}, v = 3;  1.5PN
  /* 4  */ phasing += pPhase->phi4  * powers_of_Mf->four_thirds;                      // f^{-1/3}, v = 4;  2.0PN
  /* 5  */ phasing += pPhase->phi5  * powers_of_Mf->five_thirds;                      // f^{0},    v = 5;  2.5PN; phi_initial = - LAL_PI_4
  /* 5L */ phasing += pPhase->phi5L * powers_of_Mf->five_thirds * powers_of_Mf->log;  // f^{0},    v = 5;  2.5PN Log terms.
  /* 6  */ phasing += pPhase->phi6  * powers_of_Mf->two;                              // f^{+1/3}; v = 6;  3.0PN
  /* 6L */ phasing += pPhase->phi6L * powers_of_Mf->two * powers_of_Mf->log;          // f^{+1/3}; v = 6;  3.0PN Log terms.
  /* 7  */ phasing += pPhase->phi7  * powers_of_Mf->seven_thirds;                     // f^{+2/3}: v = 7;  3.5PN
  /* 8  */ phasing += pPhase->phi8  * powers_of_Mf->eight_thirds;                     // f^{+3/3}; v = 8;  4.0PN
  /* 8  */ phasing += pPhase->phi8L * powers_of_Mf->eight_thirds * powers_of_Mf->log; // f^{+3/3}; v = 8;  4.0PN Log terms.
  /* 9  */ phasing += pPhase->phi9  * powers_of_Mf->three;                            // f^{+4/3}; v = 9;  4.5PN
  /* 9  */ phasing += pPhase->phi9L * powers_of_Mf->three * powers_of_Mf->log;        // f^{+4/3}; v = 9;  4.5PN

  // Now add in the pseudo-PN Coefficients
  phasing += (  pPhase->sigma1 * powers_of_Mf->eight_thirds
              + pPhase->sigma2 * powers_of_Mf->three
              + pPhase->sigma3 * powers_of_Mf->one_third  * powers_of_Mf->three
              + pPhase->sigma4 * powers_of_Mf->two_thirds * powers_of_Mf->three
              + pPhase->sigma5 * powers_of_Mf->itself     * powers_of_Mf->three
            );

  // This completes the TaylorF2 PN phasing series
  phasing = phasing * pPhase->phiNorm * powers_of_Mf->m_five_thirds;

  /* Add initial phasing: -pi/4 */
  //phasing += pPhase->phi_initial;

  return phasing;
}
