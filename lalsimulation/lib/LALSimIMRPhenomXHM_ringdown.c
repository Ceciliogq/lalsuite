/*
 * Copyright (C) 2019 Marta Colleoni, Cecilio Garcia Quiros
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
 *
 */

#include "LALSimIMRPhenomXHM_ringdown.h"

/* Include fits of the Ringdown related quantities and the ansatz for amplitude and phase */

/* These parameter-space fits are documented in the supplementary material of https://dcc.ligo.org/P2000011-v2.
   There are 2 Mathematica notebooks (one for amplitude and one for phase) that read the fits data and automatically generate the C-code below.
   For more information read https://git.ligo.org/waveforms/reviews/imrphenomx/blob/master/documentation/ParspaceFits/README and the documentation in the notebooks. */


/****************************************/
/*                                      */
/*              AMPLITUDE               */
/*                                      */
/****************************************/

/* Fits over parameter space for the ringdown amplitude coefficients (alambda, lambda, sigma) for each mode. */

// The spin parameter S = (m1^2*chi1 + m2^2*chi2)/(m1^2 + m2^2)


// alambda, lambda and sigma are the coefficients of the ringdown ansatz, see Eq. (6.2)


static double IMRPhenomXHM_RD_Amp_21_alambda(double eta, double chi1, double chi2, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,S2;
            eta2 = pow(eta,2);
            S2 = pow(S,2);
            double noSpin = sqrt(eta - 4.*eta2)*(0.00734983387668636 - 0.0012619735607202085*eta + 0.01042318959002753*eta2);
            double eqSpin = sqrt(eta - 4.*eta2)*S*(-0.004839645742570202 - 0.0013927779195756036*S + eta2*(-0.054621206928483663 + 0.025956604949552205*S + 0.020360826886107204*S2));
            double uneqSpin = -0.018115657394753674*(chi1 - 1.*chi2)*eta2*(1. - 10.539795474715346*eta2*delta);
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_a: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_21_lambda(double eta, double chi1, double chi2, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta), eta2;
            eta2 = pow(eta,2);
            double noSpin = 0.5566284518926176 + 0.12651770333481904*eta + 1.8084545267208734*eta2;
            double eqSpin = (0.29074922226651545 + eta2*(-2.101111399437034 - 3.4969956644617946*S) + eta*(0.059317243606471406 - 0.31924748117518226*S) + 0.27420263462336675*S)*S;
            double uneqSpin = 1.0122975748481835*(chi1 - 1.*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_lambda: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_33_alambda(double eta, double chi1, double chi2, int RDAmpFlag) {
  double total=0;
  switch (RDAmpFlag){
    case 122018:{
      double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
      double eta2,eta4,delta=sqrt(1-4*eta);
      eta2 = pow(eta,2);
      eta4 = pow(eta,4);
      double noSpin = sqrt(eta - 4.*eta2)*(0.013700854227665184 + 0.01202732427321774*eta + 0.0898095508889557*eta2);
      double eqSpin = sqrt(eta - 4.*eta2)*(0.0075858980586079065 + eta*(-0.013132320758494439 - 0.018186317026076343*S) + 0.0035617441651710473*S)*S;
      double uneqSpin = eta4*(chi2*(-0.09802218411554885 - 0.05745949361626237*S) + chi1*(0.09802218411554885 + 0.05745949361626237*S) + eta2*(chi1*(-4.2679864481479886 - 11.877399902871485*S) + chi2*(4.2679864481479886 + 11.877399902871485*S))*delta);
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_alambda: version is not valid. Recommended version is 122018.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_33_lambda(double eta, double chi1, double chi2, int RDAmpFlag) {
  double total=0;
  switch (RDAmpFlag){
    case 122018:{
      double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
      double eta2,S2,delta=sqrt(1-4*eta);
      eta2 = pow(eta,2);
      S2 = pow(S,2);
      double noSpin = 0.7435306475478924 - 0.06688558533374556*eta + 1.471989765837694*eta2;
      double eqSpin = S*(0.19457194111990656 + 0.07564220573555203*S + eta*(-0.4809350398289311 + 0.17261430318577403*S - 0.1988991467974821*S2));
      double uneqSpin = 1.8881959341735146*(chi1 - 1.*chi2)*eta2*delta;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_lambda: version is not valid. Recommended version is 122018.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_alambda(double eta, double chi1, double chi2, int RDAmpFlag) {
  double total=0;
  switch (RDAmpFlag){
    case 122018:{
      double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
      double eta2,eta3,eta4,eta5,eta6;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      eta4 = pow(eta,4);
      eta5 = pow(eta,5);
      eta6 = pow(eta,6);
      double noSpin = 0.00012587900257140724 + 0.03927886286971654*eta - 0.8109309606583066*eta2 + 8.820604907164254*eta3 - 51.43344812454074*eta4 + 141.81940900657446*eta5 - 140.0426973304466*eta6;
      double eqSpin = S*(-0.00006001471234796344 + eta4*(-0.7849112300598181 - 2.09188976953315*S) + eta2*(0.08311497969032984 - 0.15569578955822236*S) + eta*(-0.01083175709906557 + 0.00568899459837252*S) - 0.00009363591928190229*S + 1.0670798489407887*eta3*S);
      double uneqSpin = -0.04537308968659669*pow(chi1 - 1.*chi2,2)*eta2*(1. - 8.711096029480697*eta + 18.362371966229926*eta2) + (chi1 - 1.*chi2)*(-297.36978685672733 + 3103.2516759087644*eta - 10001.774055779177*eta2 + 9386.734883473799*eta3)*eta6;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_alambda: version is not valid. Recommended version is 122018.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_lambda(double eta, double chi1, double chi2, int RDAmpFlag) {
  double total=0;
  switch (RDAmpFlag){
    case 122018:{
      double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
      double eta2,eta3,delta=sqrt(1.-4*eta);
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      double noSpin = (sqrt(1. - 3.*eta)*(0.0341611244787871 - 0.3197209728114808*eta + 0.7689553234961991*eta2))/(0.048429644168112324 - 0.43758296068790314*eta + eta2);
      double eqSpin = sqrt(1. - 3.*eta)*S*(0.11057199932233873 + eta2*(25.536336676250748 - 71.18182757443142*S) + 9.790509295728649*eta*S + eta3*(-56.96407763839491 + 175.47259563543165*S));
      double uneqSpin = -5.002106168893265*pow(chi1 - 1.*chi2,2)*eta2*delta;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_lambda: version is not valid. Recommended version is 122018.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_44_alambda(double eta, double chi1, double chi2, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = sqrt(eta - 3.*eta2)*(0.007904587819112173 + 0.09558474985614368*eta - 2.663803397359775*eta2 + 28.298192768381554*eta3 - 136.10446022757958*eta4 + 233.23167528016833*eta5);
            double eqSpin = sqrt(eta - 3.*eta2)*S*(0.0049703757209330025 + 0.004122811292229324*S + eta*(-0.06166686913913691 + 0.014107365722576927*S)*S + eta2*(-0.2945455034809188 + 0.4139026619690879*S - 0.1389170612199015*S2) + eta3*(0.9225758392294605 - 0.9656098473922222*S + 0.19708289555425246*S2) + 0.000657528128497184*S2);
            double uneqSpin = 0.00659873279539475*(chi1 - 1.*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_a: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_44_lambda(double eta, double chi1, double chi2, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,eta3,eta4,eta5,eta6,eta7;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            eta7 = pow(eta,7);
            double noSpin = 0.7702864948772887 + 32.81532373698395*eta - 1625.1795901450212*eta2 + 31305.458876573215*eta3 - 297375.5347399236*eta4 + 1.4726521941846698e6*eta5 - 3.616582470072637e6*eta6 + 3.4585865843680725e6*eta7;
            double eqSpin = (-0.03011582009308575*S + 0.09034746027925727*eta*S + 1.8738784391649446*eta2*S - 5.621635317494836*eta3*S)/(-1.1340218677260014 + S);
            double uneqSpin = 0.959943270591552*pow(chi1 - 1.*chi2,2)*eta2 + 0.853573071529436*(chi1 - 1.*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_lambda: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_21_sigma(double eta, double chi1, double chi2, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2;
            eta2 = pow(eta,2);
            double noSpin = 1.2922261617161441 + 0.0019318405961363861*eta;
            double eqSpin = (0.04927982551108649 - 0.6703778360948937*eta + 2.6625014134659772*eta2)*S;
            double uneqSpin = 1.2001101665670462*(chi1 + pow(chi1,2) - 2.*chi1*chi2 + (-1. + chi2)*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_sigma: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_33_sigma(double eta, double chi1, double chi2, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double noSpin = 1.3;
            double eqSpin = 0.*eta*chi1*chi2;
            double uneqSpin = 0.;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_sigma: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_32_sigma(double eta, double chi1, double chi2, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double noSpin = 1.33;
            double eqSpin = 0.*eta*chi1*chi2;
            double uneqSpin = 0.;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_sigma: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_44_sigma(double eta, double chi1, double chi2, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double noSpin = 1.33;
            double eqSpin = 0.*eta*chi1*chi2;
            double uneqSpin = 0.;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_sigma: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_21_rdcp1(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 20211005:{
            double sqroot = sqrt(eta);
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            total = fabs(0.021343580475583254*chidiff2*(-97.10213947883408*eta1 + 993.9419016367331*eta2 - 2415.2952329189093*eta3) + 0.06885556929608708*chidiff3*(-8.730894204871106*eta1 + 227.0731668146707*eta2 - 753.9036430173783*eta3) - 1.4794519220618498*chidiff1*(5.819403790711993*eta1 + 4.073349627685885*eta2 - 45.048000086504345*eta3) + chidiff1*(-5.454189024650867*eta1 + 44.00711699964243*eta2 - 76.67467658585973*eta3)*S1 + chidiff2*(-52.18611304114328*eta1 + 535.7193397214598*eta2 - 1310.9470325164314*eta3)*S2 - 2*delta*S1*(-4.787287592131079*(8.291001621719554 - 24.823274968389548*eta1 + 34.936954517931994*eta2) - 0.31633164526706853*(-140.02377571397162 + 1583.7109216263243*eta1 - 4200.092689110534*eta2)*S1 + 0.05115393253661533*(-311.60769393707653 + 3100.9703631911066*eta1 - 7553.18583213418*eta2)*S2 + 0.21296236917654673*(-460.08984180006064 + 4933.726620936231*eta1 - 12793.825497761809*eta2)*S3) + delta*(1.8097567957155936 + 20.705083289490865*eta1 - 31.521701047602075*eta2)*sqroot);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_rdcp1: version is not valid. Recommended version is 20211005.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_21_rdcp2(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 20211005:{
            double sqroot = sqrt(eta);
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            total = fabs(0.06217100078922635*chidiff2*(-36.23814454304675*eta1 + 474.1435980389167*eta2 - 1308.6987641784285*eta3) + 0.07005497022476302*chidiff3*(-16.4283006765896*eta1 + 328.4248672870587*eta2 - 1054.1138140908247*eta3) - 1.1290041519971297*chidiff1*(5.758727553404927*eta1 + 9.092542215903991*eta2 - 66.40567292801774*eta3) + chidiff1*(-2.9660799805812914*eta1 + 28.91279278310904*eta2 - 51.588702285255735*eta3)*S1 + chidiff2*(-57.99831575665916*eta1 + 579.7806815769884*eta2 - 1395.532582061268*eta3)*S2 - 2*delta*S1*(-3.8916102062974005*(7.948095363324947 - 19.694949453296825*eta1 + 18.80334143250552*eta2) + 0.08979076698011658*(484.22620761236476 - 5274.746267308187*eta1 + 14021.902498736457*eta2)*S1 + 0.1846979967618848*(-56.610712798444666 + 590.3041758705608*eta1 - 1418.876857249019*eta2)*S2 + 0.09181644225213847*(-883.5345791636685 + 9268.6941843865*eta1 - 23726.423601674967*eta2)*S3) + delta*(1.2563504844958977 + 16.88962808267342*eta1 - 30.36941251095806*eta2)*sqroot);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_rdcp2: version is not valid. Recommended version is 20211005.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_21_rdcp3(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 20211005:{
            double sqroot = sqrt(eta);
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double S1 = S;
            double S2 = S * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            total = fabs(0.016863636703419223*chidiff3*(60.38220415807639*eta1 - 1689.4200900001954*eta2 + 11718.340178255758*eta3 - 23614.039833551666*eta4) + 0.02782712912173075*chidiff2*(56.6553038839133*eta1 - 1150.0941807912996*eta2 + 7490.92086917164*eta3 - 15172.388816670165*eta4) - 0.34667377106618386*chidiff1*(6.3361113336856025*eta1 - 6.547368783784894*eta2 + 50.662556793046576*eta3 - 262.712932979166*eta4) + chidiff1*(2.372131170615384*eta1 - 40.73712211004431*eta2 + 235.82761201471803*eta3 - 420.0591921653423*eta4)*S1 - 2*delta*S1*(-1.281401921665647*(-22.11070198522965 + 778.8261072843103*eta1 - 7690.759538266513*eta2 + 32163.58378849671*eta3 - 49056.724441802864*eta4) + 0.12767199915348384*(109.96070663847016 - 2633.5345130439096*eta1 + 23050.508735520925*eta2 - 87986.2682530337*eta3 + 125900.75672954152*eta4)*S1 + 0.09752841056632801*(-646.8575300205428 + 17239.313661874556*eta1 - 168379.18074462915*eta2 + 712720.3740487956*eta3 - 1.1013984334251098e6*eta4)*S2) + delta*(0.3731239199170631 + 5.655281608189427*eta1 - 11.830784164505003*eta2)*sqroot);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_rdcp3: version is not valid. Recommended version is 20211005.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_33_rdcp1(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 20211004:{
            double sqroot = sqrt(eta);
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2. + ((chi1 - chi2)*delta)/(1 + delta*delta);
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double eta5 = eta * eta4;
            double eta6 = eta * eta5;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            total = 0. + chidiff1*(0.2524655582890894 + (0.0000290277043860363 - 0.030906480381509324*chidiff1)*chidiff1) + (-0.01951546993047934*chidiff1 + 0.0011682409298398493*chidiff2)*(1.0499887340728007*S1 + 0.2431074545433419*S2) - 2*delta*S1*(0.21236819635185444*(1108.9346339713297 - 32962.866282257666*eta1 + 384746.0602986751*eta2 - 2.196869441850867e6*eta3 + 6.148960990371435e6*eta4 - 6.764860170957791e6*eta5) + 0.10835388942776275*(1826.033317116125 - 56373.153895140145*eta1 + 691726.0015073321*eta2 - 4.1836344170188503e6*eta3 + 1.2427564706888165e7*eta4 - 1.4483700980919447e7*eta5)*S1 + 0.05094257271759571*(-14991.00895618971 + 457376.64823641774*eta1 - 5.431239187854356e6*eta2 + 3.1472664546135046e7*eta3 - 8.92116525172971e7*eta4 + 9.917056995544848e7*eta5)*S2 + 0.08934985391310894*(-6830.438762038396 + 208976.4158706997*eta1 - 2.5069956553574502e6*eta2 + 1.475735884957653e7*eta3 - 4.264131946788113e7*eta4 + 4.8418392805496916e7*eta5)*S3) + delta*(1.8400858236629902 + 21.61154631206381*eta1 - 14.584158310276576*eta2)*sqroot + delta*(chidiff1*(507.4168681750715*eta1 - 17174.845204200545*eta2 + 220846.01029788316*eta3 - 1.3778924375555837e6*eta4 + 4.196721013647862e6*eta5 - 5.009252072292809e6*eta6) + chidiff2*(-1033.895211170711*eta1 + 32395.681145404924*eta2 - 396202.6393629291*eta3 + 2.369580450262266e6*eta4 - 6.944353450288207e6*eta5 + 7.993387792788644e6*eta6))*sqroot + chidiff1*delta*(1025.6441471783724*eta1 - 32020.43792360211*eta2 + 389947.68146473507*eta3 - 2.321110920478478e6*eta4 + 6.765893542715989e6*eta5 - 7.744103334665835e6*eta6)*S1*sqroot;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_rdcp1: version is not valid. Recommended version is 20211004.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_33_rdcp2(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 20211004:{
            double sqroot = sqrt(eta);
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2. + ((chi1 - chi2)*delta)/(1 + delta*delta);
            double chidiff = (chi1 - chi2)/2.;
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
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            total = 0. + chidiff1*(0.15396367473933462 + (0.000013955035631259617 - 0.031417978742236716*chidiff1)*chidiff1) + (-0.6316642216730622*chidiff1 + 0.043533435104519284*chidiff2)*(0.005339064093990394*S1 - 0.09353589687729917*S2) + delta*(0.9854964663227015 + 22.629543248314114*eta1 - 52.786371267707636*eta2 + 95.12291231077296*eta3)*sqroot + delta*(chidiff1*(3377.824989547245*eta1 - 155896.35784518704*eta2 + 3.024314970609635e6*eta3 - 3.2061195695426248e7*eta4 + 2.0079123446360356e8*eta5 - 7.433697547907646e8*eta6 + 1.5072394709272633e9*eta7 - 1.291883211046093e9*eta8) + chidiff2*(-9920.619325846208*eta1 + 435567.6418477763*eta2 - 8.055972546865637e6*eta3 + 8.142490980385853e7*eta4 - 4.861230700152241e8*eta5 + 1.7157400493115883e9*eta6 - 3.3175525662680845e9*eta7 + 2.713310974489422e9*eta8))*sqroot + chidiff1*delta*(12548.359569607217*eta1 - 555179.8222092787*eta2 + 1.0346380607039602e7*eta3 - 1.0534734632466829e8*eta4 + 6.334173592256904e8*eta5 - 2.2508941928327947e9*eta6 + 4.380956256272274e9*eta7 - 3.605767996583121e9*eta8)*S1*sqroot + delta*S1*(0.05401016903161101*eta1*(42610.37891214464 - 1.802216417634322e6*eta1 + 3.204851212871152e7*eta2 - 3.107792457136508e8*eta3 + 1.7766467548191123e9*eta4 - 5.994249336696646e9*eta5 + 1.106421088303287e10*eta6 - 8.628808324226036e9*eta7) + 0.05749015729925959*eta1*(83404.02279460595 - 3.5830150142831686e6*eta1 + 6.4786476257475086e7*eta2 - 6.398276480643461e8*eta3 + 3.7316561168440638e9*eta4 - 1.2866606551284636e10*eta5 + 2.4308645778823795e10*eta6 - 1.9431100443602737e10*eta7)*S1)*sqroot;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_rdcp2: version is not valid. Recommended version is 20211004.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_33_rdcp3(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 20211004:{
            double sqroot = sqrt(eta);
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double eta5 = eta * eta4;
            double eta6 = eta * eta5;
            double eta7 = eta * eta6;
            double S1 = S;
            double S2 = S * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            total = 0. + 0.03229832594453838*chidiff1 - 2*chidiff1*delta*(-0.24835616071544364*chidiff1*(268.8696046383242 - 7771.69658444856*eta1 + 80815.02331134354*eta2 - 336225.1002673399*eta3 + 123652.91014861298*eta4 + 2.828282937658085e6*eta5 - 5.586240132035846e6*eta6) - 0.09968224294665423*chidiff2*(1887.3534473722264 - 68445.02307956877*eta1 + 1.0327649969848702e6*eta2 - 8.298312129929935e6*eta3 + 3.7162933298529305e7*eta4 - 8.719414171019551e7*eta5 + 8.317527867576158e7*eta6) - 0.23734249008862335*(3451.2553425416577 - 133539.4739742386*eta1 + 2.106068930280985e6*eta2 - 1.7340014503134035e7*eta3 + 7.877232925014251e7*eta4 - 1.876011920011646e8*eta5 + 1.833514664929403e8*eta6)) - 0.4713669642603374*chidiff1*(-0.07566568410151524*S1 - 0.12044349808351801*S2) + delta*(0.42930396489230793 + 4.238925841192109*eta1)*sqroot + chidiff1*delta*(-1932.8080017631773*eta1 + 75392.24957488505*eta2 - 1.1979935413098074e6*eta3 + 9.937063202990264e6*eta4 - 4.545380168446339e7*eta5 + 1.0886858293793398e8*eta6 - 1.0682238998861058e8*eta7)*S1*sqroot + delta*S1*(-0.029171650156916855*eta1*(9568.55640912175 - 373434.1511165578*eta1 + 5.894745452817182e6*eta2 - 4.839462768313137e7*eta3 + 2.1857374172437274e8*eta4 - 5.160418715575672e8*eta5 + 4.984898967894294e8*eta6) - 0.0350989498033451*eta1*(6146.018549412825 - 241763.64316530153*eta1 + 3.8559824903990524e6*eta2 - 3.205512106135446e7*eta3 + 1.467822767118959e8*eta4 - 3.51541750967596e8*eta5 + 3.4455004973366094e8*eta6)*S1)*sqroot;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_rdcp3: version is not valid. Recommended version is 20211004.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rd1(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 202109302:{
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
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
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            double chidiff4 = chidiff * chidiff3;
            total = -0.1739466482698482 + 7.537678569978139*eta1 - 129.42727418189276*eta2 + 1152.1678779290485*eta3 - 5596.252695952866*eta4 + 14017.5301926682*eta5 - 14129.114641268448*eta6 + 0.000032716123144995364*chidiff3*(1.7514248158198886e6*eta1 - 7.825410685009141e7*eta2 + 1.4611023217836463e9*eta3 - 1.4829330841844187e10*eta4 + 8.855501547863487e10*eta5 - 3.115981111593222e11*eta6 + 5.988603396371404e11*eta7 - 4.8541926795013525e11*eta8) - 0.0011504323702108442*chidiff4*(89421.6032922539*eta1 - 3.9487058437301917e6*eta2 + 7.311932009098087e7*eta3 - 7.35337880171298e8*eta4 + 4.339662381699521e9*eta5 - 1.504758312610315e10*eta6 + 2.84287894654536e10*eta7 - 2.261106587253412e10*eta8) + 0.0026257865625151427*chidiff2*(19489.797200772107*eta1 - 869442.4106905428*eta2 + 1.6217329772540698e7*eta3 - 1.6395738885778844e8*eta4 + 9.711262524740118e8*eta5 - 3.3744195142689323e9*eta6 + 6.37953481369706e9*eta7 - 5.071009072314243e9*eta8) + 0.0005533741226681383*chidiff1*(-25471.83494236916*eta1 + 1.2050448218541564e6*eta2 - 2.3732453963716105e7*eta3 + 2.5275415580652657e8*eta4 - 1.5764953271678002e9*eta5 + 5.769481872542867e9*eta6 - 1.148537727171976e10*eta7 + 9.603940838488998e9*eta8) - 0.004123291099903633*(7467.260269877974*eta1 - 320776.13181768835*eta2 + 5.803568331907248e6*eta3 - 5.7328979589147605e7*eta4 + 3.343198918683303e8*eta5 - 1.1524656680446343e9*eta6 + 2.1769626223021135e9*eta7 - 1.7398870974551086e9*eta8)*S1 + chidiff1*(-25.40918083624006*eta1 + 1072.6660833824446*eta2 - 18734.56965899836*eta3 + 175584.44390866853*eta4 - 954645.9732908192*eta5 + 3.0131821343661095e6*eta6 - 5.112974877273543e6*eta7 + 3.595572602032916e6*eta8)*S1 + 0.0002218810086627597*(-29120.10723030688*eta1 + 1.3090040821205066e6*eta2 - 2.45088881975362e7*eta3 + 2.4758152746985695e8*eta4 - 1.4598029440735455e9*eta5 + 5.034032942554581e9*eta6 - 9.41638851714836e9*eta7 + 7.3801345429165e9*eta8)*S2 + 0.00013495436007594094*(570390.9574890276*eta1 - 2.446143478307294e7*eta2 + 4.415031811313711e8*eta3 - 4.350180843245395e9*eta4 + 2.5291375612468384e10*eta5 - 8.683903085006102e10*eta6 + 1.6319461266875494e11*eta7 - 1.2960614091975647e11*eta8)*S3 + 0.0008797432969120685*(-11704.866296895168*eta1 + 481280.4209533073*eta2 - 8.345583292347357e6*eta3 + 7.967348582015552e7*eta4 - 4.5426225624814487e8*eta5 + 1.5501476873964584e9*eta6 - 2.9329951983381686e9*eta7 + 2.372380577478077e9*eta8)*S4;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rd1: version is not valid. Recommended version is 202109302.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rd2(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 202109302:{
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double eta5 = eta * eta4;
            double S1 = S;
            double S2 = S * S1;
            double chidiff1 = chidiff;
            total = 7.81375597690496 - 258.65714426246353*eta1 + 3767.7760969102014*eta2 - 26452.84191675964*eta3 + 89239.33699339391*eta4 - 115612.66243797209*eta5 + 0.52998103554632*chidiff1*(467.53708731217074*eta1 - 11598.704847566307*eta2 + 103373.48646655097*eta3 - 392604.7402387375*eta4 + 539094.2581523317*eta5) + 0.3647065841610074*(94.91880994621141*eta1 - 1940.4940890647033*eta2 + 13767.008956271022*eta3 - 37252.69514989471*eta4 + 28983.867168291123*eta5)*S1 + chidiff1*(292.0590115406781*eta1 - 7076.162615718454*eta2 + 62610.165560211804*eta3 - 238356.53017706363*eta4 + 329776.44607064215*eta5)*S1 + 0.5123842751736173*(261.4202559048511*eta1 - 6426.164633817362*eta2 + 57499.62330895333*eta3 - 219351.59027355633*eta4 + 302006.83156850166*eta5)*S2;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rd2: version is not valid. Recommended version is 202109302.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rd3(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 202109302:{
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2. + ((chi1 - chi2)*delta)/(1 + delta*delta);
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double eta5 = eta * eta4;
            double eta6 = eta * eta5;
            double eta7 = eta * eta6;
            double eta8 = eta * eta7;
            double eta9 = eta * eta8;
            double S1 = S;
            double S2 = S * S1;
            double chidiff1 = chidiff;
            total = 0.48617782742174875 + 26.916693044800414*eta1 - 302.0413282620642*eta2 + 1427.2655792675407*eta3 - 2386.187638662855*eta4 + chidiff1*(-144961.8787858238*eta1 + 1.379656363337113e12*eta10 + 8.2665958651494e6*eta2 - 2.062214465573237e8*eta3 + 2.954099729303266e9*eta4 - 2.6785574987600887e10*eta5 + 1.5947855969584122e11*eta6 - 6.23752451145218e11*eta7 + 1.5461416362290193e12*eta8 - 2.2052168027849614e12*eta9) + chidiff1*(-353272.87070252787*eta1 + 3.3243223242432656e12*eta10 + 2.007989123603646e7*eta2 - 4.9952012670656335e8*eta3 + 7.139167243151436e9*eta4 - 6.461787557352506e10*eta5 + 3.8424830749208215e11*eta6 - 1.5017813048502031e12*eta7 + 3.721767700104094e12*eta8 - 5.309675698201082e12*eta9)*S1 + 0.3104462975762972*(156396.45275680657*eta1 - 1.595876620195859e12*eta10 - 9.171576994507995e6*eta2 + 2.337669026117386e8*eta3 - 3.401484343332711e9*eta4 + 3.1166514726007515e10*eta5 - 1.866628252292533e11*eta6 + 7.314858159519032e11*eta7 - 1.810337746377304e12*eta8 + 2.5699683625630913e12*eta9)*S1 + 0.3124850515168781*(228516.64405544935*eta1 - 2.469327742200185e12*eta10 - 1.3508643388277967e7*eta2 + 3.468055801307971e8*eta3 - 5.079991587824657e9*eta4 + 4.683865440595272e10*eta5 - 2.8221617977192847e11*eta6 + 1.1124169415556743e12*eta7 - 2.769005826116588e12*eta8 + 3.953569416167334e12*eta9)*S2;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rd3: version is not valid. Recommended version is 202109302.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rdaux1(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 202109302:{
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2. + ((chi1 - chi2)*delta)/(1 + delta*delta);
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double eta5 = eta * eta4;
            double eta6 = eta * eta5;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double S4 = S * S3;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            total = 0.2783873165821764 + 3.652393673687489*eta1 - 261.0162506101187*eta2 + 4656.703615652947*eta3 + chidiff1*(19.54109863660863*eta1 - 371.0282661670557*eta2 + 2174.2353069067017*eta3 - 4008.503026782427*eta4) - 36472.82351655576*eta4 + chidiff2*(-13.368638942143532*eta1 + 212.61001435565476*eta2 - 1076.601616881712*eta3 + 1777.4742472828732*eta4) + 129522.47727536858*eta5 - 169546.85338957966*eta6 + 0.0037930129337354133*(-1010.1040446456974*eta1 + 2530.4630516602615*eta2 + 49441.510244218756*eta3 - 181272.87932676548*eta4)*S1 + chidiff1*(-15.240239669443447*eta1 + 262.71457835137227*eta2 - 1427.7251992662223*eta3 + 2494.842010666985*eta4)*S1 + chidiff2*(87.0270216368149*eta1 - 1349.9508521769994*eta2 + 6717.145823185167*eta3 - 10842.04817972617*eta4)*S2 + 0.01115535749678476*(-973.4873814340503*eta1 + 17036.12800917856*eta2 - 89611.16620811584*eta3 + 147387.7558838333*eta4)*S2 - 0.0012517354331563185*(2632.745081157317*eta1 - 53348.294873378305*eta2 + 315217.23112188344*eta3 - 577673.4018830252*eta4)*S3 + chidiff3*(31.757198193332766*eta1 - 998.8868060258987*eta2 + 8014.343811421408*eta3 - 18751.4758438168*eta4)*S3 - 0.023902267955219568*(-4.493240627753124*eta1 + 258.46282436993675*eta2 - 724.5325731051831*eta3 - 1153.9665623138912*eta4)*S4;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rdaux1: version is not valid. Recommended version is 202109302.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rdaux2(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 202109302:{
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2. + ((chi1 - chi2)*delta)/(1 + delta*delta);
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double eta5 = eta * eta4;
            double eta6 = eta * eta5;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double S4 = S * S3;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            total = -5.072342758846551 + 228.9146665593112*eta1 - 4063.3257679872577*eta2 + 37578.394837325526*eta3 - 190192.81563875917*eta4 + chidiff1*(-57.636073394827456*eta1 + 1468.3990693500923*eta2 - 13633.914196217638*eta3 + 54240.86479550858*eta4 - 78038.45196262475*eta5) + 496407.58110565646*eta5 + chidiff2*(45.98894463535864*eta1 - 1061.1097079000988*eta2 + 8762.880861572887*eta3 - 30862.21281696972*eta4 + 39459.77367817098*eta5) + chidiff3*(45.10188659892899*eta1 - 1079.8868634962898*eta2 + 9345.839177508953*eta3 - 34664.876296339*eta4 + 46696.16319785507*eta5) - 520199.1781702057*eta6 - 0.010013895381547334*(-1511.9943500701688*eta1 + 42129.00385580998*eta2 - 376847.8036312572*eta3 + 1.3628049108208485e6*eta4 - 1.7125431386277918e6*eta5)*S1 + chidiff1*(-20.799402899185985*eta1 + 437.9449035140799*eta2 - 3246.2743948613584*eta3 + 9986.418004154457*eta4 - 10709.946754695258*eta5)*S1 + 0.043303488328785666*(-12.623968469483795*eta1 - 290.7836647563115*eta2 + 5955.446310745487*eta3 - 27573.32731044381*eta4 + 36300.351871964245*eta5)*S2 - 0.007224662452434095*(-205.14492266045667*eta1 + 1302.66806151529*eta2 + 5252.374839355036*eta3 - 43639.08076775458*eta4 + 56016.408724754525*eta5)*S3 - 0.03826948517545464*(-358.03267489459506*eta1 + 7325.425181804928*eta2 - 55018.050741618754*eta3 + 184179.63565019946*eta4 - 234633.13099987642*eta5)*S4;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rdaux2: version is not valid. Recommended version is 202109302.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_44_rdcp1(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 20211005:{
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double eta5 = eta * eta4;
            double eta6 = eta * eta5;
            double eta7 = eta * eta6;
            double S1 = S;
            double S2 = S * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            total = -0.52095121699195 + 27.42194590011417*eta1 - 384.2536024626785*eta2 + 2922.635377940521*eta3 - 11438.793140756907*eta4 + 17345.43514217596*eta5 + 0.009481234525653326*chidiff2*(-15524.393069698435*eta1 + 581774.3819353866*eta2 - 8.853869756545806e6*eta3 + 7.020494902488686e7*eta4 - 3.0628992955013716e8*eta5 + 6.979902656542205e8*eta6 - 6.499602051522098e8*eta7) - 0.00612865219736664*chidiff1*(5420.414875733676*eta1 - 209160.45819239676*eta2 + 3.2522689751582462e6*eta3 - 2.632666907701782e7*eta4 + 1.1723210513635881e8*eta5 - 2.723375445677972e8*eta6 + 2.5794037853203726e8*eta7) + chidiff1*(-6.366100477331467*eta1 + 103.71144335659596*eta2 + 303.2593390162942*eta3 - 14537.494074600547*eta4 + 106725.47105224151*eta5 - 330031.4910258338*eta6 + 385149.34959807544*eta7)*S1 - 0.09931510990047862*(1688.1088419484909*eta1 - 63994.0902671211*eta2 + 986285.660380263*eta3 - 7.92705330613304e6*eta4 + 3.507633564255217e7*eta5 - 8.103429147055472e7*eta6 + 7.637928542051452e7*eta7)*S1 - 0.007105414565836703*(41588.7389987755*eta1 - 1.5620802118043897e6*eta2 + 2.38797669731504e7*eta3 - 1.9051633859958e8*eta4 + 8.373734711594706e8*eta5 - 1.9242303107415984e9*eta6 + 1.8079995693160067e9*eta7)*S2;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_rdcp1: version is not valid. Recommended version is 20211005.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_44_rdcp2(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 20211005:{
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double eta5 = eta * eta4;
            double eta6 = eta * eta5;
            double eta7 = eta * eta6;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double S4 = S * S3;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            total = -0.8921739566519082 + 39.19379994645401*eta1 - 556.4162926720528*eta2 + 4054.1083830728735*eta3 - 14813.592901146494*eta4 + 21013.956560307623*eta5 + 0.006110380472895679*chidiff2*(-21805.738224018303*eta1 + 812861.0966083927*eta2 - 1.228168856519143e7*eta3 + 9.656191258303785e7*eta4 - 4.17310096615611e8*eta5 + 9.412641957546971e8*eta6 - 8.669540252837677e8*eta7) + 0.01891676433608453*chidiff1*(1055.8580624803935*eta1 - 49250.520253771*eta2 + 903412.1563472042*eta3 - 8.366109985745795e6*eta4 + 4.163025432262069e7*eta5 - 1.0627603445266993e8*eta6 + 1.0928666073114336e8*eta7) - 0.007083805558059964*chidiff3*(14170.042249201197*eta1 - 586515.9107090512*eta2 + 9.764512135837745e6*eta3 - 8.407745091446699e7*eta4 + 3.960168112384019e8*eta5 - 9.694266415205941e8*eta6 + 9.652363862135582e8*eta7) + chidiff1*(53.025164514922395*eta1 - 2184.5740799026084*eta2 + 35707.02369091528*eta3 - 297471.4895130072*eta4 + 1.3419716256351192e6*eta5 - 3.1284513455984164e6*eta6 + 2.9592474456660985e6*eta7)*S1 - 0.04944793967591508*(1886.9480113920245*eta1 - 69274.01583127167*eta2 + 1.0317037995547823e6*eta3 - 7.997483257662831e6*eta4 + 3.406950546241736e7*eta5 - 7.56370825569486e7*eta6 + 6.83656438797118e7*eta7)*S1 + 0.005043803638745514*(-46105.71566891879*eta1 + 1.628633008707851e6*eta2 - 2.3407638860121876e7*eta3 + 1.7540200556402352e8*eta4 - 7.224896700585026e8*eta5 + 1.5508871808849366e9*eta6 - 1.3558868725072336e9*eta7)*S2 + 0.0056766356748289785*(10272.45823588376*eta1 - 355375.7052892479*eta2 + 4.982579315066213e6*eta3 - 3.615404189277855e7*eta4 + 1.4334132527639598e8*eta5 - 2.949290496426164e8*eta6 + 2.4639243067756915e8*eta7)*S3 + 0.00567212750252267*(40863.9445166239*eta1 - 1.4051105619774864e6*eta2 + 1.971963078789975e7*eta3 - 1.445422284086013e8*eta4 + 5.8300964780643e8*eta5 - 1.2260569388042226e9*eta6 + 1.0499772436697292e9*eta7)*S4;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_rdcp2: version is not valid. Recommended version is 20211005.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_44_rdcp3(double eta, double chi1, double chi2, int RDAmpFlag)
	double total=0;
	switch (RDAmpFlag){
        case 20211005:{
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double eta5 = eta * eta4;
            double eta6 = eta * eta5;
            double eta7 = eta * eta6;
            double S1 = S;
            double S2 = S * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            total = 2.683018365152368 - 114.10242699898323*eta1 + 2002.9406560823122*eta2 - 18127.552659334895*eta3 + 89704.28559863062*eta4 - 231168.3079499347*eta5 + 242784.40116496125*eta6 + 0.00037297780455571646*chidiff1*(-179903.33076495328*eta1 + 6.636154225901366e6*eta2 - 9.9526430350374e7*eta3 + 7.801287504454938e8*eta4 - 3.3787609428612604e9*eta5 + 7.677775735838752e9*eta6 - 7.158915761518967e9*eta7) + 0.0023222484372248053*chidiff3*(42983.27846422164*eta1 - 1.509284011539547e6*eta2 + 2.147543492454528e7*eta3 - 1.5898490678560492e8*eta4 + 6.4758563255312e8*eta5 - 1.3792147851103504e9*eta6 + 1.2022543550908108e9*eta7) - 0.0012012266582370902*chidiff2*(17790.137418972667*eta1 - 794560.4828125148*eta2 + 1.4116476301443396e7*eta3 - 1.285674496848263e8*eta4 + 6.34895272669848e8*eta5 - 1.615830877187729e9*eta6 + 1.6601807612560244e9*eta7) + chidiff1*(75.91210831639803*eta1 - 2738.553293655603*eta2 + 39835.52733198384*eta3 - 298682.741250514*eta4 + 1.2176959493423172e6*eta5 - 2.561495847731494e6*eta6 + 2.17397364484616e6*eta7)*S1 - 0.011497210363361817*(4357.044152560465*eta1 - 162151.86828655298*eta2 + 2.4450811280249325e6*eta3 - 1.918307699163099e7*eta4 + 8.275358575436088e7*eta5 - 1.86340632494114e8*eta6 + 1.7126784118567467e8*eta7)*S1 + 0.004389111594983586*(390.26483617033887*eta1 + 6380.682311913337*eta2 - 419663.76034454996*eta3 + 5.967830385558099e6*eta4 - 3.77092550381064e7*eta5 + 1.1236394121557645e8*eta6 - 1.2873294832443959e8*eta7)*S2;
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_rdcp3: version is not valid. Recommended version is 20211005.");}
  }
  return total;
}

/* End of Parameter Space Fits */


/************** Amplitude Ringdown Ansatz *************/

// For the modes with mixing this is the ansatz of the spheroidal part.
static double IMRPhenomXHM_RD_Amp_Ansatz(double ff, IMRPhenomXHMWaveformStruct *pWF, IMRPhenomXHMAmpCoefficients *pAmp){

    int RDAmpFlag = pWF->IMRPhenomXHMRingdownAmpVersion;
    double frd = pWF->fRING;
    double fda = pWF->fDAMP;
    double dfr = ff - frd;
    double dfd = fda * pAmp->sigma;
    double lc  = pAmp->lc;
    double ampRD;

    switch ( RDAmpFlag )
    {
        case 0: /* Canonical, 3 fitted coefficients + fring, fdamp, lc that are fixed. sigma is also fixed except for the 21 mode. */
        {
            ampRD = (fda *fabs(pAmp->alambda) * pAmp->sigma)*exp(- dfr * pAmp->lambda / dfd )/ (dfr*dfr + dfd*dfd)*pow(ff,-lc);
            break;
        }
        default:
        {
            XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Ringdown_Amp_lm_Ansatz: IMRPhenomXHMRingdownAmpFitsVersion is not valid. Use version 0. \n");
        }
    }

    return ampRD;
}


/*** Derivative of the RD Ansatz for modes without mixing ***/

static double IMRPhenomXHM_RD_Amp_DAnsatz(double ff, IMRPhenomXHMWaveformStruct *pWF, IMRPhenomXHMAmpCoefficients *pAmp){

    int RDAmpFlag = pWF->IMRPhenomXHMRingdownAmpVersion;
    double frd = pWF->fRING;
    double fda = pWF->fDAMP;
    double dfd = fda * pAmp->sigma;
    double lc  = pAmp->lc;
    double lambda = pAmp->lambda;
    double DampRD;
    double numerator,denom,expon;

    switch ( RDAmpFlag )
    {
        case 0:  /* Canonical, 3 fitted coefficients + fring, fdamp, lc that are fixed. sigma is also fixed except for the 21 mode. */
        {
            numerator = fabs(pAmp->alambda)*pow(ff,-1.-lc)*( dfd*lc*(frd*frd + dfd*dfd)
                                                           + ff*( frd*frd*lambda - 2.*dfd*frd*(1.+lc) + dfd*dfd*lambda)
                                                           + ff*ff*( -2.*lambda*frd + dfd*(2.+lc) )
                                                           + ff*ff*ff*( lambda )
                                                           );
            denom = frd*frd + dfd*dfd - ff*2.*frd + ff*ff;
            expon = exp(-(ff-frd)*lambda/dfd);
            DampRD = -numerator*expon/(denom*denom);
            break;
        }
        default:
        {
            XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Ringdown_Amp_lm_Ansatz: IMRPhenomXHMRingdownAmpFitsVersion is not valid. Use version 0. \n");
        }
    }

    return DampRD;
}

/*** Derivative of the RD Ansatz for modes with mixing ***/
// It can not be obtained analytically, so we use finite difference of 4th order
static double IMRPhenomXHM_RD_Amp_NDAnsatz(double ff, IMRPhenomXHMAmpCoefficients *pAmp,  IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXAmpCoefficients *pAmp22,  IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22){

    double df = 10e-10;
    double Nder;
    double fun2R = ff + 2*df;
    double funR  = ff + df;
    double funL  = ff - df;
    double fun2L = ff - 2*df;
    IMRPhenomX_UsefulPowers powers_of_fun;

    IMRPhenomX_Initialize_Powers(&powers_of_fun,fun2R);
    fun2R = cabs(SpheroidalToSpherical(fun2R, &powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    IMRPhenomX_Initialize_Powers(&powers_of_fun,funR);
    funR  = cabs(SpheroidalToSpherical(funR, &powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    IMRPhenomX_Initialize_Powers(&powers_of_fun,funL);
    funL  = cabs(SpheroidalToSpherical(funL, &powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    IMRPhenomX_Initialize_Powers(&powers_of_fun,fun2L);
    fun2L = cabs(SpheroidalToSpherical(fun2L, &powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    Nder = (-fun2R + 8*funR - 8*funL + fun2L )/(12*df);

    return Nder;

}

/* There are cases for some modes where the ringdown amplitude is very low compared to the inspiral and the intermediate reconstruction fails to connect them.
   For those cases we remove the two collocation points and reconstruct with a third order polynomial or with a linear one for the 32 mode. */
void IMRPhenomXHM_Ringdown_Amplitude_Veto(double *V2, double *V3, double V4, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){

    double threshold;
    if(pWFHM->modeTag == 32){
      threshold = 0.1/(pWF22->ampNorm);
      if(1./V4 < threshold){
          *V2 = 1.;
          *V3 = 1.;
          pWFHM->IMRPhenomXHMIntermediateAmpVersion = 101;
      }
    }
    else{
      threshold = 0.01/(pWF22->ampNorm);
      if(1./V4 < threshold){
          *V2 = 1.;
          *V3 = 1.;
          pWFHM->IMRPhenomXHMIntermediateAmpVersion = 1032;
      }
    }
}



/****************************************/
/*                                      */
/*              PHASE                   */
/*                                      */
/****************************************/

/* Fits over parameter space for the ringdown phase quantities. */

// The spin parameter S = (m1^2*chi1 + m2^2*chi2)/(m1^2 + m2^2)

static double IMRPhenomXHM_RD_Phase_22_alpha2(double eta, double chi1, double chi2, int RDPhaseFlag) {
    double total=0;
    switch (RDPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3,eta4,delta=sqrt(1-4*eta);
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            double noSpin = 0.2088669311744758 - 0.37138987533788487*eta + 6.510807976353186*eta2 - 31.330215053905395*eta3 + 55.45508989446867*eta4;
            double eqSpin = ((0.2393965714370633 + 1.6966740823756759*eta - 16.874355161681766*eta2 + 38.61300158832203*eta3)*S)/(1. - 0.633218538432246*S);
            double uneqSpin = (chi1 - 1.*chi2)*(0.9088578269496244*pow(eta,2.5) + 15.619592332008951*(chi1 - 1.*chi2)*pow(eta,3.5))*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_22_alpha2: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Phase_22_alphaL(double eta, double chi1, double chi2, int RDPhaseFlag) {
    double total=0;
    switch (RDPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double delta=sqrt(1.- 4.*eta),eta2,eta3,eta4,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            double noSpin = eta*(-1.1926122248825484 + 2.5400257699690143*eta - 16.504334734464244*eta2 + 27.623649807617376*eta3);
            double eqSpin = eta3*S*(35.803988443700824 + 9.700178927988006*S - 77.2346297158916*S2) + eta*S*(0.1034526554654983 - 0.21477847929548569*S - 0.06417449517826644*S2) + eta2*S*(-4.7282481007397825 + 0.8743576195364632*S + 8.170616575493503*S2) + eta4*S*(-72.50310678862684 - 39.83460092417137*S + 180.8345521274853*S2);
            double uneqSpin = (-0.7428134042821221*chi1*pow(eta,3.5) + 0.7428134042821221*chi2*pow(eta,3.5) + 17.588573345324154*pow(chi1,2)*pow(eta,4.5) - 35.17714669064831*chi1*chi2*pow(eta,4.5) + 17.588573345324154*pow(chi2,2)*pow(eta,4.5))*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_22_alphaL: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

/**************** 32 specific fits ***************/
static double IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift(double eta, double chi1, double chi2, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3,eta4,eta5,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 11.851438981981772 + 167.95086712701223*eta - 4565.033758777737*eta2 + 61559.132976189896*eta3 - 364129.24735853914*eta4 + 739270.8814129328*eta5;
            double eqSpin = (9.506768471271634 + 434.31707030999445*eta - 8046.364492927503*eta2 + 26929.677144312944*eta3)*S + (-5.949655484033632 - 307.67253970367034*eta + 1334.1062451631644*eta2 + 3575.347142399199*eta3)*S2 + (3.4881615575084797 - 2244.4613237912527*eta + 24145.932943269272*eta2 - 60929.87465551446*eta3)*S3 + (15.585154698977842 - 2292.778112523392*eta + 24793.809334683185*eta2 - 65993.84497923202*eta3)*S4;
            double uneqSpin = 465.7904934097202*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift(double eta, double chi1, double chi2, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3,eta4,eta5,eta6,eta7,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            eta7 = pow(eta,7);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = -1.3328895897490733 - 22.209549522908667*eta + 1056.2426481245027*eta2 - 21256.376324666326*eta3 + 246313.12887984765*eta4 - 1.6312968467540336e6*eta5 + 5.614617173188322e6*eta6 - 7.612233821752137e6*eta7;
            double eqSpin = (S*(-1.622727240110213 + 0.9960210841611344*S - 1.1239505323267036*S2 - 1.9586085340429995*S3 + eta2*(196.7055281997748 + 135.25216875394943*S + 1086.7504825459278*S2 + 546.6246807461155*S3 - 312.1010566468068*S4) + 0.7638287749489343*S4 + eta*(-47.475568056234245 - 35.074072557604445*S - 97.16014978329918*S2 - 34.498125910065156*S3 + 24.02858084544326*S4) + eta3*(62.632493533037625 - 22.59781899512552*S - 2683.947280170815*S2 - 1493.177074873678*S3 + 805.0266029288334*S4)))/(-2.950271397057221 + 1.*S);
            double uneqSpin = (sqrt(1. - 4.*eta)*(chi2*pow(eta,2.5)*(88.56162028006072 - 30.01812659282717*S) + chi2*eta2*(43.126266433486435 - 14.617728550838805*S) + chi1*eta2*(-43.126266433486435 + 14.617728550838805*S) + chi1*pow(eta,2.5)*(-88.56162028006072 + 30.01812659282717*S)))/(-2.950271397057221 + 1.*S);
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

// collocation points
static double IMRPhenomXHM_Ringdown_Phase_32_p1(double eta, double chi1, double chi2, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3,eta4,eta5,S2,S3,S4,S5;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            S5 = pow(S,5);
            double noSpin = 3169.372056189274 + 426.8372805022653*eta - 12569.748101922158*eta2 + 149846.7281073725*eta3 - 817182.2896823225*eta4 + 1.5674053633767858e6*eta5;
            double eqSpin = (19.23408352151287 - 1762.6573670619173*eta + 7855.316419853637*eta2 - 3785.49764771212*eta3)*S + (-42.88446003698396 + 336.8340966473415*eta - 5615.908682338113*eta2 + 20497.5021807654*eta3)*S2 + (13.918237996338371 + 10145.53174542332*eta - 91664.12621864353*eta2 + 201204.5096556517*eta3)*S3 + (-24.72321125342808 - 4901.068176970293*eta + 53893.9479532688*eta2 - 139322.02687945773*eta3)*S4 + (-61.01931672442576 - 16556.65370439302*eta + 162941.8009556697*eta2 - 384336.57477596396*eta3)*S5;
            double uneqSpin = (chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2*(641.2473192044652 - 1600.240100295189*chi1*eta + 1600.240100295189*chi2*eta + 13275.623692212472*eta*S);
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Ringdown_Phase_32_p1: version is not valid. Recommended version is 122019.");}
    }
    return total;
}
// collocation points
static double IMRPhenomXHM_Ringdown_Phase_32_p2(double eta, double chi1, double chi2, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3131.0260952676376 + 206.09687819102305*eta - 2636.4344627081873*eta2 + 7475.062269742079*eta3;
            double eqSpin = (49.90874152040307 - 691.9815135740145*eta - 434.60154548208334*eta2 + 10514.68111669422*eta3)*S + (97.3078084654917 - 3458.2579971189534*eta + 26748.805404989867*eta2 - 56142.13736008524*eta3)*S2 + (-132.49105074500454 + 429.0787542102207*eta + 7269.262546204149*eta2 - 27654.067482558712*eta3)*S3 + (-227.8023564332453 + 5119.138772157134*eta - 34444.2579678986*eta2 + 69666.01833764123*eta3)*S4;
            double uneqSpin = 477.51566939885424*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Ringdown_Phase_32_p2: version is not valid. Recommended version is 122019.");}
    }
    return total;
}
// collocation points
static double IMRPhenomXHM_Ringdown_Phase_32_p3(double eta, double chi1, double chi2, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3082.803556599222 + 76.94679795837645*eta - 586.2469821978381*eta2 + 977.6115755788503*eta3;
            double eqSpin = (45.08944710349874 - 807.7353772747749*eta + 1775.4343704616288*eta2 + 2472.6476419567534*eta3)*S + (95.57355060136699 - 2224.9613131172046*eta + 13821.251641893134*eta2 - 25583.314298758105*eta3)*S2 + (-144.96370424517866 + 2268.4693587493093*eta - 10971.864789147161*eta2 + 16259.911572457446*eta3)*S3 + (-227.8023564332453 + 5119.138772157134*eta - 34444.2579678986*eta2 + 69666.01833764123*eta3)*S4;
            double uneqSpin = 378.2359918274837*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Ringdown_Phase_32_p3: version is not valid. Recommended version is 122019.");}
    }
    return total;
}
// collocation points
static double IMRPhenomXHM_Ringdown_Phase_32_p4(double eta, double chi1, double chi2, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,S2,S3,S4;
            eta2 = pow(eta,2);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3077.0657367004565 + 64.99844502520415*eta - 357.38692756785395*eta2;
            double eqSpin = (34.793450080444714 - 986.7751755509875*eta - 9490.641676924794*pow(eta,3) + 5700.682624203565*eta2)*S + (57.38106384558743 - 1644.6690499868596*eta - 19906.416384606226*pow(eta,3) + 11008.881935880598*eta2)*S2 + (-126.02362949830213 + 3169.3397351803583*eta + 62863.79877094988*pow(eta,3) - 26766.730897942085*eta2)*S3 + (-169.30909412804587 + 4900.706039920717*eta + 95314.99988114933*pow(eta,3) - 41414.05689348732*eta2)*S4;
            double uneqSpin = 390.5443469721231*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Ringdown_Phase_32_p4: version is not valid. Recommended version is 122019.");}
    }
    return total;
}


/**************  ANSATZ PHASE DERIVATIVE **************/

static double IMRPhenomXHM_RD_Phase_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase){

    double invf2 = powers_of_f->m_two;
    double frd   = pWFHM->fRING;
    double fda   = pWFHM->fDAMP;
    double dphaseRD;

    switch ( pWFHM->MixingOn )
    {
        case 0:
        {
            /*  rescaling of the 22 ansatz -- used for (21),(33),(44) */
            /*ansatz:
             alpha0 + ((fRDlm^2) alpha2)/(f^2)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)*/
            dphaseRD = ( pPhase->alpha0 +  frd*frd*(pPhase->alpha2)*invf2 + ( (pPhase->alphaL)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
            break;

        }
        case 1:
        {
            /*  calibration of spheroidal ringdown waveform for (32) */
            /* ansatz: alpha0 + (alpha2)/(f^2)+ (alpha4)/(f^4)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)*/
            double invf4 = powers_of_f->m_four;
            dphaseRD = ( pPhase->alpha0_S +  (pPhase->alpha2_S)*invf2 + (pPhase->alpha4_S)*invf4 +( (pPhase->alphaL_S)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
            break;
        }
        default:
        {XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_RD_Phase_Ansatz: version is not valid. Use version 0 for modes (2,1),(3,3),(4,4) and 1 for (3,2).\n");}
    }
    return dphaseRD;
}

/**************  ANSATZ INTEGRATED PHASE **************/

static double IMRPhenomXHM_RD_Phase_AnsatzInt(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase){

    double invf   = powers_of_f->m_one;
    double phaseRD;
    double frd=pWFHM->fRING;
    double fda=pWFHM->fDAMP;

    switch ( pWFHM->MixingOn )
    {
        case 0:
        {
            /*  rescaling of the 22 ansatz -- used for (21),(33),(44) */
            /*ansatz:
             alpha0 f - fRDlm^2*alpha2)/f  + alphaL*ArcTan[(f - fRDlm)/fdamplm]*/
            phaseRD = pPhase->alpha0*ff -frd*frd *(pPhase->alpha2)*invf +  (pPhase->alphaL)* atan((ff-frd)/fda);
            break;
        }
        case 1:
        {
            /*  calibration of spheroidal ringdown waveform for (32) */
            /* ansatz: f alpha0 - (alpha4)/(3 f^3) - (alpha2)/f + alphaL ArcTan[(f - fRDlm)/fdamplm]*/
            double invf3 = powers_of_f->m_three;
            phaseRD = pPhase->phi0_S+pPhase->alpha0_S*ff -(pPhase->alpha2_S)*invf -1./3.*(pPhase->alpha4_S)*invf3 +(pPhase->alphaL_S)* atan((ff-frd)/fda);
            break;
        }
        default:
        {XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_RD_Phase_AnsatzInt: version is not valid. Use version 0 for modes (2,1),(3,3),(4,4) and 1 for (3,2).\n");}
    }
    return phaseRD;
}
