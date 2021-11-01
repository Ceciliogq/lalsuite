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
 */
 //
//  LALSimIMRPhenomXHM_inspiral.c
//
//
//  Created by Marta on 06/02/2019.
//

#include "LALSimIMRPhenomXHM_inspiral.h"

/*************************************/
/*                                   */
/*            AMPLITUDE              */
/*                                   */
/*************************************/

/* Fits of the collocation points in the inspiral across parameter space. 3 for each mode: iv1, iv2, iv3
 [iv1,iv2,iv3]=[0.5*f^Ins_lm,0.75*f^Ins_lm,f^Ins_lm]
 For more details about the reconstruction procedure, see Sec. IV.A
 */

 /* These parameter-space fits are documented in the supplementary material of https://dcc.ligo.org/P2000011-v2.
    There are 2 Mathematica notebooks (one for amplitude and one for phase) that read the fits data and automatically generate the C-code below.
    For more information read https://git.ligo.org/waveforms/reviews/imrphenomx/blob/master/documentation/ParspaceFits/README and the documentation in the notebooks. */

// Spin parameter S = chi_PN = chi_eff - 38/113 * eta * (chi1 + chi2)
// chi_eff = (m1*chi1 + m2*chi2)/(m1 + m2)



static double IMRPhenomXHM_Insp_Amp_21_iv1(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = sqrt(1. - 4.*eta)*(0.037868557189995156 + 0.10740090317702103*eta + 1.963812986867654*eta2 - 16.706455229589558*eta3 + 69.75910808095745*eta4 - 98.3062466823662*eta5);
            double eqSpin = sqrt(1. - 4.*eta)*S*(-0.007963757232702219 + 0.10627108779259965*eta - 0.008044970210401218*S + eta2*(-0.4735861262934258 - 0.5985436493302649*S - 0.08217216660522082*S2));
            double uneqSpin = -0.257787704938017*(chi1 - 1.*chi2)*eta2*(1. + 8.75928187268504*eta2) - 0.2597503605427412*(chi1 - 1.*chi2)*eta2*S;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 20211005:{
            double sqroot = sqrt(eta);
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double S1 = S;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            total = fabs(0.09618927141811567*chidiff3*(51.28678251481265 - 936.9600668368321*eta1 + 5375.783716041888*eta2 - 9799.968425404946*eta3) + 0.2869008784518925*chidiff2*(10.594431595445412 - 182.91407634736888*eta1 + 1171.3908669139817*eta2 - 2441.541434213661*eta3) - 4.178501993823551*chidiff1*(1.0048009016862354 - 1.9574345147433532*eta1 + 25.42833395476959*eta2 - 75.22693603906846*eta3) + chidiff1*(1.825328803570667 - 21.394684272995*eta1 + 146.3486729478616*eta2 - 298.3480010605445*eta3)*S1 + delta*(20.335650259456585 - 4.368413695663519*eta1 + 3.913451522666885*eta2)*sqroot + delta*S1*(-16.536382404976695*(1.0681605787881077 - 1.9002932308566096*eta1 + 12.835261670662922*eta2 - 25.488872335203112*eta3) + 1.824623400650354*(-1.476577171602304 + 46.84318703646066*eta1 - 286.64649084684925*eta2 + 567.0127660847531*eta3)*S1)*sqroot);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_21_iv1: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_21_iv2(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,eta3,eta4,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            double noSpin = sqrt(1. - 4.*eta)*(0.05511628628738656 - 0.12579599745414977*eta + 2.831411618302815*eta2 - 14.27268643447161*eta3 + 28.3307320191161*eta4);
            double eqSpin = sqrt(1. - 4.*eta)*S*(-0.008692738851491525 + eta*(0.09512553997347649 + 0.116470975986383*S) - 0.009520793625590234*S + eta2*(-0.3409769288480959 - 0.8321002363767336*S - 0.13099477081654226*S2) - 0.006383232900211555*S2);
            double uneqSpin = -0.2962753588645467*(chi1 - 1.*chi2)*eta2*(1. + 1.3993978458830476*eta2) - 0.17100612756133535*(chi1 - 1.*chi2)*eta2*S*(1. + 18.974303741922743*eta2*delta);
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
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
            double chidiff4 = chidiff * chidiff3;
            total = fabs(0.03320820321337354*chidiff4*(186.04907601572523 - 2962.7993051419835*eta1 + 15256.056010309754*eta2 - 25518.163169201693*eta3) + 0.04121452943443768*chidiff3*(109.78128883103665 - 1801.840754573055*eta1 + 9496.192077089061*eta2 - 16158.284011524627*eta3) + 0.09415738550484183*chidiff2*(-6.1242303133396145 + 81.15411802232259*eta1 - 192.1850757620191*eta2 - 129.3425617198742*eta3) - 2.9392787904093005*chidiff1*(0.747519780351714 + 0.6795726347887833*eta1 + 11.964083108531323*eta2 - 44.59245982741375*eta3) + chidiff1*(0.2827151630503135 - 0.38805556359477933*eta1 + 10.345922166447316*eta2 - 19.400992727701283*eta3)*S1 + delta*(12.331458520468804 + 4.70759595852107*eta1 - 4.5247922031129155*eta2)*sqroot + delta*S1*(-10.751647771786997*(0.9870891733934575 - 2.169433875805002*eta1 + 20.160064111201418*eta2 - 45.612753440833906*eta3) + 0.9676910374370745*(-0.9098423902709749 + 45.191683460044366*eta1 - 322.813012719325*eta2 + 711.5613549845066*eta3)*S1 - 0.10048430896719716*(-8.944059244043336 + 242.3555828340873*eta1 - 1835.2582168980891*eta2 + 4257.113822264582*eta3)*S2 - 0.1734779903213141*(4.106169368377872 + 100.084476561933*eta1 - 1340.4288908114488*eta2 + 3722.1696745789427*eta3)*S3)*sqroot);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_21_iv2: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_21_iv3(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,S2;
            eta2 = pow(eta,2);
            S2 = pow(S,2);
            double noSpin = sqrt(1. - 4.*eta)*(0.059110044024271766 - 0.0024538774422098405*eta + 0.2428578654261086*eta2);
            double eqSpin = sqrt(1. - 4.*eta)*S*(-0.007044339356171243 - 0.006952154764487417*S + eta2*(-0.016643018304732624 - 0.12702579620537421*S + 0.004623467175906347*S2) - 0.007685497720848461*S2);
            double uneqSpin = -0.3172310538516028*(chi1 - 1.*chi2)*(1. - 2.9155919835488024*eta2)*eta2 - 0.11975485688200693*(chi1 - 1.*chi2)*eta2*S*(1. + 17.27626751837825*eta2*delta);
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
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
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            double chidiff4 = chidiff * chidiff3;
            total = fabs(-0.018177418073618133*chidiff4*(25.946752764599115 - 726.1447914232601*eta1 + 5343.876388980924*eta2 - 11616.284569767238*eta3) - 2.1647287364246477*chidiff1*(0.17641101512929638 + 6.487675576000363*eta1 - 12.762792675197355*eta2 + 1.106682870723617*eta3) + 0.0044191142462998035*chidiff3*(37.109918412415716 - 357.3737444492902*eta1 + 859.1370929216105*eta2 + 80.1850916194976*eta3) - 0.009838956636928556*chidiff2*(-59.95523144992317 + 1517.3089354542278*eta1 - 10021.658931664717*eta2 + 19926.504715417363*eta3) + chidiff1*(-0.03768957318623445 - 3.2309644745095123*eta1 + 26.00184618973167*eta2 - 41.40286253552052*eta3)*S1 + delta*(7.142091514368739 + 16.990644229645053*eta1 - 22.790661406974447*eta2)*sqroot + delta*S1*(-7.129105192514667*(0.38872631299582927 + 6.127406154470924*eta1 - 22.25974184077449*eta2 + 31.79117838370405*eta3) - 0.058359125701488046*(5.609225151837126 - 39.789320791123885*eta1 + 202.96379867262016*eta2 - 527.6189899676901*eta3)*S1 - 0.2782544009729164*(-1.447678955057067 + 21.68945100454068*eta1 - 138.85720355511737*eta2 + 408.69077793565464*eta3)*S2)*sqroot);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_21_iv3: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_33_iv1(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double eta2,S2;
            eta2 = pow(eta,2);
            S2 = pow(S,2);
            double noSpin = (sqrt(1. - 4.*eta)*(-0.056586690934283326 - 0.14374841547279146*eta + 0.5584776628959615*eta2))/(-0.3996185676368123 + eta);
            double eqSpin = sqrt(1. - 4.*eta)*S*((0.056042044149691175 + 0.12482426029674777*S)*S + eta*(2.1108074577110343 - 1.7827773156978863*S2) + eta2*(-7.657635515668849 - 0.07646730296478217*S + 5.343277927456605*S2));
            double uneqSpin = 0.45866449225302536*(chi1 - 1.*chi2)*(1. - 9.603750707244906*eta2)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
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
            double S1 = S;
            double S2 = S * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            total = 0.0024267979273393334*chidiff2*(-1347.4738797665993*eta1 + 53156.52198359776*eta2 - 565742.4998524422*eta3 + 2.3439328146050153e6*eta4 - 3.3809063203941695e6*eta5) - 0.4149986591903826*chidiff1*(116.61932404124704*eta1 - 1064.1301078818083*eta2 + 4067.3539268594977*eta3 - 8331.534333694231*eta4 + 6120.913248785117*eta5) + chidiff1*(22.363751444775524*eta1 - 443.49927745657635*eta2 + 3456.682699403338*eta3 - 12157.679774497858*eta4 + 16009.245888692432*eta5)*S1 + delta*(23.043576100563094 - 9.362533181914708*eta1 + 46.155587153990446*eta2 - 82.90080892085865*eta3)*sqroot + delta*eta1*S1*(-0.004923965617227187*(-315.3797195474725 + 7069.9086514657565*eta1 - 57489.257115876826*eta2 + 205949.18978557832*eta3 - 275688.6458064782*eta4) - 9.767567865244697*(30.14678718484687 - 364.89941569793604*eta1 + 2162.9638637005564*eta2 - 6203.842964455104*eta3 + 6851.948685599858*eta4) + 0.21710270374100293*(-61.67924175183982 + 2295.367117244805*eta1 - 22269.203439238547*eta2 + 87179.45997325864*eta3 - 123383.90347315231*eta4)*S1 - 0.2030055801411031*(-315.68478259857017 + 7178.32640644418*eta1 - 59871.58346979412*eta2 + 214790.4216896625*eta3 - 276986.5537130501*eta4)*S2)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_33_iv1: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_33_iv2(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double eta2,eta3,eta4,eta5,eta6,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            double noSpin = sqrt(1. - 4.*eta)*(0.2137734510411439 - 0.7692194209223682*eta + 26.10570221351058*eta2 - 316.0643979123107*eta3 + 2090.9063511488234*eta4 - 6897.3285171507105*eta5 + 8968.893362362503*eta6);
            double eqSpin = sqrt(1. - 4.*eta)*S*(0.018546836505210842 + 0.05924304311104228*S + eta*(1.6484440612224325 - 0.4683932646001618*S - 2.110311135456494*S2) + 0.10701786057882816*S2 + eta2*(-6.51575737684721 + 1.6692205620001157*S + 8.351789152096782*S2));
            double uneqSpin = 0.3929315188124088*(chi1 - 1.*chi2)*(1. - 11.289452844364227*eta2)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
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
            double S1 = S;
            double S2 = S * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            double chidiff3 = chidiff * chidiff2;
            total = chidiff1*(-0.8807843425190279*eta1 + 71.91536569294988*eta2 - 447.98728399486777*eta3 + 1142.4296192866352*eta4 - 1617.7496275708356*eta5) + chidiff3*(0.535156492326376*eta1 - 4.918611496801541*eta2 - 34.72561668518005*eta3 + 437.8414792978095*eta4 - 1020.9753643711776*eta5) + chidiff2*(1.2509879894003024*eta1 - 38.53669796626464*eta2 + 419.8092433712076*eta3 - 1914.6316515537737*eta4 + 3087.6753029984084*eta5) + chidiff1*(-0.7280929454571974*eta1 + 29.573257930221125*eta2 - 408.54982862391313*eta3 + 2114.225532364241*eta4 - 3589.842210439146*eta5)*S1 + delta*(13.258985792886566 + 6.710713920515631*eta1 - 2.2564912763507476*eta2)*sqroot + delta*eta1*S1*(-0.002283406076482348*(-440.2184482145835 + 13640.299646910553*eta1 - 148614.83932691297*eta2 + 682425.8478295282*eta3 - 1.1220757655166206e6*eta4) - 5.244962119966554*(24.736997888518 - 283.06805561001045*eta1 + 1668.2498946537862*eta2 - 4823.178614618638*eta3 + 5418.386019846486*eta4) + 0.11925207859450397*(-191.51461256280646 + 6802.885442609262*eta1 - 74724.93914018388*eta2 + 337982.9038238593*eta3 - 547277.7435123762*eta4)*S1 - 0.10889881177569258*(325.7241660741732 - 8320.864863984167*eta1 + 80001.91956153851*eta2 - 337109.0777453007*eta3 + 523468.8048624649*eta4)*S2)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_33_iv2: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_33_iv3(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double eta2,eta3,eta4,eta5,eta6,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            double noSpin = sqrt(1. - 4.*eta)*(0.2363760327127446 + 0.2855410252403732*eta - 10.159877125359897*eta2 + 162.65372389693505*eta3 - 1154.7315106095564*eta4 + 3952.61320206691*eta5 - 5207.67472857814*eta6);
            double eqSpin = sqrt(1. - 4.*eta)*S*(0.04573095188775319 + 0.048249943132325494*S + eta*(0.15922377052827502 - 0.1837289613228469*S - 0.2834348500565196*S2) + 0.052963737236081304*S2);
            double uneqSpin = 0.25187274502769835*(chi1 - 1.*chi2)*(1. - 12.172961866410864*eta2)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 20211004:{
            double sqroot = sqrt(eta);
            double delta = sqrt(1.-4*eta);
            double S = (chi1 + chi2)/2. + ((chi1 - chi2)*delta)/(1 + delta*delta);
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double S1 = S;
            double chidiff1 = chidiff;
            total = chidiff1*(0.3822941361891989*eta1 + 23.394725970792088*eta2 - 89.02968826588216*eta3) + chidiff1*(0.6888482819246591*eta1 - 6.152942932618743*eta2 + 16.762103483576162*eta3)*S1 + delta*(7.08666488892205 + 28.898258925160643*eta1 - 84.12055070513236*eta2 + 131.9798563985693*eta3)*sqroot - 2.772652493812293*delta*eta1*(8.824582778492038 - 26.44016956450907*eta1 + 31.407961775827033*eta2)*S1*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_33_iv3: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_32_iv1(double eta, double chi1, double chi2, int InspAmpFlag) {
  double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,eta3,eta4,eta5,eta6,eta8,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            eta8 = pow(eta,8);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = sqrt(1. - 3.*eta)*(0.019069933430190773 - 0.19396651989685837*eta + 11.95224600241255*eta2 - 158.90113442757382*eta3 + 1046.65239329071*eta4 - 3476.940285294999*eta5 + 4707.249209858949*eta6);
            double eqSpin = sqrt(1. - 3.*eta)*S*(0.0046910348789512895 + 0.40231360805609434*eta - 0.0038263656140933152*S + 0.018963579407636953*S2 + eta2*(-1.955352354930108 + 2.3753413452420133*S - 0.9085620866763245*S3) + 0.02738043801805805*S3 + eta3*(7.977057990568723 - 7.9259853291789515*S + 0.49784942656123987*S2 + 5.2255665027119145*S3));
            double uneqSpin = 0.058560321425018165*pow(chi1 - 1.*chi2,2)*(1. - 19.936477485971217*eta2)*eta2 + 1635.4240644598524*(chi1 - 1.*chi2)*eta8*delta + 0.2735219358839411*(chi1 - 1.*chi2)*eta2*S*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 202109302:{
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
            total = 0.4069095467494331 + 7.6463652286543295*eta1 - 47.58638982589533*eta2 + 0.0033897148952675015*chidiff2*(-34.41541801301956*eta1 + 559.3086910265454*eta2 - 1852.3485360675372*eta3) + 96.77253600595304*eta3 - 0.04260150130176263*chidiff1*(63.22305559530258*eta1 - 476.6886977694978*eta2 + 895.6987326524188*eta3) - 97.91513270775462*eta4 + 0.18525923875120004*(-11.714278168387857*eta1 + 116.63177104985033*eta2 - 137.1262739946879*eta3)*S1 + chidiff1*(-0.1631961479221324*eta1 + 0.8758711426146618*eta2 - 0.8718841849857611*eta3)*S1 - 0.024067112467842263*(4.047640308654183*eta1 + 17.53141472533972*eta2 - 56.93005472204932*eta3)*S2;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_32_iv1: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_32_iv2(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,eta3,eta4,eta8,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta8 = pow(eta,8);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = sqrt(1. - 3.*eta)*(0.024621376891809633 - 0.09692699636236377*eta + 2.7200998230836158*eta2 - 16.160563094841066*eta3 + 32.930430889650836*eta4);
            double eqSpin = sqrt(1. - 3.*eta)*S*(0.008522695567479373 - 1.1104639098529456*eta2 - 0.00362963820787208*S + 0.016978054142418417*S2 + eta*(0.24280554040831698 + 0.15878436411950506*S - 0.1470288177047577*S3) + 0.029465887557447824*S3 + eta3*(4.649438233164449 - 0.7550771176087877*S + 0.3381436950547799*S2 + 2.5663386135613093*S3));
            double uneqSpin = -0.007061187955941243*pow(chi1 - 1.*chi2,2)*(1. - 2.024701925508361*eta2)*eta2 + 215.06940561269835*(chi1 - 1.*chi2)*eta8*delta + 0.1465612311350642*(chi1 - 1.*chi2)*eta2*S*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 202109302:{
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double S4 = S * S3;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            total = 0.3008854906650421 + 4.79432244822206*eta1 - 23.153246494084026*eta2 + 0.007783904371774833*chidiff1*(-181.35755037527343*eta1 + 1775.5144589231934*eta2 - 4159.091977656661*eta3) + 0.003803955185577416*chidiff2*(22.039759280696874*eta1 + 53.36288966295568*eta2 - 759.8070829274377*eta3) + 21.170745266656343*eta3 + 0.19454219211046542*(-3.883622985627817*eta1 + 56.04312382527554*eta2 - 36.72545843990798*eta3)*S1 + chidiff1*(-0.0842819829481105*eta1 + 4.118506271695715*eta2 - 14.988761908562413*eta3)*S1 + 0.012309632443087231*(42.01567432049474*eta1 - 320.53645408071924*eta2 + 666.1189447956892*eta3)*S2 + 0.014017824614215165*(-16.99550693587934*eta1 + 241.59873005742887*eta2 - 653.3951037597933*eta3)*S3 + 0.005213026697840646*(-101.14975208681948*eta1 + 895.1325551732767*eta2 - 1745.1880867509071*eta3)*S4;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_32_iv2: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_32_iv3(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,eta3,eta8,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta8 = pow(eta,8);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = (sqrt(1. - 3.*eta)*(-0.006150151041614737 + 0.017454430190035*eta + 0.02620962593739105*eta2 - 0.019043090896351363*eta3))/(-0.2655505633361449 + eta);
            double eqSpin = sqrt(1. - 3.*eta)*S*(0.011073381681404716 + 0.00347699923233349*S + eta*S*(0.05592992411391443 - 0.15666140197050316*S2) + 0.012079324401547036*S2 + eta2*(0.5440307361144313 - 0.008730335213434078*S + 0.04615964369925028*S2 + 0.6703688097531089*S3) + 0.016323101357296865*S3);
            double uneqSpin = -0.020140175824954427*pow(chi1 - 1.*chi2,2)*(1. - 12.675522774051249*eta2)*eta2 - 417.3604094454253*(chi1 - 1.*chi2)*eta8*delta + 0.10464021067936538*(chi1 - 1.*chi2)*eta2*S*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 202109302:{
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
            double chidiff4 = chidiff * chidiff3;
            total = 0.14066540287008636 + 4.807597003793498*eta1 - 21.214001631235412*eta2 + 0.004126312632112899*chidiff3*(-534.0437600435655*eta1 + 5154.281071470763*eta2 - 12147.217280678944*eta3) + 0.04373945550372854*chidiff1*(9.599190348713053*eta1 + 14.1913689000324*eta2 - 197.02723368208484*eta3) + 20.427725665433723*eta3 + 0.014196753602525978*chidiff2*(107.08641555302717*eta1 - 870.229885982198*eta2 + 1746.9161643226535*eta3) - 0.011810992850764899*chidiff4*(353.5333407112912*eta1 - 3232.9566920980556*eta2 + 7338.079925772161*eta3) + chidiff1*(-0.31978252407346885*eta1 + 9.23980920135081*eta2 - 31.483250322052058*eta3)*S1 + 0.2251940005638608*(3.157047284081725*eta1 + 2.538558495639996*eta2 + 48.364268979743535*eta3)*S1 + chidiff2*(-5.561447567624153*eta1 + 58.658642741390274*eta2 - 157.22158385358009*eta3)*S2 + 0.037851592504054385*(15.454353806705347*eta1 - 110.67078509657317*eta2 + 302.3573826058085*eta3)*S2 + 0.011757923978727484*(-37.3963726826167*eta1 + 450.59484812078193*eta2 - 1187.7497744127027*eta3)*S3;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_32_iv3: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_44_iv1(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double eta2,eta3,eta4,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            double noSpin = sqrt(1. - 3.*eta)*(0.06190013067931406 + 0.1928897813606222*eta + 1.9024723168424225*eta2 - 15.988716302668415*eta3 + 35.21461767354364*eta4);
            double eqSpin = sqrt(1. - 3.*eta)*S*(0.011454874900772544 + 0.044702230915643903*S + eta*(0.6600413908621988 + 0.12149520289658673*S - 0.4482406547006759*S2) + 0.07327810908370004*S2 + eta2*(-2.1705970511116486 - 0.6512813450832168*S + 1.1237234702682313*S2));
            double uneqSpin = 0.4766851579723911*(chi1 - 1.*chi2)*(1. - 15.950025762198988*eta2)*eta2 + 0.127900699645338*pow(chi1 - 1.*chi2,2)*(1. - 15.79329306044842*eta2)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 20211005:{
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double eta4 = eta * eta3;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double S4 = S * S3;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            total = 0.5494362280069187 + 16.03267002586217*eta1 - 120.54748420070754*eta2 - 0.026614311102952082*chidiff2*(-31.834811683204887*eta1 + 417.7035136017306*eta2 - 1170.524836516587*eta3) + 370.23728365458703*eta3 - 0.13092341449661188*chidiff1*(38.73191289424052*eta1 - 250.28517626276223*eta2 + 384.8502325578163*eta3) - 525.5854217161606*eta4 + chidiff1*(-0.22834314965196695*eta1 - 0.8410747803180845*eta2 + 6.948032445473613*eta3)*S1 - 0.29923593493253414*(20.18377586808974*eta1 - 107.31322068116003*eta2 + 149.92458931089345*eta3)*S1 - 0.06301028169570495*(-18.645447322881797*eta1 + 275.7505496529843*eta2 - 788.6065039280193*eta3)*S2 - 0.0017886550996462738*(-675.6809286860356*eta1 + 6084.185191830418*eta2 - 13077.915347877959*eta3)*S3 + 0.03762905088957379*(-14.691432843801525*eta1 + 264.3051807401833*eta2 - 840.0342815758083*eta3)*S4;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_44_iv1: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_44_iv2(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,eta3,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 0.08406011695496626 - 0.1469952725049322*eta + 0.2997223283799925*eta2 - 1.2910560244510723*eta3;
            double eqSpin = (0.023924074703897662 + 0.26110236039648027*eta - 1.1536009170220438*eta2)*S + (0.04479727299752669 - 0.1439868858871802*eta + 0.05736387085230215*eta2)*S2 + (0.06028104440131858 - 0.4759412992529712*eta + 1.1090751649419717*eta2)*S3;
            double uneqSpin = 0.10346324686812074*pow(chi1 - 1.*chi2,2)*(1. - 16.135903382018213*eta2)*eta2 + 0.2648241309154185*(chi1 - 1.*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 20211005:{
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double S4 = S * S3;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            total = 0.3381544232049437 + 10.06008429623452*eta1 - 53.681257236298705*eta2 + 0.04295369903619075*chidiff2*(-24.782780358758952*eta1 + 349.21124008765423*eta2 - 1005.8238243368642*eta3) + 69.35931628999721*eta3 - 0.05563960189253639*chidiff1*(47.17257907206492*eta1 - 326.7877848639625*eta2 + 552.7740849283236*eta3) + chidiff1*(0.6896400934055237*eta1 - 5.987778010071329*eta2 + 12.826137446833702*eta3)*S1 - 0.17675187457986588*(17.393952522394144*eta1 - 83.35003387484079*eta2 + 101.47590151031*eta3)*S1 + 0.06336946652883135*(-37.07553727399484*eta1 + 469.6840045480143*eta2 - 1298.6845936476304*eta3)*S2 + chidiff2*(4.390737638519823*eta1 - 65.69733914675669*eta2 + 193.35200320705968*eta3)*S2 - 0.0024278468881925976*(-81.86999210994276*eta1 + 734.0284496244753*eta2 - 1434.3682520551265*eta3)*S3 - 0.07459672942277429*(-37.34101772582944*eta1 + 461.7445664640425*eta2 - 1248.7894799031521*eta3)*S4;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_44_iv2: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_44_iv3(double eta, double chi1, double chi2, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double S = XLALSimIMRPhenomXchiPNHat(eta,chi1,chi2);
            double delta=sqrt(1. - 4.*eta),eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = 0.08212436946985402 - 0.025332770704783136*eta - 3.2466088293309885*eta2 + 28.404235115663706*eta3 - 111.36325359782991*eta4 + 157.05954559045156*eta5;
            double eqSpin = S*(0.03488890057062679 + 0.039491331923244756*S + eta*(-0.08968833480313292 - 0.12754920943544915*S - 0.11199012099701576*S2) + 0.034468577523793176*S2);
            double uneqSpin = 0.2062291124580944*(chi1 - 1.*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 20211005:{
            double S = (chi1 + chi2)/2.;
            double chidiff = (chi1 - chi2)/2.;
            double eta1 = eta;
            double eta2 = eta * eta1;
            double eta3 = eta * eta2;
            double S1 = S;
            double S2 = S * S1;
            double S3 = S * S2;
            double S4 = S * S3;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff * chidiff1;
            total = 0.2876827874457881 + 5.466538383405843*eta1 - 18.679351622365438*eta2 - 0.028349922071010156*chidiff2*(-35.88619456124589*eta1 + 461.2943776115792*eta2 - 1287.770109371502*eta3) - 0.0008396201927794044*chidiff1*(28.98138903346638*eta1 - 15.398743378180084*eta2 - 590.5291752312685*eta3) - 0.07796322375825818*(-1.7422987766739346*eta1 + 88.03615861551013*eta2 - 268.26375482029175*eta3)*S1 + chidiff1*(0.14125942169694206*eta1 - 2.1136003675110677*eta2 + 6.324825434451538*eta3)*S1 - 0.0650220088381942*(-24.842806748756406*eta1 + 336.41007298898455*eta2 - 935.8435158284761*eta3)*S2 + chidiff2*(2.48622985311047*eta1 - 7.206595615272273*eta2 - 12.171414344611799*eta3)*S2 - 0.012129918434263021*(-31.007170686661848*eta1 + 376.1550710902337*eta2 - 972.4749314680503*eta3)*S3 + 0.04693318643539241*(-10.852895526339623*eta1 + 219.45391706541372*eta2 - 710.0126034515091*eta3)*S4;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_44_iv3: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

/* End of Parameter Space Fits */


/***********************************************/
/*     Pseudo PN Amplitude Coefficients    */
/***********************************************/

/* The whole inspiral ansatz is given by

    PNAnsatz(f) + pseudo-PN(f),  with pseudo-PN(f) = rho1 *(f/fcutInsp)^(7/3) + rho2(f/fcutInsp)^(8/3) + rho3(f/fcutInsp)^(9/3).

  The coefficients are computed demanding that the pseudo-PN part at the three collocation points frequencies returns the actual collocation points: iv1, iv2, iv3.

*/

/* The input arguments for these rho coefficients are the value of the collocation points: v1, v2, v3, useful powers of the frequencies of these points f1, f2, f3 and the cutting frequency for the
inspiral fcutInsp. The values v1, v2, v3 are given by the fit of the collocation points minus the value of the PNAnsatz at that frequency. */

/* There are 3 functions, one for each pseudoPN coefficient rho1,2,3.
Every function contains three possible results, depending on if you reconstruct with 3, 2 or 1 collocation points.
With 2 colloc points rho3 = 0, and with 1 colloc point rho2=rho3=0 */

/* Get Pseudo PN Amp coefficient rho1 at f^(7/3) */
static double IMRPhenomXHM_Inspiral_Amp_rho1(double v1, double v2, double v3, IMRPhenomX_UsefulPowers *powers_of_fcutInsp, IMRPhenomX_UsefulPowers *powers_of_f1, IMRPhenomX_UsefulPowers *powers_of_f2, IMRPhenomX_UsefulPowers *powers_of_f3, IMRPhenomXHMWaveformStruct *pWFHM){
  double retVal = 0.;
  switch(pWFHM->IMRPhenomXHMInspiralAmpVersion){
    case 3: // 3 inspiral collocation points
    {
      // default version
      retVal = (powers_of_fcutInsp->seven_thirds*(-(powers_of_f1->three*powers_of_f3->eight_thirds*v2) + powers_of_f1->eight_thirds*powers_of_f3->three*v2 + powers_of_f2->three*(powers_of_f3->eight_thirds*v1 - powers_of_f1->eight_thirds*v3) + powers_of_f2->eight_thirds*(-(powers_of_f3->three*v1) + powers_of_f1->three*v3)))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
      break;
    }
    case 2:  // 2 inspiral collocation points
    {
      retVal=(powers_of_fcutInsp->seven_thirds*(-(powers_of_f2->eight_thirds*v1) + powers_of_f1->eight_thirds*v2))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds);
      break;
    }
    case 1:  // 1 inspiral collocation points
    {
      retVal=(powers_of_fcutInsp->seven_thirds*v1)/(powers_of_f1->seven_thirds);
      break;
    }
    case 0:  // No collocation points, reconstruct with PN only
    {
      retVal=0;
      break;
    }
    default: {
      XLALPrintError("Error in IMRPhenomXHM_Inspiral_Amp_rho1: version is not valid. Versions available are 0,1,2,3.\n");
      XLAL_ERROR(XLAL_EDOM, "IMRPhenomXHM_Inspiral_Amp_rho1 version is not valid.  Aborting. Versions available are 0,1,2,3.\n");
    }
  }
  return retVal;

}

/* Get Pseudo PN Amp coefficient rho2 at f^(8/3) */
static double IMRPhenomXHM_Inspiral_Amp_rho2(double v1, double v2, double v3, IMRPhenomX_UsefulPowers *powers_of_fcutInsp, IMRPhenomX_UsefulPowers *powers_of_f1, IMRPhenomX_UsefulPowers *powers_of_f2, IMRPhenomX_UsefulPowers *powers_of_f3, IMRPhenomXHMWaveformStruct *pWFHM){

  double retVal = 0;
  switch(pWFHM->IMRPhenomXHMInspiralAmpVersion){
    case 3:  // 3 inspiral collocation points
    {
      // default version
      retVal=(-(powers_of_fcutInsp->eight_thirds*(-(powers_of_f1->three*powers_of_f3->seven_thirds*v2) + powers_of_f1->seven_thirds*powers_of_f3->three*v2 + powers_of_f2->three*(powers_of_f3->seven_thirds*v1 - powers_of_f1->seven_thirds*v3) + powers_of_f2->seven_thirds*(-(powers_of_f3->three*v1) + powers_of_f1->three*v3))))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
      break;
    }
    case 2: // 2 inspiral collocation points
    {

      retVal=(-(powers_of_fcutInsp->eight_thirds*(-(powers_of_f2->seven_thirds*v1) + powers_of_f1->seven_thirds*v2)))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds);
      break;
    }
    case 1:  // 1 inspiral collocation points
    {
      retVal = 0;
      break;
    }
    case 0:  // No collocation points, reconstruct with PN only
    {
      retVal=0;
      break;
    }
    default: {XLALPrintError("Error in IMRPhenomXHM_Inspiral_Amp_rho2: version is not valid. Versions avilable are 0,1,2,3.\n");XLAL_ERROR(XLAL_EDOM, "IMRPhenomXHM_Inspiral_Amp_rho2 version is not valid. Aborting. Versions available are 0,1,2,3.\n");}
  }
  return retVal;

}

/* Get Pseudo PN Amp coefficient rho3 at f^(9/3) */
static double IMRPhenomXHM_Inspiral_Amp_rho3(double v1, double v2, double v3, IMRPhenomX_UsefulPowers *powers_of_fcutInsp, IMRPhenomX_UsefulPowers *powers_of_f1, IMRPhenomX_UsefulPowers *powers_of_f2, IMRPhenomX_UsefulPowers *powers_of_f3, IMRPhenomXHMWaveformStruct *pWFHM){

  double retVal = 0;
  switch(pWFHM->IMRPhenomXHMInspiralAmpVersion){
    case 3:   // 3 inspiral collocation points
    {
      // default version
      retVal=(powers_of_fcutInsp->three*(powers_of_f1->seven_thirds*(-powers_of_f1->one_third + powers_of_f3->one_third)*powers_of_f3->seven_thirds*v2 + powers_of_f2->seven_thirds*(-(powers_of_f3->eight_thirds*v1) + powers_of_f1->eight_thirds*v3) + powers_of_f2->eight_thirds*(powers_of_f3->seven_thirds*v1 - powers_of_f1->seven_thirds*v3)))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
      break;
    }
    case 2:   // 2 inspiral collocation points
    {
      retVal = 0;
      break;
    }
    case 1:   // 1 inspiral collocation point
    {
      retVal = 0;
      break;
    }
    case 0:  // No collocation points, reconstruct with PN only
    {
      retVal=0;
      break;
    }
      default: {XLALPrintError("Error in IMRPhenomXHM_Inspiral_Amp_rho3: version is not valid.\n");XLAL_ERROR(XLAL_EDOM, "IMRPhenomXHM_Inspiral_Amp_rho3 version is not valid.  Aborting. Versions available: 0,1,2,3.\n");}
  }
  return retVal;
}


/************* INSPIRAL AMPLITUDE ANSATZ ******************/

/*
   The ansatz is built by performing the Stationary Phase Approximation of the Time-Domain Amplitude up to 3PN.
   The result is rexpanded in power series up to 3PN, except for the 21 mode, which is better behaved without the reexpansion.
*/

//Return the Fourier Domain Post-Newtonian ansatz up to 3PN without the pseudoPN terms for a particular frequency
static double IMRPhenomXHM_Inspiral_PNAmp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){

  // The 21 mode is special, is not a power series
  if(pWFHM->useFAmpPN==1){
    return IMRPhenomXHM_Inspiral_PNAmp_21Ansatz(powers_of_Mf, pWFHM, pAmp);
  }

   //This returns the amplitude strain rescaled with the prefactor of the 22 mode: divided by sqrt(2*eta/3.)/pi^(1/6)
  double complex CpnAmp;
  double pnAmp;
  /*int InsAmpFlag = pWFHM->IMRPhenomXHMInspiralAmpFitsVersion;
  switch(InsAmpFlag)
  {
    case 122018:
    case 102021://FIXME
    {*/
      CpnAmp = (pAmp->pnInitial
        + powers_of_Mf->one_third    * pAmp->pnOneThird
        + powers_of_Mf->two_thirds   * pAmp->pnTwoThirds
        + powers_of_Mf->itself       * pAmp->pnThreeThirds
        + powers_of_Mf->four_thirds  * pAmp->pnFourThirds
        + powers_of_Mf->five_thirds  * pAmp->pnFiveThirds
        + powers_of_Mf->two          * pAmp->pnSixThirds
      );
      pnAmp = (pAmp->PNglobalfactor)*cabs(CpnAmp);
    /*  break;
    }
    default :
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Inspiral_PNAmp_Ansatz: IMRPhenomXInspiralAmpVersion is not valid. Recommended version is 122018.\n");
    }
}*/

  return pnAmp;
}

//Return the 21 mode Fourier Domain Post-Newtonian ansatz up to 3PN without the pseudoPN terms for a particular frequency
static double IMRPhenomXHM_Inspiral_PNAmp_21Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){

  double complex CpnAmpTD; //This is not the real strain of the lm mode. It is the strain rescaled with the prefactor of the 22 mode: divided by sqrt(2*eta/3.)/pi^(1/6)
  double pnAmp, XdotT4, x_to_m_one_four, two_to_m_one_sixths = 0.8908987181403393, three_to_m_one_second = 0.5773502691896257;
  x_to_m_one_four = two_to_m_one_sixths * powers_of_lalpiHM.m_one_sixth * powers_of_Mf->m_one_sixth;
  /*int InsAmpFlag = pWFHM->IMRPhenomXHMInspiralAmpFitsVersion;
  switch(InsAmpFlag)
  {
    case 122018:
    case 102021://FIXME
    {*/
      // Complex time-domain Post-Newtonina amplitude, power series
      CpnAmpTD = (
         powers_of_Mf->one_third    * pAmp->x05
        + powers_of_Mf->two_thirds  * pAmp->x1
        + powers_of_Mf->itself      * pAmp->x15
        + powers_of_Mf->four_thirds * pAmp->x2
        + powers_of_Mf->five_thirds * pAmp->x25
        + powers_of_Mf->two         * pAmp->x3
      );
      CpnAmpTD = CpnAmpTD*powers_of_Mf->two_thirds*pAmp->PNTDfactor;

      // Value of the the derivative of the PN expansion parameter X given by TaylorT4
      XdotT4 = powers_of_Mf->five_thirds * powers_of_Mf->five_thirds * pAmp->xdot5
      + powers_of_Mf->four          * pAmp->xdot6
      + powers_of_Mf->eight_thirds*powers_of_Mf->five_thirds * pAmp->xdot65
      + powers_of_Mf->seven_thirds*powers_of_Mf->seven_thirds * pAmp->xdot7
      + powers_of_Mf->eight_thirds*powers_of_Mf->seven_thirds * pAmp->xdot75
      + powers_of_Mf->eight_thirds * powers_of_Mf->eight_thirds * pAmp->xdot8
      + (powers_of_Mf->log*2./3. + pAmp->log2pi_two_thirds )*powers_of_Mf->eight_thirds * powers_of_Mf->eight_thirds * pAmp->xdot8Log
      + powers_of_Mf->eight_thirds * powers_of_Mf->eight_thirds * powers_of_Mf->one_third  * pAmp->xdot85;

      // Perform the SPA, multiply time-domain by the phasing factor
      pnAmp = 2. * powers_of_lalpiHM.sqrt  * three_to_m_one_second  * cabs(CpnAmpTD) * x_to_m_one_four / sqrt(XdotT4) /powers_of_Mf->m_seven_sixths/pWFHM->ampNorm;
     /* break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Inspiral_PNAmp_Ansatz: IMRPhenomXInspiralAmpVersion is not valid. Recommended version is 122018.\n");}
}*/
  return pnAmp;
}


//This is the complete Inspiral Amplitude Ansatz: PN ansatz + pseudoPN terms
static double IMRPhenomXHM_Inspiral_Amp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp)
{
    double InspAmp; //This is the amplitude strain rescaled with the prefactor of the 22 mode: divided by [sqrt(2*eta/3.)/pi^(1/6) * f^(-7/6)]
    /*int InsAmpFlag = pWFHM->IMRPhenomXHMInspiralAmpFitsVersion;
    switch(InsAmpFlag)
    {
        case 122018:
        case 102021://FIXME
        {*/
            InspAmp = IMRPhenomXHM_Inspiral_PNAmp_Ansatz(powers_of_Mf, pWFHM, pAmp)
            + powers_of_Mf->seven_thirds / pAmp->fcutInsp_seven_thirds * pAmp->rho1
            + powers_of_Mf->eight_thirds / pAmp->fcutInsp_eight_thirds * pAmp->rho2
            + powers_of_Mf->three        / pAmp->fcutInsp_three * pAmp->rho3
            ;/*
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Inspiral_Amp_Ansatz: IMRPhenomXInspiralAmpVersion is not valid. Recommended version is 2018. \n");}
    }*/
    if(rescalefactor>0){
        InspAmp*= RescaleFactor(powers_of_Mf, pAmp, rescalefactor);
    }

    return InspAmp;
}

/* Numerical derivative of the inspiral (4th order finite differences)
   It is used for reconstructing the intermediate region */
static double IMRPhenomXHM_Inspiral_Amp_NDAnsatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){

    double df = 10e-10;
    double Nder;
    double fun2R, funR, funL, fun2L;
    double centralfreq = powers_of_Mf->itself;

    IMRPhenomX_UsefulPowers powers_of_Mf2R, powers_of_MfR, powers_of_MfL, powers_of_Mf2L;

    IMRPhenomX_Initialize_Powers(&powers_of_Mf2R, centralfreq + 2*df);
    IMRPhenomX_Initialize_Powers(&powers_of_MfR, centralfreq + df);
    IMRPhenomX_Initialize_Powers(&powers_of_MfL, centralfreq  - df);
    IMRPhenomX_Initialize_Powers(&powers_of_Mf2L, centralfreq - 2*df);

    fun2R = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_Mf2R, pWFHM, pAmp);
    funR  = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_MfR, pWFHM, pAmp);
    funL  = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_MfL, pWFHM, pAmp);
    fun2L = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_Mf2L, pWFHM, pAmp);
    Nder = (-fun2R + 8*funR - 8*funL + fun2L )/(12*df);

    return Nder;
}

/* VETO function:  when extrapolating the model outside the calibration region the collocation points can be bad behaved.
   With this function we select those that will be used in the reconstruction. */
void IMRPhenomXHM_Inspiral_Amplitude_Veto(
  double *iv1, double *iv2, double *iv3,
  IMRPhenomX_UsefulPowers *powers_of_f1,
  IMRPhenomX_UsefulPowers *powers_of_f2,
  IMRPhenomX_UsefulPowers *powers_of_f3,
  IMRPhenomXHMAmpCoefficients *pAmp,
  IMRPhenomXHMWaveformStruct *pWFHM
)
{
    double threshold = 0.2/(pWFHM->ampNorm);
    #if DEBUG == 1
    printf("\n\nf1, f2, f3 = %.16f %.16f %.16f\n\n",powers_of_f1->itself,powers_of_f2->itself,powers_of_f3->itself);
    printf("\n\nf1, f2, f3 = %.16f %.16f %.16f\n\n",threshold*powers_of_f1->seven_sixths,threshold*powers_of_f2->seven_sixths,threshold*powers_of_f3->seven_sixths);
    printf("\nInspiral Veto: AmpVersion = %i",pWFHM->IMRPhenomXHMInspiralAmpVersion);
    #endif
    // Remove too low collocation points (heuristic).
    if(pAmp->CollocationPointsValuesAmplitudeInsp[0] < threshold*powers_of_f1->seven_sixths){
        *iv1 = 0;
        pWFHM->IMRPhenomXHMInspiralAmpVersion = 2;
    }
    if(pAmp->CollocationPointsValuesAmplitudeInsp[1] < threshold*powers_of_f2->seven_sixths){
        *iv2 = 0;
        pWFHM->IMRPhenomXHMInspiralAmpVersion = pWFHM->IMRPhenomXHMInspiralAmpVersion -1;
    }
    if(pAmp->CollocationPointsValuesAmplitudeInsp[2] < threshold*powers_of_f3->seven_sixths){
        *iv3 = 0;
        pWFHM->IMRPhenomXHMInspiralAmpVersion = pWFHM->IMRPhenomXHMInspiralAmpVersion -1;
    }
}

// Check if the three collocation points are wavy
int WavyPoints(double p1, double p2, double p3){
    if((p1>p2 && p2<p3) || (p1<p2 && p2>p3)){
        return 1;
    }else{
        return 0;
    }
}


double RescaleFactor(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMAmpCoefficients *pAmp, INT4 rescalefactor){
    double factor = 0.;
    switch(rescalefactor){
        case 1:{
            factor = powers_of_Mf->m_seven_sixths;
        }
        case 2:{
            factor = pAmp->PNdominant * powers_of_Mf->m_seven_sixths; // = Pi * Sqrt(2 eta/3) (2Pi Mf / m)^(-7/6)
        }
    }
    return factor;
}


/*************************************/
/*                                   */
/*              PHASE                */
/*                                   */
/*************************************/

// Spin parameter S = (m1^2*chi1 + m2^2*chi2)/(m1^2 + m2^2)

// Below we give paramater-space fits for the weighted difference between each mode's phase and the 22-phase: phi_lm-m/2 phi_22(2/m f), see Eqs. (4.10-4.12)

static double IMRPhenomXHM_Insp_Phase_21_lambda(double eta, double chi1, double chi2, int InspPhaseFlag) {
    double total;
    switch (InspPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = 13.664473636545068 - 170.08866400251395*eta + 3535.657736681598*eta2 - 26847.690494515424*eta3 + 96463.68163125668*eta4 - 133820.89317471132*eta5;
            double eqSpin = (S*(18.52571430563905 - 41.55066592130464*S + eta3*(83493.24265292779 + 16501.749243703132*S - 149700.4915210766*S2) + eta*(3642.5891077598003 + 1198.4163078715173*S - 6961.484805326852*S2) + 33.8697137964237*S2 + eta2*(-35031.361998480075 - 7233.191207000735*S + 62149.00902591944*S2)))/(6.880288191574696 + 1.*S);
            double uneqSpin = -134.27742343186577*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_Insp_Phase_21_lambda: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Phase_33_lambda(double eta, double chi1, double chi2, int InspPhaseFlag) {
    double total=0;
    double delta=sqrt(1. - 4.*eta);
    switch (InspPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3;
            eta2 = pow(eta,2);
            eta3 = eta2*eta;
            double noSpin = 4.1138398568400705 + 9.772510519809892*eta - 103.92956504520747*eta2 + 242.3428625556764*eta3;
            double eqSpin = ((-0.13253553909611435 + 26.644159828590055*eta - 105.09339163109497*eta2)*S)/(1. + 0.11322426762297967*S);
            double uneqSpin = -19.705359163581168*(chi1 - 1.*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_Insp_Phase_32_lambda: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Phase_32_lambda(double eta, double chi1, double chi2, int InspPhaseFlag) {
    double total;
    switch (InspPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3,eta4,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = (9.913819875501506 + 18.424900617803107*eta - 574.8672384388947*eta2 + 2671.7813055097877*eta3 - 6244.001932443913*eta4)/(1. - 0.9103118343073325*eta);
            double eqSpin = (-4.367632806613781 + 245.06757304950986*eta - 2233.9319708029775*eta2 + 5894.355429022858*eta3)*S + (-1.375112297530783 - 1876.760129419146*eta + 17608.172965575013*eta2 - 40928.07304790013*eta3)*S2 + (-1.28324755577382 - 138.36970336658558*eta + 708.1455154504333*eta2 - 273.23750933544176*eta3)*S3 + (1.8403161863444328 + 2009.7361967331492*eta - 18636.271414571278*eta2 + 42379.205045791656*eta3)*S4;
            double uneqSpin = (chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2*(-105.34550407768225 - 1566.1242344157668*chi1*eta + 1566.1242344157668*chi2*eta + 2155.472229664981*eta*S);
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_Insp_Phase_32_lambda: version is not valid. Recommended version is 122019.");}
    }
    return total;
}


static double IMRPhenomXHM_Insp_Phase_44_lambda(double eta, double chi1, double chi2, int InspPhaseFlag) {
    double total;
    switch (InspPhaseFlag){
        case 122019:{
            double S = XLALSimIMRPhenomXSTotR(eta,chi1,chi2);
            double eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = 5.254484747463392 - 21.277760168559862*eta + 160.43721442910618*eta2 - 1162.954360723399*eta3 + 1685.5912722190276*eta4 - 1538.6661348106031*eta5;
            double eqSpin = (0.007067861615983771 - 10.945895160727437*eta + 246.8787141453734*eta2 - 810.7773268493444*eta3)*S + (0.17447830920234977 + 4.530539154777984*eta - 176.4987316167203*eta2 + 621.6920322846844*eta3)*S2;
            double uneqSpin = -8.384066369867833*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_Insp_Phase_44_lambda: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

// Returns linear-in-f term coming from the complex part of each mode's amplitude that will be added to the orbital phase, this is an expansion of Eq. (4.9) truncated at linear order
static double IMRPhenomXHM_Insp_Phase_LambdaPN(double eta, int modeInt){

    double output;
    const double PI = powers_of_lalpiHM.itself;

    switch(modeInt){
        case 21:
        {//2 f \[Pi] (-(1/2) - Log[16]/2)
            output=(2.*PI*(-0.5-2.*log(2)));
            break;
        }
        case 33:
        {//2/3 f \[Pi] (-(21/5) + 6 Log[3/2])
            output=2./3.*PI*(-21./5.+6.*log(1.5));
            break;
        }
        case 32:
        { //-((2376 f \[Pi] (-5 + 22 \[Eta]))/(-3960 + 11880 \[Eta]))
            output=-((2376.*PI*(-5.+22.*eta))/(-3960.+11880*eta));
            break;
        }
        case 44:
        { //(45045 f \[powers_of_lalpiHM.itself] (336 - 1193 \[Eta] +320 (-1 + 3 \[Eta]) Log[2]))/(2 (-1801800 + 5405400 \[Eta]))
            output=45045.*PI*(336.-1193.*eta+320.*(-1.+3.*eta)*log(2))/(2.*(-1801800. + 5405400.*eta));
            break;
        }
        default: output=0.;
    }

    return -1.*output;
}


/************* INSPIRAL PHASE ANSATZ *************/

/*
 The function loads the rescaled phenX coefficients for the phase and correct the result with a linear-in-f term coming from the complex part of the PN amplitude

 The rescaling of the 22-coefficients is explained in App. D
*/
static double IMRPhenomXHM_Inspiral_Phase_AnsatzInt(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMPhaseCoefficients *pPhase)
{
    //compute the orbital phase by laoding the rescaled phenX coefficients
    double philm=0.;
    double freqs[]={powers_of_Mf->m_five_thirds,powers_of_Mf->m_four_thirds,powers_of_Mf->m_one,powers_of_Mf->m_two_thirds,powers_of_Mf->m_one_third,1.,powers_of_Mf->one_third,powers_of_Mf->two_thirds, Mf, powers_of_Mf->four_thirds, powers_of_Mf->five_thirds, powers_of_Mf->two,powers_of_Mf->seven_thirds};
    double logMf=powers_of_Mf->log;
    for(int i=0; i<N_MAX_COEFFICIENTS_PHASE_INS; i++)
        philm+=(pPhase->phi[i]+pPhase->phiL[i]*(logMf))*freqs[i];
    return philm;
}

static double IMRPhenomXHM_Inspiral_Phase_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf,IMRPhenomXHMPhaseCoefficients *pPhase)
{
    //compute the orbital phase by laoding the rescaled phenX coefficients
    double dphilm=0.;
    double coeffs[]={-5./3,-4./3,-1.,-2./3,-1./3,0.,1./3, 2./3, 1., 4./3, 5./3, 2., 7./3};
    double freqs[]={powers_of_Mf->m_eight_thirds,powers_of_Mf->m_seven_thirds,powers_of_Mf->m_two,powers_of_Mf->m_five_thirds,powers_of_Mf->m_four_thirds,powers_of_Mf->m_one,powers_of_Mf->m_two_thirds,powers_of_Mf->m_one_third,1.,powers_of_Mf->one_third, powers_of_Mf->two_thirds, Mf,powers_of_Mf->four_thirds};
    double logMf=powers_of_Mf->log;

    for(int i=0; i<N_MAX_COEFFICIENTS_PHASE_INS; i++)
        dphilm+=((pPhase->phi[i]+pPhase->phiL[i]*(logMf))*coeffs[i]+pPhase->phiL[i])*freqs[i];
    return dphilm;
}
