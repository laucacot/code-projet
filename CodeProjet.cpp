#include <iostream>
#include <string>
#include <cmath>
//#include <boost/serialization/array_wrapper.hpp>

#include <vector>
// Optimization
#define BOOST_UBLAS_NDEBUG
#include <fstream>
#include <exception>
#include <string>
#include <utility>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <time.h>
#include <stdlib.h>

using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::math::tools;

// type definitions
typedef double value_type;// or typedef float value_type;
typedef boost::numeric::ublas::vector< value_type > state_type;
typedef boost::numeric::ublas::matrix< value_type > matrix_type;
typedef rosenbrock4< value_type > stepper_type;

// constants
const value_type pression = 0.1/760; //atm soit 0.1 torr et 13 pascal 
const value_type L = 3.e-2; //distance netre deux plaques en m
const value_type pi = M_PI;
const value_type diff = pow((pi/L), 2.);//facteur pour la diffusion 
const value_type n_Ar =  pression*2.69e25; //densite d'argon en m-3
const value_type n_SiH4_ini = n_Ar/80.; //densite de SiH4 initiale
const int Nbr_espece=21;
const value_type DP = 1.e23;//eV/s.m3 puissance totale du systeme par unite de volume imposee
const float C=1.e20;// m-3/s taux d'injection du SiH4 dans le réacteur
const value_type Tg =0.02758 ; //eV soit 320 K
const float D=1.;// m-3/s taux d'injection du SiH4 dans le réacteur

//calcul des diffusions
state_type DL(Nbr_espece, 0.0); //vecteur de diffusion libre en m2/s
state_type mu(Nbr_espece, 0.0); //vecteur de mobilite en m2/(V.s)
state_type DA(Nbr_espece, 0.0); //vecteur de diffusion ambipolaire en m2/s


//calcul des K dependant de Te

value_type k1 (value_type Te) //K1 Ar + e -> Ar+ + 2e
{
  value_type K1;
  K1=7.06E-17*pow((Te),0.6)*exp(-(16.14)/(Te));
  return K1;
}

value_type k2 (value_type Te) //K2 Ar + e -> Ar* + e
{
  value_type K2;
  K2=11.69E-15*exp(-(12.31)/(Te));
  return K2;
}

value_type k3 (value_type Te) //K3 Ar* + e -> Ar+ + 2e
{
  value_type K3;
  K3=124.92E-15*exp(-(5.39)/(Te));
  return K3;
}

value_type k4 (value_type Tg) //K4 Ar* + Ar* -> Ar + Ar+ + e
{
  value_type K4;
  K4=6.144e-16;
  return K4;
}

value_type k5 (value_type Te) //K5 Ar* + e -> Ar + e
{
  value_type K5;
  K5=431.89E-18*pow((Te),0.74);
  return K5;
}

value_type k6 (value_type Te) //K6 SiH4 + e -> SiH3 + H + e
{
    value_type K6;
    K6=1.83E-9*pow((Te),-1.)*exp(-(10.68)/(Te))/D;
    return K6;

}

value_type k7 (value_type Te) //K7 SiH4 + e -> SiH2 + 2H + e
{
    value_type K7;
    K7=8.97E-9*pow((Te),-1.)*exp(-(10.68)/(Te));
    return K7;

}

value_type k8 (value_type Te) //K8 SiH4 + e -> SiH3- + H
{
    value_type K8;
    K8=3.77E-9*pow((Te),-1.63)*exp(-(8.29)/(Te));
    return K8;

}

value_type k9 (value_type Te) //K9 SiH4 + e -> SiH2- + 2H
{
    value_type K9;
    K9=3.77E-9*pow((Te),-1.63)*exp(-(8.29)/(Te));
    return K9;

}

value_type k10 (value_type Te) //K10 SiH4 + e -> SiH3+ + H + 2e
{

    value_type K10;
    K10= 2.5e2*pow((Te),-2.93)*exp(-(24.1)/(Te))/D;
    return K10;

}

value_type k11 (value_type Te) //K11 SiH3 + e -> SiH2- + H
{

    value_type K11;
    K11= 5.71E-9*pow((Te),-0.5)*exp(-(1.94)/(Te));
    return K11;

}

value_type k12 (value_type Te) //K12 SiH3 + e -> SiH3+  + 2e
{
    value_type K12;
    K12= 2.26E-16*pow((Te),0.5)*exp(-(1.30)/(Te))/D;
    return K12;

}

value_type k13 (value_type Te) //K13 SiH3- + e -> SiH3 + 2e
{
    value_type K13;
    K13=3.15E-16*pow((Te),0.5)*exp(-(1.16)/(Te))/D;
    return K13;

}

value_type k14 (value_type Te) //K14 SiH2- + e -> SiH2  + 2e
{

    value_type K14;
    K14= 3.15E-16*pow((Te),0.5)*exp(-(1.16)/(Te));
    return K14;

}

value_type k15 (value_type Te) //K15 SiH2 + e -> SiH2-
{
    value_type K15;
    K15=5.71E-16*pow((Te),-0.5);
    return K15;

}

value_type k16 (value_type Te) //K16 H2 + e ->  2H + e
{
    value_type K16;
    K16=4.73E-14*pow((Te),-0.23)*exp(-(10.09)/(Te));
    return K16;

}

value_type k17 (value_type Te) //K17 H2 + e ->  H2+ + 2e
{
    value_type K17;
    K17=1.1E-14*pow((Te),0.42)*exp(-(16.05)/(Te));
    return K17;

}

value_type k18(value_type Tg) //K18 SiH4 + Ar* -> SiH3 + H + Ar
{
    value_type K18;
    K18 = 1.400e-16/D;
    return K18;
}

value_type k19(value_type Tg)
{
    value_type K19;
    K19=2.591e-16;
    return K19;
}

value_type k20(value_type Tg)
{
    value_type K20;
    K20= 99.67e-18;
    return K20;
}

value_type k21(value_type Tg)
{
    value_type K21;
    K21= 9.963e-17;
    return K21;
}

value_type k22(value_type Tg)
{
    value_type K22;
    K22= 6.974e-17;
    return K22;
}

value_type k23 (value_type Tg) //k23(Tg)%K23 SiH3 + SiH3 -> SiH2 + SiH4
{
    value_type K23;
    K23= 2.99e-17;
    return K23;
}

value_type k24 (value_type Tg) //K24 SiH4 + SIH3 -> Si2H5 + H2
{

    value_type K24;
    K24=2.94e-18*exp(-0.1908/Tg);
    return K24;
}

value_type k25 (value_type Tg) //K25 SiH2 + H2 -> SiH4
{
    value_type K25;
    K25=2.e-19 ;
    return K25;
}

value_type k26 (value_type Tg) //K26 SiH2  -> Si + H2
{
    value_type K26;
    K26=1.51E-9*pow((Tg),1.658)*exp(-(1.66)/(Tg));
    return K26;
}

value_type k27 (value_type Tg) //K27 SiH4 + H -> H2 + SiH3
{
    value_type K27;
    K27=2.44E-22*pow((Tg),1.9)*exp(-(0.09)/(Tg))/D;
    return K27;
}

value_type k28 (value_type Tg) //K28 SiH2 + SiH2 -> Si2H2 + H2
{
    value_type K28;
    K28=1.08e-15 ;
    return K28;
}

value_type k29 (value_type Tg) //K29 SiH2 + H -> SiH + H2
{
    value_type K29;
    K29=2.31e-17;
    return K29;
}

value_type k30 (value_type Tg) //K30 SiH3-> SiH + H2
{
    value_type K30;
    K30=328.9E-6*pow((Tg),-3.1)*exp(-(1.94)/(Tg));
    return K30;
}

value_type k31 (value_type Tg) //K31 SiH3 + H -> SiH2 + H2
{
    value_type K31;
    K31=2.49E-17*exp(-(0.1084)/(Tg));
    return K31;
}

value_type k32 (value_type Tg) //K32 SiH2- + H2+ -> SiH2 + H2
{
    value_type K32;
    K32= 5.55E-12*pow((Tg),-0.5);
    return K32;
}

value_type k33 (value_type Tg) //K33 SiH3- + H2+ -> SiH3 + H2
{
    value_type K33;
    K33=5.55E-12*pow((Tg),-0.5)/D;
    return K33;
}

value_type k34 (value_type Tg) //K34 SiH3- +H2+ -> SiH3 + H2
{
    value_type K34;
    K34=2.11e-20;
    return K34;
}

value_type k35 (value_type Tg) //K35 SiH3- + SiH3+ -> Si2H6
{
    value_type K35;
    K35= 2.11e-20 ;
    return K35;
}

value_type k36 (value_type Tg) //K36 SiH2- + SiH4 -> Si2H4- + H2
{
    value_type K36;
    K36= 2.11e-20;
    return K36;
}

value_type k37 (value_type Tg) //K37 SiH2- + SiH3 -> SiH2 + SiH3-
{
    value_type K37;
    K37=2.11e-20 ;
    return K37;
}

value_type k38 (value_type Tg) //K38 SiH3- + SiH2 -> Si2H3- + H2
{
    value_type K38;
    K38=2.11e-20 ;
    return K38;
}

value_type k39 (value_type Tg) //K39 SiH3- + SiH4 -> Si2H5- + H2
{
    value_type K39;
    K39=2.11e-20 ;
    return K39;
}

value_type k40 (value_type Tg) //K40 SiH2- + SiH3+ -> Si2H5
{
    value_type K40;
    K40=2.11e-20 ;
    return K40;
}

value_type k41 (value_type Tg) //K41 SiH2- + Ar+ -> SiH2 + Ar
{
    value_type K41;
    K41=1.44e-12*pow((Tg),-0.5);
    return K41;
}

value_type k42 (value_type Tg) //K42 SiH3- + Ar+ -> SiH3 + Ar
{
    value_type K42;
    K42=1.44E-12*pow((Tg),-0.5)/D;
    return K42;
}

value_type k43 (value_type Te) //K43 SiH3 + e ->  SiH3-
{
    value_type K43;
    K43=5.71E-16*pow((Te),-0.5);
    return K43;

}

value_type k44 (value_type Te) //K44 SiH + e ->  SiH-
{
    value_type K44;
    K44=5.71E-15*pow((Te),-0.5);
    return K44;

}

value_type k45 (value_type Te) //K45 SiH- + e ->  SiH + 2e
{
    value_type K45;
    K45=3.16E-16*pow((Te),0.5)*exp(-(1.25)/(Te));
    return K45;

}

value_type k46 (value_type Tg) //K46 SiH2- + SiH ->  SiH2 + SiH2m
{
    value_type K46;
    K46=2.31e-17;
    return K46;
}

value_type k47 (value_type Tg) //K47 SiH- + H2p ->  SiH + H2
{
    value_type K47;
    K47= 3.21e-13;
    return K47;
}


/*value_type k58 (value_type Tg) //K58 SiH3- + SiH3+ ->  Si2H4 + H2
{
    value_type K58;
    K58=  4.87E-13*0.5*pow((Tg),-0.5);
    return K58;
} pas de Si2H4*/

struct Condition
{
  value_type tol=1.e-9;
  bool operator() (value_type min, value_type max)  {
    return abs(min - max) <= tol;
  }
};

struct nsystem
{
  void operator()(const state_type &n, state_type &dndt, const value_type &t)
  {



    /*0=e, 1=Armet, 2=SiH3-, 3=SiH2-, 4=SiH3+, 5=SiH4, 6=SiH3,
    7=H, 8=SiH2, 9=H2, 10=H2+, 11=Si2H5, 12=Si2H2, 13=Si2H4-,
    14=Si2H6, 15=Si2H3-, 16=Si2H5-, 17=SiH-, 18=SiH, 19=Si, 20=Arp*/

 value_type   dSi= n[2] + n[3] + n[4] + n[5] + n[6] + n[8] + 2.*n[11] 
      + 2.*n[12] + 2.*n[13] + 2.*n[14] + 2.*n[15] + 2.*n[16] + n[17] 
      + n[18] + n[19]; //somme de tous les atomes de si 

    dndt[0]=k1(Te)*n_Ar*n[0] +k3(Te)*n[1]*n[0] +k4(Tg)*pow(n[1],2.) -k8(Te)*n[5]*n[0]
        -k9(Te)*n[5]*n[0] +k10(Te)*n[5]*n[0] -k11(Te)*n[6]*n[0] +k12(Te)*n[0]*n[6]
        +k13(Te)*n[2]*n[0] +k14(Te)*n[3]*n[0] -k15(Te)*n[0]*n[8] +k17(Te)*n[9]*n[0]
        -k43(Te)*n[6]*n[0]-k44(Te)*n[18]*n[0] +k45(Te)*n[0]*n[17]-DA[0]*n[0];

    dndt[1]= k2(Te)*n_Ar*n[0] -k3(Te)*n[1]*n[0] -2*k4(Tg)*pow(n[1],2.) -k5(Te)*n[1]*n[0]
        -k18(Tg)*n[1]*n[5] -k19(Tg)*n[1]*n[5] -k20(Tg)*n[6]*n[1] -k21(Tg)*n[8]*n[1]
        -k22(Tg)*n[9]*n[1] -DA[1]*n[1];

    dndt[2]= k8(Te)*n[5]*n[0] -k13(Te)*n[2]*n[0] -k33(Tg)*n[2]*n[10] -k34(Tg)*n[2]*n[6]
        -k35(Tg)*n[2]*n[4] + k37(Tg)*n[3]*n[6] -k38(Tg)*n[2]*n[8] -k39(Tg)*n[2]*n[5]
        -k42(Tg)*n[2]*n[20] + k43(Te)*n[6]*n[0] - C * n[2] /dSi ;

    dndt[3]= k9(Te)*n[5]*n[0] +k11(Te)*n[6]*n[0] -k14(Te)*n[3]*n[0] +k15(Te)*n[0]*n[8]
        -k32(Tg)*n[3]*n[10] -k36(Tg)*n[3]*n[5] -k37(Tg)*n[3]*n[6] -k40(Tg)*n[3]*n[4]
        -k41(Tg)*n[3]*n[20] - k46(Tg)*n[3]*n[18] - C * n[3] /dSi;

    dndt[4]=k10(Te)*n[5]*n[0] +k12(Te)*n[0]*n[6] -k35(Tg)*n[2]*n[4] -k40(Tg)*n[3]*n[4]
	 -DA[4]*n[4] - C * n[4] /dSi;

    dndt[5]= C-k6(Te)*n[5]*n[0] -k7(Te)*n[0]*n[5] -k8(Te)*n[5]*n[0] -k9(Te)*n[5]*n[0]
        -k10(Te)*n[5]*n[0] -k18(Tg)*n[5]*n[1] -k19(Tg)*n[1]*n[5] +k23(Tg)*pow(n[6],2)
        -k24(Tg)*n[5]*n[6] +k25(Tg)*n[8]*n[9] -k27(Tg)*n[5]*n[7] -k36(Tg)*n[3]*n[5]
        -k39(Tg)*n[2]*n[5] - C * n[5] /dSi;

    dndt[6]= k6(Te)*n[5]*n[0] -k11(Te)*n[6]*n[0] -k12(Te)*n[0]*n[6] +k13(Te)*n[2]*n[0]
        +k18(Tg)*n[5]*n[1] -k20(Tg)*n[6]*n[1] -2*k23(Tg)*pow(n[6],2) -k24(Tg)*n[5]*n[6]
        +k27(Tg)*n[5]*n[7] -k30(Tg)*n[6] -k31(Tg)*n[6]*n[7] +k33(Tg)*n[2]*n[10]
	-k34(Tg)*n[2]*n[6]
        -k37(Tg)*n[3]*n[6] +k42(Tg)*n[2]*n[20]-k43(Te)*n[6]*n[0] + DA[4]*n[4]
	 - C * n[6] /dSi ; //SiH3+ + e -> SiH3 sur paroi

    dndt[7]=  k6(Te)*n[5]*n[0] +2*k7(Te)*n[0]*n[5] +k8(Te)*n[5]*n[0] +2*k9(Te)*n[5]*n[0]
        +k10(Te)*n[5]*n[0] +k11(Te)*n[6]*n[0] +2*k16(Te)*n[9]*n[0] +k18(Tg)*n[5]*n[1]
        +2*k19(Tg)*n[1]*n[5] +k20(Tg)*n[6]*n[1] +k21(Tg)*n[8]*n[1] +2*k22(Tg)*n[9]*n[1]
        -k27(Tg)*n[5]*n[7] -k29(Tg)*n[8]*n[7] -k31(Tg)*n[6]*n[7]- C * n[7] /dSi;

    dndt[8]= k7(Te)*n[0]*n[5] +k14(Te)*n[3]*n[0] -k15(Te)*n[0]*n[8] +k19(Tg)*n[1]*n[5]
        +k20(Tg)*n[6]*n[1] -k21(Tg)*n[8]*n[1] +k23(Tg)*pow(n[6],2) -k25(Tg)*n[8]*n[9]
        -k26(Tg)*n[8] -2*k28(Tg)*pow(n[8],2.) -k29(Tg)*n[8]*n[7] +k31(Tg)*n[6]*n[7]
	+k32(Tg)*n[3]*n[10]
        +k37(Tg)*n[3]*n[6] -k38(Tg)*n[2]*n[8] +k41(Tg)*n[3]*n[20] +k46(Tg)*n[3]*n[18]
	- C * n[8] /dSi ;

    dndt[9]=-k16(Te)*n[9]*n[0] -k17(Te)*n[9]*n[0] -k22(Tg)*n[9]*n[1] +k24(Tg)*n[5]*n[6]
        -k25(Tg)*n[8]*n[9] +k26(Tg)*n[8] +k27(Tg)*n[5]*n[7] +k28(Tg)*pow(n[8],2.) +k29(Tg)
	*n[8]*n[7]+k30(Tg)*n[6] +k31(Tg)*n[6]*n[7] +k32(Tg)*n[3]*n[10] +k33(Tg)*n[2]*n[10]
	+k34(Tg)*n[2]*n[6]+ k36(Tg)*n[3]*n[5] +k38(Tg)*n[2]*n[8] +k39(Tg)*n[2]*n[5]
	+k47(Tg)*n[17]*n[10] + DA[10]*n[10] - C * n[9] /dSi; // H2+ + e -> H2 sur paroi

    dndt[10]= k17(Te)*n[9]*n[0] -k32(Tg)*n[3]*n[10] -k33(Tg)*n[2]*n[10] -k47(Tg)*n[17]*n[10]
	 -DA[10]*n[10] - C * n[10] /dSi;

    dndt[11]= k24(Tg)*n[5]*n[6] +k40(Tg)*n[3]*n[4]- C * n[11] /dSi;

    dndt[12]=k28(Tg)*pow(n[8],2)- C * n[12] /dSi;

    dndt[13]= k34(Tg)*n[2]*n[6] +k36(Tg)*n[3]*n[5]- C * n[13] /dSi;

    dndt[14]= k35(Tg)*n[2]*n[4] - C * n[14] /dSi;

    dndt[15]= k38(Tg)*n[2]*n[8]- C * n[15] /dSi;

    dndt[16]= k39(Tg)*n[2]*n[5] - C * n[16] /dSi;

    dndt[17]= k44(Te)*n[18]*n[0] -k45(Te)*n[0]*n[17] +k46(Tg)*n[3]*n[18] -k47(Tg)*n[17]*n[10]
	- C * n[17] /dSi;

    dndt[18]= k21(Tg)*n[8]*n[1] +k29(Tg)*n[8]*n[7] +k30(Tg)*n[6]
        -k44(Te)*n[18]*n[0] +k45(Te)*n[0]*n[17] -k46(Tg)*n[3]*n[18] +k47(Tg)*n[17]*n[10]
	- C * n[18] /dSi;

    dndt[19]=k26(Tg)*n[8]- C * n[19] /dSi;

    dndt[20]=k1(Te)*n_Ar*n[0] +k3(Te)*n[1]*n[0] +k4(Tg)*pow(n[1],2.)
        -k41(Tg)*n[3]*n[20]-k42(Tg)*n[2]*n[20] -DA[20]*n[20];

  }

  value_type Te;


};

struct jacobian
{
  void operator()(const state_type &n, matrix_type &jacobi,
                  const value_type &t, state_type &dfdt ) const
  {

value_type   dSi= n[2] + n[3] + n[4] + n[5] + n[6] + n[8] + 2.*n[11] 
      + 2.*n[12] + 2.*n[13] + 2.*n[14] + 2.*n[15] + 2.*n[16] + n[17] 
      + n[18] + n[19];


    jacobi( 0 , 0 )=k1(Te)*n_Ar +k3(Te)*n[1]  -k8(Te)*n[5]
        -k9(Te)*n[5] +k10(Te)*n[5] -k11(Te)*n[6] +k12(Te)*n[6]
        +k13(Te)*n[2] +k14(Te)*n[3] -k15(Te)*n[8] +k17(Te)*n[9]
        -k43(Te)*n[6]-k44(Te)*n[18] +k45(Te)*n[17] -DA[0];
    jacobi( 0 , 1 )=  +k3(Te)*n[0] +k4(Tg)*n[1] ;
    jacobi( 0 , 2 )= +k13(Te)*n[0] ;
    jacobi( 0 , 3 )= +k14(Te)*n[0] ;
    jacobi( 0 , 4 )= 0.0;
    jacobi( 0 , 5 )= -k8(Te)*n[0] -k9(Te)*n[0] +k10(Te)*n[0] ;
    jacobi( 0 , 6 )=-k11(Te)*n[0] +k12(Te)*n[0]-k43(Te)*n[0];
    jacobi( 0 , 7 )=0.0;
    jacobi( 0 , 8 )=-k15(Te)*n[0] ;
    jacobi( 0 , 9 )=+k17(Te)*n[0];
    jacobi( 0 , 10 )=0.0;
    jacobi( 0 , 11 )=0.0;
    jacobi( 0 , 12 )=0.0;
    jacobi( 0 , 13 )=0.0;
    jacobi( 0 , 14 )=0.0;
    jacobi( 0 , 15 )=0.0;
    jacobi( 0 , 16 )=0.0;
    jacobi( 0 , 17 )= +k45(Te)*n[0];
    jacobi( 0 , 18 )=-k44(Te)*n[0] ;
    jacobi( 0 , 19 )=0.0;
    jacobi( 0 , 20 )=0.0;
    jacobi( 1 , 0 ) = k2(Te)*n_Ar -k3(Te)*n[1]  -k5(Te)*n[1];
    jacobi( 1 , 1 ) =  -k3(Te)*n[0] -2*k4(Tg)*n[1] -k5(Te)*n[0]
        -k18(Tg)*n[5] -k19(Tg)*n[5] -k20(Tg)*n[6] -k21(Tg)*n[8]
        -k22(Tg)*n[9]-DA[1];
    jacobi( 1 , 2 ) = 0.0;
    jacobi( 1 , 3 ) = 0.0;
    jacobi( 1 , 4 ) = 0.0;
    jacobi( 1 , 5 ) =-k18(Tg)*n[1] -k19(Tg)*n[1] ;
    jacobi( 1 , 6 ) =  -k20(Tg)*n[1] ;
    jacobi( 1 , 7 ) = 0.0;
    jacobi( 1 , 8 ) =  -k21(Tg)*n[1];
    jacobi( 1 , 9 ) = -k22(Tg)*n[1];
    jacobi( 1 , 10 ) = 0.0;
    jacobi( 1 , 11 ) =0.0;
    jacobi( 1 , 12 ) = 0.0;
    jacobi( 1 , 13 ) = 0.0;
    jacobi( 1 , 14 ) =0.0;
    jacobi( 1 , 15 ) = 0.0;
    jacobi( 1 , 16 ) =0.0;
    jacobi( 1 , 17 ) = 0.0;
    jacobi( 1 , 18 ) = 0.0;
    jacobi( 1 , 19 ) =0.0;
    jacobi( 1 , 20 ) = 0.0;
    jacobi( 2 , 0 ) = k8(Te)*n[5] -k13(Te)*n[2]  + k43(Te)*n[6];
    jacobi( 2 , 1 ) = 0.0;
    jacobi( 2 , 2 ) =  -k13(Te)*n[0] -k33(Tg)*n[10] -k34(Tg)*n[6]
        -k35(Tg)*n[4]  -k38(Tg)*n[8] -k39(Tg)*n[5]
        -k42(Tg)*n[20] -C/dSi ;
    jacobi( 2 , 3 ) = + k37(Tg)*n[6] ;
    jacobi( 2 , 4 ) = -k35(Tg)*n[2] ;
    jacobi( 2 , 5 ) = k8(Te)*n[0] -k39(Tg)*n[2];
    jacobi( 2 , 6 ) =  -k34(Tg)*n[2]+ k37(Tg)*n[3] + k43(Te)*n[0];
    jacobi( 2 , 7 ) = 0.0;
    jacobi( 2 , 8 ) =  -k38(Tg)*n[2] ;
    jacobi( 2 , 9 ) = 0.0;
    jacobi( 2 , 10 ) =  -k33(Tg)*n[2];
    jacobi( 2 , 11 ) = 0.0;
    jacobi( 2 , 12 ) = 0.0;
    jacobi( 2 , 13 ) = 0.0;
    jacobi( 2 , 14 ) = 0.0;
    jacobi( 2 , 15 ) = 0.0;
    jacobi( 2 , 16 ) = 0.0;
    jacobi( 2 , 17 ) = 0.0;
    jacobi( 2 , 18 ) = 0.0;
    jacobi( 2 , 19 ) = 0.0;
    jacobi( 2 , 20 ) = -k42(Tg)*n[2];
    jacobi( 3 , 0 ) = k9(Te)*n[5] +k11(Te)*n[6] -k14(Te)*n[3] +k15(Te)*n[8];
    jacobi( 3 , 1 ) = 0.0;
    jacobi( 3 , 2 ) = 0.0;
    jacobi( 3 , 3 ) =   -k14(Te)*n[0] 
        -k32(Tg)*n[10] -k36(Tg)*n[5] -k37(Tg)*n[6] -k40(Tg)*n[4]
        -k41(Tg)*n[20] - k46(Tg)*n[18]-C/dSi;
    jacobi( 3 , 4 ) =  -k40(Tg)*n[3];
    jacobi( 3 , 5 ) = k9(Te)*n[0]  -k36(Tg)*n[3];
    jacobi( 3 , 6 ) = +k11(Te)*n[0]  -k37(Tg)*n[3] ;
    jacobi( 3 , 7 ) = 0.0;
    jacobi( 3 , 8 ) =  +k15(Te)*n[0];
    jacobi( 3 , 9 ) = 0.0;
    jacobi( 3 , 10 ) =-k32(Tg)*n[3];
    jacobi( 3 , 11 ) = 0.0;
    jacobi( 3 , 12 ) = 0.0;
    jacobi( 3 , 13 ) = 0.0;
    jacobi( 3 , 14 ) = 0.0;
    jacobi( 3 , 15 ) = 0.0;
    jacobi( 3 , 16 ) = 0.0;
    jacobi( 3 , 17 ) = 0.0;
    jacobi( 3 , 18 ) = - k46(Tg)*n[3];
    jacobi( 3 , 19 ) = 0.0;
    jacobi( 3 , 20 ) =-k41(Tg)*n[3];
    jacobi( 4 , 0 ) =k10(Te)*n[5] +k12(Te)*n[6] ;
    jacobi( 4 , 1 ) =0.0;
    jacobi( 4 , 2 ) = -k35(Tg)*n[4] ;
    jacobi( 4 , 3 ) = -k40(Tg)*n[4];
    jacobi( 4 , 4 ) =  -k35(Tg)*n[2] -k40(Tg)*n[3] -DA[4] -C/dSi;
    jacobi( 4 , 5 ) =k10(Te)*n[0] ;
    jacobi( 4 , 6 ) =+k12(Te)*n[0];
    jacobi( 4 , 7 ) =0.0;
    jacobi( 4 , 8 ) =0.0;
    jacobi( 4 , 9 ) =0.0;
    jacobi( 4 , 10 ) =0.0;
    jacobi( 4 , 11 ) =0.0;
    jacobi( 4 , 12 ) =0.0;
    jacobi( 4 , 13 ) =0.0;
    jacobi( 4 , 14 ) =0.0;
    jacobi( 4 , 15 ) =0.0;
    jacobi( 4 , 16 ) =0.0;
    jacobi( 4 , 17 ) =0.0;
    jacobi( 4 , 18 ) =0.0;
    jacobi( 4 , 19 ) =0.0;
    jacobi( 4 , 20 ) =0.0;
    jacobi( 5 , 0 ) =-k6(Te)*n[5] -k7(Te)*n[5] -k8(Te)*n[5] -k9(Te)*n[5]-k10(Te)*n[5] ;
    jacobi( 5 , 1 ) =-k18(Tg)*n[5] -k19(Tg)*n[5];
    jacobi( 5 , 2 ) =-k39(Tg)*n[5];
    jacobi( 5 , 3 ) = -k36(Tg)*n[5];
    jacobi( 5 , 4 ) =0.0;
    jacobi( 5 , 5 ) =-k6(Te)*n[0] -k7(Te)*n[0] -k8(Te)*n[0] -k9(Te)*n[0]-k10(Te)*n[0]
    	-k18(Tg)*n[1] -k19(Tg)*n[1] -k24(Tg)*n[6] -k27(Tg)*n[7] -k36(Tg)*n[3]-k39(Tg)*n[2]-C/dSi;
    jacobi( 5 , 6 ) = +2.*k23(Tg)*n[6]-k24(Tg)*n[5] ;
    jacobi( 5 , 7 ) = -k27(Tg)*n[5];
    jacobi( 5 , 8 ) = +k25(Tg)*n[9] ;
    jacobi( 5 , 9 ) =+k25(Tg)*n[8];
    jacobi( 5 , 10 ) =0.0;
    jacobi( 5 , 11 ) =0.0;
    jacobi( 5 , 12 ) =0.0;
    jacobi( 5 , 13 ) =0.0;
    jacobi( 5 , 14 ) =0.0;
    jacobi( 5 , 15 ) =0.0;
    jacobi( 5 , 16 ) =0.0;
    jacobi( 5 , 17 ) =0.0;
    jacobi( 5 , 18 ) =0.0;
    jacobi( 5 , 19 ) =0.0;
    jacobi( 5 , 20 ) =0.0;
    jacobi( 6 , 0 ) =k6(Te)*n[5] -k11(Te)*n[6]-k12(Te)*n[6] +k13(Te)*n[2]-k43(Te)*n[6];
    jacobi( 6 , 1 ) =+k18(Tg)*n[5] -k20(Tg)*n[6] ;
    jacobi( 6 , 2 ) =+k13(Te)*n[0]+k33(Tg)*n[10] -k34(Tg)*n[6] +k42(Tg)*n[20];
    jacobi( 6 , 3 ) =-k37(Tg)*n[6] ;
    jacobi( 6 , 4 ) =+ DA[4];
    jacobi( 6 , 5 ) =k6(Te)*n[0] +k18(Tg)*n[1] -k24(Tg)*n[6]+k27(Tg)*n[7] ;
    jacobi( 6 , 6 ) = -k11(Te)*n[0] -k12(Te)*n[0]-k20(Tg)*n[1] -2*2.*k23(Tg)*n[6]-k24(Tg)*n[5]
    	-k30(Tg) -k31(Tg)*n[7]  -k34(Tg)*n[2]-k37(Tg)*n[3] -k43(Te)*n[0] - C  /dSi;
    jacobi( 6 , 7 ) =+k27(Tg)*n[5] -k31(Tg)*n[6] ;
    jacobi( 6 , 8 ) =0.0;
    jacobi( 6 , 9 ) =0.0;
    jacobi( 6 , 10 ) =+k33(Tg)*n[2] ;
    jacobi( 6 , 11 ) =0.0;
    jacobi( 6 , 12 ) =0.0;
    jacobi( 6 , 13 ) =0.0;
    jacobi( 6 , 14 ) =0.0;
    jacobi( 6 , 15 ) =0.0;
    jacobi( 6 , 16 ) =0.0;
    jacobi( 6 , 17 ) =0.0;
    jacobi( 6 , 18 ) =0.0;
    jacobi( 6 , 19 ) =0.0;
    jacobi( 6 , 20 ) =+k42(Tg)*n[2];
    jacobi( 7 , 0 ) = k6(Te)*n[5]* +2*k7(Te)*n[5] +k8(Te)*n[5]*+2*k9(Te)*n[5]
    	+k10(Te)*n[5] +k11(Te)*n[6] +2*k16(Te)*n[9] ;
    jacobi( 7 , 1 ) = +k18(Tg)*n[5]+2*k19(Tg)*n[5] +k20(Tg)*n[6] +k21(Tg)*n[8] +2*k22(Tg)*n[9];
    jacobi( 7 , 2 ) =0.0;
    jacobi( 7 , 3 ) = 0.0;
    jacobi( 7 , 4 ) = 0.0;
    jacobi( 7 , 5 ) = k6(Te)*n[0] +2*k7(Te)*n[0] +k8(Te)*n[0] +2*k9(Te)*n[0]
    	+k10(Te)*n[0] +k18(Tg)*n[1]+2*k19(Tg)*n[1] -k27(Tg)*n[7] ;
    jacobi( 7 , 6 ) =  +k11(Te)*n[0] +k20(Tg)*n[1] -k31(Tg)*n[7];
    jacobi( 7 , 7 ) = -k27(Tg)*n[5] -k29(Tg)*n[8] -k31(Tg)*n[6] - C  /dSi;
    jacobi( 7 , 8 ) = +k21(Tg)*n[1] -k29(Tg)*n[7] ;
    jacobi( 7 , 9 ) =  +2*k16(Te)*n[0] +2*k22(Tg)*n[1];
    jacobi( 7 , 10 ) = 0.0;
    jacobi( 7 , 11 ) = 0.0;
    jacobi( 7 , 12 ) = 0.0;
    jacobi( 7 , 13 ) = 0.0;
    jacobi( 7 , 14 ) = 0.0;
    jacobi( 7 , 15 ) = 0.0;
    jacobi( 7 , 16 ) = 0.0;
    jacobi( 7 , 17 ) = 0.0;
    jacobi( 7 , 18 ) = 0.0;
    jacobi( 7 , 19 ) = 0.0;
    jacobi( 7 , 20 ) = 0.0;
    jacobi( 8 , 0 ) = k7(Te)*n[5] +k14(Te)*n[3] -k15(Te)*n[8] ;
    jacobi( 8 , 1 ) =  +k19(Tg)*n[5]+k20(Tg)*n[6] -k21(Tg)*n[8] ;
    jacobi( 8 , 2 ) = -k38(Tg)*n[8] ;
    jacobi( 8 , 3 ) = +k14(Te)*n[0] +k32(Tg)*n[10] +k37(Tg)*n[6] +k41(Tg)*n[20] +k46(Tg)*n[18];
    jacobi( 8 , 4 ) =0.0;
    jacobi( 8 , 5 ) = k7(Te)*n[0]  +k19(Tg)*n[1];
    jacobi( 8 , 6 ) = +k20(Tg)*n[1]  +k23(Tg)*2.*n[6]+k31(Tg)*n[7] +k37(Tg)*n[3] ;
    jacobi( 8 , 7 ) =  -k29(Tg)*n[8]+k31(Tg)*n[6] ;
    jacobi( 8 , 8 ) =-k15(Te)*n[0]  -k21(Tg)*n[1]  -k25(Tg)*n[9]-k26(Tg)-2*k28(Tg)*2.*n[8]
    	-k29(Tg)*n[7]-k38(Tg)*n[2]- C  /dSi ;
    jacobi( 8 , 9 ) =  -k25(Tg)*n[8];
    jacobi( 8 , 10 ) = +k32(Tg)*n[3];
    jacobi( 8 , 11 ) =0.0;
    jacobi( 8 , 12 ) =0.0;
    jacobi( 8 , 13 ) =0.0;
    jacobi( 8 , 14 ) =0.0;
    jacobi( 8 , 15 ) =0.0;
    jacobi( 8 , 16 ) =0.0;
    jacobi( 8 , 17 ) =0.0;
    jacobi( 8 , 18 ) =  +k46(Tg)*n[3];
    jacobi( 8 , 19 ) =0.0;
    jacobi( 8 , 20 ) = +k41(Tg)*n[3];
    jacobi( 9 , 0 ) =-k16(Te)*n[9] -k17(Te)*n[9];
    jacobi( 9 , 1 ) = -k22(Tg)*n[9];
    jacobi( 9 , 2 ) =+k33(Tg)*n[10] +k34(Tg)*n[6] +k38(Tg)*n[8] +k39(Tg)*n[5] ;
    jacobi( 9 , 3 ) = +k32(Tg)*n[10] + k36(Tg)*n[5];
    jacobi( 9 , 4 ) =0.0;
    jacobi( 9 , 5 ) = +k24(Tg)*n[6]+k27(Tg)*n[7] + k36(Tg)*n[3] +k39(Tg)*n[2];
    jacobi( 9 , 6 ) =+k24(Tg)*n[5]+k30(Tg)+k31(Tg)*n[7] +k34(Tg)*n[2] ;
    jacobi( 9 , 7 ) =+k27(Tg)*n[5]+k29(Tg)*n[8]+k31(Tg)*n[6];
    jacobi( 9 , 8 ) =-k25(Tg)*n[9] +k26(Tg)+k28(Tg)*2.*n[8]+k29(Tg)*n[7]+k38(Tg)*n[2];
    jacobi( 9 , 9 ) =-k16(Te)*n[0] -k17(Te)*n[0] -k22(Tg)*n[1] -k25(Tg)*n[8]- C  /dSi;
    jacobi( 9 , 10 ) = +k32(Tg)*n[3] +k33(Tg)*n[2] +k47(Tg)*n[17]+ DA[10];
    jacobi( 9 , 11 ) =0.0;
    jacobi( 9 , 12 ) =0.0;
    jacobi( 9 , 13 ) =0.0;
    jacobi( 9 , 14 ) =0.0;
    jacobi( 9 , 15 ) =0.0;
    jacobi( 9 , 16 ) =0.0;
    jacobi( 9 , 17 ) = +k47(Tg)*n[10];
    jacobi( 9 , 18 ) =0.0;
    jacobi( 9 , 19 ) =0.0;
    jacobi( 9 , 20 ) =0.0;
    jacobi( 10 , 0 ) =k17(Te)*n[9] ;
    jacobi( 10 , 1 ) =0.0;
    jacobi( 10 , 2 ) =-k33(Tg)*n[10] ;
    jacobi( 10 , 3 ) =-k32(Tg)*n[10] ;
    jacobi( 10 , 4 ) =0.0;
    jacobi( 10 , 5 ) =0.0;
    jacobi( 10 , 6 ) =0.0;
    jacobi( 10 , 7 ) =0.0;
    jacobi( 10 , 8 ) =0.0;
    jacobi( 10 , 9 ) =k17(Te)*n[0] ;
    jacobi( 10 , 10 ) = -k32(Tg)*n[3] -k33(Tg)*n[2] -k47(Tg)*n[17] -DA[10] - C  /dSi;
    jacobi( 10 , 11 ) =0.0;
    jacobi( 10 , 12 ) =0.0;
    jacobi( 10 , 13 ) =0.0;
    jacobi( 10 , 14 ) =0.0;
    jacobi( 10 , 15 ) =0.0;
    jacobi( 10 , 16 ) =0.0;
    jacobi( 10 , 17 ) =-k47(Tg)*n[10];
    jacobi( 10 , 18 ) =0.0;
    jacobi( 10 , 19 ) =0.0;
    jacobi( 10 , 20 ) =0.0;
    jacobi( 11 , 0 ) =0.0;
    jacobi( 11 , 1 ) =0.0;
    jacobi( 11 , 2 ) =0.0;
    jacobi( 11 , 3 ) =+k40(Tg)*n[4];
    jacobi( 11 , 4 ) = +k40(Tg)*n[3];
    jacobi( 11 , 5 ) =k24(Tg)*n[6] ;
    jacobi( 11 , 6 ) =k24(Tg)*n[5];
    jacobi( 11 , 7 ) =0.0;
    jacobi( 11 , 8 ) =0.0;
    jacobi( 11 , 9 ) =0.0;
    jacobi( 11 , 10 ) =0.0;
    jacobi( 11 , 11 ) =- C  /dSi;
    jacobi( 11 , 12 ) =0.0;
    jacobi( 11 , 13 ) =0.0;
    jacobi( 11 , 14 ) =0.0;
    jacobi( 11 , 15 ) =0.0;
    jacobi( 11 , 16 ) =0.0;
    jacobi( 11 , 17 ) =0.0;
    jacobi( 11 , 18 ) =0.0;
    jacobi( 11 , 19 ) =0.0;
    jacobi( 11 , 20 ) =0.0;
    jacobi( 12 , 0 ) =0.0;
    jacobi( 12 , 1 ) =0.0;
    jacobi( 12 , 2 ) =0.0;
    jacobi( 12 , 3 ) =0.0;
    jacobi( 12 , 4 ) =0.0;
    jacobi( 12 , 5 ) =0.0;
    jacobi( 12 , 6 ) =0.0;
    jacobi( 12 , 7 ) =0.0;
    jacobi( 12 , 8 ) =k28(Tg)*2.*n[8];
    jacobi( 12 , 9 ) =0.0;
    jacobi( 12 , 10 ) =0.0;
    jacobi( 12 , 11 ) =0.0;
    jacobi( 12 , 12 ) =- C  /dSi;
    jacobi( 12 , 13 ) =0.0;
    jacobi( 12 , 14 ) =0.0;
    jacobi( 12 , 15 ) =0.0;
    jacobi( 12 , 16 ) =0.0;
    jacobi( 12 , 17 ) =0.0;
    jacobi( 12 , 18 ) =0.0;
    jacobi( 12 , 19 ) =0.0;
    jacobi( 12 , 20 ) =0.0;
    jacobi( 13 , 0 ) =0.0;
    jacobi( 13 , 1 ) =0.0;
    jacobi( 13 , 2 ) =k34(Tg)*n[6];
    jacobi( 13 , 3 ) = +k36(Tg)*n[5];
    jacobi( 13 , 4 ) =0.0;
    jacobi( 13 , 5 ) = +k36(Tg)*n[3] ;
    jacobi( 13 , 6 ) =k34(Tg)*n[2];
    jacobi( 13 , 7 ) =0.0;
    jacobi( 13 , 8 ) =0.0;
    jacobi( 13 , 9 ) =0.0;
    jacobi( 13 , 10 ) =0.0;
    jacobi( 13 , 11 ) =0.0;
    jacobi( 13 , 12 ) =0.0;
    jacobi( 13 , 13 ) =- C  /dSi ;
    jacobi( 13 , 14 ) =0.0;
    jacobi( 13 , 15 ) =0.0;
    jacobi( 13 , 16 ) =0.0;
    jacobi( 13 , 17 ) =0.0;
    jacobi( 13 , 18 ) =0.0;
    jacobi( 13 , 19 ) =0.0;
    jacobi( 13 , 20 ) =0.0;
    jacobi( 14 , 0 ) =0.0;
    jacobi( 14 , 1 ) =0.0;
    jacobi( 14 , 2 ) =k35(Tg)*n[4];
    jacobi( 14 , 3 ) =0.0;
    jacobi( 14 , 4 ) =k35(Tg)*n[2];
    jacobi( 14 , 5 ) =0.0;
    jacobi( 14 , 6 ) =0.0;
    jacobi( 14 , 7 ) =0.0;
    jacobi( 14 , 8 ) =0.0;
    jacobi( 14 , 9 ) =0.0;
    jacobi( 14 , 10 ) =0.0;
    jacobi( 14 , 11 ) =0.0;
    jacobi( 14 , 12 ) =0.0;
    jacobi( 14 , 13 ) =0.0;
    jacobi( 14 , 14 ) =- C  /dSi;
    jacobi( 14 , 15 ) =0.0;
    jacobi( 14 , 16 ) =0.0;
    jacobi( 14 , 17 ) =0.0;
    jacobi( 14 , 18 ) =0.0;
    jacobi( 14 , 19 ) =0.0;
    jacobi( 14 , 20 ) =0.0;
    jacobi( 15 , 0 ) =0.0;
    jacobi( 15 , 1 ) =0.0;
    jacobi( 15 , 2 ) =k38(Tg)*n[8];
    jacobi( 15 , 3 ) =0.0;
    jacobi( 15 , 4 ) =0.0;
    jacobi( 15 , 5 ) =0.0;
    jacobi( 15 , 6 ) =0.0;
    jacobi( 15 , 7 ) =0.0;
    jacobi( 15 , 8 ) =k38(Tg)*n[2];
    jacobi( 15 , 9 ) =0.0;
    jacobi( 15 , 10 ) =0.0;
    jacobi( 15 , 11 ) =0.0;
    jacobi( 15 , 12 ) =0.0;
    jacobi( 15 , 13 ) =0.0;
    jacobi( 15 , 14 ) =0.0;
    jacobi( 15 , 15 ) =- C  /dSi;
    jacobi( 15 , 16 ) =0.0;
    jacobi( 15 , 17 ) =0.0;
    jacobi( 15 , 18 ) =0.0;
    jacobi( 15 , 19 ) =0.0;
    jacobi( 15 , 20 ) =0.0;
    jacobi( 16 , 0 ) =0.0;
    jacobi( 16 , 1 ) =0.0;
    jacobi( 16 , 2 ) =k39(Tg)*n[5];
    jacobi( 16 , 3 ) =0.0;
    jacobi( 16 , 4 ) =0.0;
    jacobi( 16 , 5 ) =k39(Tg)*n[2];
    jacobi( 16 , 6 ) =0.0;
    jacobi( 16 , 7 ) =0.0;
    jacobi( 16 , 8 ) =0.0;
    jacobi( 16 , 9 ) =0.0;
    jacobi( 16 , 10 ) =0.0;
    jacobi( 16 , 11 ) =0.0;
    jacobi( 16 , 12 ) =0.0;
    jacobi( 16 , 13 ) =0.0;
    jacobi( 16 , 14 ) =0.0;
    jacobi( 16 , 15 ) =0.0;
    jacobi( 16 , 16 ) =- C  /dSi;
    jacobi( 16 , 17 ) =0.0;
    jacobi( 16 , 18 ) =0.0;
    jacobi( 16 , 19 ) =0.0;
    jacobi( 16 , 20 ) =0.0;
    jacobi( 17 , 0 ) =k44(Te)*n[18] -k45(Te)*n[17] ;
    jacobi( 17 , 1 ) =0.0;
    jacobi( 17 , 2 ) =0.0;
    jacobi( 17 , 3 ) = +k46(Tg)*n[18] ;
    jacobi( 17 , 4 ) =0.0;
    jacobi( 17 , 5 ) =0.0;
    jacobi( 17 , 6 ) =0.0;
    jacobi( 17 , 7 ) =0.0;
    jacobi( 17 , 8 ) =0.0;
    jacobi( 17 , 9 ) =0.0;
    jacobi( 17 , 10 ) = -k47(Tg)*n[17];
    jacobi( 17 , 11 ) =0.0;
    jacobi( 17 , 12 ) =0.0;
    jacobi( 17 , 13 ) =0.0;
    jacobi( 17 , 14 ) =0.0;
    jacobi( 17 , 15 ) =0.0;
    jacobi( 17 , 16 ) =0.0;
    jacobi( 17 , 17 ) = -k45(Te)*n[0]  -k47(Tg)*n[10]- C  /dSi; 
    jacobi( 17 , 18 ) =k44(Te)*n[0]  +k46(Tg)*n[3];
    jacobi( 17 , 19 ) =0.0;
    jacobi( 17 , 20 ) =0.0;
    jacobi( 18 , 0 ) = -k44(Te)*n[18] +k45(Te)*n[17] ;
    jacobi( 18 , 1 ) = k21(Tg)*n[8] ;
    jacobi( 18 , 2 ) = 0.0;
    jacobi( 18 , 3 ) = -k46(Tg)*n[18];
    jacobi( 18 , 4 ) = 0.0;
    jacobi( 18 , 5 ) = 0.0;
    jacobi( 18 , 6 ) = +k30(Tg) ;
    jacobi( 18 , 7 ) =  +k29(Tg)*n[8];
    jacobi( 18 , 8 ) = k21(Tg)*n[1] +k29(Tg)*n[7] ;
    jacobi( 18 , 9 ) = 0.0;
    jacobi( 18 , 10 ) =  +k47(Tg)*n[17];
    jacobi( 18 , 11 ) = 0.0;
    jacobi( 18 , 12 ) = 0.0;
    jacobi( 18 , 13 ) = 0.0;
    jacobi( 18 , 14 ) = 0.0;
    jacobi( 18 , 15 ) = 0.0;
    jacobi( 18 , 16 ) = 0.0;
    jacobi( 18 , 17 ) = +k45(Te)*n[0] +k47(Tg)*n[10];
    jacobi( 18 , 18 ) = -k44(Te)*n[0] -k46(Tg)*n[3]- C  /dSi;
    jacobi( 18 , 19 ) = 0.0;
    jacobi( 18 , 20 ) = 0.0;
    jacobi( 19 , 0 ) = 0.0;
    jacobi( 19 , 1 ) = 0.0;
    jacobi( 19 , 2 ) = 0.0;
    jacobi( 19 , 3 ) = 0.0;
    jacobi( 19 , 4 ) = 0.0;
    jacobi( 19 , 5 ) = 0.0;
    jacobi( 19 , 6 ) = 0.0;
    jacobi( 19 , 7 ) = 0.0;
    jacobi( 19 , 8 ) = k26(Tg);
    jacobi( 19 , 9 ) = 0.0;
    jacobi( 19 , 10 ) = 0.0;
    jacobi( 19 , 11 ) = 0.0;
    jacobi( 19 , 12 ) = 0.0;
    jacobi( 19 , 13 ) = 0.0;
    jacobi( 19 , 14 ) = 0.0;
    jacobi( 19 , 15 ) = 0.0;
    jacobi( 19 , 16 ) = 0.0;
    jacobi( 19 , 17 ) = 0.0;
    jacobi( 19 , 18 ) = 0.0;
    jacobi( 19 , 19 ) = - C  /dSi;
    jacobi( 19 , 20 ) = 0.0;
    jacobi( 20 , 0 ) = k1(Te)*n_Ar +k3(Te)*n[1] ;
    jacobi( 20 , 1 ) =+k3(Te)*n[0] +k4(Tg)*2.*n[1];
    jacobi( 20 , 2 ) =-k42(Tg)*n[20];
    jacobi( 20 , 3 ) =-k41(Tg)*n[20];
    jacobi( 20 , 4 ) =0.0;
    jacobi( 20 , 5 ) =0.0;
    jacobi( 20 , 6 ) =0.0;
    jacobi( 20 , 7 ) =0.0;
    jacobi( 20 , 8 ) =0.0;
    jacobi( 20 , 9 ) =0.0;
    jacobi( 20 , 10 ) =0.0;
    jacobi( 20 , 11 ) =0.0;
    jacobi( 20 , 12 ) =0.0;
    jacobi( 20 , 13 ) =0.0;
    jacobi( 20 , 14 ) =0.0;
    jacobi( 20 , 15 ) =0.0;
    jacobi( 20 , 16 ) =0.0;
    jacobi( 20 , 17 ) =0.0;
    jacobi( 20 , 18 ) =0.0;
    jacobi( 20 , 19 ) =0.0;
    jacobi( 20 , 20 ) =-k41(Tg)*n[3]-k42(Tg)*n[2] -DA[20]- C  /dSi;
    

    dfdt( 0 ) = 0.0;
    dfdt( 1 ) = 0.0;
    dfdt( 2 ) = 0.0;
    dfdt( 3 ) = 0.0;
    dfdt( 4 ) = 0.0;
    dfdt( 5 ) = 0.0;
    dfdt( 6 ) = 0.0;
    dfdt( 7 ) = 0.0;
    dfdt( 8 ) = 0.0;
    dfdt( 9 ) = 0.0;
    dfdt( 10 ) = 0.0;
    dfdt( 11 ) = 0.0;
    dfdt( 12 ) = 0.0;
    dfdt( 13 ) = 0.0;
    dfdt( 14 ) = 0.0;
    dfdt( 15 ) = 0.0;
    dfdt( 16 ) = 0.0;
    dfdt( 17 ) = 0.0;
    dfdt( 18 ) = 0.0;
    dfdt( 19 ) = 0.0;
    dfdt( 20 ) = 0.0;


  }

  value_type Te;
};

struct etemperature
{
  value_type operator()(value_type const& Te)
  {// DP/ n[0] puissance par unité de volume DP constante
    return 
      -DP/n[0]                                               
       + k1(Te)*n_Ar*16.14 + k2(Te)*n_Ar*12.31 + k3(Te)*n[1]*5.39
       - k4(Tg)*n[1]*n[1]*8.48   
       - k5(Te)*n[1]*12.31  + k6(Te)*n[5]*10.68 + k7(Te)*n[5]*10.68
       + k8(Te)*n[5]*8.29 + k9(Te)*n[5]*8.29 + k10(Te)*n[5]*24.1
       + k11(Te)*n[6]*1.94 + k12(Te)*n[6]*1.30 + k13(Te)*n[2]*1.16
       + k14(Te)*n[3]*1.16 + k15(Te)*n[8]*1.5*Te 
       + k16(Te)*n[9]*10.09
       + k17(Te)*n[9]*16.05 + k43(Te)*n[6]*1.5*Te 
       + k44(Te)*n[18]*1.5*Te + k45(Te)*n[17]*1.25
       + DA[0]*4.5*Te;                            // pertes sur les parois  ! sqrt(M/(2pim) ~ 4.5 venant du potentiel de gaine     
  }

  state_type n;
};


void write_density( const value_type t, const value_type Te, const state_type &n)
{
  cout << t  << '\t' <<Te <<'\t' << n[0] << '\t' << n[1] << '\t'
               << n[2] << '\t' << n[3] <<'\t'<< n[4] << '\t' << n[5] << '\t'
               << n[6] << '\t' << n[7] << '\t' << n[8] << '\t' << n[9] << '\t'
               << n[10] << '\t' << n[11] << '\t' << n[12] << '\t'
               << n[13] << '\t' << n[14] << '\t' << n[15] << '\t'
               << n[16] << '\t' << n[17] << '\t'<< n[18] << '\t'
               << n[19] << '\t' << n[20]<< '\t' <<n[13]+n[16]<<endl;
}

int main(int argc, char **argv)
{ 
value_type Te=0.7;

    /*0=e, 1=Armet, 2=SiH3-, 3=SiH2-, 4=SiH3+, 5=SiH4, 6=SiH3,
    7=H, 8=SiH2, 9=H2, 10=H2+, 11=Si2H5, 12=Si2H2, 13=Si2H4-,
    14=Si2H6, 15=Si2H3-, 16=Si2H5-, 17=SiH-, 18=SiH, 19=Si, 20=Arp*/

  cout <<"t"<<'\t'<<"Te"<<'\t'<<"e"<<'\t'<<"Armet"<<'\t'<< "SiH3m"<<'\t'
               << "SiH2m"<<'\t'<< "SiH3p"<<'\t'<< "SiH4"<<'\t'<< "SiH3"<<'\t'
               <<"H"<<'\t'<< "SiH2"<<'\t'<< "H2"<<'\t'<< "H2p"<<'\t'<< "Si2H5"
               <<'\t'<< "Si2H2"<<'\t'<<"Si2H4m"<<'\t'<<"Si2H6"<<'\t'<< "Si2H3m"
               <<'\t'<< "Si2H5m"<<'\t'<< "SiHm"<<'\t'<<"SiH"<<'\t'<< "Si"<<'\t'
               << "Arp"<< '\t'<<"NP"<<endl;
//clock_t t1,t2;

// vecteur des densites initiales
state_type n_ini(Nbr_espece, 0.0); 
  n_ini[0] = 1.e16;
  n_ini[1] = 1.e10; 
  n_ini[5] = n_SiH4_ini;
  n_ini[20] =1.e16;
  n_ini[4] = n_ini[0]-n_ini[20];

state_type DL(Nbr_espece, 0.0); //vecteur de diffusion libre en m2/s
state_type mu(Nbr_espece, 0.0); //vecteur de mobilite en m2/(V.s)
state_type DA(Nbr_espece, 0.0); //vecteur de diffusion ambipolaire en m2/s

//Coefficients de diffusion libres de Chapman-Enskog
value_type D_mol=2.; // diametre de (molecule + argon)/2 en A
value_type D_e=1.; // diametre de (electron+ argon)/2 en A

value_type CE_e= 1.858e-3*(Tg*1.1604e4)*sqrt(Te*1.1604e4)/(pression*pow(D_e,2.))*1.e-4;//on met les temperatures en K et on convertis pour l'avoir en m2/s (*1.e-4)
value_type CE_mol= 1.858e-3*pow((Tg*1.1604e4),3./2.)/(pression*pow(D_mol,2.))*1.e-4;//on met les temperatures en K et on convertis pour l'avoir en m2/s (*1.e-4) 

 
/*DL[0]= CE_e*sqrt(1836.2 + 1./40.)/2.; // OMEGA = 2 et masse atomique de l'electron =1/1836.2    
mu[0]= DL[0]/Te; //Te en eV et DL en m2/s */
DL[0]=120./(pression*760.); //valeur de benjamin
mu[0]= DL[0]/Te;

DL[1]=0.075; //valeur de benjamin en m2/s
	
DL[4]= CE_mol*sqrt(1./31.+1./40.)/10.; // OMEGA = 10         
mu[4]= DL[4]/Tg;

DL[10]= CE_mol*sqrt(1./2.+1./40.)/10.; // OMEGA = 10          
mu[10]= DL[10]/Tg;

/*DL[20]= CE_mol*sqrt(2./40.)/10.; // OMEGA = 10         
mu[20]= DL[20]/Tg;*/
DL[20]=4.e-3/(pression*760.); //valeur de benjamin
mu[20]= DL[20]/Tg;

// Time variables
  value_type t = 0.0;
  value_type dt = 1.0e-8;
  value_type Tmax = 20.e-3;
  value_type NT = Tmax/dt;

  // Root finding variables
  value_type min = Tg;
  value_type max = 10.0;
  boost::uintmax_t max_iter = 500;
  eps_tolerance<value_type> tol(30);

  state_type n_new(Nbr_espece, 0.0);
n_new=n_ini;
  state_type n_err(Nbr_espece, 0.0); //error

   // declare the functor etemperature
  etemperature etemp;
  // assign initial values to functor etemp
  etemp.n = n_ini;

// Find Te first calculation
  pair<value_type, value_type> pair_Te =\
                toms748_solve(etemp, min, max, tol, max_iter);

  Te = pair_Te.first;
  cerr << "\n[ii] Initial Temperature  = " << Te << endl;

  // declare system and jacobian
  nsystem sys;
  jacobian jac;

  // declare stepper Rosenbrock
  stepper_type stepper;
//cerr<<"patate1"<<endl;
  for (int i = 1; i <= NT+1; i++)
  {
value_type n_mu= n_new[0]*mu[0] + n_new[4]*mu[4] + n_new[10]*mu[10] + n_new[20]*mu[20];
value_type n_DL= -n_new[0]*DL[0] + n_new[4]*DL[4] + n_new[10]*DL[10] + n_new[20]*DL[20];

DA[0]=(DL[0]+mu[0]*n_DL/n_mu)*diff; //s-1
DA[1]=DL[1]*diff;
DA[4]=(DL[4]-mu[4]*n_DL/n_mu)*diff;
DA[10]=(DL[10]-mu[10]*n_DL/n_mu)*diff;
DA[20]=(DL[20]-mu[20]*n_DL/n_mu)*diff;

    // update Te in system and jacobian
    sys.Te = Te;
    jac.Te = Te;
//cerr<<"patate2"<<endl;
    // Integrate at least one step dt
    stepper.do_step( std::make_pair( sys, jac ), n_new, t, dt, n_err);
//cerr<<"patate3"<<endl;
    // assign values to functor etemp
    etemp.n = n_new;
    if (i%((int)(NT/50))==0)
    {
      write_density(t, Te, n_new);
    }
   // cerr<<"patate4"<<endl;
    // Find new Te
    pair<value_type, value_type> pair_Te =\
                  toms748_solve(etemp, min, max, tol, max_iter);
   // cerr<<"patate5"<<endl;
    Te = pair_Te.first;
    t+= dt;
    n_ini=n_new;
  }

 value_type charge= (n_new[20]+n_new[4]+n_new[10]-n_new[0]-n_new[2]-n_new[3]-n_new[13]-n_new[15]-n_new[16]-n_new[17])/(n_new[20]+n_new[4]+n_new[10]);

  cerr<<"charge/dArp="<<charge<<endl;

 value_type Si=(n_new[2]+n_new[3]+n_new[4]+n_new[13]*2+2*n_new[15]+n_new[16]*2
          +n_new[17]+n_new[5]+n_new[6]+n_new[8]+n_new[18]+2*n_new[11]+n_new[19]
          +n_new[12]*2+n_new[14]*2)/n_SiH4_ini;

  cerr<<"Si="<<Si<<endl;


  value_type H=(3*n_new[2]+2*n_new[3]+3*n_new[4]+4*n_new[5]+3*n_new[6]+n_new[7]
	+2*n_new[8]+2*n_new[9]+2*n_new[10]+5* n_new[11]+2*n_new[12]
	+4*n_new[13]+6*n_new[14]+3*n_new[15]+5*n_new[16]+n_new[17]
         +n_new[18])/(4.*n_SiH4_ini);

  cerr<<"H="<<H<<endl;

//cerr<<Da_e<<endl;
  return 0;

}

