/**
 * @file testdata_expected.hh
 * @brief Auto-generated test data from Mathematica
 * DO NOT EDIT - regenerate with: math -script test/generate_testdata.wl
 */
#ifndef CLPOLY_TESTDATA_EXPECTED_HH
#define CLPOLY_TESTDATA_EXPECTED_HH

#include <clpoly/clpoly.hh>
#include <vector>
#include <utility>

namespace testdata {

using namespace clpoly;

// ==================== Arithmetic ====================

inline polynomial_ZZ arith_f1_simple() {
    variable x("x");
    return 3*pow(x,3)-2*pow(x,2)+x+polynomial_ZZ({{{},-5}});
}

inline polynomial_ZZ arith_g1_simple() {
    variable x("x");
    return pow(x,2)+2*x+1;
}

inline polynomial_ZZ arith_add_simple() {
    variable x("x");
    return 3*pow(x,3)-pow(x,2)+3*x+polynomial_ZZ({{{},-4}});
}

inline polynomial_ZZ arith_sub_simple() {
    variable x("x");
    return 3*pow(x,3)-3*pow(x,2)-x+polynomial_ZZ({{{},-6}});
}

inline polynomial_ZZ arith_mul_simple() {
    variable x("x");
    return 3*pow(x,5)+4*pow(x,4)-5*pow(x,2)-9*x+polynomial_ZZ({{{},-5}});
}

inline polynomial_ZZ arith_pow2_simple() {
    variable x("x");
    return 9*pow(x,6)-12*pow(x,5)+10*pow(x,4)-34*pow(x,3)+21*pow(x,2)-10*x+polynomial_ZZ({{{},25}});
}

inline polynomial_ZZ arith_pow3_simple() {
    variable x("x");
    return 27*pow(x,9)-54*pow(x,8)+63*pow(x,7)-179*pow(x,6)+201*pow(x,5)-156*pow(x,4)+286*pow(x,3)-165*pow(x,2)+75*x+polynomial_ZZ({{{},-125}});
}

inline polynomial_ZZ arith_f1_medium() {
    variable x("x"), y("y");
    return 2*pow(x,3)*pow(y,2)-3*x*pow(y,4)+pow(x,2)*y+polynomial_ZZ({{{},-7}});
}

inline polynomial_ZZ arith_g1_medium() {
    variable x("x"), y("y");
    return pow(x,2)*y-pow(y,3)+2*x+1;
}

inline polynomial_ZZ arith_mul_medium() {
    variable x("x"), y("y");
    return 2*pow(x,5)*pow(y,3)-5*pow(x,3)*pow(y,5)+3*x*pow(y,7)+5*pow(x,4)*pow(y,2)-7*pow(x,2)*pow(y,4)+2*pow(x,3)*pow(y,2)-3*x*pow(y,4)+2*pow(x,3)*y-6*pow(x,2)*y+7*pow(y,3)-14*x+polynomial_ZZ({{{},-7}});
}

inline polynomial_ZZ arith_f1_complex() {
    variable x("x"), y("y"), z("z");
    return pow(x,3)*y*z-2*pow(x,2)*pow(z,3)+3*pow(y,2)*pow(z,2)-x*y+polynomial_ZZ({{{},5}});
}

inline polynomial_ZZ arith_g1_complex() {
    variable x("x"), y("y"), z("z");
    return x*pow(y,2)*z+pow(z,3)-pow(x,2)+y-1;
}

inline polynomial_ZZ arith_mul_complex() {
    variable x("x"), y("y"), z("z");
    return pow(x,4)*pow(y,3)*pow(z,2)-2*pow(x,3)*pow(y,2)*pow(z,4)+pow(x,3)*y*pow(z,4)-2*pow(x,2)*pow(z,6)+3*x*pow(y,4)*pow(z,3)-pow(x,5)*y*z+2*pow(x,4)*pow(z,3)+3*pow(y,2)*pow(z,5)+pow(x,3)*pow(y,2)*z-pow(x,2)*pow(y,3)*z-3*pow(x,2)*pow(y,2)*pow(z,2)-2*pow(x,2)*y*pow(z,3)-pow(x,3)*y*z+2*pow(x,2)*pow(z,3)-x*y*pow(z,3)+3*pow(y,3)*pow(z,2)+pow(x,3)*y+5*x*pow(y,2)*z-3*pow(y,2)*pow(z,2)-x*pow(y,2)+5*pow(z,3)-5*pow(x,2)+x*y+5*y+polynomial_ZZ({{{},-5}});
}

// ==================== Prem/Pquo ====================

inline polynomial_ZZ prem_f1() {
    variable x("x");
    return pow(x,3)+2*pow(x,2)-x+polynomial_ZZ({{{},3}});
}

inline polynomial_ZZ prem_g1() {
    variable x("x");
    return pow(x,2)-x+1;
}

inline polynomial_ZZ prem_result1() {
    variable x("x");
    return x;
}

inline polynomial_ZZ pquo_result1() {
    variable x("x");
    return x+polynomial_ZZ({{{},3}});
}

inline polynomial_ZZ prem_f2() {
    variable x("x"), y("y");
    return pow(x,3)+pow(x,2)*y+x*pow(y,2)+pow(y,3)-x;
}

inline polynomial_ZZ prem_g2() {
    variable x("x"), y("y");
    return pow(x,2)+x*y+pow(y,2);
}

inline polynomial_ZZ prem_result2() {
    variable x("x"), y("y");
    return pow(y,3)-x;
}

inline polynomial_ZZ pquo_result2() {
    variable x("x"), y("y");
    return x;
}

// ==================== GCD ====================

inline polynomial_ZZ gcd_f1() {
    variable x("x");
    return pow(x,4)-2*pow(x,3)+4*pow(x,2)-2*x+polynomial_ZZ({{{},3}});
}

inline polynomial_ZZ gcd_g1() {
    variable x("x");
    return pow(x,3)+5*pow(x,2)+x+polynomial_ZZ({{{},5}});
}

inline polynomial_ZZ gcd_result1() {
    variable x("x");
    return pow(x,2)+1;
}

inline polynomial_ZZ gcd_f2_coprime() {
    variable x("x");
    return pow(x,3)+x+1;
}

inline polynomial_ZZ gcd_g2_coprime() {
    variable x("x");
    return pow(x,2)+x+1;
}

inline polynomial_ZZ gcd_f3_mv() {
    variable x("x"), y("y");
    return pow(x,3)*y-x*pow(y,2)+pow(x,2)-y;
}

inline polynomial_ZZ gcd_g3_mv() {
    variable x("x"), y("y");
    return x*pow(y,3)+pow(x,2)*y+pow(y,2)+x;
}

inline polynomial_ZZ gcd_result3_mv() {
    variable x("x"), y("y");
    return x*y+1;
}

// ==================== Resultant / Discriminant ====================

inline polynomial_ZZ res_f1() {
    variable x("x");
    return pow(x,2)+x+1;
}

inline polynomial_ZZ res_g1() {
    variable x("x");
    return pow(x,3)-1;
}

inline ZZ res_result1_val() { return ZZ("0"); }

inline polynomial_ZZ res_f2() {
    variable x("x"), y("y");
    return pow(x,2)+x*y+pow(y,2);
}

inline polynomial_ZZ res_g2() {
    variable x("x"), y("y");
    return pow(x,3)-pow(y,3);
}

inline polynomial_ZZ res_result2() {
    variable y("y");
    return polynomial_ZZ();
}

inline polynomial_ZZ res_f3_common() {
    variable x("x"), y("y");
    return pow(x,2)-pow(y,2);
}

inline polynomial_ZZ res_g3_common() {
    variable x("x"), y("y");
    return pow(x,2)+x*y+x+y;
}

inline polynomial_ZZ disc_f1_repeated() {
    variable x("x");
    return pow(x,3)-3*x+polynomial_ZZ({{{},2}});
}

inline ZZ disc_result1_val() { return ZZ("0"); }

inline polynomial_ZZ disc_f2() {
    variable x("x");
    return pow(x,3)+x+1;
}

inline ZZ disc_result2_val() { return ZZ("-31"); }

inline polynomial_ZZ disc_f3_mv() {
    variable x("x"), y("y");
    return pow(x,2)+x*y+pow(y,2)-1;
}

inline polynomial_ZZ disc_result3_mv() {
    variable y("y");
    return -3*pow(y,2)+polynomial_ZZ({{{},4}});
}

// ==================== Squarefree ====================

inline polynomial_ZZ sqf_f1() {
    variable x("x");
    return pow(x,3)-3*x+polynomial_ZZ({{{},-2}});
}

inline polynomial_ZZ sqf_f2() {
    variable x("x");
    return pow(x,5)-pow(x,4)-2*pow(x,3)+2*pow(x,2)+x-1;
}

inline polynomial_ZZ sqf_f3_already() {
    variable x("x");
    return pow(x,3)+x+1;
}

inline polynomial_ZZ sqf_f4_mv() {
    variable x("x"), y("y");
    return pow(x,3)+pow(x,2)*y-x*pow(y,2)-pow(y,3)+pow(x,2)+2*x*y+pow(y,2);
}

// ==================== Complex (random with fixed seed) ====================

inline polynomial_ZZ complex_f1() {
    variable x("x"), y("y"), z("z");
    return 4*pow(x,5)-3*pow(x,4)*z-7*pow(x,2)*pow(y,2)*z-9*x*pow(y,3)*z-7*pow(y,5)+5*pow(x,2)*y*z;
}

inline polynomial_ZZ complex_g1() {
    variable x("x"), y("y"), z("z");
    return -4*pow(x,5)-9*pow(x,4)*z+9*pow(x,3)*pow(y,2)-6*pow(x,3)*y*z+3*x*pow(y,3)*z-pow(y,3)*pow(z,2);
}

inline polynomial_ZZ complex_mul1() {
    variable x("x"), y("y"), z("z");
    return -16*pow(x,10)-24*pow(x,9)*z+36*pow(x,8)*pow(y,2)-24*pow(x,8)*y*z+27*pow(x,8)*pow(z,2)+pow(x,7)*pow(y,2)*z+18*pow(x,7)*y*pow(z,2)+48*pow(x,6)*pow(y,3)*z+63*pow(x,6)*pow(y,2)*pow(z,2)+28*pow(x,5)*pow(y,5)-63*pow(x,5)*pow(y,4)*z+110*pow(x,5)*pow(y,3)*pow(z,2)-18*pow(x,4)*pow(y,5)*z+54*pow(x,4)*pow(y,4)*pow(z,2)+3*pow(x,4)*pow(y,3)*pow(z,3)-63*pow(x,3)*pow(y,7)+42*pow(x,3)*pow(y,6)*z-21*pow(x,3)*pow(y,5)*pow(z,2)-27*pow(x,2)*pow(y,6)*pow(z,2)+7*pow(x,2)*pow(y,5)*pow(z,3)-21*x*pow(y,8)*z+9*x*pow(y,6)*pow(z,3)+7*pow(y,8)*pow(z,2)-20*pow(x,7)*y*z-45*pow(x,6)*y*pow(z,2)+45*pow(x,5)*pow(y,3)*z-30*pow(x,5)*pow(y,2)*pow(z,2)+15*pow(x,3)*pow(y,4)*pow(z,2)-5*pow(x,2)*pow(y,4)*pow(z,3);
}

inline polynomial_ZZ complex_gcd_f() {
    variable x("x"), y("y"), z("z");
    return 4*pow(x,6)*y-3*pow(x,5)*y*z-7*pow(x,3)*pow(y,3)*z-9*pow(x,2)*pow(y,4)*z-7*x*pow(y,6)+4*pow(x,5)*z-3*pow(x,4)*pow(z,2)+5*pow(x,3)*pow(y,2)*z-7*pow(x,2)*pow(y,2)*pow(z,2)-9*x*pow(y,3)*pow(z,2)-7*pow(y,5)*z+4*pow(x,5)-3*pow(x,4)*z-7*pow(x,2)*pow(y,2)*z+5*pow(x,2)*y*pow(z,2)-9*x*pow(y,3)*z-7*pow(y,5)+5*pow(x,2)*y*z;
}

inline polynomial_ZZ complex_gcd_g() {
    variable x("x"), y("y"), z("z");
    return -4*pow(x,6)*y-9*pow(x,5)*y*z+9*pow(x,4)*pow(y,3)-6*pow(x,4)*pow(y,2)*z+3*pow(x,2)*pow(y,4)*z-x*pow(y,4)*pow(z,2)-4*pow(x,5)*z-9*pow(x,4)*pow(z,2)+9*pow(x,3)*pow(y,2)*z-6*pow(x,3)*y*pow(z,2)+3*x*pow(y,3)*pow(z,2)-pow(y,3)*pow(z,3)-4*pow(x,5)-9*pow(x,4)*z+9*pow(x,3)*pow(y,2)-6*pow(x,3)*y*z+3*x*pow(y,3)*z-pow(y,3)*pow(z,2);
}

inline polynomial_ZZ complex_gcd_common() {
    variable x("x"), y("y"), z("z");
    return x*y+z+1;
}

inline polynomial_ZZ complex_res_result() {
    variable y("y"), z("z");
    return 306850544*pow(y,25)-1840444704*pow(y,24)*z+786836512*pow(y,23)*pow(z,2)-756499296*pow(y,22)*pow(z,3)-461159440*pow(y,21)*pow(z,4)-237115392*pow(y,20)*pow(z,5)-450705664*pow(y,19)*pow(z,6)+1254669648*pow(y,18)*pow(z,7)-1003179008*pow(y,17)*pow(z,8)+273473744*pow(y,16)*pow(z,9)-255636064*pow(y,15)*pow(z,10)-53161776*pow(y,14)*pow(z,11)-173664*pow(y,13)*pow(z,12)-3888*pow(y,12)*pow(z,13)+55319040*pow(y,23)*z+902027280*pow(y,22)*pow(z,2)-1915116000*pow(y,21)*pow(z,3)+810552960*pow(y,20)*pow(z,4)-558284160*pow(y,19)*pow(z,5)-838933200*pow(y,18)*pow(z,6)+748560*pow(y,17)*pow(z,7)-331603440*pow(y,16)*pow(z,8)+17228560*pow(y,15)*pow(z,9)-57210480*pow(y,14)*pow(z,10)+61466400*pow(y,13)*pow(z,11)+544320*pow(y,12)*pow(z,12)-44452800*pow(y,21)*pow(z,2)+72441600*pow(y,20)*pow(z,3)+401864400*pow(y,19)*pow(z,4)-460773600*pow(y,18)*pow(z,5)+368863200*pow(y,17)*pow(z,6)-255241200*pow(y,16)*pow(z,7)-41914800*pow(y,15)*pow(z,8)-90070000*pow(y,14)*pow(z,9)-15393600*pow(y,13)*pow(z,10)-33825600*pow(y,12)*pow(z,11)-194400*pow(y,11)*pow(z,12)-21168000*pow(y,18)*pow(z,4)+14112000*pow(y,17)*pow(z,5)+44982000*pow(y,16)*pow(z,6)-17334000*pow(y,15)*pow(z,7)+54924000*pow(y,14)*pow(z,8)+8942000*pow(y,13)*pow(z,9)-246000*pow(y,12)*pow(z,10)+14256000*pow(y,11)*pow(z,11)-2520000*pow(y,15)*pow(z,6)-1080000*pow(y,13)*pow(z,8)-4370000*pow(y,12)*pow(z,9)+1980000*pow(y,11)*pow(z,10)-2430000*pow(y,10)*pow(z,11)+200000*pow(y,11)*pow(z,9);
}

} // namespace testdata

#endif
