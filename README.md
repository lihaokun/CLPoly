# CLPoly

CLPoly 是一个开发中的多项式c++库，目标是一个高效易且的多项式库。

CLPoly is a c++ library of polynomial.

## Example

```
#include "clpoly.hh"

int main(int argc, char const *argv[])
{
    clpoly::variable x("x");
    clpoly::variable y("y");
    clpoly::variable z("z");
    clpoly::polynomial_ZZ f=-6*pow(x,4)-7*pow(x,2)*y-x*pow(y,3)*z-3*x*y-9;
    clpoly::polynomial_ZZ g=10*pow(x,4)*y+2*pow(x,4)-7*pow(x,3)*pow(y,2)+9*pow(x,2)+9*pow(z,4)-9*y+2*z;
    std::cout<<"f="<<f<<std::endl;
    std::cout<<"g="<<g<<std::endl;
    std::cout<<"prem(f,g,x)="<<clpoly::prem(f,g,x)<<std::endl;
    std::cout<<"resultant(f,g,x)="<<clpoly::resultant(f,g,x)<<std::endl;
}
```
编译指令
```
g++ example.cc -I./ -O3 -lgmpxx -lgmp  
```
结果
```
f=-x*y^3*z-6*x^4-7*x^2*y-3*x*y-9
g=10*x^4*y-7*x^3*y^2+2*x^4+9*z^4+9*x^2-9*y+2*z
prem(f,g,x)=-10*x*y^4*z-42*x^3*y^2-2*x*y^3*z-70*x^2*y^2+54*z^4-14*x^2*y-30*x*y^2+54*x^2-6*x*y-144*y+12*z-18
resultant(f,g,x)=9000*y^15*z^8+18522*y^15*z^7+5400*y^14*z^8+44100*y^14*z^7+1080*y^13*z^8-9000*y^16*z^4+2000*y^15*z^5+125640*y^13*z^7+72*y^12*z^8+428652*y^10*z^10-18522*y^16*z^3-1284*y^15*z^4+1200*y^14*z^5+166698*y^13*z^6-35496*y^12*z^7+1360800*y^9*z^10-44100*y^15*z^3+8720*y^14*z^4+240*y^13*z^5+396900*y^12*z^6-7452*y^11*z^7+544320*y^8*z^10+3306744*y^5*z^13-125640*y^14*z^3+27848*y^13*z^4+16*y^12*z^5-283986*y^11*z^6+191376*y^10*z^7+904932*y^8*z^9-732888*y^7*z^10-166698*y^14*z^2+72540*y^13*z^3-7888*y^12*z^4+500094*y^11*z^5-5251392*y^10*z^6+604800*y^9*z^7+7831404*y^7*z^9-157464*y^6*z^10+8503056*z^16-396900*y^13*z^2+95652*y^12*z^3-1656*y^11*z^4-2309958*y^10*z^5-2416068*y^9*z^6+241920*y^8*z^7+8168202*y^7*z^8-17125668*y^6*z^9+2204496*y^5*z^10+35639352*y^3*z^12-144666*y^12*z^2-63972*y^11*z^3+21360*y^10*z^4-9125136*y^9*z^5+1878012*y^8*z^6-325728*y^7*z^7+18305028*y^6*z^8-6234408*y^5*z^9+5143824*y^2*z^12-500094*y^12*z+4001724*y^11*z^2-1166976*y^10*z^3+567294*y^9*z^4-32633928*y^8*z^5+4004856*y^7*z^6-69984*y^6*z^7-7978824*y^5*z^8-944784*y^4*z^9-110539728*y*z^12+7558272*z^13+2309958*y^11*z+1358424*y^10*z^2-536904*y^9*z^3-25620798*y^8*z^4+40189176*y^7*z^5-7340220*y^6*z^6+489888*y^5*z^7-248391360*y^4*z^8+23759568*y^3*z^9-11337408*z^12+8220204*y^10*z-2041740*y^9*z^2+372648*y^8*z^3-61088646*y^7*z^4+25303032*y^6*z^5-2770848*y^5*z^6-67831992*y^3*z^8+3429216*y^2*z^9-500094*y^10+26413938*y^9*z-7181352*y^8*z^2+503232*y^7*z^3-4589460*y^6*z^4-1931652*y^5*z^5-419904*y^4*z^6+462113100*y^2*z^8-73693152*y*z^9+2519424*z^10+17479476*y^9-19764324*y^8*z+8343852*y^7*z^2-785448*y^6*z^3+465676128*y^5*z^4-108769032*y^4*z^5+5279904*y^3*z^6+97312752*y*z^8-7558272*z^9+53249400*y^8-15635628*y^7*z+4724776*y^6*z^2-307872*y^5*z^3+203411412*y^4*z^4-30147552*y^3*z^5+762048*y^2*z^6+31177872*z^8+42066054*y^7-6575076*y^6*z-35240*y^5*z^2-46656*y^4*z^3-647949780*y^3*z^4+205383600*y^2*z^5-16376256*y*z^6+373248*z^7-216860652*y^6+100760724*y^5*z-11904656*y^4*z^2+391104*y^3*z^3-214603020*y^2*z^4+43250112*y*z^5-1679616*z^6-183344958*y^5+45237528*y^4*z-3349728*y^3*z^2+56448*y^2*z^3-193838184*y*z^4+13856832*z^5+259748532*y^4-143988840*y^3*z+22820400*y^2*z^2-1213056*y*z^3-18245088*z^4+134812512*y^3-47689560*y^2*z+4805568*y*z^2-124416*z^3+189980316*y^2-43075152*y*z+1539648*z^2+38053800*y-4059072*z+22071204
```
