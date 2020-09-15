#include <clpoly.hh>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
int main(int argc, char const *argv[])
{
    clpoly::variable x("x");
    clpoly::variable y("y");
    clpoly::variable z("z");
    clpoly::polynomial_ZZ f1,f2,f3,f4,f,g,G;

    
    time_t t;
    double s1=0,s2=0;
    for (int i=0;i<100;++i)
    {
        //std::cout<<"test "<<i<<":\n";
        f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.02,10,-10);
        g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.02,10,-10);
        std::cout<<"f["<<i<<"]="<<f<<";"<<std::endl;
        std::cout<<"g["<<i<<"]="<<g<<";"<<std::endl; 
        t=clock();
        f1=resultant(f,g,x);
        s1+=double(clock()-t)/CLOCKS_PER_SEC;
        std::cout<<"o["<<i<<"]="<<f1<<";"<<std::endl; 
        std::cout<<"t["<<i<<"]="<<double(clock()-t)/CLOCKS_PER_SEC <<";"<<std::endl; 
        //std::cout<<"resultant    time:"<<double(clock()-t)/CLOCKS_PER_SEC<<std::endl;
        // t=clock();
        // f2=resultant_v1(f,g,x);
        // s2+=double(clock()-t)/CLOCKS_PER_SEC;
        // std::cout<<"resultant_v1 time:"<<double(clock()-t)/CLOCKS_PER_SEC<<std::endl;
        // if (f1==f2)
        //     std::cout<<"一致\n";
        // else
        // {
        //     std::cout<<"不一致\n";
        //     std::cout<<"f1="<<f1<<std::endl;
        //     std::cout<<"f2="<<f2<<std::endl;
             
             
        // }   
        // assert(f1==f2);
        
    }
    std::cout<<"(*resultant    total time:"<<s1<<"*)\n";
    // std::cout<<"resultant_v1 total time:"<<s2<<std::endl;
    return 0;

       // f=-pow(x,2)*pow(z,3) - pow(x,4) - pow(z,4) + pow(x,2) + 2*pow(z,2) - 1;
    // g=-pow(r,2)*pow(x,2) + pow(x,4) + pow(x,2)*pow(z,2) + pow(z,4) - 2*pow(z,2) + 1;
    //g=-4*pow(x,10)+8*pow(x,9)*z-7*pow(x,8)*y*z-4*pow(x,7)*pow(z,3)+10*pow(x,4)*pow(y,2)*pow(z,4)-pow(x,4)*pow(z,6)-7*pow(x,3)*pow(y,6)*z+2*pow(x,3)*y*pow(z,6)-2*x*pow(y,4)*pow(z,5)-4*pow(y,8)*pow(z,2)-6*pow(y,7)*pow(z,3)+2*y*pow(z,9)-6*pow(x,7)*y*z+7*pow(x,5)*pow(y,3)*z-2*pow(x,3)*pow(y,2)*pow(z,4)-7*pow(x,2)*pow(y,5)*pow(z,2)-pow(x,2)*pow(y,2)*pow(z,5)-8*x*pow(y,4)*pow(z,4)-8*x*y*pow(z,7)-10*pow(y,8)*z+8*pow(x,8)+5*pow(x,7)*y+6*pow(x,7)*z-5*pow(x,5)*pow(y,3)-5*pow(x,4)*pow(y,4)+9*pow(x,3)*pow(y,3)*pow(z,2)-9*x*pow(y,6)*z+x*pow(y,4)*pow(z,3)-7*y*pow(z,7)-10*pow(x,5)*pow(z,2)+9*pow(x,3)*pow(y,4)+3*pow(x,3)*pow(y,2)*pow(z,2)+6*pow(x,2)*pow(y,3)*pow(z,2)-8*x*pow(y,4)*pow(z,2)+7*pow(y,2)*pow(z,5)-4*y*pow(z,6)+8*pow(x,4)*pow(y,2)-pow(x,4)*y*z-pow(x,2)*y*pow(z,3)+8*pow(x,2)*pow(z,4)-2*x*pow(y,5)-8*pow(y,6)+2*y*pow(z,5)+pow(x,2)*pow(y,2)*z+pow(x,2)*y*pow(z,2)+4*x*y*pow(z,3)+3*pow(y,4)*z-9*pow(x,4)+7*pow(x,3)*z+3*x*pow(z,3)+4*pow(y,2)*pow(z,2)-4*pow(z,4)-6*pow(x,2)*y+8*x*pow(y,2)+7*x*y*z-8*pow(y,3)-2*x*y-3*x*z-8*pow(z,2)+1;
    //f=-8*pow(x,9)*z-pow(x,8)*pow(y,2)-2*pow(x,6)*pow(y,4)-6*pow(x,5)*pow(z,5)+2*pow(x,3)*pow(y,6)*z-3*pow(x,3)*pow(y,3)*pow(z,4)+2*pow(x,2)*pow(y,4)*pow(z,4)+9*pow(x,2)*pow(z,8)+3*x*pow(y,9)+3*x*pow(y,7)*pow(z,2)+5*pow(y,10)-5*pow(x,9)+5*pow(x,7)*pow(z,2)-6*pow(x,5)*pow(z,4)-3*pow(x,3)*pow(y,6)+pow(x,3)*pow(y,4)*pow(z,2)-2*pow(x,2)*pow(y,3)*pow(z,4)-4*pow(x,2)*pow(y,2)*pow(z,5)-8*pow(x,2)*y*pow(z,6)-3*x*pow(y,6)*pow(z,2)-9*x*pow(y,3)*pow(z,5)+2*pow(x,3)*pow(y,5)-pow(x,2)*pow(y,6)+4*x*pow(y,5)*pow(z,2)+2*pow(y,4)*pow(z,4)-3*pow(x,2)*pow(y,5)+pow(x,2)*pow(y,4)*z-8*pow(x,2)*pow(y,2)*pow(z,3)+8*pow(x,2)*pow(z,5)+5*x*pow(z,6)+7*pow(x,4)*pow(z,2)+5*pow(x,3)*y*pow(z,2)-7*pow(x,2)*pow(y,2)*pow(z,2)-5*x*pow(y,3)*pow(z,2)-4*x*y*pow(z,4)+3*pow(x,4)*y-3*pow(x,4)*z+pow(x,3)*y*z+6*pow(x,3)*pow(z,2)+9*pow(x,2)*pow(y,2)*z-2*x*pow(y,4)-9*x*pow(z,4)-6*pow(x,2)*y*z+5*pow(x,2)*pow(z,2)+9*x*y*pow(z,2)-8*pow(y,3)*z-pow(y,2)*pow(z,2)-9*x*pow(z,2)+y*z+z;
    //g=-6*pow(x,10)+5*pow(x,9)*y+2*pow(x,4)*pow(y,2)*pow(z,4)+10*pow(x,3)*pow(y,6)*z+2*pow(x,3)*pow(z,7)-5*pow(x,2)*pow(y,5)*pow(z,3)-7*x*pow(y,6)*pow(z,3)+7*pow(y,10)+9*pow(y,8)*pow(z,2)-3*pow(y,6)*pow(z,4)+8*pow(y,2)*pow(z,8)+3*pow(x,8)*y-2*pow(x,3)*pow(y,6)+5*pow(x,3)*pow(y,4)*pow(z,2)-3*pow(x,2)*pow(y,7)+3*x*pow(y,8)-2*x*pow(y,3)*pow(z,5)+8*pow(y,8)*z+pow(x,7)*y+3*pow(x,4)*pow(y,4)+9*pow(x,4)*pow(y,2)*pow(z,2)+6*pow(x,3)*pow(y,4)*z-7*pow(x,2)*pow(y,2)*pow(z,4)+5*x*pow(y,7)+4*x*pow(y,3)*pow(z,4)+9*pow(y,7)*z+3*pow(y,6)*pow(z,2)-4*pow(y,3)*pow(z,5)-4*pow(x,5)*pow(y,2)-6*pow(x,3)*pow(y,3)*z+6*pow(x,2)*pow(y,3)*pow(z,2)-pow(x,2)*pow(y,2)*pow(z,3)+7*x*pow(y,4)*pow(z,2)-9*x*y*pow(z,5)+7*pow(y,5)*pow(z,2)+pow(x,2)*pow(z,4)-6*pow(x,3)*pow(y,2)+2*pow(x,2)*y*pow(z,2)-10*x*pow(y,3)*z-5*pow(y,2)*pow(z,3)-6*pow(z,5)-10*pow(x,2)*pow(z,2)-4*x*y*pow(z,2)-3*x*pow(z,3)-10*pow(z,4)-7*pow(x,3)-5*pow(x,2)*y-3*pow(y,3)+y*pow(z,2)+4*pow(z,3)-7*x*y+x*z-9;
    //f=3*pow(x,9)*z+4*pow(x,7)*pow(y,2)*z+6*pow(x,6)*pow(y,4)-8*pow(x,4)*pow(y,4)*pow(z,2)-8*pow(x,3)*y*pow(z,6)+pow(x,2)*pow(y,6)*pow(z,2)+9*x*pow(y,9)-x*pow(y,8)*z+4*x*pow(y,5)*pow(z,4)+6*pow(x,3)*pow(y,4)*pow(z,2)-9*pow(x,3)*pow(y,2)*pow(z,4)-3*pow(x,3)*y*pow(z,5)-7*pow(x,3)*pow(z,6)-4*pow(x,2)*pow(y,7)+9*pow(x,2)*pow(y,6)*z-4*x*pow(y,7)*z-2*pow(y,6)*pow(z,3)+6*pow(y,3)*pow(z,6)-5*pow(y,2)*pow(z,7)+7*pow(x,4)*y*pow(z,3)-pow(x,3)*y*pow(z,4)-4*pow(x,2)*pow(y,2)*pow(z,4)+4*x*pow(y,6)*z+x*pow(y,5)*pow(z,2)-10*x*y*pow(z,6)-3*pow(y,7)*z+2*pow(x,6)*y+2*pow(x,4)*y*pow(z,2)-7*pow(x,2)*y*pow(z,4)+3*x*pow(y,2)*pow(z,4)+2*x*pow(z,6)+8*pow(y,6)*z+pow(y,5)*pow(z,2)-4*pow(x,4)*pow(z,2)+pow(y,5)*z+5*pow(y,3)*pow(z,3)-10*pow(y,2)*pow(z,4)-8*y*pow(z,5)-4*pow(x,3)*y*z+7*x*pow(y,4)-2*pow(y,4)*z-3*pow(y,3)*pow(z,2)-pow(y,2)*pow(z,3)-7*x*y*pow(z,2)+10*pow(y,3)*z-5*pow(z,4)+6*x*pow(y,2)+6*pow(z,3)-3*pow(z,2)-10*y+2*z-7;
    // g=4*x*pow(y,3)*z-4*x*pow(y,2)*pow(z,2)-8*x*pow(z,4)-8*pow(y,4)*z-3*pow(x,2)*y*z+2*pow(x,2)*pow(z,2)-8*x*pow(y,3)-x*pow(z,3)+2*pow(x,2)*z-5*x*pow(y,2)-8*x*pow(z,2)+9*y*pow(z,2)-4*pow(z,3)+5*x*y+3;
    // f=7*pow(x,2)*pow(y,3)+x*pow(y,4)-8*x*pow(z,4)+9*pow(y,5)-9*y*pow(z,4)+2*pow(x,4)-7*pow(x,2)*pow(y,2)-7*x*pow(y,3)+6*pow(y,4)-9*pow(y,2)*pow(z,2)-3*y*pow(z,3)+2*pow(z,4)+pow(x,2)*z+5*x*pow(y,2)-3*x*pow(z,2)-6*pow(y,3)+9*y*z+4;
    // f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.2,10,-10);
    // g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.2,10,-10);
    //g=-6*pow(x,3)*pow(z,2)+4*pow(x,2)*pow(z,3)-7*pow(z,5)+pow(x,3)*y+10*x*pow(y,3)+3*pow(y,4)+y*pow(z,3)+2*pow(z,2)+x-3*z-3;
    //f=-6*pow(x,5)-10*pow(x,3)*pow(y,2)-9*pow(x,3)*y*z-8*pow(x,3)*pow(z,2)+4*pow(x,2)*pow(y,3)+5*pow(x,2)*y*pow(z,2)-7*x*y*pow(z,3)-3*pow(y,5)-9*pow(x,2)*pow(z,2)+9*y*pow(z,2)-3*pow(y,2)-7;
    // g=5*pow(x,5)+9*pow(x,4)*y+pow(x,4)*z-10*pow(x,3)*y*z+10*pow(x,3)*pow(z,2)-8*pow(x,2)*y*pow(z,2)-9*x*pow(z,4)+2*pow(y,3)*pow(z,2)-5*pow(x,4)+7*pow(x,3)*y-pow(x,2)*z+6*x*pow(z,2)+6*pow(y,3)-8;
    // f=-2*pow(x,3)*pow(z,2)+7*pow(x,2)*pow(y,2)*z+x*y*pow(z,3)+2*pow(y,5)-10*pow(z,5)-3*pow(x,3)*z-3*pow(y,3)*z+6*pow(x,3)-5*pow(y,3)-3*z+1;
    // g=5*pow(x,5)*pow(y,5)+7*pow(x,6)*pow(y,2)+8*pow(x,5)*pow(y,3)-6*pow(y,8)-2*pow(x,4)*y+3*x*pow(y,2)-6*pow(y,3)-3;
    // f=-4*pow(y,6)+7*pow(x,4)*y-10*pow(x,4)-4;
    // g=-5*pow(x,7)*y-3*pow(x,6)*y-4*pow(x,6)+9*x*pow(y,3);
    // f=3*pow(x,6)*pow(y,4)+pow(x,6)*pow(y,2)+3*pow(x,5)*pow(y,2)+7;
      // f=x*y*pow(z,5)-pow(d,6)*r;
    // g=pow(x,5)*pow(y,3)*pow(z,17)-pow(d,22)*pow(r,3);
    // std::cout<<"g="<<g<<";"<<std::endl; 
    // std::cout<<"f="<<f<<";"<<std::endl;
    // std::cout<<clpoly::resultant(f,g,z)<<std::endl;

}
