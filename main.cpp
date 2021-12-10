#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

double alpha_k(double k, double l);
double r1k(double k, double l, double B, double L, double sigma, double mu);
double r2k(double k, double l, double B, double L, double sigma, double mu);
double N(double k, double l, double B, double L, double sigma, double mu);

double V2(double k, double l, double B, double L, double sigma, double mu, double xi, double eta);
double V3(double k, double l, double B, double L, double sigma, double mu, double xi, double eta);

double past_sum(double k, double l, double B, double L, double sigma, double mu, double xi, double eta);
double sum(double l, double B, double L, double sigma, double mu, double xi, double eta);

int main()
{
    ofstream myfile;

    myfile.open ("MHD_Analytical.txt");


double B =0.00052552761542017795;
double L =1.0;
double mu =1.2566156992899439E-06;//Permeability
double sigma =4.55e6;//Electrical conductivity

double Ha = B*L*sqrt(sigma/mu);

cout<<"Hartmann Number = "<<Ha<<endl;

double pressure_gradient = 376.84/20.0;

double kinematic_viscosity = 1.322401481089659E-07;
double density = 11343;//

double asqr=1.0*1.0;

for(double dist=0.0;dist<=1.01;dist+=0.01){

//    cout<<dist<<"\t"<<(1.0/(kinematic_viscosity*density))*pressure_gradient*asqr*sum(1.0,B,L,sigma,mu,dist,0.0)<<endl;
    myfile<<dist<<"\t"<<(1.0/(kinematic_viscosity*density))*pressure_gradient*asqr*sum(1.0,B,L,sigma,mu,dist,0.0)<<endl;
}





myfile.close();
    return 0;
}

double alpha_k(double k, double l){

double alpha_k_out = (k+0.5)*(M_PI/l);


return alpha_k_out;
}

double r1k(double k, double l, double B, double L, double sigma, double mu){

double Ha = B*L*sqrt(sigma/mu);

//double r1k_out = 0.5*(Ha+((Ha*Ha)+(4.0*pow(alpha_k(k,l),2.0))));
double r1k_out = 0.5*(Ha+N(k,l,B,L,sigma,mu));


return r1k_out;
}

double r2k(double k, double l, double B, double L, double sigma, double mu){

double Ha = B*L*sqrt(sigma/mu);

//double r2k_out = 0.5*(-Ha+((Ha*Ha)+(4.0*pow(alpha_k(k,l),2.0))));
double r2k_out = 0.5*(-Ha+N(k,l,B,L,sigma,mu));


return r2k_out;
}

double N(double k, double l, double B, double L, double sigma, double mu){

double Ha = B*L*sqrt(sigma/mu);

double N_out = sqrt((Ha*Ha)+(4.0*(pow(alpha_k(k,l),2.0))));//0.5*(-Ha+((Ha*Ha)+(4.0*alpha_k(k,l))));


return N_out;
}

double V2(double k, double l, double B, double L, double sigma, double mu, double xi, double eta){

    double db=1e99;

    double N_=N(k,l,B,L,sigma,mu);
    double r1k_=r1k(k,l,B,L,sigma,mu);
    double r2k_=r2k(k,l,B,L,sigma,mu);

    double numerator=((db*r2k_)+((1.0-exp(-2.0*r2k_))/(1.0+exp(-2.0*r2k_))))*(0.5*(exp(-r1k_*(1.0-eta))+exp(-r1k_*(1.0+eta))));
    double denominator=(((1.0+exp(-2.0*r1k_))/2.0)*(db*N_))+((1.0+exp(-2.0*(r1k_+r2k_)))/(1.0+exp(-2.0*r2k_)));

double V2_out=(numerator)/(denominator);

return V2_out;
}

double V3(double k, double l, double B, double L, double sigma, double mu, double xi, double eta){

    double db=1e99;

    double N_=N(k,l,B,L,sigma,mu);
    double r1k_=r1k(k,l,B,L,sigma,mu);
    double r2k_=r2k(k,l,B,L,sigma,mu);

    double numerator=((db*r1k_)+((1.0-exp(-2.0*r1k_))/(1.0+exp(-2.0*r1k_))))*(0.5*(exp(-r2k_*(1.0-eta))+exp(-r2k_*(1.0+eta))));
    double denominator=(((1.0+exp(-2.0*r2k_))/2.0)*(db*N_))+((1.0+exp(-2.0*(r1k_+r2k_)))/(1.0+exp(-2.0*r1k_)));

double V3_out=(numerator)/(denominator);

return V3_out;
}

double past_sum(double k, double l, double B, double L, double sigma, double mu, double xi, double eta){

double alpha_k_=alpha_k(k,l);

double V2_ = V2(k,l,B,L,sigma,mu,xi,eta);
double V3_ = V3(k,l,B,L,sigma,mu,xi,eta);

double value = ((2.0*(pow(-1.0,k))*cos(alpha_k_*xi))/(l*pow(alpha_k_,3.0)))*(1.0-V2_-V3_);

return value;

}

double sum(double l, double B, double L, double sigma, double mu, double xi, double eta){



double valuesum = 0.0;

for(double k=0.0;k<50000.0;k+=1.0){

    valuesum+=past_sum(k,l,B,L,sigma,mu,xi,eta);

}

return valuesum;

}
