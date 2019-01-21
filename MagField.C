#include "MagField.h" 

using namespace std;

MagField::MagField() {
    B0= 3.07e-5;
}

double fz(double x,double z) {
    double f = 0.;

    double s = sqrt(5.*z*z*x*x + 4.*z*z*z*z + x*x*x*x);
    f = (2.*z*z-x*x)/s;
    return f;
}

double fx(double x,double z) {
    double f = 0.;

    double s = sqrt(5.*z*z*x*x + 4.*z*z*z*z + x*x*x*x);
    f = (3.*z*x)/s;
    
    return f;
}

TPolyLine3D* MagField::GetLine_Num_Cart(double s, double x0, double y0, double z0, int f) {
    
    vector<double> z;
    vector<double> y;
    vector<double> x;
    x.reserve(f);
    y.reserve(f);
    z.reserve(f);
    
    double k1,k2,k3,k4;
    
    x.push_back(sqrt(x0*x0+y0*y0));
    z.push_back(z0);
    
    double h=s;
    double al= atan2(y0,x0);
    
    for (int i = 1; i != f; ++i) {
        k1 = h*fx(x[i-1], z[i-1]);
        k2 = h*fx(x[i-1]+k1/2., z[i-1] + k1/2.);
        k3 = h*fx(x[i-1]+k2/2., z[i-1] + k2/2.);
        k4 = h*fx(x[i-1]+k3, z[i-1]+k3);
        
        x.push_back(x[i-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
        
        k1 = h*fz(x[i-1], z[i-1]);
        k2 = h*fz(x[i-1]+k1/2., z[i-1] + k1/2.);
        k3 = h*fz(x[i-1]+k2/2., z[i-1] + k2/2.);
        k4 = h*fz(x[i-1]+k3, z[i-1]+k3);
        
        z.push_back(z[i-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
        
        
        if (fabs(z[i]-z[i-1])<1e-10) {
            h = -s;
            
            x[i]=sqrt(x0*x0+y0*y0);
            z[i]=z0;
            
            for (int j = i+1; j != f; ++j) {
                
                k1 = h*fx(x[j-1], z[j-1]);
                k2 = h*fx(x[j-1]+k1/2., z[j-1] + k1/2.);
                k3 = h*fx(x[j-1]+k2/2., z[j-1] + k2/2.);
                k4 = h*fx(x[j-1]+k3, z[j-1]+k3);
                
                x.push_back(x[j-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
                
                k1 = h*fz(x[j-1], z[j-1]);
                k2 = h*fz(x[j-1]+k1/2., z[j-1] + k1/2.);
                k3 = h*fz(x[j-1]+k2/2., z[j-1] + k2/2.);
                k4 = h*fz(x[j-1]+k3, z[j-1]+k3);
                
                z.push_back(z[j-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
                
                if (fabs(z[j]-z[j-1])<1e-10) {
                    double sn= sin(al);
                    double csn= cos(al);
                    
                    TPolyLine3D* line = new TPolyLine3D();
                    for (int g=0; g != f; ++g){
                        line->SetPoint(g,x[g]*csn,x[g]*sn,z[g]);
                    }
        
                    return line;
                }
            }
            double sn= sin(al);
            double csn= cos(al);
            
            TPolyLine3D* line = new TPolyLine3D();
            for (int g=0; g != f; ++g){
                line->SetPoint(g,x[g]*csn,x[g]*sn,z[g]);
            }
            
            return line;
        }
    }
    
    double sn= sin(al);
    double csn= cos(al);
    
    TPolyLine3D* line = new TPolyLine3D();
    for (int g=0; g != f; ++g){
        line->SetPoint(g,x[g]*csn,x[g]*sn,z[g]);
    }
    
    return line;
}


TPolyLine3D* MagField::GetLShell_Pol(double ls, double ph) {
    TPolyLine3D* line = new TPolyLine3D();
    double l,r,c=M_PI/50.;
    for (int i=0; i <=50;++i){
        l=-0.5*M_PI-c*i;
        r=ls*cos(l)*cos(l);
        line->SetPoint(i,r*cos(l)*cos(ph),r*cos(l)*sin(ph),r*sin(l));
    }
    return line;
}

TPolyLine3D* MagField::GetLShell_Cart(double x, double y) {
    double ls = sqrt(x*x+y*y);
    double ph = atan2(y,x);
    
    return GetLShell_Pol(ls,ph);
}





