#include "Trajectory.h"

#define C 2.99792458e8
#define B 3.07e-5
#define RE 6371000.0

using namespace std;

Particle::Particle(string s, double c, double m) {
    name = s;
    charge =c;
    mass = m;
}
Particle::Particle() {
    name = "000";
    charge =1;
    mass = 1;
}
Particle::Particle(string s) {
    if (s=="p+") {
        name = "p+";
        charge = 1.60217653e-19;
        mass = 1.672621777e-27;
    }
    else if (s=="pi+"){
        name = "pi+";
        charge =1.60217653e-19;
        mass = 2.488064e-28;
    }
    else {
        name = "000";
        charge =1;
        mass = 1;
    }
}
vector<double> Particle::Mom_Conversion (int u, double p1,double p2,double p3) {
    double v, vx, vy, vz, inv, t2, t3;
    
    switch (u) {
        case 0:  //Unit: Gev c-1   Polar(deg) to Cartesian
            inv = 5.34428576e-19;
            v = sqrt(1/((mass*mass)/(p1*p1*inv*inv)+1/(C*C)))/C;
      
            if (p1<0) {
                t2 = p2*M_PI/180.;
                t3 = p3*M_PI/180.;
                vx = -v*sin(t2)*cos(t3);
                vy = -v*sin(t2)*sin(t3);
                vz = -v*cos(t2);
            }
            else {
                t2 = p2*M_PI/180.;
                t3 = p3*M_PI/180.;
                vx = v*sin(t2)*cos(t3);
                vy = v*sin(t2)*sin(t3);
                vz = v*cos(t2);
            }
            break;
        default:
            break;
    }
    
    vector<double> res;
    res.push_back(vz);
    res.push_back(vy);
    res.push_back(vx);
    
    return res;
}



Trajectory::Trajectory(Particle* p) {
    CT = -(p->charge)*B/((p->mass));
    h=1;
    f=100;
    traj = new TPolyLine3D();
}
double Trajectory::acex(double x, double y, double z, double vx, double vy, double vz) {
    
    double r=sqrt(x*x + y*y + z*z);
    double v2=(vx*vx + vy*vy + vz*vz);
    double c=CT*sqrt(1-v2)/(r*r*r*r*r);
    
    double a= ((2.*z*z-x*x-y*y)*vy - 3.*y*z*vz);
    return a*c;
}
double Trajectory::acey(double x, double y, double z, double vx, double vy, double vz) {
    
    double r=sqrt(x*x + y*y + z*z);
    double v2=(vx*vx + vy*vy + vz*vz);
    double c=CT*sqrt(1-v2)/(r*r*r*r*r);
    
    double a= (3.*x*z*vz - (2.*z*z-x*x-y*y)*vx);
    return a*c;
}
double Trajectory::acez(double x, double y, double z, double vx, double vy, double vz) {
    
    double r=sqrt(x*x + y*y + z*z);
    double v2=(vx*vx + vy*vy + vz*vz);
    double c=CT*sqrt(1-v2)/(r*r*r*r*r);
    
    double a= (3.*z*(y*vx-x*vy));
    return a*c;
}
Trajectory::Trajectory(Particle* p,double s, double x0, double y0, double z0, double vx0, double vy0, double vz0,  int fi) {
    
    CT = -(p->charge)*B/((p->mass));
    h = s;
    f = fi;
    x.reserve(f);
    y.reserve(f);
    z.reserve(f);
    vx.reserve(f);
    vy.reserve(f);
    vz.reserve(f);
    double k1,k2,k3,k4;
    
    x.push_back(x0);
    y.push_back(y0);
    z.push_back(z0);
    
    vx.push_back(vx0);
    vy.push_back(vy0);
    vz.push_back(vz0);
    
    double r;
    
    for (int i = 1; i != f; ++i) {
        
        k1 = h*vx[i-1]*C/RE;
        k2 = h*(vx[i-1]+k1/2.)*C/RE;
        k3 = h*(vx[i-1]+k2/2.)*C/RE;
        k4 = h*(vx[i-1]+k3)*C/RE;
        
        x.push_back(x[i-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
        
        k1 = h*vy[i-1]*C/RE;
        k2 = h*(vy[i-1]+k1/2.)*C/RE;
        k3 = h*(vy[i-1]+k2/2.)*C/RE;
        k4 = h*(vy[i-1]+k3)*C/RE;
        
        y.push_back(y[i-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
        
        k1 = h*vz[i-1]*C/RE;
        k2 = h*(vz[i-1]+k1/2.)*C/RE;
        k3 = h*(vz[i-1]+k2/2.)*C/RE;
        k4 = h*(vz[i-1]+k3)*C/RE;
        
        z.push_back(z[i-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
        
        k1 = h*acex(x[i-1], y[i-1], z[i-1], vx[i-1], vy[i-1], vz[i-1]);
        k2 = h*acex(x[i-1]+k1/2., y[i-1]+k1/2., z[i-1]+k1/2., vx[i-1]+k1/2., vy[i-1]+k1/2., vz[i-1]+k1/2.);
        k3 = h*acex(x[i-1]+k2/2., y[i-1]+k2/2., z[i-1]+k2/2., vx[i-1]+k2/2., vy[i-1]+k2/2., vz[i-1]+k2/2.);
        k4 = h*acex(x[i-1]+k3, y[i-1]+k3, z[i-1]+k3, vx[i-1]+k3, vy[i-1]+k3, vz[i-1]+k3);
        
        vx.push_back(vx[i-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
        
        k1 = h*acey(x[i-1], y[i-1], z[i-1], vx[i-1], vy[i-1], vz[i-1]);
        k2 = h*acey(x[i-1]+k1/2., y[i-1]+k1/2., z[i-1]+k1/2., vx[i-1]+k1/2., vy[i-1]+k1/2., vz[i-1]+k1/2.);
        k3 = h*acey(x[i-1]+k2/2., y[i-1]+k2/2., z[i-1]+k2/2., vx[i-1]+k2/2., vy[i-1]+k2/2., vz[i-1]+k2/2.);
        k4 = h*acey(x[i-1]+k3, y[i-1]+k3, z[i-1]+k3, vx[i-1]+k3, vy[i-1]+k3, vz[i-1]+k3);
        
        vy.push_back(vy[i-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
        
        k1 = h*acez(x[i-1], y[i-1], z[i-1], vx[i-1], vy[i-1], vz[i-1]);
        k2 = h*acez(x[i-1]+k1/2., y[i-1]+k1/2., z[i-1]+k1/2., vx[i-1]+k1/2., vy[i-1]+k1/2., vz[i-1]+k1/2.);
        k3 = h*acez(x[i-1]+k2/2., y[i-1]+k2/2., z[i-1]+k2/2., vx[i-1]+k2/2., vy[i-1]+k2/2., vz[i-1]+k2/2.);
        k4 = h*acez(x[i-1]+k3, y[i-1]+k3, z[i-1]+k3, vx[i-1]+k3, vy[i-1]+k3, vz[i-1]+k3);
        
        vz.push_back(vz[i-1]+(1/6.)*(k1+2.*k2+2.*k3+k4));
        
        r=x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
        if (r<1. || r>625.) {
            cout << "BUM!!!" << endl;
            break;
        }
    }
    
    traj = new TPolyLine3D();
}
TPolyLine3D* Trajectory::GetLine(int i) {
    int nP = x.size();
    double posx, posy, posz;
    for (int N=0; N<nP; N=N+i) {
        posx = x[N];
        posy = y[N];
        posz = z[N];
        traj->SetPoint(N/i,posx,posy,posz);
    }
    
    return traj;
}


