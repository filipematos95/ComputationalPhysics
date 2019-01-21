
#ifndef __Particle__
#define __Particle__

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <limits>
#include <iomanip>
#include <TF1.h> 
#include <TPolyLine3D.h>

#include "cFCgraphics.h"


using namespace std;

class Particle {
public:
    Particle();
    Particle(string s);
    Particle(string s, double c, double m);
    vector<double> Mom_Conversion (int, double,double,double); 
    
    string name;
    double charge;
    double mass;
};

class Trajectory {
public:
    Trajectory(Particle *);
    Trajectory(Particle*, double, double, double, double, double, double, double,int);
    TPolyLine3D* GetLine(int);
    
private:

    double acex(double, double, double, double, double, double);
    double acey(double, double, double, double, double, double);
    double acez(double, double, double, double, double, double);

    double CT;
    double h;
    int f;
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> vx;
    vector<double> vy;
    vector<double> vz;
    TPolyLine3D *traj;
};

#endif