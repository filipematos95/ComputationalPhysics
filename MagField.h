
#ifndef __MagField__
#define __MagField__

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


class MagField {
public:
    MagField();
    MagField(double b) {B0 = b;}
    TPolyLine3D* GetLine_Num_Cart(double, double, double, double, int);
    TPolyLine3D* GetLShell_Cart(double, double);
    TPolyLine3D* GetLShell_Pol(double, double);
    
private:
    double B0;
};

#endif