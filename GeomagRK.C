//
//  main.cpp
//  Projeto2parte
// 
//  Created by Joao Ferreira on 21/12/14.
//  Copyright (c) 2014 Joao Ferreira & Filipe Matos. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>
#include <limits>
#include <iomanip>
#include <TStyle.h>
#include <TPolyLine3D.h>
#include <TSPHE.h>
#include <TMaterial.h>
#include <TObjArray.h>

#include "cFCgraphics.h"
#include "Trajectory.h"
#include "MagField.h"

using namespace std;
int main(int argc, const char * argv[]) {

    Particle* pr = new Particle("p+");
    Particle* pi = new Particle("pi+");
   
    vector<double> vel1= pr->Mom_Conversion(0,-0.1,45,0);
    Trajectory* tr1 = new Trajectory(pr,-0.000001, 0.7626,0, 0.7626, vel1[0], vel1[1], vel1[2], 3000000);
    TPolyLine3D *g1 = new TPolyLine3D();
    g1= tr1->GetLine(200);
    g1->SetLineColor(2);
    g1->SetLineWidth(4);
    
    vector<double> vel2= pr->Mom_Conversion(0,-10,55,0);
    Trajectory* tr2 = new Trajectory(pr,-0.000001, 0.7626,0,0.7626, vel2[0], vel2[1],vel2[2], 3000001);
    TPolyLine3D *g2 = new TPolyLine3D();
    g2= tr2->GetLine(200);
    g2->SetLineColor(3);
    g2->SetLineWidth(4);
   /*
    Trajectory* tr3 = new Trajectory(pr,-0.000001, 0 ,-0.7626, -0.7626, vel1[0], vel1[1], vel1[2], 3000000);
    TPolyLine3D *g3 = new TPolyLine3D();
    g3= tr3->GetLine(2000);
    g3->SetLineColor(5);
    g3->SetLineWidth(4);
    */
    MagField* mf1 = new MagField();
    TPolyLine3D *tp;
    
    cFCgraphics can;
    TPad *pad = can.CreatePad("Pad1");
    gStyle->SetPalette(1,0);
    
    TSPHE *sph = new TSPHE("Esfera","algo","void",1);
    
    can.AddObject(sph, "Pad1");
    
    can.AddObject(g1, "Pad1","same");
    can.AddObject(g2, "Pad1","same");
//    can.AddObject(g3, "Pad1","same");
    
    for (int i=0; i!=5; ++i) {
        tp = mf1->GetLShell_Pol(3,-0.3*M_PI*i);
        tp->SetLineColor(9);
        tp->SetLineWidth(1);
        can.AddObject(tp, "Pad1","same");
    }
    for (int i=0; i!=5; ++i) {
        tp = mf1->GetLShell_Pol(6,-0.3*M_PI*i);
        tp->SetLineColor(51);
        tp->SetLineWidth(1);
        can.AddObject(tp, "Pad1","same");
    }
    can.DumpPad("Pad1");
    can.AddObject(pad);
    can.Draw();
   
    return 0;
}
