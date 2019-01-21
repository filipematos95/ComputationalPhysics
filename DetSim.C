
//
//  DetSim.C
//  Trab2_Ex1
//
//  Created by Joao Ferreira on 29/11/14.
//  Copyright (c) 2014 Joao Ferreira & Filipe Matos. All rights reserved.
//

#include <iostream>
#include "PixelDet.h"
#include "TH1F.h"
#include <string>
#include "cFCgraphics.h"
#include "TH2F.h"
#include <ctgmath>
#include <ctime>
#include <vector>

using namespace std;

int main(int argc, const char * argv[]) {
    
    // Alinea a)
    srand(time(NULL));
    
    PixelDet a(80);
    
    //TH1F *grafico = a.GetNoise(0.5,10000,"y");
    //gr.AddObject(grafico);
    
    cFCgraphics pr;
    double* pos = new double[2];
    pos[0] = 190.;
    pos[1] = 230.;
    
    double* pos2 = new double[2];
    double distancia;
    TH2C *h = a.DrawEvent(pos,0.85,0.5);
    
    
    
    a.Event_Rec();
    cout << "Ponto de impacto: " << a.Get_ponto_impacto() << endl;
    cout << "raio simulado: " << 40*sqrt(1.3*1.3*0.85*0.85-1) << endl;
    cout << "raio reconstruido: " << a.Get_raio_Rec() << endl;
    cout << "Velocidade reconsturuida: " << a.Get_veloc_Rec() << endl;
    
    h->SetOption("colz");
    pr.AddObject(h);
    pr.Draw();
    
    int j=0;
    TH1F *h_raios = new TH1F("hist_raio ","Raios; mm", 20,0,40);
    TH1F *h_velo = new TH1F("hist_velocidade ","Velocidades", 200,0,0.26);
    TH1F *h_ponto = new TH1F("hist_ponto ","Ponto de impacto", 3,0,3);
    TH1F *h_dist = new TH1F("hist_dist ","Distancia ao ponto de impacto;mm", 30,0,20);
    TH1F *h_tempo = new TH1F("hist_tempo ","Tempo gasto na reconstrucao de um aconteciemnto; s", 20,0,0.0005);
    TH1F *h_phi = new TH1F("hist_phi ","Angulo azimutal;rad", 40,0,7);
    TH1F *h_tracos = new TH1F("hist_tracos","Comprimento dos tracos;mm", 50,0,30);
    TH1F *h_noise = new TH1F("hist_noise","Distribuicao dos pixeis ruidosos;ID",1600,0,6400);
    TH1F *h_noise_100 = new TH1F("hist_noise_100","Distribuicao dos pixeis ruidosos 100 canais;ID",100,0,6400);
    
    pos2 = a.ID_CU(a.Get_ponto_impacto());
    
    clock_t t,tempo;
    vector<double> tracos;
    vector<double> angulos;
    double angulo_med=0;
    int angulo_it=0;
    vector<int> noise;
    
    for (  double i = 0.7; i < 1; i = i + 0.00003){
        pos[0] = (((double)rand())/RAND_MAX)*200+100;
        pos[1] = (((double)rand())/RAND_MAX)*200+100;
        
        a.SimulateEvent(pos,i,0.5);
        t = clock();
        a.Event_Rec();
        tempo = clock()-t;
        
        noise = a.EventNoise(0.5);
        for( int it = 0; it < noise.size() ; it++){
            h_noise->Fill(noise[it]);
            h_noise_100->Fill(noise[it]);
        }
        
        double raio_simu = 40*sqrt(1.3*1.3*i*i-1);
        int ponto_simu =  a.CU_ID((pos[0]),(pos[1]));
        
        pos2 = a.ID_CU(a.Get_ponto_impacto());
        
        
        distancia  = sqrt((pos[0]-pos2[0])*(pos[0]-pos2[0]) + (pos[1]-pos2[1])*(pos[1]-pos2[1]));
        
        tracos = a.Get_tracos();
        angulos = a.Get_angulo();
        
        for( int it = 0; it < tracos.size(); it++ ){
            h_tracos->Fill(tracos[it]);
            h_phi->Fill(angulos[i]);
            angulo_med += angulos[i];
            angulo_it++;
        }
        
        h_raios->Fill(fabs(raio_simu - a.Get_raio_Rec()));
        h_velo->Fill(fabs(i-a.Get_veloc_Rec()));
        h_ponto->Fill(abs(a.Get_ponto_impacto()-ponto_simu));
        h_dist->Fill(distancia);
        h_tempo->Fill(((double)tempo)/CLOCKS_PER_SEC);
        
        j++;
        
        a.tracos_clear();
        a.angulo_clear();
        
    }
    
    cFCgraphics gr;
    cout << "Ângulo médio: " << angulo_med/angulo_it << endl;
    gr.AddObject(h_ponto);
    gr.AddObject(h_raios);
    gr.AddObject(h_velo);
    gr.AddObject(h_dist);
    gr.AddObject(h_tempo);
    gr.AddObject(h_phi);
    gr.AddObject(h_tracos);
    gr.AddObject(h_noise);
    gr.AddObject(h_noise_100);
    
    gr.Draw();
    
    cFCgraphics pdf;
    
    pdf.AddObject(h_velo);
    
    pdf.Draw();
  
  
    cFCgraphics pdf2;
  
    pdf2.AddObject(h_dist);
  
    pdf2.Draw();
  
    cFCgraphics pdf3;
  
    pdf3.AddObject(h_phi);
  
    pdf3.Draw();
  
    cFCgraphics pdf4;
  
    pdf4.AddObject(h_noise100);
  
    pdf8.Draw();
  
    cFCgraphics pdf8;
  
    pdf4.AddObject(h_tempo);
  
    pdf8.Draw();
  
  
    cFCgraphics pdf5;
  
    pdf5.AddObject(h_tracos);
  
    pdf5.Draw();
  
  
    cFCgraphics pdf6;
  
    pdf6.AddObject(h_raio);
  
    pdf6.Draw();
  
    cFCgraphics pdf7;
  
    pdf7.AddObject(h_ponto);
  
    pdf7.Draw();
  
  
  
  
  
    return 0;
}
