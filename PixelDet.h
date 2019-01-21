//
//  PixelDet.h
//  Trab2_Ex1
//
//  Created by Joao Ferreira on 29/11/14.
//  Copyright (c) 2014 Joao Ferreira & Filipe Matos. All rights reserved.
//

#ifndef __Trab2_Ex1__PixelDet__
#define __Trab2_Ex1__PixelDet__

#include <stdio.h>
#include "TH1F.h"
#include <string>
#include <iostream>
#include "TH2F.h"
#include <vector>
using namespace std;

class PixelDet{
public:
  PixelDet(int);
  
  double XY_theta(double , double);
  vector<int> EventNoise(double); //Funciona
  TH1F* GetNoise( double , int , string ); //Funciona
  vector<int> EventSignal1(double*); // Está a funcionar bem
  vector<int> EventSignal2(double,double*); // está a funcionar
 
  double* ID_CU(int); // Converte ID para coordenadas
  int **ID_COO(int*); // Método que recebe um aray de ID e devolve um array com as ocordenadas
  int CU_ID(double , double ); // Converte Coordenadas para ID
  
  TH2C* DrawEvent(double* , double , double ); //Simula um evento e desenha-o
  
  vector<double> drandshuflle(vector<double>,int); //
  vector<int> irandshuflle(vector<int>,int);
  double* theta_coo(double ,double);
  
  void Event_Rec();
  void SimulateEvent(double*,double,double);
  vector<double> Get_tracos(){return tracos;}

  double Get_raio_Rec(){return raio_rec;}
  double Get_veloc_Rec(){return veloc_rec;}
  int Get_ponto_impacto(){return ponto_impacto;}
  vector<int> Get_PixID(){return PixID;}
  vector<double> Get_angulo(){return angulo;}
  
  void tracos_clear(){tracos.clear();}
  void angulo_clear(){angulo.clear();}
  
private:
  int dim;
  double raio_rec;
  double  veloc_rec;
  int ponto_impacto;
  vector<double> tracos;
  vector<int> PixID;
  vector<double> angulo;
};


#endif /* defined(__Trab2_Ex1__PixelDet__) */
