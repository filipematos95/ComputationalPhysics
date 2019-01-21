//
//  PixelDet.C
//  Trab2_Ex1
//
//  Created by Joao Ferreira on 29/11/14.
//  Copyright (c) 2014 Joao Ferreira & Filipe Matos. All rights reserved.
//

#include <stdio.h>
#include "TH1F.h"
#include <string>
#include <iostream>
#include "PixelDet.h"
#include <ctime>
#include <cmath>
#include "TH2F.h"
#include <string>
#include <vector>


using namespace std;

PixelDet::PixelDet(int dimensao){
    dim = dimensao;
}

//Funcção random shuffle

vector<double> PixelDet::drandshuflle(vector<double> vec,int len){
    
    int j;
    double temp;
    for (int i=len; i < 1; i--){
        j=rand()%i;
        if (j != i ){
            temp = vec[j];
            vec[j]= vec[i];
            vec[i]= temp;
        }
    }
    
    return vec;
}


vector<int> PixelDet::irandshuflle(vector<int> vec,int len){
    
    int j;
    int temp;
    for (int i=len; i > 1; i--){
        j=rand()%i;
        if (j != i ){
            temp = vec[j];
            vec[j]= vec[i];
            vec[i]= temp;
        }
    }
    
    return vec;
}


double* PixelDet::theta_coo(double R,double theta){
    
    double* xy = new double[2];
    double x,y;
    x = R*cos(theta);
    y = R*sin(theta);
    
    xy[0] = x;
    xy[1] = y;
    
    return xy;
    
}


// Gera um evento de ruido - Retorna um array com os ID dos pixeis ruidosos

vector<int> PixelDet::EventNoise(double probnoise){
    vector<int> PixID_noise;
    PixID_noise.reserve(6400);
    double n;
    for (int i = 0; i < (dim*dim); i++){
        
        n = (rand()/double(RAND_MAX))*100.;
        
        if ( n < probnoise){
            PixID_noise.push_back(i);
        }
    }
    
    
    return PixID_noise;
}


//Gera N eventos ruidosos - Retorna um histograma com a sua projecção

/*TH1F* PixelDet::GetNoise(double probnoise, int N, string proj){
 
 double *Pix_noise_x = new double[100000];
 double *Pix_noise_y = new double[100000];
 int n = 0;
 double *data = new double[100000];
 
 
 //Recber os ID dos pixeis ruidosos;
 int *Pix_ID = new int[dim*dim];
 int it;
 double *xy = new double[2];
 
 int lo = 0;
 
 for ( int i = 0; i < N; i++){
 Pix_ID = EventNoise(0.5);
 for ( int j = 0 ; Pix_ID[j] != 0 ;j++){
 xy = ID_CU(Pix_ID[j]);
 Pix_noise_x[lo] = xy[0];
 Pix_noise_y[lo] = xy[1];
 lo++;
 }
 }
 
 if ( proj == 'x')
 {
 for (int i = 0 ; i < lo ; i++)
 data[i] = Pix_noise_x[i];
 }
 
 else if (proj == 'y')
 {
 for( int i = 0; i < lo; i++)
 data[i] = Pix_noise_y[i];
 }
 
 else{
 cout << "Erro: Argumento de projeção deverá ser 'x' ou 'y'  -> Returnado um NULL"<< endl;
 }
 
 
 
 TH1F* proje = new TH1F("projecao", "Projeção", dim, 0,dim*5);
 
 for(int i = 0; data[i] != 0 ; i++){
 proje->Fill(data[i]);
 }
 
 return proje;
 }*/

/* Converte um array de ID de pixel em coordenadas ...
 
 int** PixelDet::ID_COO(int* Pix_ID_C){
 
 int* ID_X = new int[dim];
 int* ID_Y = new int[dim];
 int k;
 
 for (int i = 0; i < (dim); i++){
 k = Pix_ID_C[i]/dim;
 ID_Y[i] = k+1;
 ID_X[i] = (((double(Pix_ID_C[i])/double(dim))-double(k)) * dim);
 }
 
 
 int **ID = new int* [6400];
 
 ID[0]  = ID_X;
 ID[1] = ID_Y;
 
 return ID;
 }
 
 */

//Função converte ID em Coordenadas
double* PixelDet::ID_CU(int Pixel){
    
    double x;
    double y;
    int k;
    
    k = Pixel/dim;
    y = double(k)*5. + 2.5;
    x = ((double(Pixel)/dim)-double(k))*dim*5. + 2.5;
    
    double* pixel_xy = new  double[2];
    
    pixel_xy[0] = x;
    pixel_xy[1] = y;
    
    return pixel_xy;
}


int PixelDet::CU_ID(double x, double y){
    
    int ID = (int(y)/5)*dim + (int(x)/5+1) - 1;
    
    return ID;
}

vector<int> PixelDet::EventSignal1(double a[2]){
    double x,y;
    
    x = a[0];
    y = a[1];
    
    
    vector<int> ID_list(9);
    
    //Quadrado de colisão
    int ID_col = CU_ID(x,y);
    
    int mit = 0;
    for (int i = ID_col-dim ; i < (ID_col+160) ; i += dim){
        for( int j = -1;   j < 2 ; j++){
            ID_list[mit] = i+j;
            mit++;
        }
    }
    
    //Vamos lá baralhar isto
    int j,temp;
    
    
    ID_list = irandshuflle(ID_list,8);
    
    
    
    vector<int> ID_chosen(ID_list.begin()+1,ID_list.begin()+6);
    
    
    
    return ID_chosen;
    
}


double PixelDet::XY_theta(double x, double y){
    double theta;
    
    if(x > 0 && y > 0){
        
        
        theta = atan(y/x);
        
    }
    else if(x < 0)
        theta = atan(y/x)+M_PI;
    else if(x > 0 && y < 0)
        theta = 2*M_PI+atan(y/x);
    else if(x == 0 && y > 0)
        theta = M_PI/2;
    else if(x == 0 && y < 0)
        theta = 3*M_PI/2;
    else{
        theta = 0;
    }
    return theta;
}

vector<int> PixelDet::EventSignal2(double beta,double a[2]){
    
    int N0 = int(20.*(1-(1/((1.30*beta)*(1.30*beta)))));
    
    double R;
    if(1.3*beta >= 1 && N0 > 0)
        R = 40*sqrt((1.3*beta)*(1.3*beta)-1);
    else{
        cout << "NULL" << endl;
    }
    
    //Bora lá descobrir os pontos em y
    
    double z;
    double x,y,ym,yn,xm,xn;
    
    int ID = CU_ID(a[0],a[1]);
    double *xy = new double[2];
    xy = ID_CU(ID);
    
    double x0 = xy[0]-2.5;
    double y0 = xy[1]-2.5;
    
    
    delete[] xy;
    
    vector<double> Pontosy;
    vector<double> Pontosx;
    
    int conta = 0;
    
    double lim = double(int(R)/5+1)*5;
    
    
    for ( double x = x0-lim; x <= x0+lim; x = x + 5 ){
        
        z = (x-x0)*(x-x0)-R*R+y0*y0;
        
        if ( (4.*y0*y0 - 4.*z ) >= 0. ){
            
            ym = (2.*y0 + sqrt(4.*y0*y0 - 4.*z))/(2.);
            yn = ( 2.*y0 - sqrt(4.*y0*y0 - 4.*z))/(2.);
            
            if ( ym > 0 ){
                Pontosy.push_back(ym);
                Pontosx.push_back(x);
                conta++;
            }
            if ( yn > 0 ){
                Pontosy.push_back(yn);
                Pontosx.push_back(x);
                conta++;
            }
        }
    }
    
    
    for ( double y = y0-lim; y <= y0+lim; y= y+5 ){
        
        z = (y-y0)*(y-y0) -R*R + x0*x0;
        
        if (( 4.*x0*x0 - 4.*z) >= 0.){
            
            xm = (2.*x0 + sqrt(4.*x0*x0 - 4.*z))/(2.);
            xn = (2.*x0 - sqrt(4.*x0*x0 - 4.*z))/(2.);
            
            if ( xm > 0 ){
                
                Pontosx.push_back(xm);
                Pontosy.push_back(y);
                conta++;
            }
            if ( xn > 0 ){
                
                Pontosx.push_back(xn);
                Pontosy.push_back(y);
                conta++;
            }
            
        }
    }
    
    
    //Indice conta-1 corresponde ao ultimo elemento da lista
    
    //Conversão em polar
    vector<double> thetas(conta);
    for ( int i = 0; i < conta; i++){
        thetas[i] = (XY_theta(Pontosx[i]-x0,Pontosy[i]-y0));
    }
    
    
    double aux;
    for(int i=conta-1; i >= 1; i--) {
        for( int j=0; j < i ; j++) {
            if(thetas[j]>thetas[j+1]) {
                aux = thetas[j];
                thetas[j] = thetas[j+1];
                thetas[j+1] = aux;
            }
        }
    }
    
    vector<pair<int,double>> ID_THETA(conta);
    
    double theta_av;
    
    for ( int i = 0; i < conta; i++){
        theta_av = (thetas[i+1]+thetas[i])/2;
        ID_THETA[i] = make_pair(CU_ID(theta_coo(R,theta_av)[0]+x0,theta_coo(R,theta_av)[1]+y0),thetas[i+1]-thetas[i]);
    }
    
    ID_THETA[conta-1].first = CU_ID(theta_coo(R,(2*M_PI+thetas[conta-1])/2)[0],theta_coo(R,(2*M_PI+thetas[conta-1])/2)[0]);
    ID_THETA[conta-1].second = 2*M_PI-thetas[conta-1];
    
    
    
    //Comprimentos dos traços:
    
    double len_s;
    for ( int i = 0; i < conta-2; i++){
        ID_THETA[i].second = R*ID_THETA[i].second;
    }
    
    ID_THETA[conta-1].second = R*(2*M_PI-thetas[conta-1]);
    
    
    //Escolha aleatória dos elementos
    
    int j;
    
    pair<int,double> temp;
    
    for (int i=conta-2; i >= 1; i--){
        j=rand()%i;
        if (j != i ){
            temp = ID_THETA[j];
            ID_THETA[j]= ID_THETA[i];
            ID_THETA[i]= temp;
        }
    }
    
    vector<int> ID_chosen(N0);
    j = 0;
    for (int i = 0; i < N0; i++){
        if( ID_THETA[i].second > 1.5){
            ID_chosen[j] = ID_THETA[i].first;
            tracos.push_back(ID_THETA[i].second);
            angulo.push_back(XY_theta(ID_CU(ID_chosen[j])[0]-x0,ID_CU(ID_chosen[j])[1]-y0));
            j++;
        }
        else{
            N0++;
        }
    }
    
    
    return ID_chosen;
    
}




// Função para mostrar o evento num array bidimensional

TH2C* PixelDet::DrawEvent(double a[2], double beta, double probnoise){
    
    int N0 = int(20.*(1-(1/((1.30*beta)*(1.30*beta)))));
    
    if (N0  <= 0)
        N0 = 0;
    
    vector<int> ID_noise;
    vector<int> ID_signal1;
    vector<int> ID_signal2;
    
    ID_noise = EventNoise(probnoise);
    
    ID_signal1 = EventSignal1(a);
    
    TH2C *event = new TH2C("evento","Evento",dim,0,400,dim,0,400);
    
    vector<int> ID_total;
    int ito = 0;
    
    double* xy = new double[2];
    
    for (int j = 0 ; j < ID_noise.size() ;j++){
        ID_total.push_back(ID_noise[j]);
        xy = ID_CU(ID_noise[j]);
        event->Fill(xy[0],xy[1]);
        
    }
    
    
    for ( int j = 0 ; j < 5 ;j++){
        ID_total.push_back(ID_signal1[j]);
        xy = ID_CU(ID_signal1[j]);
        event->Fill(xy[0],xy[1],2);
        
        
    }
    
    if (N0 > 0){
        ID_signal2 = EventSignal2(beta,a);
        
        for ( int j = 0 ; j < N0 ;j++){
            ID_total.push_back(ID_signal2[j]);
            xy = ID_CU(ID_signal2[j]);
            event->Fill(xy[0],xy[1],3);
            
        }
    }
    
    PixID = ID_total;
    
    return event;
}

//Simulação
void PixelDet::SimulateEvent(double a[2], double beta, double probnoise){
    
    int N0 = int(20.*(1-(1/((1.30*beta)*(1.30*beta)))));
    
    if (N0  <= 0)
        N0 = 0;
    
    vector<int> ID_noise;
    vector<int> ID_signal1;
    vector<int> ID_signal2;
    
    ID_noise = EventNoise(probnoise);
    
    ID_signal1 = EventSignal1(a);
    
    vector<int> ID_total;
    int ito = 0;
    
    double* xy = new double[2];
    
    for (int j = 0 ; j < ID_noise.size() ;j++){
        ID_total.push_back(ID_noise[j]);
    }
    
    
    for ( int j = 0 ; j < 5 ;j++){
        ID_total.push_back(ID_signal1[j]);
    }
    
    if (N0 > 0){
        ID_signal2 = EventSignal2(beta,a);
        for ( int j = 0 ; j < N0 ;j++){
            ID_total.push_back(ID_signal2[j]);
        }
    }
    
    PixID = ID_total;
    
}


void PixelDet::Event_Rec(){
    
    vector<int> activos = PixID;
    vector<int> ponto(2);
    
    int i = 0, j=0, k=0;
    
    int tamanho = activos.size();
    
    //Ordenação dos pixeis activos
    
    int m,aux;
    
    for(m=tamanho-1; m >= 1; m--) {
        for( int n=0; n < m ; n++) {
            if(activos[n]>activos[n+1]) {
                aux = activos[n];
                activos[n] = activos[n+1];
                activos[n+1] = aux;
            }
        }
    }
    
    
    int ID_ponto = 0;
    int conta = 0;
    int s;
    
    //Lista com pontos entre x(100 e 300) e y(100 e 300)
    
    double* p_xy = new double[2];
    
    
    vector<int> activos_r;
    
    for (i = 0; i < tamanho ;i++){
        p_xy = ID_CU(activos[i]);
        if (p_xy[0] >= 100 && p_xy[0] <= 300 & p_xy[1] >= 100 & p_xy[1] <= 300){
            activos_r.push_back(activos[i]);
        }
    }
    
    int tam = activos_r.size();
    
    delete[] p_xy;
    
    
    
    
    //Reconhecimento do ponto de impacto
    
    for (i = 0; i < tam; i++){
        conta = 0;
        for ( j = 0; j < 240 ; j = j + dim ){
            for ( k = 0; k < 3; k++){
                for( int mn = 0; mn < tam; mn++){
                    if(activos_r[i]+j+k == activos_r[mn]){
                        conta++;
                    }
                }
            }
        }
        if (conta >= 5){
            
            
            ID_ponto = activos_r[i]+81;
            break;
        }
    }
    
    
    
    if (ID_ponto == 0){
        for (i = 0; i < tam; i++){
            conta = 0;
            for ( j = 0; j < 240; j = j + dim ){
                for ( k = 0; k > -3; k--){
                    for( int mn = 0; mn < tam; mn++){
                        if(activos_r[i]+j+k == activos_r[mn]){
                            conta++;
                        }
                    }
                }
            }
            if (conta >= 5){
              
                ID_ponto = activos_r[i]+79;
                break;
            }
        }
    }
    
    if (ID_ponto == 0){
        for (i = 0; i < tam; i++){
            conta = 0;
            for ( j = 0; j > -240; j = j - dim ){
                for ( k = 0; k < 3; k++){
                    for( int mn = 0; mn < tam; mn++){
                        if(activos_r[i]+j+k == activos_r[mn]){
                            conta++;
                        }
                    }
                }
            }
            if (conta >= 5){
                ID_ponto = activos_r[i]-79;
              
                
                break;
            }
        }
    }
    
    if (ID_ponto == 0){
        for (i = 0; i < tam; i++){
            conta = 0;
            for ( j = 0; j > -240; j = j - dim ){
                for ( k = 0; k > -3; k--){
                    for( int mn = 0; mn < tam; mn++){
                        if(activos_r[i]+j+k == activos_r[mn]){
                            conta++;
                        }
                    }
                }
            }
            if (conta >= 5){
                ID_ponto = activos_r[i]-81;
              
                
                break;
            }
        }
    }
    
    
    
    
    
    
    if (ID_ponto == 0){
        vector<int> ID(6400);
        k = 0;
        
        for (i = 0; i < (dim*dim) ; i++ ){
            if( activos_r[k] == i)
                ID[i] = i;
            else
                ID[i] =0;
        }
        
        for (i = 0; i < tam; i++){
            if(activos_r[i]+dim == ID[activos_r[i]+dim]){
                if (activos_r[i]+79 == ID[activos_r[i]+79]) {
                    if (activos_r[i]+81 == ID[activos_r[i]+81]) {
                        if(activos_r[i]+160 == ID[activos_r[i]+160]){
                            ID_ponto = activos_r[i]+dim;
                          
                        }
                    }
                }
            }
        }
    }
    
    
    
    
    //Pontos candidatos a pixeis cerenkov
    
    vector<int> pontos_cerenkov;
    int ponto_c;
    m = 0;
    
    for (i = 0; i < tamanho ; i++){
        ponto_c = activos[i];
        if(ponto_c != ID_ponto && ponto_c != ID_ponto+dim && ponto_c != ID_ponto-dim){
            if(ponto_c != ID_ponto+1 && ponto_c != ID_ponto-1){
                if(ponto_c != ID_ponto+81 && ponto_c != ID_ponto-81){
                    if(ponto_c != ID_ponto+79 && ponto_c != ID_ponto-79){
                        pontos_cerenkov.push_back(ponto_c);
                        m++;
                    }
                }
            }
        }
    }
    
    
    
    
    //Reconhecer o raio
    
    vector<double> raios;
    
    double raio;
    
    double* xy_r = new double[2];
    double* xy_p = new double[2];
    
    xy_p = ID_CU(ID_ponto);
    
    for ( i = 0; i < pontos_cerenkov.size() ; i++){
        xy_r = ID_CU(pontos_cerenkov[i]);
        raios.push_back(sqrt((xy_r[0]-xy_p[0])*(xy_r[0]-xy_p[0])+(xy_r[1]-xy_p[1])*(xy_r[1]-xy_p[1])));
    }
    
    
    delete[] xy_r;
    delete[] xy_p;
    
    
    // ordenaçao dos raios
    int aux_r;
    
    for(i=raios.size()-1; i >= 1; i--) {
        for( j=0; j < i ; j++) {
            if(raios[j]>raios[j+1]) {
                aux_r = raios[j];
                raios[j] = raios[j+1];
                raios[j+1] = aux_r;
            }
        }
    }
    
    int n = 0;
    
    
    vector<double> ini;
    int verified = 0;
    for (i = 0; i < raios.size()-1; i++){
        if (raios[i] < 39 && raios[i+1] < 39){
            if((raios[i+1]-raios[i]) < 6){
                if (verified == 0)
                    ini.push_back(raios[i]);
                ini.push_back(raios[i+1]);
                n++;
                verified = 1;
            }
        }
    }
    
    
    
    double soma = 0;
    
    for( i = 0; i < ini.size(); i++){
        soma = soma + ini[i] ;
    }
    
    
    
    
    double raio_final;
    
    
    if ( n == 0 ){
        if (raios[0] < 16){
            raio_final = raios[0];
            n++;
        }
        else
            raio_final = 0;
    }
    else{
        n++;
        raio_final = soma/n;
    }
    
    
    double beta_rec;
    
    if( n > 0 )
        beta_rec = sqrt((raio_final/40.)*(raio_final/40.) + 1.)/1.3;
    else
        beta_rec = 0.745;
    
    
    veloc_rec = beta_rec;
    raio_rec = raio_final;
    ponto_impacto = ID_ponto;
    
}




