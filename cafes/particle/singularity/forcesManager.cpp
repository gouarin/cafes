#ifndef FORCESMANAGER_CPP
#define FORCESMANAGER_CPP

#include <petsc.h>

#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "DtoN/UandPNormal.h"
#include "DtoN/UandPTang.h"

namespace DtoN{

  template<std::size_t Dimensions>
  class ForcesManager{

    static double mu;

    PetscReal _H1;
    PetscReal _H2;
    PetscReal _a;
    PetscReal _UN;
    PetscReal _UT;
    PetscReal _l;
    PetscReal _eps;
    double    _cutoffdist; 
    double    _theta;
    double    _A;
    int       _N; 

    std::vector<std::vector< double> > _Spoints;
    std::vector< std::vector<double> > _n;

    std::vector<double> _FN;
    std::vector<double> _FT;


    //Initialize --------------------------------------------------------------------------------
    template<int T>
    void initialize(){}
    //-------------------------------------------------------------------------------------------

    //Initialize normal forces--------------------------------------------------------------------
    template<int T>
    void initForcesNormal(){}
    //-------------------------------------------------------------------------------------------
    
    //Initialize tangential forces----------------------------------------------------------------
    template<int T>
    void initForcesTang(){}
    //-------------------------------------------------------------------------------------------

  public:

    ForcesManager(PetscReal H1, PetscReal H2, PetscReal a, PetscReal UN, PetscReal UT, PetscReal l, PetscReal eps, double cutoffdist);

    void setN(int N);

    void getForcesNormal(std::vector<double>& FN);
    void printForcesNormal();
    void clearForcesNormal();

    void getForcesTang(std::vector<double>& FT);
    void printForcesTang();
    void clearForcesTang();

    void getForces(std::vector<double>& F);	 
    void printForces();
    void clearForces();
  };

  template<std::size_t Dimensions> double ForcesManager<Dimensions>::mu = 1.0;





  //Initialize normal forces--------------------------------------------------------------------
  template<>
  template<>
  void ForcesManager<2>::initialize<2>(){

    _Spoints.clear();
    _n.clear();

    std::vector<double> tmpSpoints;
    std::vector<double> tmpn;
    tmpSpoints.resize(2);
    tmpn.resize(2);

    for(int i=-_N;i<=_N;i++)
      {
	//position points sur sphere (A AJOUTER : S'assurer que H1 rayon de la sphere et pas du mur.)
	tmpSpoints[0] = (cos(i*_theta/_N) - 1)/_H1;
	tmpSpoints[1] =  sin(i*_theta/_N)/_H1;
	_Spoints.push_back(tmpSpoints);

	//normal sortante
	tmpn[0]     = cos(i*_theta/_N);
	tmpn[1]     = sin(i*_theta/_N);
	_n.push_back(tmpn);
      }
  }

  template<>
  template<>
  void ForcesManager<3>::initialize<3>(){
    _A = 2*M_PI*_theta/_N/_N/_H1/_H1;

    _Spoints.clear();
    _n.clear();

    std::vector<double> tmpSpoints;
    std::vector<double> tmpn;
    tmpSpoints.resize(3);
    tmpn.resize(3);


    for(int i=1;i<=_N;i++){
      for(int j=1;j<=_N;j++){
	//position points sur sphere (A AJOUTER : S'assurer que H1 rayon de la sphere et pas du mur.)
	tmpSpoints[0] =  sin(i*_theta/_N)*cos(2*M_PI*j/_N)/_H1;
	tmpSpoints[1] =  sin(i*_theta/_N)*sin(2*M_PI*j/_N)/_H1;
	tmpSpoints[2] = (cos(i*_theta/_N)-1)/_H1;
	_Spoints.push_back(tmpSpoints);

	//std::cout << "_N = " << _N << " ; size = " << _Spoints.size() << std::endl;

	//normal sortante
	tmpn[0]     = sin(i*_theta/_N)*cos(2*M_PI*j/_N);
	tmpn[1]     = sin(i*_theta/_N)*sin(2*M_PI*j/_N);
	tmpn[2]     = cos(i*_theta/_N);
	_n.push_back(tmpn);
      }
    }
  }
  //-------------------------------------------------------------------------------------------


  //Initialize normal forces--------------------------------------------------------------------
  template<>
  template<>
  void ForcesManager<2>::initForcesNormal<2>(){
    std::vector< std::vector<double> > sigma(2);
    for(int i=0;i<2;i++)
      sigma[i].resize(2);

    _FN.clear();
    _FN.resize(2);

    for(int i=0;i<=2*_N;i++)
      {
	double dxux = dxux_sing_normalMvt2D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);
	double dxuz = dxuz_sing_normalMvt2D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);

	double dzux = dzux_sing_normalMvt2D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);
	double dzuz = dzuz_sing_normalMvt2D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);

	double P = p_sing_withT_normalMvt2D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);


	sigma[0][0] = 2*mu* dxux - P;
	sigma[0][1] =   mu*(dzux + dxuz);
	sigma[1][0] = sigma[0][1];
	sigma[1][1] = 2*mu* dzuz - P;

	for(int j=0;j<2;j++)
	  for(int k=0;k<2;k++)
	    _FN[j] += _theta/_N/_H1*sigma[j][k]*_n[i][k];
      }
  }

  template<>
  template<>
  void ForcesManager<3>::initForcesNormal<3>(){
    std::vector< std::vector<double> > sigma(3);
    for(int i=0;i<3;i++)
      sigma[i].resize(3);

    _FN.clear();
    _FN.resize(3);

    for(int i=0;i<_N*_N;i++)
      {
	double dxux = dxux_sing_normalMvt3D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);
	//dxuy = dyux
	double dxuz = dxuz_sing_normalMvt3D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);

	double dyux = dyux_sing_normalMvt3D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);
	double dyuy = dyuy_sing_normalMvt3D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);
	double dyuz = dyuz_sing_normalMvt3D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);

	double dzux = dzux_sing_normalMvt3D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);
	double dzuy = dzuy_sing_normalMvt3D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);
	double dzuz = dzuz_sing_normalMvt3D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);

	double P = p_sing_withT_normalMvt3D(_Spoints[i],_H1,_H2,_a,_UN,_l,_eps,NULL);


	sigma[0][0] = 2*mu* dxux - P;
	sigma[0][1] = 2*mu* dyux;
	sigma[0][2] =   mu*(dzux + dxuz);

	sigma[1][0] = sigma[0][1];
	sigma[1][1] = 2*mu* dyuy - P;
	sigma[1][2] =   mu*(dzuy + dyuz) ;

	sigma[2][0] = sigma[0][2];
	sigma[2][1] = sigma[1][2];
	sigma[2][2] = 2*mu* dzuz - P;

	for(int j=0;j<3;j++)
	  for(int k=0;k<3;k++)
	    _FN[j] += _A*sigma[j][k]*_n[i][k]*sin(ceil(i/_N)*_theta/_N);
      }
  }
  //-------------------------------------------------------------------------------------------



  //Initialize tangential forces----------------------------------------------------------------
  template<>
  template<>
  void ForcesManager<2>::initForcesTang<2>(){
    std::vector< std::vector<double> > sigma(2);
    for(int i=0;i<2;i++)
      sigma[i].resize(2);

    _FT.clear();
    _FT.resize(2);

    for(int i=0;i<=2*_N;i++)
      {
	double dxux = dxux_sing_tangMvt2D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);
	double dxuz = dxuz_sing_tangMvt2D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);

	double dzux = dzux_sing_tangMvt2D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);
	double dzuz = dzuz_sing_tangMvt2D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);

	double P = p_sing_withT_tangMvt2D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);


	sigma[0][0] = 2*mu* dxux - P;
	sigma[0][1] =   mu*(dzux + dxuz);
	sigma[1][0] = sigma[0][1];
	sigma[1][1] = 2*mu* dzuz - P;

	for(int j=0;j<2;j++)
	  for(int k=0;k<2;k++)
	    _FT[j] += _theta/_N/_H1*sigma[j][k]*_n[i][k];
      }
  }

  template<>
  template<>
  void ForcesManager<3>::initForcesTang<3>(){
    std::vector< std::vector<double> > sigma(3);
    for(int i=0;i<3;i++)
      sigma[i].resize(3);

    _FT.clear();
    _FT.resize(3);

    for(int i=0;i<_N*_N;i++)
      {
	double dxux = dxux_sing_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);
	double dxuy = dxuy_sing_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);
	double dxuz = dxuz_sing_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);

	double dyux = dyux_sing_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);
	double dyuy = dyuy_sing_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);
	double dyuz = dyuz_sing_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);

	double dzux = dzux_sing_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);
	double dzuy = dzuy_sing_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);
	double dzuz = dzuz_sing_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);

	double P = p_sing_withT_tangMvt3D(_Spoints[i],_H1,_H2,_a,_UT,_l,_eps,NULL);

     
	sigma[0][0] = 2*mu* dxux - P;
	sigma[0][1] =   mu*(dyux + dxuy);
	sigma[0][2] =   mu*(dzux + dxuz);

	sigma[1][0] = sigma[0][1];
	sigma[1][1] = 2*mu* dyuy - P;
	sigma[1][2] =   mu*(dzuy + dyuz) ;

	sigma[2][0] = sigma[0][2];
	sigma[2][1] = sigma[1][2];
	sigma[2][2] = 2*mu* dzuz - P;

	for(int j=0;j<3;j++)
	  for(int k=0;k<3;k++)
	    _FT[j] += _A*sigma[j][k]*_n[i][k]*sin(ceil(i/_N)*_theta/_N);
      }
  }
  //-------------------------------------------------------------------------------------------

  template<std::size_t Dimensions>
  ForcesManager<Dimensions>::ForcesManager(PetscReal H1, PetscReal H2, PetscReal a, PetscReal UN, PetscReal UT, PetscReal l, PetscReal eps, double cutoffdist) : _H1(H1), _H2(H2), _a(a), _UN(UN), _UT(UT), _l(l), _eps(eps), _cutoffdist(cutoffdist), _N(100){
    _theta = asin(cutoffdist*H1);
    initialize<Dimensions>();
  }
  
  template<std::size_t Dimensions>
  void ForcesManager<Dimensions>::setN(int N){_N=N;initialize<Dimensions>();}

  template<std::size_t Dimensions>
  void ForcesManager<Dimensions>::getForcesNormal(std::vector<double>& FN){
    if(_FN.empty())
      initForcesNormal<Dimensions>();
    FN = _FN;
  }
    
  template<std::size_t Dimensions>
  void ForcesManager<Dimensions>::printForcesNormal(){
    if(_FN.empty())
      initForcesNormal<Dimensions>();
    std::cout << "FNx = " << _FN[0] << "; FNy = " << _FN[1] ;
    if(Dimensions==3)
      std::cout << "; FNz = " << _FN[2] ;
    std::cout << std::endl;
  }

   template<std::size_t Dimensions>
   void ForcesManager<Dimensions>:: clearForcesNormal(){
     _FN.clear();
   }

  template<std::size_t Dimensions>
  void ForcesManager<Dimensions>::getForcesTang(std::vector<double>& FT){
    if(_FT.empty())
      initForcesTang<Dimensions>();
    FT = _FT;
  }  

  template<std::size_t Dimensions>
  void ForcesManager<Dimensions>::printForcesTang(){
    if(_FT.empty())
      initForcesTang<Dimensions>();
    std::cout << "FTx = " << _FT[0] << "; FTy = " << _FT[1];
    if(Dimensions==3)
      std::cout << "; FTz = " << _FT[2] ;
    std::cout << std::endl;
  }

   template<std::size_t Dimensions>
   void ForcesManager<Dimensions>:: clearForcesTang(){
     _FT.clear();
   }

  template<std::size_t Dimensions>
  void ForcesManager<Dimensions>::getForces(std::vector<double>& F){
    if(_FN.empty())
      initForcesNormal<Dimensions>();
    
    if(_FT.empty())
      initForcesTang<Dimensions>();

    F    = _FN;
    for(int i; i<Dimensions; i++)
      F[i] += _FT[i];
  }	 

  template<std::size_t Dimensions>
  void ForcesManager<Dimensions>::printForces(){
    if(_FN.empty())
      initForcesNormal<Dimensions>();

    if(_FT.empty())
      initForcesTang<Dimensions>();

    std::cout << "Fx = " << _FT[0]+_FN[0] << "; Fy = " << _FT[1]+_FN[1];
    if(Dimensions==3)
      std::cout << "; Fz = " << _FT[2]+_FN[2] ;
    std::cout << std::endl;
  }

   template<std::size_t Dimensions>
   void ForcesManager<Dimensions>:: clearForces(){
     _FN.clear();
     _FT.clear();
   }
};

#endif
