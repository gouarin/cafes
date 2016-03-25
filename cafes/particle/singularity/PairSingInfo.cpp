#ifndef PAIRSINGINFO_CPP
#define PAIRSINGINFO_CPP

#include <petsc.h>

#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "SEM/cafesparticle.h"

template <typename T> int sgn(T val) {
  return (T(0) <= val) - (val < T(0));
}

namespace DtoN{

  template<typename PTYPE, std::size_t Dimensions>
  struct InitializeOXYZ_impl{};

  template<typename PTYPE, std::size_t Dimensions>
  class PairSingInfo{

  protected:


    std::vector<double> _center_p1;
    std::vector<double> _vel_p1;
    double              _H1;


    std::vector<double> _center_p2;
    std::vector<double> _vel_p2;
    double              _H2;

    std::vector<double> _posDiff;
    PetscScalar         _distC1C2;
    std::vector<double> _vdiff;

    std::vector< std::vector<double> > _base;
    std::vector<double>                _posOrigin;

    PetscScalar         _UX;
    PetscScalar         _UY;
    PetscScalar         _UZ;


    friend struct InitializeOXYZ_impl<PTYPE,Dimensions>;
    void initializeOXYZ(){return InitializeOXYZ_impl<PTYPE,Dimensions>::initializeOXYZ(*this);}
 
  public:

    PairSingInfo(PTYPE& part1, PTYPE& part2);

    //Simple accessor
    double getH1(){return _H1;}
    double getH2(){return _H2;}
    

    void getVelocity(std::vector<double>& U);
    void printVelocity();
    void clearVelocity();

    void getBase(std::vector<double>& base);
    void clearBase();


  };

  template<typename PTYPE>
  struct InitializeOXYZ_impl<PTYPE,2>{
    static void initializeOXYZ(PairSingInfo<PTYPE,2>& p){

      p._base.clear();
      p._base.resize(2);

      p._base[0] = p._posDiff;
      for(int i=0; i<2; i++)
        p._base[0][i]   /=  p._distC1C2;

      p._UX = std::inner_product(p._vdiff.begin(), p._vdiff.end(), p._base[0].begin(), 0.);

      p._base[1].resize(2);
      p._base[1][0]   -= p._base[0][1];
      p._base[1][1]    = p._base[0][0];

      p._UZ = std::inner_product(p._vdiff.begin(), p._vdiff.end(), p._base[1].begin(), 0.);

      double w = 1./p._H1;
      auto posOrigin_comp = [w](double x, double y){return x + w*y;};
      p._posOrigin.resize(2);
      std::transform(p._center_p1.begin(), p._center_p1.end(), p._base[0].begin(), p._posOrigin.begin(), posOrigin_comp);
    }
  };


  template<typename PTYPE>
  struct InitializeOXYZ_impl<PTYPE,3>{
    static void initializeOXYZ(PairSingInfo<PTYPE,3>& p){

      p._base.resize(3);

      p._base[2] = p._posDiff;
      for(int i=0; i<3; i++)
        p._base[2][i]   /=  p._distC1C2;

      p._UZ            = std::inner_product(p._vdiff.begin(), p._vdiff.end(), p._base[2].begin(), 0.);
      double UNorm   = std::inner_product(p._vdiff.begin(), p._vdiff.end(), p._vdiff  .begin(), 0.);

      p._base[0].resize(3);
      if(p._UZ*p._UZ == UNorm){
        p._base[0][0]    = -sgn(p._base[2][2])*p._base[2][2]-sgn(p._base[2][1])*p._base[2][1];
        p._base[0][1]    =  sgn(p._base[2][1])*p._base[2][0];
        p._base[0][2]    =  sgn(p._base[2][2])*p._base[2][0];
      }
      else{
        p._base[0][0]    =  p._vdiff[0]-p._vdiff[0]*p._base[2][0];
        p._base[0][1]    =  p._vdiff[1]-p._vdiff[1]*p._base[2][1];
        p._base[0][2]    =  p._vdiff[2]-p._vdiff[2]*p._base[2][2];
      }

      double tmp = sqrt(p._base[0][0]*p._base[0][0]+p._base[0][1]*p._base[0][1]+p._base[0][2]*p._base[0][2]);
      for(int i=0; i<3; i++)
        p._base[0][i] /= tmp;
   
      p._base[1].resize(3);
      p._base[1][0]    = p._base[2][1]*p._base[0][2]-p._base[0][1]*p._base[2][2];
      p._base[1][1]    = p._base[2][2]*p._base[0][0]-p._base[0][2]*p._base[2][0];
      p._base[1][2]    = p._base[2][0]*p._base[0][1]-p._base[0][0]*p._base[2][1];
      //deja norme

      p._UX = std::inner_product(p._vdiff.begin(), p._vdiff.end(), p._base[0].begin(), 0.);
      p._UY = std::inner_product(p._vdiff.begin(), p._vdiff.end(), p._base[1].begin(), 0.);

      double w = 1./p._H1;
      auto posOrigin_comp = [w](double x, double y){return x + w*y;};
      p._posOrigin.resize(3);
      std::transform(p._center_p1.begin(), p._center_p1.end(), p._base[2].begin(), p._posOrigin.begin(), posOrigin_comp);
    }
  };


  template<class PTYPE, std::size_t Dimensions>
  PairSingInfo<PTYPE,Dimensions>::PairSingInfo(PTYPE& part1, PTYPE& part2) :   _center_p1(part1->_center), _center_p2(part2->_center), _vel_p1(part1->_velocity), _vel_p2(part2->_velocity), _H1(1./part1->_radius), _H2(1./part2->_radius){

    auto diff = [](double x, double y){return x-y;};
    _posDiff.resize(Dimensions);
    std::transform(_center_p2.begin(), _center_p2.end(), _center_p1.begin(), _posDiff.begin(), diff);
    _vdiff.resize(Dimensions);
    std::transform(_vel_p2.cbegin()  , _vel_p2.cend()  , _vel_p1.cbegin()  , _vdiff.begin()  , diff); 

    _distC1C2 = std::inner_product(_posDiff.begin(), _posDiff.end(), _posDiff.begin(), 0.);
    _distC1C2 = sqrt(_distC1C2);
  }

  template<class PTYPE, std::size_t Dimensions>
  void PairSingInfo<PTYPE,Dimensions>::getVelocity(std::vector<double>& U){

    U.clear();

    if(_base.empty())
      initializeOXYZ();

    U.push_back(_UX);
    if(Dimensions==2)
      U.push_back(_UZ);
    else{
      U.push_back(_UY);
      U.push_back(_UZ);
    }
  }

  template<class PTYPE, std::size_t Dimensions>
  void PairSingInfo<PTYPE,Dimensions>::printVelocity(){

    if(_base.empty())
      initializeOXYZ();

    std::cout << "UX = " << _UX;
    if(Dimensions==3)
      std::cout << "; UY = " << _UY;
    std::cout << "; UZ = " << _UZ << std::endl;
  }

  template<class PTYPE, std::size_t Dimensions>
  void PairSingInfo<PTYPE,Dimensions>::clearVelocity(){
    _UX = 0;
    _UY = 0;
    _UZ = 0;
  }

  template<class PTYPE, std::size_t Dimensions>
  void PairSingInfo<PTYPE,Dimensions>::getBase(std::vector<double>& base){
    base = _base;
  }

  template<class PTYPE, std::size_t Dimensions>
  void PairSingInfo<PTYPE,Dimensions>::clearBase(){
    _base.clear();
  }

};


#endif
