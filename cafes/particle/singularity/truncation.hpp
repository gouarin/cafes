// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
// OF SUCH DAMAGE.

#ifndef CAFES_PARTICLE_SINGULARITY_TRUNCATION_HPP_INCLUDED
#define CAFES_PARTICLE_SINGULARITY_TRUNCATION_HPP_INCLUDED

#include <petsc.h>
#include <cassert>
#include <iostream>

namespace cafes
{
  namespace singularity
  {
    PetscScalar alphaTrunc  (PetscReal X){

      PetscReal X2 = 0;
      PetscReal X3 = 0;
      PetscReal X5 = 0;
      PetscReal X7 = 0;
      PetscReal X9 = 0;

      X2 = X*X;
      X3 = X*X2;
      X5 = X3*X2;
      X7 = X5*X2;
      X9 = X7*X2;
      return 315/128.0*(X9/9.0 - 4.0/7.0*X7 + 6.0/5.0*X5 - 4.0/3.0*X3 + X);
        
    }

    PetscScalar dalphaTrunc (PetscReal X){

      PetscReal X2Mun  = 0;
      PetscReal X2Mun2 = 0;
      PetscReal X2Mun4 = 0;

      X2Mun  = X*X-1;
      X2Mun2 = X2Mun*X2Mun;
      X2Mun4 = X2Mun2*X2Mun2;

      return 315/128.0*X2Mun4;
    }

    PetscScalar d2alphaTrunc(PetscReal X){

      PetscReal X2Mun  = 0;
      PetscReal X2Mun3 = 0;

      X2Mun  = X*X-1;
      X2Mun3 = X2Mun*X2Mun*X2Mun;

      return 2520/128.0*X2Mun3*X;
    }
    PetscScalar d3alphaTrunc(PetscReal X){

      PetscReal X2Mun     = 0;
      PetscReal X2Mun2    = 0;
      PetscReal septX2Mun = 0;

      X2Mun     = X*X-1;
      X2Mun2    = X2Mun*X2Mun;
      septX2Mun = 7*X*X-1;

      return 2520/128.0*X2Mun2*septX2Mun;
    }


    PetscScalar d4alphaTrunc(PetscReal X){

      PetscReal X2 = 0;
      PetscReal X4 = 0;

      X2 = X*X;
      X4 = X2*X2;

      return 15120/128.0*(7*X4-10*X2+3);
    }

   PetscScalar chiTrunc(PetscReal X, PetscReal l, PetscReal eps)
   {
     assert(l>=eps);

     double X2 = X*X;

     if(0<=X2<l-eps)
       return 1;
     else 
       if(l-eps<=X2 && X2<l+eps)
         return 0.5*( 1-alphaTrunc((X2-l)/eps) );
       else
         return 0;
   }

   PetscScalar dchiTrunc(PetscReal X, PetscReal l, PetscReal eps)
   {
     assert(l>=eps);
      
     double X2 = X*X;
     X = sqrt(X2); //on veut la norme

     if(0<=X2<l-eps)
       return 1;
     else if(l-eps<=X2 && X2<l+eps)
       //return -X/eps*dalphaTrunc((X2-l)/eps) ;
       return - dalphaTrunc((X2-l)/eps)/eps/2.;
     else
       return 0;
   }

   PetscScalar d2chiTrunc(PetscReal X, PetscReal l, PetscReal eps)
   {
     assert(l>=eps);

     double X2 = X*X;
     //X = sqrt(X2); //on veut la norme

     if(0<=X<l-eps)
       return 1;
     else if(l-eps<=X2 && X2<l+eps)
       //return -2*X2/eps/eps*d2alphaTrunc((X2-l)/eps) - dalphaTrunc((X2-l)/eps)/eps;
       return - d2alphaTrunc((X2-l)/eps)/eps/eps/2.;
     else
       return 0;
   }

   PetscScalar d3chiTrunc(PetscReal X, PetscReal l, PetscReal eps)
   {
     assert(l>=eps);

     double X2 = X*X;
     //X = sqrt(X2); //on veut la norme

     if(0<=X<l-eps)
       return 1;
     else if(l-eps<=X2 && X2<l+eps)
       //return -4*X2*X/eps/eps/eps*d3alphaTrunc((X2-l)/eps) - 6*X/eps/eps*d2alphaTrunc((X2-l)/eps);
       return - d3alphaTrunc((X2-l)/eps)/eps/eps/eps/2.;
     else
       return 0;
   }

   PetscScalar d4chiTrunc(PetscReal X, PetscReal l, PetscReal eps)
   {
     assert(l>=eps);
     double X2 = X*X;
     //X = sqrt(X2); //on veut la norme

     if(0<=X<l-eps)
       return 1;
     else if(l-eps<=X2 && X2<l+eps)
       //return -4*X2*X/eps/eps/eps*d3alphaTrunc((X2-l)/eps) - 6*X/eps/eps*d2alphaTrunc((X2-l)/eps);
       return - d4alphaTrunc((X2-l)/eps)/eps/eps/eps/eps/2.;
     else
       return 0;  
   }


     // PetscScalar chiTrunc(PetscReal X, PetscReal l, PetscReal eps)
     // {
     //   return 1;
     // }

     // PetscScalar dchiTrunc(PetscReal X, PetscReal l, PetscReal eps)
     // {
     //   return 0;
     // }

     // PetscScalar d2chiTrunc(PetscReal X, PetscReal l, PetscReal eps)
     // {
     //   return 0;
     // }

     // PetscScalar d3chiTrunc(PetscReal X, PetscReal l, PetscReal eps)
     // {
     //   return 0;
     // }

     // PetscScalar d4chiTrunc(PetscReal X, PetscReal l, PetscReal eps)
     // {
     //   return 0;  
     // }
  }
}

#endif