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

#ifndef CAFES_PARTICLE_SINGULARITY_UANDPNORMAL_HPP_INCLUDED
#define CAFES_PARTICLE_SINGULARITY_UANDPNORMAL_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/singularity/truncation.hpp>

#include <petsc.h>

#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace cafes
{
  namespace singularity
  {

    using position2d = geometry::position<double, 2>;
    using position3d = geometry::position<double, 3>;

    ///////////////////
    // PRESSURE
    ///////////////////

    PetscScalar p_sing_withT_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r = X[1];
      double z = X[0];
      double mu = 1;
      double t4 = r*r;
      double t8 = pow(2 * a + M * t4 + H * t4, 2);
      double t12 = chiTrunc(r, l, eps); // * chiTrunc(z,l,eps);
      return (-24 / (H + M) / t8 * U * mu * t12) * (z > -1./H) * (z< a + 1./M);
    }

    ///////////////////
    // Ur
    ///////////////////

    PetscScalar ur_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double r  = X[1];
      double z  = X[0];

      double t6 = 2 * z;
      double t7 = r * r;
      double t8 = M * t7;
      double t9 = 2 * a;
      double t10 = -t6 + t8 + t9;
      double t12 = H * t7;
      double t16 = (t12 + t6) / (H + M);
      double t17 = t9 + t8 + t12;
      double t18 = t17 * t17;
      double t20 = 1 / t18 / t17;
      double t21 = chiTrunc(r,l, eps);
      double t25 = dchiTrunc(r, l, eps);
      return  (-3 * U * (4 * H * r + 4 * M * r) * t10 * t16 * t20 * t21 - 3 * U * (-2 * t25 * a - t25 * M * t7 - t25 * H * t7) * t10 * t16 * t20 ) * (z > -1./H) * (z < a + 1./M);
    }

    PetscScalar drur_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double rr  = X[1];
      double zz  = X[0];
      double r1b = 1./H;
      double r2b = 1./M;
      double t1 = rr * rr;
      double t2 = H * t1;
      double t3 = zz * zz;
      double t7 = H * H;
      double t8 = H * t7;
      double t9 = t1 * t1;
      double t10 = t8 * t9;
      double t14 = t9 * t1;
      double t17 = M * M;
      double t21 = t17 * M;
      double t28 = t7 * t1;
      double t29 = a * a;
      double t38 = zz * a;
      double t42 = t9 * M;
      double t45 = -160 * t2 * t3 * M + 24 * t10 * a + 4 * t8 * M * t14 + 8 * t7 * t17 * t14 + 4 * H * t21 * t14 + 24 * t21 * zz * t9 - 48 * t28 * t29 - 32 * H * t29 * zz - 32 * t29 * M * zz + 128 * t28 * t38 - 24 * t7 * zz * t42;
      double t46 = H * t9;
      double t50 = H * a;
      double t51 = t17 * t9;
      double t54 = t1 * t17;
      double t60 = H * M;
      double t72 = a * M;
      double t75 = t72 * zz;
      double t78 = 24 * t46 * t17 * zz - 16 * t50 * t51 + 32 * t54 * t38 + 8 * t7 * a * t42 - 48 * t60 * t1 * t29 - 24 * t10 * zz + 32 * t50 * t3 - 80 * t54 * t3 - 80 * t28 * t3 + 32 * t72 * t3 + 160 * t2 * t75;
      double t84 = M * t1;
      double t86 = pow((2 * a + t84 + t2),2);
      double t87 = t86 * t86;
      double t89 = 1 / (H + M) / t87;
      double t90 = chiTrunc(rr, l, eps);
      double t93 = d2chiTrunc(rr, l, eps);
      double t94 = t29 * a;
      double t95 = t93 * t94;
      double t98 = t93 * t29;
      double t101 = t1 * rr;
      double t102 = t7 * t101;
      double t103 = dchiTrunc(rr, l, eps);
      double t104 = t103 * a;
      double t105 = t104 * zz;
      double t108 = t9 * rr;
      double t110 = t103 * t17;
      double t114 = t7 * t108;
      double t115 = t103 * M;
      double t120 = t103 * t29;
      double t121 = t120 * zz;
      double t136 = t93 * a;
      double t141 = 16 * t95 * zz - 16 * t98 * t3 - 16 * t105 * t102 - 12 * H * t108 * t110 * zz + 12 * t114 * t115 * zz - 48 * M * rr * t121 - 48 * t17 * t101 * t105 + 8 * t98 * t2 * zz + 12 * t98 * t42 * H + 24 * t98 * t84 * zz + 6 * t136 * t17 * t14 * H;
      double t145 = t93 * t7;
      double t150 = t93 * t17;
      double t155 = t14 * M;
      double t170 = a * rr * t3;
      double t173 = H * t103;
      double t181 = t9 * t101;
      double t185 = 12 * t136 * t51 * zz - 4 * t145 * t9 * a * zz + 2 * t150 * t14 * H * zz - 2 * t145 * t155 * zz - 24 * t114 * t104 * M + 8 * t136 * t155 * t7 - 12 * t17 * t108 * t104 * H + 64 * t115 * t170 + 64 * t173 * t170 - 16 * t102 * t120 - 12 * t8 * t108 * t104 - 8 * t7 * t181 * t110;
      double t193 = t93 * t8;
      double t197 = t9 * t9;
      double t212 = t93 * t21;
      double t215 = t14 * zz;
      double t218 = t9 * t3;
      double t221 = t101 * t3;
      double t224 = -4 * t8 * t181 * t115 + 8 * t98 * t7 * t9 + 2 * t193 * t14 * a + 2 * t150 * t197 * t7 + t193 * t197 * M + 8 * t95 * t2 - 4 * t21 * t181 * t173 - 12 * t21 * t108 * t103 * zz + t212 * t197 * H + 2 * t212 * t215 - 4 * t145 * t218 + 32 * t110 * t221;
      double t252 = a * t1 * t3;
      double t266 = 32 * t7 * t103 * t221 - 4 * t150 * t218 + 16 * t173 * t94 * rr + 12 * t8 * t103 * t108 * zz - 2 * t193 * t215 - 8 * t60 * t93 * t9 * t3 - 80 * H * rr * t121 + 64 * t60 * t103 * t101 * t3 - 16 * H * t93 * t252 - 16 * M * t93 * t252 - 64 * H * t101 * t103 * t75 + 8 * t136 * M * t46 * zz;
      return (3 * U * (t45 + t78) * t89 * t90 + 3 * U * (t141 + t185 + t224 + t266) * t89)  * (zz > -r1b) * (zz < a + r2b);
    }

    PetscScalar dzur_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double rr  = X[1];
      double zz  = X[0];

      double r1b = 1./H;
      double r2b = 1./M;
      double t6 = rr * rr;
      double t7 = H * t6;
      double t9 = M * t6;
      double t10 = 2 * a;
      double t11 = t7 + 4 * zz - t9 - t10;
      double t14 = 1 / (H + M);
      double t15 = t10 + t9 + t7;
      double t16 = t15 * t15;
      double t18 = 1 / t16 / t15;
      double t20 = chiTrunc(rr, l, eps);
      double t23 = dchiTrunc(rr, l, eps);
      return  (6 * U * (4 * H * rr + 4 * M * rr) * t11 * t14 * t18 * t20 + 6 * U * (-2 * t23 * a - t23 * M * t6 - t23 * H * t6) * t11 * t14 * t18 ) * (zz > -r1b) * (zz < a + r2b);
    }

    ///////////////////
    // Uz
    ///////////////////

    PetscScalar uz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double rr  = X[1];
      double zz  = X[0];

      double r1b = 1./H;
      double r2b = 1./M;
      double t1 = rr * rr;
      double t2 = H * t1;
      double t6 = H * H;
      double t8 = t1 * t1;
      double t9 = t8 * t1;
      double t21 = M * M;
      double t24 = H * t8;
      double t28 = a * M;
      double t31 = a * a;
      double t37 = zz * zz;
      double t53 = t6 * H * t9 + 4 * M * t9 * t6 + 8 * a * t6 * t8 - 2 * t6 * t8 * zz + 3 * H * t9 * t21 - 20 * t24 * M * zz + 36 * t24 * t28 + 60 * t2 * t31 - 88 * t2 * a * zz + 40 * t2 * t37 - 18 * t21 * t8 * zz - 24 * t28 * t1 * zz + 40 * t1 * t37 * M - 16 * a * t37 + 24 * t31 * zz;
      double t54 = chiTrunc(rr, l, eps);
      double t59 = pow((2 * a + M * t1 + t2), 2);
      double t60 = t59 * t59;
      return (U * (t2 + 2 * zz) * t53 * t54 / t60)  * (zz > -r1b) * (zz < a + r2b);
    }

    PetscScalar druz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double rr  = X[1];
      double zz  = X[0];
      double r1b = 1./H;
      double r2b = 1./M;
      double t1 = rr * rr;
      double t2 = t1 * rr;
      double t3 = zz * zz;
      double t4 = t3 * zz;
      double t5 = t2 * t4;
      double t6 = H * H;
      double t9 = M * M;
      double t12 = t1 * t1;
      double t13 = t12 * rr;
      double t14 = t13 * t3;
      double t15 = t6 * H;
      double t18 = t9 * M;
      double t21 = t13 * t15;
      double t22 = a * a;
      double t25 = t2 * t6;
      double t26 = t22 * a;
      double t29 = rr * t4;
      double t30 = H * a;
      double t33 = rr * t3;
      double t34 = M * t22;
      double t37 = a * M;
      double t44 = zz * t26;
      double t50 = t2 * t3;
      double t54 = -480 * t5 * t6 - 480 * t5 * t9 - 144 * t14 * t15 + 144 * t14 * t18 - 144 * t21 * t22 + 480 * t25 * t26 + 576 * t29 * t30 - 576 * t33 * t34 + 576 * t29 * t37 - 1152 * t33 * H * t22 + 576 * rr * H * t44 - 960 * t5 * H * M + 1440 * t50 * t6 * a;
      double t61 = zz * t22;
      double t64 = zz * a;
      double t67 = t12 * t2;
      double t68 = t67 * t6;
      double t69 = zz * t9;
      double t72 = t67 * t15;
      double t73 = M * zz;
      double t77 = zz * t18;
      double t82 = t13 * t6;
      double t85 = a * t9;
      double t92 = t61 * M;
      double t96 = t64 * t9;
      double t99 = t64 * M;
      double t102 = -144 * t14 * t6 * M + 144 * t14 * H * t9 - 1440 * t25 * t61 + 288 * t21 * t64 + 48 * t68 * t69 + 24 * t72 * t73 + 24 * t67 * H * t77 - 24 * t72 * t37 + 192 * t82 * t34 - 24 * t68 * t85 + 1440 * t50 * t30 * M - 480 * t2 * H * t92 - 336 * t13 * H * t96 - 48 * t82 * t99;
      double t106 = M * t1;
      double t107 = H * t1;
      double t108 = 2 * a + t106 + t107;
      double t109 = t108 * t108;
      double t110 = t109 * t109;
      double t112 = 1 / t110 / t108;
      double t113 = chiTrunc(rr, l, eps);
      double t116 = t6 * t6;
      double t118 = t12 * t12;
      double t119 = t118 * t1;
      double t131 = t6 * t12;
      double t141 = t12 * t1;
      double t144 = t15 * t141;
      double t156 = a * t4;
      double t161 = t22 * t3;
      double t166 = t116 * H * t119 - 64 * t22 * t4 + 96 * t26 * t3 + 10 * t116 * t118 * a + 5 * t116 * t119 * M + 80 * t131 * t4 + 80 * t9 * t12 * t4 + 7 * t9 * t119 * t15 + 76 * t22 * t15 * t141 + 36 * t144 * t3 - 36 * t18 * t141 * t3 + 3 * t18 * t119 * t6 + 120 * t26 * t6 * t12 + 128 * t156 * t106 + 128 * t156 * t107 - 48 * t161 * t106 - 336 * t161 * t107;
      double t167 = a * t3;
      double t173 = t141 * t6;
      double t176 = t118 * t6;
      double t179 = t118 * t15;
      double t193 = H * t12;
      double t214 = -120 * t131 * t167 - 120 * t85 * t12 * t3 + 132 * t34 * t173 + 42 * t85 * t176 + 52 * t37 * t179 - 72 * t144 * t64 - 24 * t176 * t69 - 12 * t179 * t73 - 12 * H * t118 * t77 + 288 * t107 * t44 + 160 * t193 * t4 * M + 36 * t173 * t3 * M - 36 * t9 * t141 * t3 * H + 240 * t193 * t92 + 24 * H * t141 * t96 - 48 * t173 * t99 - 240 * t193 * t167 * M;
      double t217 = dchiTrunc(rr, l, eps);
      return (U * (t54 + t102) * t112 * t113 + U * (t166 + t214) * t112 * t217 ) * (zz > -r1b) * (zz < a + r2b);
    }

    PetscScalar dzuz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double rr  = X[1];
      double zz  = X[0];
      double r1b = 1./H;
      double r2b = 1./M;

      double t1 = rr * rr;
      double t2 = H * t1;
      double t6 = t1 * t1;
      double t8 = a * M;
      double t14 = H * H;
      double t18 = t6 * t1;
      double t20 = M * M;
      double t24 = a * a;
      double t33 = zz * zz;
      double t43 = 32 * t2 * a * zz - 4 * H * t6 * t8 + 8 * t8 * t1 * zz - 6 * t14 * t6 * zz + H * t18 * t20 + M * t18 * t14 - 12 * t2 * t24 + 6 * a * t14 * t6 + 6 * t20 * t6 * zz - 20 * t1 * t33 * M - 20 * t2 * t33 - 8 * t24 * zz + 8 * a * t33;
      double t48 = pow((2 * a + M * t1 + t2), 2);
      double t49 = t48 * t48;
      double t51 = chiTrunc(rr, l, eps);
      return (-12 * U * t43 / t49 * t51)  * (zz > -r1b) * (zz < a + r2b);
    }



    // !!!!!!!!!!   OLD FUNCTIONS   !!!!!!!!!!

    // PetscScalar uz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    // {
    //   double mu = 1;
    //   double r  = X[1];
    //   double z  = X[0];

    //   double t1 = r * r;
    //   double t2 = H * t1;
    //   double t3 = 2 * z;
    //   double t5 = M * t1;
    //   double t6 = 2 * a;
    //   double t10 = t6 + t5 + t2;
    //   double t11 = t10 * t10;
    //   double t15 = chiTrunc(r, l, eps);
    //   return -12 * (t2 + t3) * (-t3 + t5 + t6) * U * r / t11 / t10 * t15;

    // }


    //Attention : r cordonnee radial, donc z en 2D et donc z = x .....
    PetscScalar ux_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double r  = X[1];
      double z  = X[0];

      double t1 = r * r;
      double t2 = H * t1;
      double t3 = M * a;
      double t4 = z * z;
      double t5 = t3 * t4;
      double t8 = a * a;
      double t9 = t8 * M;
      double t10 = t9 * z;
      double t13 = t1 * t1;
      double t14 = H * t13;
      double t15 = M * M;
      double t16 = t15 * a;
      double t17 = t16 * z;
      double t20 = H * H;
      double t21 = t20 * t13;
      double t22 = t3 * z;
      double t31 = t15 * t1;
      double t32 = a * t4;
      double t35 = t13 * t1;
      double t36 = H * t35;
      double t37 = t15 * M;
      double t41 = t20 * t1;
      double t45 = t20 * H;
      double t46 = t45 * t13;
      double t51 = t45 * t35;
      double t55 = t20 * t35;
      double t59 = t15 * t15;
      double t61 = t13 * t13;
      double t64 = t8 * t45;
      double t67 = t20 * t20;
      double t68 = a * t67;
      double t71 = t8 * a;
      double t74 = t15 * t61;
      double t89 = -24 * t51 * M * z - 48 * t55 * t15 * z - 5 * H * t59 * t61 + 96 * t64 * t13 + 8 * t68 * t35 - 32 * t41 * t71 + 4 * t74 * t45 - 4 * t37 * t61 * t20 + 5 * M * t61 * t67 - 32 * t71 * t15 * t1 - 24 * t8 * t37 * t13;
      double t94 = t4 * z;
      double t99 = H * a;
      double t100 = t99 * t94;
      double t102 = t3 * t94;
      double t109 = H * t8;
      double t110 = t109 * t4;
      double t112 = t9 * t4;
      double t114 = t8 * t8;
      double t118 = t67 * H;
      double t126 = a * t20;
      double t133 = M * t1;
      double t136 = t15 * t13;
      double t139 = t37 * t35;
      double t147 = t118 * t61 - 16 * M * t114 - t59 * M * t61 + 56 * t3 * t51 + 24 * t126 * t35 * t15 + 48 * t21 * t9 - 64 * H * t71 * t133 - 72 * t109 * t136 - 32 * t99 * t139 + 320 * t2 * M * t94 - 384 * t41 * t32;
      double t155 = pow((2 * a + t133 + t2), 2);
      double t156 = t155 * t155;
      double t158 = 1. / (H + M) / t156;
      double t159 = chiTrunc(r, l, eps);
      double t162 = d2chiTrunc(r, l, eps);
      double t163 = t162 * t13;
      double t166 = t61 * t1;
      double t167 = t162 * t166;
      double t168 = t3 * t45;
      double t173 = t162 * t35;
      double t174 = H * M;
      double t175 = t174 * t94;
      double t186 = t162 * t61;
      double t212 = t45 * t61;
      double t226 = 128 * t163 * t102 - 56 * t167 * t168 + 128 * t163 * t100 + 64 * t173 * t175 - 144 * t173 * t16 * t4 - 36 * t167 * t16 * t20 - 288 * t163 * t112 - 72 * t186 * t9 * t20 + 48 * t173 * t126 * t4 - 96 * t163 * t110 + 24 * t186 * t20 * M * t4 - 24 * t186 * H * t15 * t4 - 192 * t14 * t162 * t71 * z - 24 * H * t166 * t162 * t37 * z - 48 * t212 * t162 * a * z - 192 * t55 * t162 * t8 * z - 24 * t45 * t166 * t162 * M * z;
      double t232 = dchiTrunc(r, l, eps);
      double t233 = t232 * t1;
      double t236 = t232 * t13;
      double t241 = t232 * t61;
      double t244 = t232 * t8;
      double t245 = t244 * t4;
      double t250 = t232 * t166;
      double t251 = t67 * M;
      double t256 = t20 * t94;
      double t259 = t15 * t94;
      double t262 = t15 * t45;
      double t265 = t232 * t4;
      double t274 = t232 * t71;
      double t280 = t162 * t61 * t13;
      double t283 = -48 * t20 * t166 * t162 * t15 * z - 192 * t233 * t100 - 224 * t236 * t175 - 192 * t233 * t102 - 12 * t241 * t168 + 432 * t2 * t245 + 6 * t241 * t68 + 5 * t250 * t251 - 52 * t244 * t51 - 112 * t236 * t256 - 112 * t236 * t259 + 7 * t250 * t262 - 60 * t51 * t265 + 60 * t139 * t265 + 3 * t37 * t166 * t232 * t20 - 120 * t274 * t21 - 20 * t167 * t68 - 10 * t280 * t251;
      double t285 = t162 * t1;
      double t323 = H * t61;
      double t330 = 128 * t285 * t8 * t94 - 56 * t186 * t64 + 32 * t173 * t259 - 14 * t280 * t262 + 32 * t173 * t256 - 192 * t285 * t71 * t4 - 48 * t173 * t71 * t20 - 24 * t186 * t37 * t4 - 6 * t280 * t37 * t20 + 24 * t186 * t45 * t4 + 64 * t244 * t94 - 96 * t274 * t4 + t250 * t118 - 2 * t280 * t118 - 96 * t173 * a * t174 * t4 - 144 * t323 * t162 * t17 - 288 * t36 * t162 * t10;
      double t331 = t20 * t61;
      double t339 = t244 * z;
      double t343 = t232 * a;
      double t344 = t343 * z;
      double t350 = t343 * t4;
      double t353 = t232 * t15;
      double t357 = t232 * M;
      double t389 = -192 * t331 * t162 * t22 + 336 * t14 * t232 * t5 - 144 * t14 * M * t339 + 72 * t36 * t15 * t344 + 192 * t55 * t232 * t22 + 120 * t21 * t350 + 60 * t36 * t353 * t4 - 60 * t55 * t357 * t4 + 144 * t133 * t245 - 108 * M * t35 * t244 * t20 + 216 * t136 * t350 - 18 * t74 * t343 * t20 + 36 * t323 * t37 * t232 * z + 96 * t21 * t339 + 120 * t51 * t344 + 72 * t331 * t353 * z + 36 * t212 * t357 * z - 288 * t2 * t274 * z;
      return U * (-480 * t2 * t5 + 288 * t2 * t10 + 96 * t14 * t17 - 48 * t21 * t22 + 72 * t21 * M * t4 - 72 * t14 * t15 * t4 - 96 * t31 * t32 - 24 * t36 * t37 * z + 288 * t41 * t8 * z - 144 * t46 * a * z + t89 - 8 * a * t59 * t35 + 160 * t41 * t94 + 160 * t31 * t94 - 64 * t100 - 64 * t102 + 72 * t46 * t4 - 72 * t37 * t13 * t4 + 96 * t110 + 96 * t112 - 16 * H * t114 + t147) * t158 * t159 / 2 + U * (t226 + t283 + t330 + t389) * t158 / 2;
    }

    // PetscScalar dzuz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    // {
    //   double mu  = 1;
    //   double z   = X[0];
    //   double r   = X[1];

    //   double t1 = r * r;
    //   double t2 = t1 * t1;
    //   double t3 = z * z;
    //   double t4 = t2 * t3;
    //   double t14 = t2 * t1;
    //   double t15 = t14 * H;
    //   double t16 = a * M;
    //   double t22 = H * H;
    //   double t25 = M * M;
    //   double t29 = a * a;
    //   double t33 = H * t2;
    //   double t39 = t2 * t2;
    //   double t46 = 8 * M * t4 + 16 * t1 * a * t3 - 16 * t2 * a * M * z - 8 * t15 * t16 + 8 * H * t4 + 4 * t14 * z * t22 - 4 * t14 * t25 * z - 16 * t1 * t29 * z - 8 * t33 * t29 - 4 * t14 * a * t22 - 2 * t39 * H * t25 - 2 * t39 * M * t22;
    //   double t50 = H * t1;
    //   double t52 = pow(2 * a + M * t1 + t50, 2);
    //   double t53 = t52 * t52;
    //   double t54 = 0.1e1 / t53;
    //   double t55 = dchiTrunc(r, l, eps);
    //   double t90 = -20 * t1 * t3 * M + 8 * t16 * t1 * z + 8 * a * t3 - 8 * t29 * z + 32 * t50 * a * z - 4 * t33 * t16 - 6 * z * t22 * t2 + t15 * t25 + M * t14 * t22 - 12 * t50 * t29 + 6 * a * t22 * t2 + 6 * t25 * t2 * z - 20 * t3 * H * t1;
    //   double t92 = chiTrunc(r, l, eps);
    //   return 12 * U * t46 * t54 * t55 + 12 * U * t90 * t54 * t92;

    // }


    PetscScalar dxux_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[0];
      double r   = X[1];

      double t1 = r * r;
      double t2 = t1 * t1;
      double t3 = t2 * t1;
      double t4 = t3 * a;
      double t8 = t2 * M;
      double t9 = a * z;
      double t13 = H * H;
      double t16 = t2 * t2;
      double t17 = M * M;
      double t25 = a * a;
      double t30 = z * z;
      double t31 = t2 * t30;
      double t34 = t3 * t17;
      double t45 = 8 * t4 * H * M + 16 * t8 * t9 - 4 * t3 * z * t13 + 2 * t16 * t17 * H + 2 * t16 * M * t13 + 8 * t2 * H * t25 + 4 * t4 * t13 - 8 * t31 * H + 4 * t34 * z + 16 * t1 * t25 * z - 16 * t1 * a * t30 - 8 * t31 * M;
      double t48 = M * t1;
      double t49 = H * t1;
      double t51 = pow(2 * a + t48 + t49, 2);
      double t52 = t51 * t51;
      double t53 = 0.1e1 / t52;
      double t54 = dchiTrunc(r, l, eps);
      double t88 = -32 * t49 * t9 + 4 * a * H * t8 + 6 * z * t13 * t2 - t34 * H - M * t3 * t13 + 12 * t49 * t25 - 6 * a * t13 * t2 - 6 * t17 * t2 * z + 20 * t30 * H * t1 + 8 * t25 * z - 8 * a * t30 + 20 * t30 * M * t1 - 8 * t48 * t9;
      double t90 = chiTrunc(r, l, eps);
      return 12 * U * t45 * t53 * t54 + 12 * U * t88 * t53 * t90;


    }

    PetscScalar dxuz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double r   = X[1];
      double z   = X[0];

      double t2 = r * r;
      double t3 = chiTrunc(r, l, eps);
      double t5 = M * t2;
      double t6 = 2 * a;
      double t7 = H * t2;
      double t10 = t6 + t5 + t7;
      double t11 = t10 * t10;
      return 24 * U * r * t3 * (4 * z - t5 - t6 + t7) / t11 / t10;
    }

    PetscScalar dzux_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double r   = X[1];
      double z   = X[0];

      double t1 = r * U;
      double t2 = r * r;
      double t3 = t2 * M;
      double t4 = z * z;
      double t5 = t4 * z;
      double t6 = a * t5;
      double t9 = H * t2;
      double t12 = t2 * t2;
      double t13 = H * t12;
      double t14 = M * t5;
      double t17 = a * a;
      double t18 = t17 * t4;
      double t21 = H * H;
      double t22 = a * t21;
      double t23 = t12 * t4;
      double t26 = t12 * t2;
      double t33 = M * M;
      double t34 = a * t33;
      double t41 = t17 * a;
      double t45 = t12 * t12;
      double t46 = H * t45;
      double t47 = t33 * M;
      double t51 = 384 * t3 * t6 + 384 * t9 * t6 + 576 * t13 * t14 - 1248 * t9 * t18 - 528 * t22 * t23 + 120 * M * t26 * t21 * t4 + 96 * t3 * t18 - 336 * t34 * t23 - 120 * t33 * t26 * H * t4 + 1344 * t9 * t41 * z - 24 * t46 * t47 * z;
      double t52 = t21 * t45;
      double t53 = t33 * z;
      double t56 = t21 * H;
      double t57 = t56 * t45;
      double t58 = M * z;
      double t61 = t56 * t26;
      double t62 = a * z;
      double t65 = t21 * t12;
      double t66 = t17 * z;
      double t74 = a * M;
      double t78 = t26 * H;
      double t81 = t41 * H;
      double t82 = t12 * M;
      double t85 = t17 * M;
      double t86 = t26 * t21;
      double t89 = a * H;
      double t93 = -48 * t52 * t53 - 24 * t57 * t58 - 240 * t61 * t62 + 192 * t65 * t66 - 40 * a * t47 * t46 + 144 * t34 * t52 + 224 * t74 * t57 - 120 * t17 * t33 * t78 - 160 * t81 * t82 + 480 * t85 * t86 - 864 * t89 * t82 * t4;
      double t95 = t85 * z;
      double t98 = t34 * z;
      double t101 = t17 * t17;
      double t109 = t41 * t21;
      double t115 = t17 * t56;
      double t118 = t33 * t33;
      double t122 = t21 * t21;
      double t123 = a * t122;
      double t126 = t45 * t2;
      double t130 = 1248 * t13 * t95 + 240 * t78 * t98 - 80 * t3 * t101 - 80 * t9 * t101 - 80 * t41 * t33 * t12 + 448 * t109 * t12 - 40 * t17 * t47 * t26 + 320 * t115 * t26 - 10 * a * t118 * t45 + 50 * t123 * t45 - 5 * t118 * t126 * H;
      double t140 = t33 * t12;
      double t152 = t122 * H;
      double t155 = t17 * t5;
      double t157 = t41 * t4;
      double t161 = 8 * t47 * t126 * t21 + 32 * t33 * t126 * t56 + 25 * M * t126 * t122 + 288 * t140 * t5 + 288 * t65 * t5 + 120 * t61 * t4 - 120 * t47 * t26 * t4 - t118 * M * t126 + 5 * t152 * t126 - 384 * t155 + 576 * t157 - 32 * t101 * a;
      double t165 = 2 * a + t3 + t9;
      double t166 = t165 * t165;
      double t167 = t166 * t166;
      double t169 = 0.1e1 / t167 / t165;
      double t171 = dchiTrunc(r, l, eps);
      double t178 = t21 * t2;
      double t194 = t56 * t12;
      double t208 = -960 * t9 * t14 - 24 * t74 * t61 + 1440 * t178 * a * t4 - 144 * t65 * M * t4 + 144 * t13 * t33 * t4 - 24 * t86 * t34 + 192 * t85 * t65 - 1440 * t178 * t66 + 288 * t194 * t62 + 48 * t86 * t53 + 24 * t61 * t58 + 24 * H * t47 * t26 * z + 1440 * t9 * t74 * t4;
      double t209 = H * M;
      double t217 = t74 * z;
      double t220 = t89 * t5;
      double t227 = t74 * t5;
      double t237 = H * t17 * t4;
      double t239 = t85 * t4;
      double t245 = -480 * t209 * t2 * t17 * z - 336 * t89 * t140 * z - 48 * t65 * t217 + 576 * t220 - 480 * t178 * t5 - 480 * t33 * t2 * t5 + 576 * t227 - 144 * t194 * t17 - 144 * t194 * t4 + 144 * t47 * t12 * t4 - 1152 * t237 - 576 * t239 + 480 * t178 * t41 + 576 * t81 * z;
      double t248 = chiTrunc(r, l, eps);
      double t251 = d2chiTrunc(r, l, eps);
      double t252 = t2 * t251;
      double t255 = t45 * t251;
      double t258 = t26 * t251;
      double t263 = t45 * t12 * t251;
      double t270 = t126 * t251;
      double t302 = t12 * t251;
      double t305 = -256 * t252 * t155 + 112 * t255 * t115 - 64 * t258 * t33 * t5 + 28 * t263 * t33 * t56 - 64 * t258 * t21 * t5 + 40 * t270 * t123 + 20 * t263 * t122 * M - 48 * t255 * t56 * t4 + 48 * t255 * t47 * t4 + 12 * t263 * t47 * t21 + 384 * t252 * t157 + 96 * t258 * t109 + 192 * t258 * a * t209 * t4 + 576 * t78 * t251 * t95 + 288 * t46 * t251 * t98 + 384 * t52 * t251 * t217 - 256 * t302 * t227;
      double t367 = 112 * t270 * t74 * t56 - 256 * t302 * t220 - 128 * t258 * t209 * t5 + 192 * t302 * t237 - 96 * t258 * t22 * t4 + 48 * t255 * t33 * H * t4 - 48 * t255 * M * t21 * t4 + 576 * t302 * t239 + 144 * t255 * t85 * t21 + 288 * t258 * t34 * t4 + 72 * t270 * t34 * t21 + 384 * t13 * t251 * t41 * z + 384 * t86 * t251 * t17 * z + 96 * t57 * t251 * a * z + 48 * H * t126 * t251 * t47 * z + 96 * t21 * t126 * t251 * t33 * z + 48 * t56 * t126 * t251 * M * z + 4 * t263 * t152;
      return t1 * (t51 + t93 + t130 + t161) * t169 * t171 + t1 * (t208 + t245) * t169 * t248 + t1 * (t305 + t367) * t169;

    }
    PetscScalar DIVCART_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      return dxux_sing_normalMvt2D(X, H,  M, d, U, l, eps, param)
        +    dzuz_sing_normalMvt2D(X, H,  M, d, U, l, eps, param);
    }













    PetscScalar p_sing_withT_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t2 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t2);

      double t3 = chiTrunc(r, l, eps);
      double t11 = pow(t2 * H + t2 * M + 2 * a, 2);//std::cout << "U = " << U << std::endl;
      return  -12 / t11 / (H + M) * t3 * U * mu;

    }

      
    PetscScalar ur_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t2 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t2);
      double z  = X[2];

      double t1 = U * r;
      double t3 = H * t2;
      double t4 = M * t2;
      double t5 = 2 * a;
      double t6 = -t3 - t4 - t5;
      double t7 = 2 * z;
      double t8 = t4 + t5 - t7;
      double t11 = t3 + t7;
      double t12 = dchiTrunc(r, l, eps);
      double t14 = H + M;
      double t15 = 0.1e1 / t14;
      double t16 = -t6;
      double t17 = t16 * t16;
      double t19 = 0.1e1 / t17 / t16;
      double t27 = chiTrunc(r, l, eps);
      return  -3 * t1 * t11 * t12 * t15 * t19 * t6 * t8 - 6 * t1 * t11 * t14 * t15 * t19 * t27 * t8;

    }

    PetscScalar uz_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[2];
      double t4  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(t4);

      double t1 = H * H;
      double t2 = t1 * t1;
      double t3 = t2 * H;
      double t5 = t4 * t4;
      double t6 = t5 * t5;
      double t7 = t6 * t5;
      double t9 = t2 * M;
      double t12 = t1 * H;
      double t13 = M * M;
      double t14 = t12 * t13;
      double t17 = t13 * M;
      double t18 = t1 * t17;
      double t21 = t2 * a;
      double t22 = t6 * t4;
      double t25 = t12 * M;
      double t26 = a * t22;
      double t29 = t22 * z;
      double t32 = t1 * t13;
      double t37 = H * t17;
      double t40 = a * a;
      double t41 = t12 * t40;
      double t44 = t12 * a;
      double t45 = t6 * z;
      double t49 = z * z;
      double t52 = t1 * M;
      double t56 = a * t6;
      double t57 = t56 * z;
      double t60 = t6 * t49;
      double t63 = H * t13;
      double t66 = 12 * t12 * t6 * t49 - 36 * t52 * t40 * t6 - 7 * t14 * t7 - 3 * t18 * t7 - 10 * t21 * t22 - 28 * t25 * t26 - 12 * t25 * t29 - 18 * t32 * t26 - 24 * t32 * t29 - 12 * t37 * t29 - t3 * t7 - 28 * t41 * t6 - 24 * t44 * t45 - 96 * t52 * t57 + 12 * t52 * t60 - 72 * t63 * t57 - 5 * t9 * t7;
      double t72 = t40 * a;
      double t73 = t1 * t72;
      double t74 = t5 * t4;
      double t77 = t1 * t40;
      double t78 = t74 * z;
      double t81 = t1 * a;
      double t82 = t74 * t49;
      double t86 = t49 * z;
      double t89 = H * M;
      double t90 = t40 * t74;
      double t94 = a * t74;
      double t101 = t13 * a;
      double t107 = H * t72;
      double t108 = t5 * z;
      double t111 = H * t40;
      double t112 = t5 * t49;
      double t115 = H * a;
      double t116 = t5 * t86;
      double t119 = M * t40;
      double t122 = M * a;
      double t125 = t72 * t4;
      double t128 = t40 * t4;
      double t131 = 16 * t1 * t74 * t86 + 16 * t13 * t74 * t86 - 12 * t17 * t49 * t6 - 48 * t49 * t89 * t94 + 32 * t74 * t86 * t89 - 144 * t89 * t90 * z - 72 * t101 * t82 - 96 * t107 * t108 - 48 * t111 * t112 - 144 * t112 * t119 + 64 * t115 * t116 + 64 * t116 * t122 - 96 * t125 * t49 + 64 * t128 * t86 - 12 * t60 * t63 - 24 * t73 * t74 - 96 * t77 * t78 + 24 * t81 * t82;
      double t134 = d2chiTrunc(r, l, eps);
      double t136 = 0.1e1 / (H + M);
      double t142 = pow(H * t4 + M * t4 + 2 * a,2);
      double t143 = t142 * t142;
      double t144 = 0.1e1 / t143;
      double t178 = -24 * t12 * t49 * t74 + 24 * t17 * t49 * t74 + 48 * t52 * t94 * z - 2 * t21 * t6 + 12 * t25 * t45 - 20 * t25 * t56 + 24 * t32 * t45 - 18 * t32 * t56 + 12 * t37 * t45 - 40 * t41 * t74 + 48 * t44 * t78 - 24 * t52 * t82 - 72 * t52 * t90 + 24 * t63 * t82;
      double t186 = t40 * t5;
      double t190 = a * t5;
      double t201 = t4 * z;
      double t204 = t4 * t49;
      double t207 = t4 * t86;
      double t216 = -48 * t1 * t5 * t86 - 48 * t13 * t5 * t86 - 144 * t186 * t89 * z + 144 * t190 * t49 * t89 + 72 * t101 * t112 - 192 * t107 * t201 + 192 * t111 * t204 + 72 * t112 * t81 - 64 * t115 * t207 - 96 * t116 * t89 - 64 * t122 * t207 + 64 * t40 * t86 - 96 * t49 * t72 - 72 * t5 * t73;
      double t219 = dchiTrunc(r, l, eps);
      double t230 = t13 * t13;
      double t252 = t190 * z;
      double t263 = -5 * H * t230 * t6 - M * t230 * t6 - 8 * a * t230 * t74 + 24 * t12 * t49 * t5 - 48 * t108 * t44 + 24 * t112 * t52 - 24 * t112 * t63 + 4 * t14 * t6 - 4 * t18 * t6 - 72 * t186 * t63 + 8 * t21 * t74 + 32 * t25 * t94 + 48 * t252 * t52 + 96 * t252 * t63 + t3 * t6 - 32 * t37 * t94 + 48 * t41 * t5 + 5 * t6 * t9;
      double t296 = t40 * t40;
      double t309 = -192 * a * t4 * t49 * t89 + 64 * t1 * t4 * t86 + 192 * t128 * t89 * z - 32 * t13 * t4 * t72 + 64 * t13 * t4 * t86 - 24 * t17 * t40 * t5 - 24 * t17 * t49 * t5 - 16 * H * t296 - 16 * M * t296 + 96 * t111 * t49 - 64 * t115 * t86 + 96 * t119 * t49 - 64 * t122 * t86 - 64 * t125 * t89 + 192 * t201 * t77 - 192 * t204 * t81 + 128 * t207 * t89 - 32 * t4 * t73;
      double t313 = chiTrunc(r, l, eps);
      return U * (t131 + t66) * t134 * t136 * t144 / 2 + U * (t216 + t178) * t219 * t136 * t144 / 2 + U * (t309 + t263) * t136 * t144 * t313 / 2;

    }

    PetscScalar drur_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[2];
      double t4  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(t4);

      double t1 = H * H;
      double t2 = t1 * H;
      double t3 = t2 * M;
      double t5 = t4 * t4;
      double t6 = t5 * t5;
      double t7 = t6 * t4;
      double t10 = M * M;
      double t11 = t1 * t10;
      double t14 = t10 * M;
      double t15 = H * t14;
      double t18 = t2 * a;
      double t24 = t1 * M;
      double t25 = a * t6;
      double t28 = t6 * z;
      double t31 = H * t10;
      double t39 = a * a;
      double t40 = t1 * t39;
      double t41 = t5 * t4;
      double t44 = t1 * a;
      double t45 = t41 * z;
      double t48 = 4 * t14 * t6 * z - 4 * t2 * t6 * z + 4 * t11 * t7 + 2 * t15 * t7 + 4 * t18 * t6 + 16 * t24 * t25 - 4 * t24 * t28 + 12 * t31 * t25 + 4 * t31 * t28 + 2 * t3 * t7 + 16 * t40 * t41 - 8 * t44 * t45;
      double t50 = z * z;
      double t53 = H * M;
      double t57 = a * t41;
      double t64 = t10 * a;
      double t70 = t39 * a;
      double t71 = H * t70;
      double t74 = H * t39;
      double t75 = t5 * z;
      double t78 = H * a;
      double t79 = t5 * t50;
      double t82 = M * t39;
      double t85 = M * a;
      double t91 = t39 * t4;
      double t94 = -8 * t1 * t41 * t50 - 8 * t10 * t41 * t50 + 24 * t53 * t39 * t41 + 32 * t70 * t4 * z - 16 * t53 * t41 * t50 + 16 * t53 * t57 * z + 24 * t64 * t45 + 16 * t71 * t5 - 32 * t91 * t50 + 16 * t74 * t75 + 48 * t82 * t75 - 32 * t78 * t79 - 32 * t85 * t79;
      double t97 = d2chiTrunc(r, l, eps);
      double t99 = 0.1e1 / (H + M);
      double t105 = pow(H * t4 + M * t4 + 2 * a, 2);
      double t106 = t105 * t105;
      double t107 = 0.1e1 / t106;
      double t136 = -10 * t14 * t41 * z + 10 * t2 * t41 * z - 6 * t11 * t6 - 3 * t15 * t6 - 10 * t18 * t41 + 10 * t24 * t45 - 16 * t24 * t57 - 3 * t3 * t6 - 10 * t31 * t45 - 6 * t31 * t57 - 8 * t40 * t5 - 20 * t44 * t75;
      double t143 = a * t5;
      double t156 = t4 * z;
      double t159 = t4 * t50;
      double t170 = 28 * t1 * t5 * t50 + 28 * t10 * t5 * t50 - 56 * t53 * t143 * z + 12 * t53 * t39 * t5 - 72 * t74 * t156 - 24 * t82 * t156 + 48 * t78 * t159 + 48 * t85 * t159 - 16 * t39 * t50 + 24 * t71 * t4 + 56 * t53 * t79 - 36 * t64 * t75 + 16 * t70 * z;
      double t173 = dchiTrunc(r, l, eps);
      double t201 = 12 * t14 * t5 * z - 12 * t2 * t5 * z + 4 * t11 * t41 + 4 * t24 * t143 - 8 * t31 * t143 + 2 * t15 * t41 + 12 * t18 * t5 - 12 * t24 * t75 + 2 * t3 * t41 + 12 * t31 * t75 - 24 * t40 * t4;
      double t228 = 80 * t53 * a * t4 * z - 40 * t1 * t4 * t50 - 40 * t10 * t4 * t50 + 64 * t44 * t156 + 16 * t64 * t156 - 80 * t53 * t159 + 16 * t78 * t50 + 16 * t85 * t50 - 24 * t53 * t91 - 16 * t74 * z - 16 * t82 * z;
      double t232 = chiTrunc(r, l, eps);
      return  3 * U * (t94 + t48) * t97 * t99 * t107 + 3 * U * (t170 + t136) * t173 * t99 * t107 + 3 * U * (t228 + t201) * t99 * t107 * t232;
    }

    PetscScalar dzur_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[2];
      double t2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(t2);

      double t3 = H * t2;
      double t4 = M * t2;
      double t5 = 2 * a;
      double t9 = H + M;
      double t11 = t3 + t4 + t5;
      double t12 = t11 * t11;
      double t17 = chiTrunc(r, l, eps);
      double t20 = dchiTrunc(r, l, eps);
      return 6 * r * U * (t3 - t4 - t5 + 4 * z) / t9 / t12 / t11 * (-t11 * t20 + 2 * t17 * t9);
    }

    PetscScalar druz_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[2];
      double t5  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(t5);

      double t1 = r * U;
      double t2 = H * H;
      double t3 = t2 * t2;
      double t4 = t3 * t2;
      double t6 = t5 * t5;
      double t7 = t6 * t5;
      double t8 = t6 * t6;
      double t9 = t8 * t7;
      double t11 = t3 * H;
      double t12 = t11 * M;
      double t15 = M * M;
      double t16 = t3 * t15;
      double t19 = t2 * H;
      double t20 = t15 * M;
      double t21 = t19 * t20;
      double t24 = t15 * t15;
      double t25 = t2 * t24;
      double t28 = t11 * a;
      double t29 = t8 * t6;
      double t32 = t3 * M;
      double t33 = a * t29;
      double t36 = t29 * z;
      double t39 = t19 * t15;
      double t44 = t2 * t20;
      double t49 = 6 * t12 * t9 + 12 * t16 * t9 + 10 * t21 * t9 + 3 * t25 * t9 + 12 * t28 * t29 + 48 * t32 * t33 + 12 * t32 * t36 + 60 * t39 * t33 + 24 * t44 * t33 + 36 * t39 * t36 + 36 * t44 * t36 + t4 * t9;
      double t50 = H * t24;
      double t53 = a * a;
      double t54 = t3 * t53;
      double t55 = t8 * t5;
      double t58 = t3 * a;
      double t63 = z * z;
      double t66 = t19 * M;
      double t67 = t53 * t55;
      double t70 = a * t55;
      double t71 = t70 * z;
      double t74 = t55 * t63;
      double t77 = t2 * t15;
      double t82 = H * t20;
      double t90 = t53 * a;
      double t91 = t19 * t90;
      double t94 = 12 * t24 * t55 * t63 - 12 * t3 * t55 * t63 + 24 * t58 * t55 * z + 12 * t50 * t36 + 48 * t54 * t55 + 120 * t66 * t67 + 144 * t66 * t71 - 24 * t66 * t74 + 72 * t77 * t67 + 216 * t77 * t71 + 96 * t82 * t71 + 24 * t82 * t74 + 80 * t91 * t8;
      double t96 = t19 * t53;
      double t97 = t8 * z;
      double t100 = t19 * a;
      double t101 = t8 * t63;
      double t105 = t63 * z;
      double t108 = t2 * M;
      double t112 = t53 * t8;
      double t113 = t112 * z;
      double t116 = t8 * t105;
      double t119 = H * t15;
      double t122 = a * t8;
      double t128 = t20 * a;
      double t134 = t53 * t53;
      double t135 = t2 * t134;
      double t138 = t2 * t90;
      double t139 = t7 * z;
      double t142 = -16 * t19 * t8 * t105 - 16 * t20 * t8 * t105 + 96 * t108 * t90 * t8 + 144 * t119 * t122 * t63 - 48 * t100 * t101 + 96 * t128 * t101 + 432 * t108 * t113 - 48 * t108 * t116 + 288 * t119 * t113 - 48 * t119 * t116 + 48 * t135 * t7 + 288 * t138 * t139 + 144 * t96 * t97;
      double t143 = t2 * a;
      double t144 = t7 * t105;
      double t147 = H * M;
      double t148 = t90 * t7;
      double t152 = t53 * t7;
      double t156 = a * t7;
      double t160 = t15 * t53;
      double t161 = t7 * t63;
      double t164 = t15 * a;
      double t167 = H * t134;
      double t168 = t6 * z;
      double t171 = H * t90;
      double t172 = t6 * t63;
      double t175 = H * t53;
      double t176 = t6 * t105;
      double t179 = M * t90;
      double t182 = M * t53;
      double t185 = t134 * t5;
      double t188 = t90 * t5;
      double t191 = -192 * t147 * t156 * t105 + 384 * t147 * t148 * z + 288 * t147 * t152 * t63 - 128 * t188 * t105 - 96 * t143 * t144 - 96 * t164 * t144 + 288 * t160 * t161 + 192 * t167 * t168 + 192 * t171 * t172 + 384 * t179 * t172 - 192 * t175 * t176 - 192 * t182 * t176 + 192 * t185 * t63;
      double t195 = d3chiTrunc(r, l, eps);
      double t197 = 0.1e1 / (H + M);
      double t202 = H * t5 + M * t5 + 2 * a;
      double t203 = t202 * t202;
      double t204 = t203 * t203;
      double t206 = 0.1e1 / t204 / t202;
      double t231 = 12 * t12 * t29 + 24 * t16 * t29 + 20 * t21 * t29 + 6 * t25 * t29 + 24 * t28 * t55 + 2 * t4 * t29 + 120 * t32 * t70 + 168 * t39 * t70 + 72 * t44 * t70 + 144 * t54 * t8 - 48 * t58 * t97;
      double t241 = t122 * z;
      double t258 = 64 * t19 * t7 * t105 - 24 * t24 * t8 * t63 + 24 * t3 * t8 * t63 - 96 * t100 * t161 + 48 * t66 * t101 - 48 * t82 * t101 + 432 * t66 * t112 + 288 * t77 * t112 + 144 * t77 * t241 + 96 * t82 * t241 + 352 * t91 * t7;
      double t262 = t152 * z;
      double t265 = t156 * t63;
      double t285 = 64 * t20 * t7 * t105 + 192 * t108 * t144 + 480 * t108 * t148 + 576 * t108 * t262 - 288 * t108 * t265 + 192 * t119 * t144 + 576 * t119 * t262 - 288 * t119 * t265 - 96 * t128 * t161 + 288 * t135 * t6 + 576 * t138 * t168;
      double t286 = t2 * t53;
      double t291 = t90 * t6;
      double t295 = t53 * t6;
      double t299 = a * t6;
      double t305 = t5 * z;
      double t308 = t5 * t63;
      double t317 = 384 * t147 * t299 * t105 + 1152 * t147 * t291 * z - 576 * t147 * t295 * t63 - 256 * t90 * t105 + 384 * t134 * t63 + 192 * t143 * t176 + 192 * t164 * t176 + 768 * t167 * t305 - 384 * t171 * t308 - 576 * t286 * t172 + 384 * t179 * t308;
      double t321 = d2chiTrunc(r, l, eps);
      double t332 = t24 * M;
      double t355 = 6 * H * t332 * t55 + 10 * t332 * a * t8 + t24 * t15 * t55 - 6 * t12 * t55 - 50 * t32 * t122 - 40 * t39 * t122 + 40 * t44 * t122 + 50 * t50 * t122 + 96 * t58 * t139 - 9 * t16 * t55 + 9 * t25 * t55 - 10 * t28 * t8 - t4 * t55 - 88 * t54 * t7;
      double t365 = t156 * z;
      double t389 = -160 * t19 * t6 * t105 + 40 * t24 * t53 * t7 + 48 * t24 * t7 * t63 - 48 * t3 * t7 * t63 + 432 * t100 * t172 - 64 * t66 * t152 + 144 * t77 * t152 + 160 * t82 * t152 - 96 * t66 * t161 + 96 * t82 * t161 - 384 * t96 * t168 - 288 * t77 * t365 - 192 * t82 * t365 + 32 * t91 * t6;
      double t393 = t295 * z;
      double t396 = t299 * t63;
      double t423 = -160 * t20 * t6 * t105 + 80 * t20 * t90 * t6 - 480 * t108 * t176 + 384 * t108 * t291 - 1056 * t108 * t393 + 912 * t108 * t396 - 480 * t119 * t176 + 240 * t119 * t291 - 672 * t119 * t393 + 528 * t119 * t396 + 48 * t128 * t172 + 368 * t135 * t5 - 960 * t138 * t305 + 576 * t286 * t308;
      double t424 = t5 * t105;
      double t432 = t53 * t5;
      double t436 = a * t5;
      double t447 = t134 * a;
      double t462 = -128 * t147 * t436 * t105 + 80 * t15 * t134 * t5 - 384 * t147 * t188 * z + 192 * t147 * t432 * t63 + 32 * H * t447 + 32 * M * t447 + 512 * t175 * t105 + 512 * t182 * t105 - 64 * t143 * t424 + 160 * t147 * t185 - 384 * t160 * t308 - 64 * t164 * t424 + 384 * t167 * z - 960 * t171 * t63 - 576 * t179 * t63;
      double t466 = dchiTrunc(r, l, eps);
      double t483 = t299 * z;
      double t504 = 192 * t19 * t5 * t105 - 48 * t24 * t6 * t63 + 48 * t3 * t6 * t63 - 672 * t100 * t308 - 288 * t108 * t188 - 96 * t58 * t168 + 96 * t66 * t172 - 96 * t82 * t172 - 96 * t66 * t295 - 144 * t77 * t295 + 768 * t96 * t305 + 288 * t77 * t483 + 192 * t82 * t483 - 288 * t91 * t5 + 48 * t54 * t6;
      double t505 = t432 * z;
      double t508 = t436 * t63;
      double t543 = -768 * t147 * a * t105 + 192 * t20 * t5 * t105 + 1152 * t147 * t53 * t63 - 384 * t147 * t90 * z - 384 * t143 * t105 - 384 * t164 * t105 + 576 * t108 * t424 + 960 * t108 * t505 - 1248 * t108 * t508 + 576 * t119 * t424 + 192 * t119 * t505 - 480 * t119 * t508 + 96 * t128 * t308 - 384 * t138 * z + 384 * t160 * t63 + 768 * t286 * t63;
      double t547 = chiTrunc(r, l, eps);
      return -t1 * (t191 + t142 + t94 + t49) * t195 * t197 * t206 - t1 * (t317 + t285 + t258 + t231) * t321 * t197 * t206 - t1 * (t462 + t423 + t389 + t355) * t466 * t197 * t206 - t1 * (t543 + t504) * t197 * t206 * t547;

    }









    PetscScalar ux_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
        return ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]/r;
      else
        return ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param);
    }

    PetscScalar uy_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
        return ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[1]/r;
      else
        return 0.;
    }









    PetscScalar dxux_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
	return drur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[0]/r2 
	  +      ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[1]*X[1]/r2/r;
      else
        return 0;
    }

    PetscScalar dyux_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
	return drur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[1]/r2 
	  -      ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[1]/r2/r;
      else
        return 0;
    }

    PetscScalar dzux_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
	return dzur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]/r;
      else 
	return dzur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param);
    }

    PetscScalar dxuy_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
	return drur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[1]/r2 
	  -      ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[1]/r2/r;
      else
        return 0;
    }

    PetscScalar dyuy_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
	return drur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[1]*X[1]/r2 
	  +      ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[0]/r2/r;
      else
        return 0;
    }

    PetscScalar dzuy_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
	return dzur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[1]/r;
      else
        return 0;
    }

    PetscScalar dxuz_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
	return druz_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]/r;
      else
        return 0;
    }


    PetscScalar dyuz_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);
        
      if(r2!=0)
	return druz_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[1]/r;
      else
        return 0;
    }

    PetscScalar dzuz_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[2];
      double t4  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(t4);

      double t1 = H * H;
      double t2 = t1 * H;
      double t3 = t2 * M;
      double t5 = t4 * t4;
      double t6 = t5 * t5;
      double t7 = t6 * t4;
      double t9 = M * M;
      double t10 = t1 * t9;
      double t13 = t9 * M;
      double t14 = H * t13;
      double t16 = t2 * a;
      double t22 = t1 * M;
      double t23 = a * t6;
      double t26 = t6 * z;
      double t29 = H * t9;
      double t37 = a * a;
      double t38 = t1 * t37;
      double t39 = t5 * t4;
      double t42 = t1 * a;
      double t43 = t39 * z;
      double t46 = 2 * t13 * t6 * z - 2 * t2 * t6 * z + 2 * t10 * t7 + t14 * t7 + 2 * t16 * t6 + 8 * t22 * t23 - 2 * t22 * t26 + 6 * t29 * t23 + 2 * t29 * t26 + t3 * t7 + 8 * t38 * t39 - 4 * t42 * t43;
      double t48 = z * z;
      double t51 = H * M;
      double t55 = a * t39;
      double t62 = t9 * a;
      double t68 = t37 * a;
      double t69 = H * t68;
      double t72 = H * t37;
      double t73 = t5 * z;
      double t76 = H * a;
      double t77 = t5 * t48;
      double t80 = M * t37;
      double t83 = M * a;
      double t89 = t37 * t4;
      double t92 = -4 * t1 * t39 * t48 + 12 * t51 * t37 * t39 - 8 * t51 * t39 * t48 - 4 * t9 * t39 * t48 + 16 * t68 * t4 * z + 8 * t51 * t55 * z + 12 * t62 * t43 - 16 * t89 * t48 + 8 * t69 * t5 + 8 * t72 * t73 + 24 * t80 * t73 - 16 * t76 * t77 - 16 * t83 * t77;
      double t95 = d2chiTrunc(r, l, eps);
      double t97 = 0.1e1 / (H + M);
      double t103 = pow( H * t4 + M * t4 + 2 * a, 2);
      double t104 = t103 * t103;
      double t105 = 0.1e1 / t104;
      double t131 = 12 * t1 * t5 * t48 - 4 * t13 * t39 * z + 4 * t2 * t39 * z - 2 * t10 * t6 - t14 * t6 - 4 * t16 * t39 + 4 * t22 * t43 - 4 * t22 * t55 - 4 * t29 * t43 - t3 * t6 - 12 * t42 * t73;
      double t135 = a * t5;
      double t148 = t4 * z;
      double t151 = t4 * t48;
      double t160 = -24 * t51 * t135 * z + 12 * t51 * t37 * t5 + 12 * t9 * t5 * t48 - 32 * t72 * t148 + 16 * t76 * t151 + 16 * t83 * t151 - 16 * t37 * t48 + 16 * t69 * t4 + 24 * t51 * t77 - 12 * t62 * t73 + 16 * t68 * z;
      double t163 = dchiTrunc(r, l, eps);
      double t209 = 32 * t51 * a * t4 * z - 16 * t1 * t4 * t48 + 4 * t13 * t5 * z - 4 * t2 * t5 * z - 16 * t9 * t4 * t48 - 4 * t22 * t135 - 8 * t29 * t135 + 32 * t42 * t148 - 32 * t51 * t151 + 4 * t16 * t5 - 4 * t22 * t73 + 4 * t29 * t73 - 16 * t38 * t4 + 16 * t76 * t48 + 16 * t83 * t48 - 16 * t51 * t89 - 16 * t72 * z - 16 * t80 * z;
      double t212 = chiTrunc(r, l, eps);
      return -6 * U * (t92 + t46) * t95 * t97 * t105 - 6 * U * (t160 + t131) * t163 * t97 * t105 - 6 * U * t209 * t97 * t105 * t212;
    }

    PetscScalar DIVCYLIND_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      return drur_sing_normalMvt3D(X, H,  M, d, U, l, eps, param)
        +      ur_sing_normalMvt3D(X, H,  M, d, U, l, eps, param)/r
        +    dzuz_sing_normalMvt3D(X, H,  M, d, U, l, eps, param);
    }

    PetscScalar DIVCART_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      return dxux_sing_normalMvt3D(X, H,  M, d, U, l, eps, param)
        +    dyuy_sing_normalMvt3D(X, H,  M, d, U, l, eps, param)
        +    dzuz_sing_normalMvt3D(X, H,  M, d, U, l, eps, param);
    }
  }
}
#endif
