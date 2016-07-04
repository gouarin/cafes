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

#ifndef CAFES_PARTICLE_SINGULARITY_UANDPTANG_HPP_INCLUDED
#define CAFES_PARTICLE_SINGULARITY_UANDPTANG_HPP_INCLUDED

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


    PetscScalar p_sing_withT_tangMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t4 = X[1]*X[1];
      double r  = X[1];

      double t9 = pow(H * t4 + M * t4 + 2 * d, 2);
      return  -8 * mu * U * (H - M) * r / t9 / (H + M) * chiTrunc(r, l, eps);
    }

    PetscScalar ux_sing_tangMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t1 = X[1]*X[1];
      double r  = X[1];
      double z  = X[0];

      double t2 = t1 * H;
      double t6 = (t2 + 2 * z) * U * r;
      double t7 = H * H;
      double t8 = t7 * H;
      double t9 = M * M;
      double t10 = t8 * t9;
      double t11 = t1 * t1;
      double t12 = t11 * t1;
      double t15 = t9 * M;
      double t16 = t7 * t15;
      double t19 = t9 * t9;
      double t20 = H * t19;
      double t23 = t7 * t9;
      double t24 = d * t11;
      double t27 = t11 * z;
      double t30 = H * t15;
      double t38 = d * d;
      double t39 = t8 * t38;
      double t42 = t8 * d;
      double t43 = t1 * z;
      double t47 = z * z;
      double t50 = -48 * t8 * t1 * t47 - 24 * t19 * t11 * z - 48 * t39 * t1 + 12 * t10 * t12 + 24 * t16 * t12 + 12 * t20 * t12 + 96 * t23 * t24 - 24 * t23 * t27 + 96 * t30 * t24 - 48 * t30 * t27 + 96 * t42 * t43;
      double t51 = t7 * M;
      double t52 = t38 * t1;
      double t56 = d * t1 * z;
      double t59 = t1 * t47;
      double t62 = H * t9;
      double t72 = t38 * d;
      double t73 = t72 * t7;
      double t75 = t7 * t38;
      double t78 = t7 * d;
      double t81 = t9 * t38;
      double t84 = t9 * d;
      double t87 = 48 * t15 * t1 * t47 + 96 * t78 * t47 - 96 * t84 * t47 + 144 * t51 * t52 - 96 * t51 * t56 - 48 * t51 * t59 + 144 * t62 * t52 - 192 * t62 * t56 + 48 * t62 * t59 - 192 * t75 * z + 96 * t81 * z + 96 * t73;
      double t89 = chiTrunc(r, l, eps);
      double t94 = pow(M * t1 + 2 * d + t2, 2);
      double t95 = t94 * t94;
      double t99 = 0.1e1 / t95 / (H + M);
      double t102 = t7 * t7;
      double t103 = t102 * H;
      double t104 = t11 * t11;
      double t106 = t102 * M;
      double t117 = t102 * d;
      double t123 = t8 * M;
      double t124 = d * t12;
      double t127 = t12 * z;
      double t151 = t38 * t11;
      double t154 = -3 * t19 * M * t104 - 24 * t19 * d * t12 + 4 * t102 * t12 * z + 40 * t8 * t11 * t47 + 24 * t19 * t12 * z - 10 * t10 * t104 + t103 * t104 + t106 * t104 - 22 * t16 * t104 - 15 * t20 * t104 + 72 * t39 * t11 + 8 * t117 * t12 + 12 * t123 * t124 + 4 * t123 * t127 - 104 * t23 * t124 - 132 * t30 * t124 + 20 * t23 * t127 + 44 * t30 * t127 - 72 * t51 * t151 - 64 * t42 * t27;
      double t155 = t24 * z;
      double t158 = t11 * t47;
      double t170 = t15 * d;
      double t183 = H * M;
      double t195 = t38 * t38;
      double t198 = H * t72;
      double t201 = H * t38;
      double t206 = M * t72;
      double t209 = M * t38;
      double t212 = -336 * t183 * t72 * t1 - 96 * t9 * t72 * t1 + 336 * t183 * t52 * z - 48 * H * t195 - 48 * M * t195 + 192 * t198 * z - 96 * t201 * t47 - 96 * t206 * z + 96 * t209 * t47 + 32 * t78 * t59 - 32 * t84 * t59;
      double t215 = dchiTrunc(r, l, eps);
      double t219 = t104 * t1;
      double t235 = d * t104;
      double t238 = t104 * z;
      double t259 = t38 * t12;
      double t262 = -4 * t102 * t104 * z - 12 * t19 * t104 * z - 16 * t8 * t12 * t47 + 4 * t10 * t219 + 2 * t103 * t219 + 20 * t117 * t104 + 8 * t106 * t219 + 56 * t39 * t12 + 36 * t123 * t235 + 8 * t123 * t238 + 8 * t42 * t127 - 8 * t16 * t219 - 6 * t20 * t219 - 20 * t23 * t235 + 16 * t23 * t238 - 36 * t30 * t235 - 8 * t30 * t238 + 16 * t51 * t259;
      double t263 = t124 * z;
      double t266 = t12 * t47;
      double t304 = -48 * t183 * t72 * t11 + 16 * t15 * t12 * t47 + 64 * t183 * t151 * z + 48 * t73 * t11 - 72 * t170 * t127 - 64 * t78 * t158 + 64 * t84 * t158 + 96 * t198 * t43 - 64 * t201 * t59 - 96 * t206 * t43 + 64 * t209 * t59 - 72 * t62 * t259 + 72 * t51 * t263 - 8 * t62 * t263 - 16 * t51 * t266 + 16 * t62 * t266 + 80 * t75 * t27 - 144 * t81 * t27;
      double t306 = d2chiTrunc(r, l, eps);
      return  -t6 * (t87 + t50) * t89 * t99 / 6. - t6 * (-72 * t15 * t38 * t11 - 40 * t15 * t11 * t47 + 96 * t73 * t1 - 360 * t62 * t151 + 72 * t51 * t155 + 208 * t62 * t155 + 40 * t51 * t158 - 40 * t62 * t158 + 72 * t170 * t27 - 48 * t75 * t43 + t154 + t212) * t215 * t99 / 6. - t6 * (t304 + t262) * t306 * t99 / 6.;
    }

    PetscScalar uz_sing_tangMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t3 = X[1]*X[1];
      double r  = X[1];
      double z  = X[0];

      double t1 = H * H;
      double t2 = t1 * t1;
      double t4 = t3 * t3;
      double t5 = t4 * t3;
      double t7 = t1 * H;
      double t8 = t7 * M;
      double t11 = M * M;
      double t12 = t11 * M;
      double t13 = H * t12;
      double t16 = t11 * t11;
      double t18 = t7 * d;
      double t24 = t1 * M;
      double t25 = d * t4;
      double t28 = H * t11;
      double t37 = d * d;
      double t38 = t1 * t37;
      double t41 = -6 * d * t12 * t4 + 16 * t12 * t4 * z + 16 * t4 * t7 * z + 4 * t13 * t5 - t16 * t5 - 10 * t18 * t4 + t2 * t5 + 2 * t24 * t25 - 2 * t25 * t28 + 4 * t3 * t38 - 4 * t5 * t8;
      double t43 = t3 * z;
      double t47 = z * z;
      double t50 = H * M;
      double t61 = t11 * d;
      double t67 = t37 * d;
      double t70 = H * t37;
      double t73 = d * H;
      double t78 = d * M;
      double t81 = 48 * d * t3 * t50 * z - 16 * d * t1 * t43 + 24 * t1 * t3 * t47 - 12 * t11 * t3 * t37 - 24 * t11 * t3 * t47 - 24 * t3 * t37 * t50 - 8 * H * t67 - 8 * M * t67 + 32 * t43 * t61 - 16 * t47 * t73 + 16 * t47 * t78 + 32 * t70 * z;
      double t87 = H * t3 + M * t3 + 2 * d;
      double t88 = t87 * t87;
      double t90 = 0.1e1 / t88 / t87;
      double t92 = 0.1e1 / (H + M);
      double t94 = chiTrunc(r, l, eps);
      double t97 = t4 * t4;
      double t107 = d * t5;
      double t110 = t5 * z;
      double t139 = t3 * t47;
      double t147 = -32 * M * t37 * t43 - 16 * t1 * t4 * t47 + 16 * t11 * t4 * t47 - 8 * t12 * t5 * z + 32 * t25 * t50 * z - 16 * t37 * t4 * t50 - 32 * t4 * t61 * z - 8 * t5 * t7 * z + 8 * t107 * t24 - 16 * t107 * t28 + 8 * t110 * t24 + 8 * t110 * t28 - 4 * t13 * t97 - 32 * t139 * t73 + 32 * t139 * t78 + 8 * t18 * t5 + 16 * t38 * t4 + 32 * t43 * t70 + 4 * t8 * t97;
      double t149 =  dchiTrunc(r, l, eps);
      return  U * (t81 + t41) * t90 * t92 * t94 / 2 + U * t147 * t149 * t90 * t92 / 2;

    }

    PetscScalar dxux_sing_tangMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t5 = X[1]*X[1];
      double r  = X[1];
      double z  = X[0];

      double t1 = U * r;
      double t2 = H * H;
      double t3 = t2 * t2;
      double t4 = t3 * d;
      double t6 = t5 * t5;
      double t12 = t2 * H;
      double t13 = t12 * M;
      double t14 = d * t6;
      double t17 = t6 * z;
      double t20 = M * M;
      double t21 = t20 * M;
      double t22 = H * t21;
      double t27 = t20 * t20;
      double t31 = d * d;
      double t32 = t12 * t31;
      double t35 = t12 * d;
      double t36 = t5 * z;
      double t40 = z * z;
      double t43 = t2 * M;
      double t44 = t31 * t5;
      double t47 = -48 * t12 * t5 * t40 - 16 * t27 * t6 * z - 16 * t3 * t6 * z - 16 * t13 * t14 - 16 * t13 * t17 + 32 * t22 * t14 - 16 * t22 * t17 - 48 * t32 * t5 + 96 * t35 * t36 + 16 * t4 * t6 + 48 * t43 * t44;
      double t49 = d * t5 * z;
      double t52 = t5 * t40;
      double t55 = H * t20;
      double t65 = t31 * d;
      double t66 = t65 * t2;
      double t68 = t2 * t31;
      double t71 = t2 * d;
      double t74 = t20 * t31;
      double t77 = t20 * d;
      double t80 = 48 * t21 * t5 * t40 + 96 * t71 * t40 - 96 * t77 * t40 - 64 * t43 * t49 - 48 * t43 * t52 + 64 * t55 * t44 - 160 * t55 * t49 + 48 * t55 * t52 - 128 * t68 * z + 64 * t74 * z + 32 * t66;
      double t83 = chiTrunc(r, l, eps);
      double t88 = pow(H * t5 + M * t5 + 2 * d, 2);
      double t89 = t88 * t88;
      double t90 = 0.1e1 / t89;
      double t93 = 0.1e1 / (H + M);
      double t97 = t6 * t6;
      double t99 = t3 * M;
      double t101 = H * t27;
      double t105 = t6 * t5;
      double t111 = d * t105;
      double t114 = t105 * z;
      double t134 = t31 * t6;
      double t137 = t14 * z;
      double t140 = t6 * t40;
      double t143 = t3 * H * t97 - t27 * M * t97 - 8 * t27 * d * t105 + 16 * t27 * t105 * z + 16 * t3 * t105 * z + 40 * t12 * t6 * t40 - t101 * t97 - 8 * t4 * t105 + 16 * t13 * t111 - 32 * t22 * t111 + 16 * t13 * t114 + 16 * t22 * t114 + 32 * t43 * t134 + 48 * t43 * t137 + 40 * t43 * t140 - 32 * t35 * t17 + 16 * t32 * t6 + t99 * t97;
      double t153 = t21 * d;
      double t165 = H * M;
      double t177 = t31 * t31;
      double t180 = H * t65;
      double t183 = H * t31;
      double t188 = M * t65;
      double t191 = M * t31;
      double t194 = 256 * t165 * t44 * z - 128 * t165 * t65 * t5 - 32 * t20 * t65 * t5 - 24 * t21 * t31 * t6 - 40 * t21 * t6 * t40 - 16 * H * t177 - 16 * M * t177 - 120 * t55 * t134 + 128 * t55 * t137 - 40 * t55 * t140 + 48 * t153 * t17 + 128 * t180 * z - 96 * t183 * t40 - 64 * t188 * z + 96 * t191 * t40 - 64 * t68 * t36 + 64 * t66 * t5 + 32 * t71 * t52 - 32 * t77 * t52;
      double t197 = dchiTrunc(r, l, eps);
      double t201 = t97 * t5;
      double t217 = d * t97;
      double t220 = t2 * t20;
      double t238 = t31 * t105;
      double t241 = t111 * z;
      double t244 = -16 * t12 * t105 * t40 + 4 * t12 * t20 * t201 - 4 * t2 * t21 * t201 + 16 * t220 * t97 * z - 8 * t27 * t97 * z - 8 * t3 * t97 * z - 4 * t101 * t201 + 32 * t32 * t105 - 16 * t35 * t114 + 24 * t13 * t217 + 4 * t99 * t201 - 24 * t22 * t217 - 8 * t220 * t217 + 16 * t43 * t238 + 48 * t43 * t241 + 8 * t4 * t97;
      double t245 = t105 * t40;
      double t283 = 16 * t21 * t105 * t40 + 64 * t165 * t134 * z - 32 * t165 * t65 * t6 - 48 * t153 * t114 - 64 * t71 * t140 + 64 * t77 * t140 + 32 * t68 * t17 - 96 * t74 * t17 + 64 * t180 * t36 - 64 * t183 * t52 - 64 * t188 * t36 + 64 * t191 * t52 - 48 * t55 * t238 + 16 * t55 * t241 - 16 * t43 * t245 + 16 * t55 * t245 + 32 * t66 * t6;
      double t286 = d2chiTrunc(r, l, eps);
      return  -t1 * (t80 + t47) * t83 * t90 * t93 - t1 * (t194 + t143) * t197 * t90 * t93 - t1 * (t283 + t244) * t286 * t90 * t93;


    }

    PetscScalar dzux_sing_tangMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t6 = X[1]*X[1];
      double r  = X[1];
      double z  = X[0];

      double t1 = H * H;
      double t2 = t1 * t1;
      double t3 = t2 * H;
      double t4 = M * M;
      double t5 = t3 * t4;
      double t7 = t6 * t6;
      double t8 = t7 * t7;
      double t9 = t8 * t6;
      double t12 = t4 * M;
      double t13 = t2 * t12;
      double t16 = t1 * H;
      double t17 = t4 * t4;
      double t18 = t16 * t17;
      double t21 = t17 * M;
      double t22 = t1 * t21;
      double t25 = t2 * t4;
      double t26 = d * t8;
      double t29 = t16 * t12;
      double t32 = t1 * t17;
      double t35 = d * d;
      double t36 = t3 * t35;
      double t37 = t7 * t6;
      double t40 = t3 * d;
      double t41 = t37 * z;
      double t45 = z * z;
      double t48 = t2 * M;
      double t49 = t35 * t37;
      double t52 = t37 * t45;
      double t55 = t16 * t4;
      double t58 = 144 * t3 * t37 * t45 + 36 * t13 * t9 + 36 * t18 * t9 + 12 * t22 * t9 + 120 * t25 * t26 + 240 * t29 * t26 + 120 * t32 * t26 + 144 * t36 * t37 - 288 * t40 * t41 - 288 * t48 * t49 + 288 * t48 * t52 + 480 * t55 * t49 + 12 * t5 * t9;
      double t59 = d * t37;
      double t60 = t59 * z;
      double t65 = t1 * t12;
      double t72 = H * t17;
      double t80 = t35 * d;
      double t81 = t2 * t80;
      double t84 = t2 * t35;
      double t85 = t7 * z;
      double t88 = t2 * d;
      double t89 = t7 * t45;
      double t93 = t45 * z;
      double t96 = t16 * M;
      double t97 = t80 * t7;
      double t100 = 480 * t2 * t7 * t93 + 144 * t21 * t37 * t45 + 912 * t65 * t49 + 144 * t55 * t52 + 144 * t65 * t52 + 288 * t72 * t52 + 288 * t55 * t60 - 576 * t65 * t60 - 576 * t72 * t60 - 960 * t81 * t7 + 2400 * t84 * t85 - 1920 * t88 * t89 + 960 * t96 * t97;
      double t102 = t35 * t7;
      double t103 = t102 * z;
      double t106 = d * t7;
      double t107 = t106 * t45;
      double t110 = t7 * t93;
      double t113 = t1 * t4;
      double t120 = H * t12;
      double t125 = t17 * d;
      double t131 = t35 * t35;
      double t132 = t16 * t131;
      double t135 = t16 * t80;
      double t136 = t6 * z;
      double t139 = t16 * t35;
      double t140 = t6 * t45;
      double t143 = -480 * t17 * t7 * t93 - 3360 * t113 * t103 - 960 * t96 * t103 + 3360 * t113 * t107 + 1920 * t120 * t107 - 960 * t96 * t107 - 960 * t120 * t110 + 960 * t96 * t110 + 1440 * t113 * t97 - 480 * t125 * t89 + 576 * t132 * t6 - 3072 * t135 * t136 + 4416 * t139 * t140;
      double t144 = t16 * d;
      double t145 = t6 * t93;
      double t148 = t1 * M;
      double t149 = t80 * t6;
      double t150 = t149 * z;
      double t154 = t35 * t6 * t45;
      double t158 = d * t6 * t93;
      double t161 = H * t4;
      double t168 = t12 * t35;
      double t171 = t12 * d;
      double t174 = t1 * t131;
      double t177 = t1 * t80;
      double t180 = t1 * t35;
      double t183 = t4 * t80;
      double t186 = t4 * t35;
      double t189 = -1344 * t168 * t140 - 1920 * t144 * t145 + 1920 * t171 * t145 + 384 * t148 * t150 + 1536 * t148 * t154 - 1920 * t148 * t158 + 2304 * t161 * t150 - 4224 * t161 * t154 + 1920 * t161 * t158 + 384 * t174 * z - 768 * t177 * t45 + 384 * t180 * t93 + 384 * t183 * t45 - 384 * t186 * t93;
      double t196 = H * t6 + M * t6 + 2 * d;
      double t197 = t196 * t196;
      double t198 = t197 * t197;
      double t200 = 0.1e1 / t198 / t196;
      double t202 = 0.1e1 / (H + M);
      double t204 = chiTrunc(r, l, eps);
      double t207 = t2 * t16;
      double t208 = t8 * t7;
      double t211 = t2 * t1;
      double t212 = t211 * M;
      double t223 = t17 * t4;
      double t224 = H * t223;
      double t227 = t211 * d;
      double t233 = t3 * M;
      double t234 = d * t9;
      double t238 = t9 * z;
      double t251 = H * t21;
      double t261 = t8 * z;
      double t264 = -6 * t223 * t9 * z + 12 * t233 * t238 - 72 * t25 * t234 - 90 * t251 * t234 - 240 * t29 * t234 - 246 * t32 * t234 + 6 * t25 * t238 - 12 * t251 * t238 - 6 * t32 * t238 + 348 * t40 * t261 - 24 * t36 * t8;
      double t269 = t35 * t8;
      double t272 = t26 * z;
      double t275 = t8 * t45;
      double t297 = t21 * d;
      double t312 = t80 * t37;
      double t315 = t49 * z;
      double t318 = -432 * t2 * t37 * t93 - 144 * t21 * t8 * t45 - 60 * t297 * t261 - 360 * t72 * t269 + 516 * t72 * t272 - 288 * t72 * t275 + 480 * t96 * t312 + 96 * t96 * t315 + 720 * t81 * t37 - 1152 * t84 * t41 + 1344 * t88 * t52;
      double t321 = t59 * t45;
      double t324 = t37 * t93;
      double t341 = t17 * t35;
      double t358 = t131 * t7;
      double t361 = t97 * z;
      double t364 = t102 * t45;
      double t367 = t106 * t93;
      double t372 = 432 * t17 * t37 * t93 + 480 * t144 * t110 + 48 * t125 * t52 + 1584 * t132 * t7 - 1728 * t135 * t85 + 384 * t139 * t89 - 2880 * t148 * t358 + 5376 * t148 * t361 - 2976 * t148 * t364 + 480 * t148 * t367 - 720 * t161 * t358;
      double t380 = t12 * t80;
      double t387 = t131 * d;
      double t388 = t1 * t387;
      double t397 = H * M;
      double t401 = -288 * t397 * t387 * t6 - 480 * t171 * t110 + 3744 * t174 * t136 - 5376 * t177 * t140 + 2496 * t180 * t145 - 864 * t161 * t361 - 1536 * t161 * t364 - 480 * t161 * t367 + 1824 * t168 * t89 - 480 * t380 * t85 - 288 * t388 * t6;
      double t409 = t4 * t131;
      double t416 = H * t387;
      double t419 = H * t131;
      double t422 = H * t80;
      double t425 = M * t387;
      double t428 = M * t131;
      double t431 = M * t80;
      double t434 = -3264 * t397 * t131 * t6 * z + 3264 * t397 * t149 * t45 - 480 * t409 * t136 + 2112 * t183 * t140 - 2496 * t186 * t145 - 192 * t416 * z + 768 * t419 * t45 - 384 * t422 * t93 - 192 * t425 * z - 384 * t428 * t45 + 384 * t431 * t93;
      double t439 = dchiTrunc(r, l, eps);
      double t443 = t8 * t37;
      double t463 = d * t208;
      double t467 = t208 * z;
      double t491 = -12 * t223 * t208 * z + 96 * t233 * t467 - 24 * t40 * t238 - 72 * t25 * t463 + 156 * t25 * t467 - 60 * t251 * t463 - 96 * t251 * t467 - 936 * t29 * t463 - 636 * t32 * t463 - 156 * t32 * t467 + 672 * t36 * t9;
      double t496 = t35 * t9;
      double t499 = t234 * z;
      double t502 = t9 * t45;
      double t538 = t80 * t8;
      double t541 = t269 * z;
      double t544 = 192 * t2 * t8 * t93 + 72 * t21 * t9 * t45 - 120 * t297 * t238 + 672 * t84 * t261 - 384 * t88 * t275 - 240 * t72 * t496 - 1176 * t72 * t499 + 168 * t72 * t502 - 480 * t96 * t538 + 3264 * t96 * t541 + 1440 * t81 * t8;
      double t547 = t26 * t45;
      double t550 = t8 * t93;
      double t583 = t131 * t37;
      double t586 = t312 * z;
      double t589 = t49 * t45;
      double t592 = t59 * t93;
      double t597 = -192 * t17 * t8 * t93 + 192 * t125 * t275 + 960 * t132 * t37 + 3264 * t135 * t41 - 1728 * t139 * t52 + 384 * t144 * t324 - 2208 * t148 * t583 + 1344 * t148 * t586 + 2880 * t148 * t589 + 384 * t148 * t592 - 480 * t161 * t583;
      double t645 = -5376 * t397 * t358 * z + 5376 * t397 * t97 * t45 + 768 * t186 * t110 - 384 * t416 * t136 - 384 * t425 * t136 + 2688 * t419 * t140 - 1920 * t428 * t140 - 1536 * t422 * t145 + 1536 * t431 * t145 - 2304 * t183 * t89 - 960 * t409 * t85;
      double t650 = d2chiTrunc(r, l, eps);
      double t654 = t8 * t8;
      double t669 = d * t443;
      double t672 = t443 * z;
      double t691 = -8 * t13 * t654 - 28 * t18 * t654 + 4 * t207 * t654 + 192 * t36 * t208 + 20 * t212 * t654 - 12 * t22 * t654 + 48 * t227 * t443 + 144 * t233 * t669 + 48 * t233 * t672 + 48 * t25 * t669 + 96 * t25 * t672 - 48 * t251 * t672 - 144 * t29 * t669 - 96 * t32 * t669 - 96 * t32 * t672 + 96 * t40 * t467 + 24 * t5 * t654;
      double t695 = t35 * t208;
      double t698 = t463 * z;
      double t701 = t208 * t45;
      double t732 = -64 * t2 * t9 * t93 - 48 * t21 * t208 * t45 - 48 * t3 * t208 * t45 + 576 * t84 * t238 + 288 * t48 * t695 + 480 * t48 * t698 - 48 * t48 * t701 - 192 * t88 * t502 - 192 * t55 * t695 + 288 * t55 * t698 + 96 * t55 * t701 - 288 * t65 * t695 - 480 * t65 * t698 + 96 * t65 * t701 - 384 * t72 * t698 - 48 * t72 * t701 + 320 * t81 * t9;
      double t734 = t80 * t9;
      double t737 = t496 * z;
      double t740 = t234 * t45;
      double t743 = t9 * t93;
      double t772 = t538 * z;
      double t775 = -192 * t148 * t131 * t8 + 64 * t17 * t9 * t93 - 384 * t113 * t734 - 576 * t113 * t737 + 576 * t113 * t740 - 1152 * t120 * t737 - 192 * t120 * t740 + 128 * t120 * t743 - 384 * t125 * t502 + 192 * t132 * t8 + 1152 * t135 * t261 - 384 * t144 * t550 + 384 * t148 * t772 + 64 * t96 * t734 + 1152 * t96 * t737 + 192 * t96 * t740 - 128 * t96 * t743;
      double t779 = t26 * t93;
      double t814 = 1152 * t148 * t269 * t45 + 768 * t397 * t312 * t45 - 768 * t397 * t583 * z - 512 * t422 * t110 + 512 * t431 * t110 - 384 * t148 * t779 - 1536 * t161 * t772 + 384 * t161 * t779 - 1152 * t168 * t275 + 384 * t171 * t550 + 768 * t174 * t41 + 768 * t177 * t52 - 768 * t180 * t324 - 1536 * t183 * t52 + 768 * t186 * t324 + 768 * t419 * t89 - 768 * t428 * t89;
      double t818 = d3chiTrunc(r, l, eps);
      return -U * (t189 + t143 + t100 + t58) * t200 * t202 * t204 / 6 - U * (42 * t233 * t234 + 3 * t207 * t208 + 6 * t212 * t208 - 3 * t5 * t208 - 24 * t13 * t208 - 39 * t18 * t208 - 30 * t22 * t208 - 9 * t224 * t208 + 30 * t227 * t9 + 528 * t96 * t321 - 864 * t96 * t324 - 3168 * t113 * t312 + 2928 * t113 * t315 - 2928 * t113 * t321 - 720 * t120 * t312 + 1440 * t120 * t315 - 2064 * t120 * t321 + 864 * t120 * t324 - 240 * t341 * t41 + 408 * t48 * t269 + 60 * t48 * t272 - 288 * t48 * t275 - 480 * t55 * t269 - 288 * t55 * t272 - 144 * t55 * t275 - 1272 * t65 * t269 + 576 * t65 * t272 - 144 * t65 * t275 + t434 + t372 + t264 + t401 + 6 * t211 * t9 * z - 144 * t3 * t8 * t45 + t318) * t439 * t200 * t202 / 6 - U * (t491 + t645 - 8256 * t161 * t586 + 4032 * t161 * t589 - 384 * t161 * t592 - 960 * t380 * t41 - 576 * t168 * t52 - 384 * t171 * t324 - 192 * t388 * t7 + 3264 * t174 * t85 - 768 * t180 * t110 + t597 + 12 * t211 * t208 * z + 72 * t3 * t9 * t45 + t544 - 480 * t120 * t538 - 4800 * t120 * t541 + 1344 * t120 * t547 - 384 * t120 * t550 - 480 * t341 * t261 + 192 * t96 * t547 + 384 * t96 * t550 - 3552 * t113 * t538 - 1728 * t113 * t541 + 1728 * t113 * t547 + 420 * t233 * t463 + 12 * t207 * t443 + 54 * t212 * t443 + 42 * t5 * t443 - 84 * t13 * t443 - 144 * t18 * t443 - 66 * t22 * t443 - 6 * t224 * t443 + 132 * t227 * t208 + 912 * t48 * t496 + 840 * t48 * t499 + 168 * t48 * t502 - 1776 * t55 * t496 + 816 * t55 * t499 + 144 * t55 * t502 - 2256 * t65 * t496 - 1104 * t65 * t499 + 144 * t65 * t502 - 192 * t397 * t387 * t7) * t650 * t200 * t202 / 6 - U * (t814 + t775 + t732 + t691) * t818 * t200 * t202 / 6;

    }

    PetscScalar dxuz_sing_tangMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t3 = X[1]*X[1];
      double r  = X[1];
      double z  = X[0];

      double t1 = H * H;
      double t2 = t1 * H;
      double t4 = t3 * t3;
      double t7 = M * M;
      double t8 = t7 * M;
      double t17 = H * M;
      double t21 = t7 * d;
      double t27 = d * d;
      double t28 = H * t27;
      double t30 = d * H;
      double t33 = d * M;
      double t41 = H * t3 + M * t3 + 2 * d;
      double t42 = t41 * t41;
      double t44 = 0.1e1 / t42 / t41;
      double t46 = 0.1e1 / (H + M);
      double t48 = chiTrunc(r, l, eps);
      double t51 = t4 * t3;
      double t71 = t3 * z;
      double t79 = -H * t7 * t51 - t1 * M * t51 + 4 * M * t27 * t3 - 4 * t17 * d * t4 + 4 * t1 * t4 * z - 4 * t7 * t4 * z + t2 * t51 + 4 * t21 * t4 - 4 * t28 * t3 + 8 * t30 * t71 - 8 * t33 * t71 + t8 * t51;
      double t81 = dchiTrunc(r, l, eps);
      return  -4 * U * (2 * t1 * d * t3 - 6 * t17 * d * t3 - 6 * t1 * t3 * z + 6 * t7 * t3 * z - 2 * t2 * t4 - 4 * t21 * t3 + 4 * t30 * z - 4 * t33 * z - 2 * t8 * t4 - 4 * t28) * t44 * t46 * t48 - 4 * U * t79 * t81 * t44 * t46;

    }

    PetscScalar dzuz_sing_tangMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t5 = X[1]*X[1];
      double r  = X[1];
      double z  = X[0];

      double t1 = U * r;
      double t2 = H * H;
      double t3 = t2 * t2;
      double t4 = t3 * d;
      double t6 = t5 * t5;
      double t12 = t2 * H;
      double t13 = t12 * M;
      double t14 = d * t6;
      double t17 = t6 * z;
      double t20 = M * M;
      double t21 = t20 * M;
      double t22 = H * t21;
      double t27 = t20 * t20;
      double t31 = d * d;
      double t32 = t12 * t31;
      double t35 = t12 * d;
      double t36 = t5 * z;
      double t40 = z * z;
      double t43 = t2 * M;
      double t44 = t31 * t5;
      double t47 = -48 * t12 * t5 * t40 - 16 * t27 * t6 * z - 16 * t3 * t6 * z - 16 * t13 * t14 - 16 * t13 * t17 + 32 * t22 * t14 - 16 * t22 * t17 - 48 * t32 * t5 + 96 * t35 * t36 + 16 * t4 * t6 + 48 * t43 * t44;
      double t49 = d * t5 * z;
      double t52 = t5 * t40;
      double t55 = H * t20;
      double t65 = t31 * d;
      double t66 = t65 * t2;
      double t68 = t2 * t31;
      double t71 = t2 * d;
      double t74 = t20 * t31;
      double t77 = t20 * d;
      double t80 = 48 * t21 * t5 * t40 + 96 * t71 * t40 - 96 * t77 * t40 - 64 * t43 * t49 - 48 * t43 * t52 + 64 * t55 * t44 - 160 * t55 * t49 + 48 * t55 * t52 - 128 * t68 * z + 64 * t74 * z + 32 * t66;
      double t87 = pow(H * t5 + M * t5 + 2 * d, 2);
      double t88 = t87 * t87;
      double t89 = 0.1e1 / t88;
      double t91 = 0.1e1 / (H + M);
      double t93 = chiTrunc(r, l, eps);
      double t97 = t6 * t6;
      double t99 = t3 * M;
      double t101 = H * t27;
      double t105 = t6 * t5;
      double t111 = d * t105;
      double t114 = t105 * z;
      double t134 = t31 * t6;
      double t137 = t14 * z;
      double t140 = t6 * t40;
      double t143 = t3 * H * t97 - t27 * M * t97 - 8 * t27 * d * t105 + 16 * t27 * t105 * z + 16 * t3 * t105 * z + 40 * t12 * t6 * t40 - t101 * t97 - 8 * t4 * t105 + 16 * t13 * t111 - 32 * t22 * t111 + 16 * t13 * t114 + 16 * t22 * t114 + 32 * t43 * t134 + 48 * t43 * t137 + 40 * t43 * t140 - 32 * t35 * t17 + 16 * t32 * t6 + t99 * t97;
      double t153 = t21 * d;
      double t165 = H * M;
      double t177 = t31 * t31;
      double t180 = H * t65;
      double t183 = H * t31;
      double t188 = M * t65;
      double t191 = M * t31;
      double t194 = 256 * t165 * t44 * z - 128 * t165 * t65 * t5 - 32 * t20 * t65 * t5 - 24 * t21 * t31 * t6 - 40 * t21 * t6 * t40 - 16 * H * t177 - 16 * M * t177 - 120 * t55 * t134 + 128 * t55 * t137 - 40 * t55 * t140 + 48 * t153 * t17 + 128 * t180 * z - 96 * t183 * t40 - 64 * t188 * z + 96 * t191 * t40 - 64 * t68 * t36 + 64 * t66 * t5 + 32 * t71 * t52 - 32 * t77 * t52;
      double t197 = dchiTrunc(r, l, eps);
      double t201 = t97 * t5;
      double t217 = d * t97;
      double t220 = t2 * t20;
      double t238 = t31 * t105;
      double t241 = t111 * z;
      double t244 = -16 * t12 * t105 * t40 + 4 * t12 * t20 * t201 - 4 * t2 * t21 * t201 + 16 * t220 * t97 * z - 8 * t27 * t97 * z - 8 * t3 * t97 * z - 4 * t101 * t201 + 32 * t32 * t105 - 16 * t35 * t114 + 24 * t13 * t217 + 4 * t99 * t201 - 24 * t22 * t217 - 8 * t220 * t217 + 16 * t43 * t238 + 48 * t43 * t241 + 8 * t4 * t97;
      double t245 = t105 * t40;
      double t283 = 16 * t21 * t105 * t40 + 64 * t165 * t134 * z - 32 * t165 * t65 * t6 - 48 * t153 * t114 - 64 * t71 * t140 + 64 * t77 * t140 + 32 * t68 * t17 - 96 * t74 * t17 + 64 * t180 * t36 - 64 * t183 * t52 - 64 * t188 * t36 + 64 * t191 * t52 - 48 * t55 * t238 + 16 * t55 * t241 - 16 * t43 * t245 + 16 * t55 * t245 + 32 * t66 * t6;
      double t286 = d2chiTrunc(r, l, eps);
      return  t1 * (t80 + t47) * t89 * t91 * t93 + t1 * (t194 + t143) * t197 * t89 * t91 + t1 * (t283 + t244) * t286 * t89 * t91;

    }

    PetscScalar DIVCART_tangMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      return dxux_sing_tangMvt2D(X, H,  M, d, U, l, eps, param)
        +    dzuz_sing_tangMvt2D(X, H,  M, d, U, l, eps, param);
    }











    PetscScalar p_sing_withT_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t3 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t3);
      double cosTheta;

      if(r!=0)
        cosTheta = X[0]/r;
      else
        cosTheta = 1.;

      double t4  = chiTrunc(r, l, eps);
      double t11 = pow(H * t3 + M * t3 + 2 * d, 2);
      double P   = (0.48e2 / 0.5e1 * r * (H - M) * t4 / t11 / (H + M));

      //std::cout << "inside p_sing_withT_tangMvt3D " << chiTrunc(r, l, eps) << " " << H << " " << M << " " << d << " " << U << " " << cosTheta << std::endl;

      return mu*U*P*cosTheta/2.;
    }



      
    PetscScalar U0_tang3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double t3 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t3);
      double z  = X[2];

      double t1 = H * H;
      double t2 = t1 * t1;
      double t4 = t3 * t3;
      double t5 = t4 * t3;
      double t8 = t1 * H;
      double t9 = t8 * M;
      double t12 = M * M;
      double t13 = t12 * M;
      double t14 = H * t13;
      double t17 = t12 * t12;
      double t26 = t1 * M;
      double t27 = d * t4;
      double t30 = t4 * z;
      double t33 = H * t12;
      double t38 = t13 * d;
      double t44 = d * d;
      double t48 = 30 * t8 * d * t4 + 60 * t1 * t44 * t3 - 56 * t13 * t4 * z - 56 * t8 * t4 * z + 8 * t14 * t5 - 5 * t17 * t5 + 5 * t2 * t5 + 26 * t26 * t27 - 24 * t26 * t30 - 2 * t33 * t27 - 24 * t33 * t30 + 26 * t38 * t4 - 8 * t9 * t5;
      double t49 = t1 * d;
      double t50 = t3 * z;
      double t54 = z * z;
      double t57 = H * M;
      double t65 = t12 * t44;
      double t74 = t44 * d;
      double t77 = H * t44;
      double t80 = H * d;
      double t85 = M * t44;
      double t88 = M * d;
      double t91 = -208 * d * t3 * t57 * z + 16 * d * t12 * t50 + 72 * t1 * t3 * t54 - 72 * t12 * t3 * t54 + 104 * t3 * t44 * t57 + 40 * H * t74 + 40 * M * t74 - 4 * t3 * t65 - 128 * t49 * t50 - 48 * t54 * t80 + 48 * t54 * t88 - 32 * t77 * z - 128 * t85 * z;
      double t94 = 0.1e1 / (H + M);
      double t99 = H * t3 + M * t3 + 2 * d;
      double t100 = t99 * t99;
      double t102 = 0.1e1 / t100 / t99;
      double t103 = chiTrunc(r, l, eps);
      double t106 = t4 * t4;
      double t114 = d * t5;
      double t117 = t5 * z;
      double t147 = t3 * t54;
      double t154 = -48 * t1 * t4 * t54 + 48 * t12 * t4 * t54 + 24 * t13 * t5 * z - 96 * t27 * t57 * z + 48 * t4 * t44 * t57 + 24 * t5 * t8 * z - 12 * t106 * t14 + 12 * t106 * t9 + 48 * t114 * t26 - 24 * t114 * t33 - 24 * t117 * t26 - 24 * t117 * t33 - 96 * t147 * t80 + 96 * t147 * t88 + 96 * t30 * t49 - 24 * t38 * t5 - 48 * t4 * t65 + 96 * t50 * t77 - 96 * t50 * t85;
      double t155 = dchiTrunc(r, l, eps);

      return -(t91 + t48) * t94 * t102 * t103 / 5 - t154 * t155 * t102 * t94 / 5;
    }

    PetscScalar V0_tang3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double t3 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t3);
      double z  = X[2];

      double t1 = H * H;
      double t4 = t3 * t3;
      double t10 = M * M;
      double t23 = H * M;
      double t36 = d * d;
      double t42 = z * z;
      double t52 = -16 * H * d * z + 5 * H * t1 * t4 - 11 * H * t10 * t4 - 64 * M * d * z + 11 * M * t1 * t4 - 5 * M * t10 * t4 + 20 * d * t1 * t3 - 12 * d * t10 * t3 + 32 * d * t23 * t3 - 8 * t1 * t3 * z - 8 * t10 * t3 * z - 64 * t23 * t3 * z + 20 * H * t36 - 24 * H * t42 + 20 * M * t36 + 24 * M * t42;
      double t53 = chiTrunc(r, l, eps);
      double t61 = pow(H * t3 + M * t3 + 2 * d, 2);

      return t52 * t53 / (H + M) / t61 / 5;

    }

    PetscScalar W0_tang3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double t1 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t1);
      double z  = X[2];

      double t2 = M * t1;
      double t4 = t2 + 2 * z;
      double t5 = H * H;
      double t6 = t5 * t5;
      double t7 = t6 * M;
      double t8 = t1 * t1;
      double t9 = t8 * t1;
      double t12 = t5 * H;
      double t13 = M * M;
      double t14 = t12 * t13;
      double t17 = t13 * M;
      double t18 = t5 * t17;
      double t24 = t12 * M;
      double t25 = d * t8;
      double t28 = t8 * z;
      double t31 = t5 * t13;
      double t36 = H * t17;
      double t41 = t12 * d;
      double t42 = t1 * z;
      double t46 = z * z;
      double t49 = t5 * M;
      double t50 = d * d;
      double t51 = t50 * t1;
      double t55 = d * t1 * z;
      double t58 = -32 * t12 * t1 * t46 + 16 * t6 * t8 * z - 40 * t14 * t9 - 20 * t18 * t9 - 136 * t24 * t25 + 56 * t24 * t28 - 160 * t31 * t25 - 24 * t36 * t25 + 64 * t31 * t28 + 24 * t36 * t28 - 32 * t41 * t42 - 208 * t49 * t51 + 224 * t49 * t55 - 20 * t7 * t9;
      double t59 = t1 * t46;
      double t62 = H * t13;
      double t69 = t17 * t50;
      double t72 = t17 * d;
      double t78 = t5 * t50;
      double t81 = t5 * d;
      double t84 = H * M;
      double t85 = t50 * d;
      double t91 = t13 * t85;
      double t93 = t13 * t50;
      double t96 = t13 * d;
      double t99 = 32 * t17 * t1 * t46 + 32 * t84 * t50 * z + 32 * t69 * t1 - 64 * t72 * t42 + 128 * t81 * t46 - 128 * t96 * t46 - 32 * t49 * t59 - 224 * t62 * t51 + 192 * t62 * t55 + 32 * t62 * t59 - 128 * t78 * z - 32 * t84 * t85 + 256 * t93 * z - 128 * t91;
      double t106 = pow(t1 * H + 2 * d + t2, 2);
      double t107 = t106 * t106;
      double t108 = 0.1e1 / t107;
      double t110 = 0.1e1 / (H + M);
      double t112 = chiTrunc(r, l, eps);
      double t116 = t8 * t8;
      double t125 = t13 * t13;
      double t126 = H * t125;
      double t129 = t125 * M;
      double t138 = d * t9;
      double t141 = t9 * z;
      double t152 = t125 * d;
      double t166 = t50 * t8;
      double t169 = 5 * t6 * H * t116 + 40 * t6 * d * t9 + 32 * t12 * t8 * t46 + 120 * t12 * t50 * t8 - 6 * t125 * t9 * z - 22 * t6 * t9 * z - 3 * t126 * t116 - 2 * t129 * t116 + 38 * t14 * t116 + 16 * t18 * t116 + 26 * t7 * t116 + 202 * t24 * t138 + 174 * t31 * t138 - 2 * t36 * t138 - 56 * t24 * t141 - 52 * t31 * t141 - 24 * t36 * t141 - 14 * t152 * t9 + 516 * t49 * t166 - 60 * t41 * t28;
      double t170 = t25 * z;
      double t173 = t8 * t46;
      double t204 = t50 * t50;
      double t207 = H * t85;
      double t210 = H * t50;
      double t215 = M * t85;
      double t218 = M * t50;
      double t221 = 160 * t5 * t85 * t1 + 472 * t84 * t85 * t1 - 32 * t17 * t8 * t46 - 464 * t84 * t51 * z + 80 * H * t204 + 80 * M * t204 - 72 * t91 * t1 + 160 * t62 * t166 - 276 * t49 * t170 - 180 * t62 * t170 + 32 * t49 * t173 - 32 * t62 * t173 + 112 * t207 * z - 128 * t210 * t46 - 272 * t215 * z + 128 * t218 * t46 + 36 * t72 * t28 + 24 * t78 * t42 - 40 * t93 * t42 - 76 * t69 * t8;
      double t225 = dchiTrunc(r, l, eps);
      double t229 = t116 * t1;
      double t243 = d * t116;
      double t246 = t116 * z;
      double t267 = t50 * t9;
      double t270 = t138 * z;
      double t273 = 4 * t125 * t116 * z + 12 * t6 * t116 * z - 16 * t12 * t9 * t46 - 20 * t152 * t116 - 8 * t126 * t229 - 2 * t129 * t229 + 8 * t14 * t229 + 72 * t41 * t141 - 4 * t18 * t229 + 6 * t7 * t229 + 36 * t24 * t243 + 8 * t24 * t246 + 20 * t31 * t243 - 36 * t36 * t243 - 16 * t31 * t246 - 8 * t36 * t246 + 72 * t49 * t267 + 8 * t49 * t270;
      double t274 = t9 * t46;
      double t314 = -64 * t84 * t166 * z + 16 * t17 * t9 * t46 + 48 * t84 * t85 * t8 - 8 * t72 * t141 - 64 * t81 * t173 + 64 * t96 * t173 + 96 * t207 * t42 - 64 * t210 * t59 - 96 * t215 * t42 + 64 * t218 * t59 - 16 * t62 * t267 - 72 * t62 * t270 - 16 * t49 * t274 + 16 * t62 * t274 + 144 * t78 * t28 - 80 * t93 * t28 - 56 * t69 * t9 - 48 * t91 * t8;
      double t318 = d2chiTrunc(r, l, eps);

      return t4 * (t99 + t58) * r * t108 * t110 * t112 / 5 + t4 * (t221 + t169) * r * t225 * t108 * t110 / 5 + t4 * (t314 + t273) * r * t318 * t108 * t110 / 5;


    }

    PetscScalar drU0_tang3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double t3 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t3);
      double z  = X[2];

      double t1 = H * H;
      double t2 = t1 * t1;
      double t4 = t3 * t3;
      double t8 = t1 * H;
      double t9 = t8 * M;
      double t10 = d * t4;
      double t13 = t4 * z;
      double t16 = M * M;
      double t17 = t1 * t16;
      double t22 = t16 * M;
      double t23 = H * t22;
      double t28 = t16 * t16;
      double t29 = t28 * d;
      double t35 = t8 * d;
      double t36 = t3 * z;
      double t40 = z * z;
      double t43 = t1 * M;
      double t44 = d * d;
      double t45 = t44 * t3;
      double t49 = d * t3 * z;
      double t52 = 56 * t2 * t4 * z + 56 * t28 * t4 * z - 144 * t8 * t3 * t40 - 24 * t17 * t10 + 24 * t23 * t10 - 104 * t9 * t10 + 48 * t17 * t13 + 80 * t23 * t13 + 80 * t9 * t13 - 56 * t29 * t4 + 32 * t35 * t36 - 224 * t43 * t45 + 576 * t43 * t49;
      double t53 = t3 * t40;
      double t56 = H * t16;
      double t63 = t22 * t44;
      double t66 = t22 * d;
      double t72 = t1 * t44;
      double t75 = t1 * d;
      double t78 = H * M;
      double t79 = t44 * d;
      double t85 = t16 * t79;
      double t87 = t16 * t44;
      double t90 = t16 * d;
      double t93 = 144 * t22 * t3 * t40 + 64 * t78 * t44 * z + 112 * t63 * t3 - 256 * t66 * t36 + 288 * t75 * t40 - 288 * t90 * t40 - 144 * t43 * t53 - 208 * t56 * t45 + 288 * t56 * t49 + 144 * t56 * t53 - 160 * t72 * z - 32 * t78 * t79 + 416 * t87 * z - 128 * t85;
      double t100 = pow(t3 * H + M * t3 + 2 * d, 2);
      double t101 = t100 * t100;
      double t102 = 0.1e1 / t101;
      double t104 = 0.1e1 / (H + M);
      double t106 = chiTrunc(r, l, eps);
      double t110 = t4 * t4;
      double t113 = t2 * M;
      double t116 = t8 * t16;
      double t119 = t1 * t22;
      double t122 = H * t28;
      double t129 = t4 * t3;
      double t135 = d * t129;
      double t138 = t129 * z;
      double t163 = t44 * t4;
      double t166 = t10 * z;
      double t169 = -56 * t28 * t129 * z + 120 * t8 * t4 * t40 + 120 * t8 * t44 * t4 + 16 * t29 * t129 - 192 * t35 * t13 + 24 * t17 * t135 - 56 * t23 * t135 - 48 * t17 * t138 - 80 * t23 * t138 + 456 * t43 * t163 - 528 * t43 * t166;
      double t171 = t4 * t40;
      double t207 = t44 * t44;
      double t210 = H * t79;
      double t213 = H * t44;
      double t218 = M * t79;
      double t221 = M * t44;
      double t224 = 480 * t78 * t79 * t3 - 960 * t78 * t45 * z + 80 * H * t207 + 80 * M * t207 + 128 * t210 * z - 288 * t213 * t40 - 448 * t218 * z + 288 * t221 * t40 - 160 * t85 * t3 + 96 * t87 * t36 - 96 * t90 * t53;
      double t228 = dchiTrunc(r, l, eps);
      double t232 = t110 * t3;
      double t244 = d * t110;
      double t264 = t44 * t129;
      double t267 = t135 * z;
      double t270 = t129 * t40;
      double t273 = -48 * t17 * t110 * z + 24 * t2 * t110 * z + 24 * t28 * t110 * z - 48 * t8 * t129 * t40 - 24 * t29 * t110 + 12 * t113 * t232 + 12 * t116 * t232 - 12 * t119 * t232 - 12 * t122 * t232 + 144 * t35 * t138 + 24 * t17 * t244 - 72 * t23 * t244 + 72 * t9 * t244 + 144 * t43 * t264 - 48 * t43 * t267 - 48 * t43 * t270;
      double t311 = 48 * t22 * t129 * t40 - 192 * t78 * t163 * z + 96 * t78 * t79 * t4 - 96 * t63 * t129 + 288 * t72 * t13 - 96 * t87 * t13 + 48 * t66 * t138 - 192 * t75 * t171 + 192 * t90 * t171 + 192 * t210 * t36 - 192 * t213 * t53 - 192 * t218 * t36 + 192 * t221 * t53 - 48 * t56 * t264 - 144 * t56 * t267 + 48 * t56 * t270 - 96 * t85 * t4;
      double t314 = d2chiTrunc(r, l, eps);
      return -0.2e1 / 0.5e1 * r * (t93 + t52) * t102 * t104 * t106 - 0.2e1 / 0.5e1 * r * (160 * t1 * t79 * t3 - 120 * t22 * t4 * t40 + 48 * t66 * t13 - 48 * t56 * t163 - 288 * t56 * t166 + 120 * t43 * t171 - 120 * t56 * t171 - 96 * t72 * t36 - 48 * t63 * t4 + 96 * t75 * t53 + t224 + 5 * t2 * H * t110 - 5 * t28 * M * t110 + 40 * t2 * d * t129 - 56 * t2 * t129 * z + 9 * t113 * t110 + 4 * t116 * t110 - 4 * t119 * t110 - 9 * t122 * t110 + 136 * t9 * t135 - 80 * t9 * t138 + t169) * t228 * t102 * t104 - 0.2e1 / 0.5e1 * r * (t311 + t273) * t314 * t102 * t104;

    }

    PetscScalar dzU0_tang3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double t3 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t3);
      double z  = X[2];

      double t1 = H * H;
      double t2 = t1 * H;
      double t4 = t3 * t3;
      double t7 = t1 * M;
      double t10 = M * M;
      double t11 = H * t10;
      double t14 = t10 * M;
      double t17 = t1 * d;
      double t23 = H * M;
      double t33 = d * d;
      double t34 = H * t33;
      double t36 = H * d;
      double t39 = M * t33;
      double t41 = M * d;
      double t44 = -2 * t10 * d * t3 + 26 * t23 * d * t3 - 18 * t1 * t3 * z + 18 * t10 * t3 * z + 3 * t11 * t4 + 7 * t14 * t4 + 16 * t17 * t3 + 7 * t2 * t4 + 12 * t36 * z + 3 * t7 * t4 - 12 * t41 * z + 4 * t34 + 16 * t39;
      double t48 = t3 * H + M * t3 + 2 * d;
      double t49 = t48 * t48;
      double t51 = 0.1e1 / t49 / t48;
      double t54 = 0.1e1 / (H + M);
      double t55 = chiTrunc(r, l, eps);
      double t58 = t4 * t3;
      double t80 = t3 * z;
      double t87 = 12 * t23 * d * t4 + 12 * t1 * t4 * z - 12 * t10 * t4 * z + 3 * t11 * t58 - 3 * t14 * t58 - 12 * t17 * t4 - 3 * t2 * t58 - 12 * t34 * t3 + 12 * t39 * t3 + 24 * t36 * t80 - 24 * t41 * t80 + 3 * t7 * t58;
      double t88 = dchiTrunc(r, l, eps);
      return  0.8e1 / 0.5e1 * t44 * t51 * t54 * t55 + 0.8e1 / 0.5e1 * t87 * t88 * t51 * t54;
    }

    PetscScalar drV0_tang3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double t3 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t3);
      double z  = X[2];

      double t1 = H * H;
      double t2 = t1 * H;
      double t7 = t1 * M;
      double t8 = d * t3;
      double t11 = t3 * z;
      double t14 = M * M;
      double t15 = H * t14;
      double t20 = t14 * M;
      double t21 = t20 * d;
      double t27 = t1 * d;
      double t30 = z * z;
      double t33 = H * M;
      double t34 = d * d;
      double t40 = t14 * t34;
      double t42 = t14 * d;
      double t47 = 32 * t33 * d * z + 8 * t2 * t3 * z + 8 * t20 * t3 * z + 48 * t1 * t30 + 72 * t15 * t11 + 72 * t7 * t11 - 48 * t14 * t30 - 64 * t15 * t8 - 8 * t21 * t3 + 16 * t27 * z - 16 * t33 * t34 + 112 * t42 * z - 8 * t7 * t8 - 64 * t40;
      double t52 = t3 * H + M * t3 + 2 * d;
      double t53 = t52 * t52;
      double t55 = 0.1e1 / t53 / t52;
      double t57 = 0.1e1 / (H + M);
      double t59 = chiTrunc(r, l, eps);
      double t62 = t1 * t1;
      double t63 = t3 * t3;
      double t64 = t63 * t3;
      double t73 = t14 * t14;
      double t82 = d * t63;
      double t85 = t63 * z;
      double t100 = -16 * H * t20 * t64 + 16 * t2 * M * t64 + 30 * t2 * d * t63 + 60 * t1 * t34 * t3 - 8 * t2 * t63 * z - 8 * t20 * t63 * z - 2 * t15 * t82 - 72 * t15 * t85 - 22 * t21 * t63 + 5 * t62 * t64 - 5 * t73 * t64 + 74 * t7 * t82 - 72 * t7 * t85;
      double t119 = t34 * d;
      double t136 = -48 * H * d * t30 - 32 * H * t34 * z + 48 * M * d * t30 - 128 * M * t34 * z - 24 * t1 * t3 * t30 + 24 * t14 * t3 * t30 + 104 * t33 * t34 * t3 - 208 * t33 * t8 * z + 40 * H * t119 + 40 * M * t119 - 32 * t27 * t11 - 80 * t42 * t11 - 4 * t40 * t3;
      double t139 = dchiTrunc(r, l, eps);

      return  0.2e1 / 0.5e1 * r * t47 * t55 * t57 * t59 + 0.2e1 / 0.5e1 * r * (t136 + t100) * t139 * t55 * t57;
    }


    PetscScalar dzV0_tang3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double t2 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t2);
      double z  = X[2];
      
      double t1 = H * H;
      double t7 = M * M;
      double t18 = chiTrunc(r, l, eps);
      double t24 = pow(t2 * H + M * t2 + 2 * d, 2);

      return -0.8e1 / 0.5e1 * (8 * H * M * t2 + 2 * H * d + 6 * H * z + 8 * M * d - 6 * M * z + t1 * t2 + t7 * t2) * t18 / t24 / (H + M);
    }


    PetscScalar drW0_tang3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double t6 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t6);
      double z  = X[2];

      double t1 = H * H;
      double t2 = t1 * t1;
      double t3 = t2 * H;
      double t4 = M * M;
      double t5 = t3 * t4;
      double t7 = t6 * t6;
      double t8 = t7 * t7;
      double t9 = t8 * t6;
      double t12 = t4 * M;
      double t13 = t2 * t12;
      double t16 = t1 * H;
      double t17 = t4 * t4;
      double t18 = t16 * t17;
      double t21 = t17 * M;
      double t22 = t1 * t21;
      double t25 = t3 * M;
      double t26 = t8 * z;
      double t29 = t2 * t4;
      double t30 = d * t8;
      double t35 = t16 * t12;
      double t38 = t1 * t17;
      double t43 = H * t21;
      double t48 = t7 * t6;
      double t50 = z * z;
      double t53 = t2 * M;
      double t54 = d * t48;
      double t55 = t54 * z;
      double t58 = t48 * t50;
      double t61 = -96 * t3 * t48 * t50 - 60 * t13 * t9 - 60 * t18 * t9 - 20 * t22 * t9 + 24 * t25 * t26 + 48 * t29 * t26 - 48 * t38 * t26 - 24 * t43 * t26 - 224 * t29 * t30 - 424 * t35 * t30 - 176 * t38 * t30 + 24 * t43 * t30 - 20 * t5 * t9 + 576 * t53 * t55 - 336 * t53 * t58;
      double t62 = t16 * t4;
      double t63 = d * d;
      double t64 = t63 * t48;
      double t71 = t1 * t12;
      double t78 = H * t17;
      double t85 = t21 * t63;
      double t88 = t21 * d;
      double t89 = t48 * z;
      double t95 = t2 * d;
      double t96 = t7 * t50;
      double t100 = t50 * z;
      double t103 = t16 * M;
      double t104 = t63 * t7;
      double t105 = t104 * z;
      double t108 = d * t7;
      double t109 = t108 * t50;
      double t112 = 320 * t2 * t7 * t100 - 96 * t21 * t48 * t50 - 320 * t103 * t105 - 1760 * t103 * t109 - 96 * t85 * t48 + 864 * t62 * t55 + 192 * t71 * t55 + 96 * t78 * t55 - 528 * t62 * t58 - 528 * t71 * t58 - 336 * t78 * t58 - 1280 * t62 * t64 - 944 * t71 * t64 + 240 * t78 * t64 + 192 * t88 * t89 + 640 * t95 * t96;
      double t114 = t7 * t100;
      double t117 = t1 * t4;
      double t118 = t63 * d;
      double t119 = t118 * t7;
      double t126 = H * t12;
      double t135 = t17 * t118;
      double t138 = t17 * t63;
      double t139 = t7 * z;
      double t142 = t17 * d;
      double t148 = t16 * t63;
      double t149 = t6 * t50;
      double t152 = t16 * d;
      double t153 = t6 * t100;
      double t156 = t1 * M;
      double t157 = t118 * t6;
      double t158 = t157 * z;
      double t162 = t63 * t6 * t50;
      double t165 = -320 * t100 * t17 * t7 + 640 * t103 * t114 + 3840 * t105 * t117 + 1920 * t105 * t126 - 3840 * t109 * t117 + 160 * t109 * t126 - 640 * t114 * t126 - 1920 * t117 * t119 - 1440 * t119 * t126 + 960 * t135 * t7 - 2240 * t138 * t139 + 1600 * t142 * t96 + 1408 * t148 * t149 - 2176 * t152 * t153 - 2816 * t156 * t158 + 4800 * t156 * t162;
      double t167 = d * t6 * t100;
      double t170 = H * t4;
      double t171 = t63 * t63;
      double t172 = t171 * t6;
      double t181 = t12 * t171;
      double t184 = t12 * t118;
      double t185 = t6 * z;
      double t188 = t12 * t63;
      double t191 = t12 * d;
      double t194 = t1 * t118;
      double t197 = t1 * t63;
      double t200 = H * M;
      double t207 = t4 * t171;
      double t210 = t4 * t118;
      double t213 = t4 * t63;
      double t216 = 128 * t118 * t200 * t50 - 128 * t171 * t200 * z + 512 * t100 * t197 - 512 * t100 * t213 - 5120 * t149 * t188 + 2176 * t153 * t191 - 2176 * t156 * t167 - 256 * t158 * t170 - 1728 * t162 * t170 + 2176 * t167 * t170 - 192 * t170 * t172 - 768 * t181 * t6 + 3712 * t184 * t185 - 512 * t194 * t50 - 512 * t207 * z + 1024 * t210 * t50;
      double t220 = 0.1e1 / (H + M);
      double t225 = H * t6 + M * t6 + 2 * d;
      double t226 = t225 * t225;
      double t227 = t226 * t226;
      double t229 = 0.1e1 / t227 / t225;
      double t230 = chiTrunc(r,l,eps);
      double t233 = t2 * t1;
      double t234 = t233 * M;
      double t235 = t8 * t7;
      double t245 = t17 * t4;
      double t246 = H * t245;
      double t249 = t17 * t12;
      double t255 = d * t9;
      double t258 = t9 * z;
      double t276 = t245 * d;
      double t282 = t3 * d;
      double t288 = -10 * t245 * t9 * z + 108 * t3 * t50 * t8 + 462 * t255 * t29 + 460 * t255 * t35 + 76 * t255 * t38 - 130 * t255 * t43 - 46 * t258 * t29 + 46 * t258 * t38 + 8 * t258 * t43 + 100 * t26 * t282 - 58 * t276 * t9;
      double t290 = t63 * t8;
      double t293 = t30 * z;
      double t296 = t8 * t50;
      double t325 = t2 * t63;
      double t333 = t118 * t48;
      double t336 = t64 * z;
      double t339 = t54 * t50;
      double t342 = -320 * t100 * t2 * t48 + 108 * t21 * t50 * t8 + 1200 * t103 * t333 - 1376 * t103 * t336 + 2240 * t103 * t339 - 316 * t26 * t88 - 108 * t293 * t78 + 348 * t296 * t78 + 400 * t325 * t89 - 256 * t58 * t95 - 112 * t8 * t85;
      double t345 = t48 * t100;
      double t372 = t16 * t118;
      double t379 = t171 * t7;
      double t382 = t119 * z;
      double t385 = t104 * t50;
      double t388 = t108 * t100;
      double t395 = 320 * t100 * t17 * t48 + 896 * t114 * t152 + 800 * t139 * t372 - 1216 * t142 * t58 - 2208 * t148 * t96 + 1200 * t156 * t379 + 1376 * t156 * t382 + 1440 * t156 * t385 + 896 * t156 * t388 + 3792 * t170 * t379 - 6944 * t170 * t382;
      double t409 = t1 * t171;
      double t416 = t171 * d;
      double t423 = 3840 * t172 * t200 * z + 480 * t200 * t416 * t6 - 896 * t114 * t191 + 928 * t139 * t184 - 2304 * t149 * t194 + 2816 * t153 * t197 + 4128 * t170 * t385 - 896 * t170 * t388 - 1632 * t181 * t7 + 800 * t185 * t409 + 480 * t188 * t96;
      double t427 = t4 * t416;
      double t436 = H * t416;
      double t439 = H * t171;
      double t442 = H * t118;
      double t445 = M * t416;
      double t448 = M * t171;
      double t451 = M * t118;
      double t454 = -3840 * t157 * t200 * t50 - 512 * t100 * t442 + 512 * t100 * t451 + 6144 * t149 * t210 - 2816 * t153 * t213 - 4640 * t185 * t207 + 480 * t427 * t6 + 320 * t436 * z + 448 * t439 * t50 + 320 * t445 * z - 1088 * t448 * t50;
      double t458 = dchiTrunc(r,l,eps);
      double t460 = t220 * t229;
      double t462 = t8 * t48;
      double t480 = d * t235;
      double t483 = t235 * z;
      double t511 = -20 * t235 * t245 * z - 64 * t3 * t50 * t9 - 148 * t235 * t276 + 200 * t258 * t282 + 852 * t29 * t480 + 244 * t29 * t483 + 1280 * t35 * t480 + 256 * t38 * t480 - 244 * t38 * t483 - 420 * t43 * t480 - 152 * t43 * t483;
      double t513 = t63 * t9;
      double t516 = t255 * z;
      double t519 = t9 * t50;
      double t555 = t118 * t8;
      double t558 = t290 * z;
      double t561 = t30 * t50;
      double t564 = 160 * t100 * t2 * t8 - 64 * t21 * t50 * t9 + 800 * t103 * t555 + 6272 * t103 * t558 - 1760 * t103 * t561 - 72 * t258 * t88 + 800 * t26 * t325 - 128 * t296 * t95 - 1192 * t516 * t78 - 224 * t519 * t78 - 704 * t85 * t9;
      double t567 = t8 * t100;
      double t600 = t171 * t48;
      double t603 = t333 * z;
      double t606 = t64 * t50;
      double t609 = t54 * t100;
      double t616 = -160 * t100 * t17 * t8 + 352 * t142 * t296 + 768 * t148 * t58 + 192 * t152 * t345 + 800 * t156 * t600 + 10432 * t156 * t603 - 5184 * t156 * t606 + 192 * t156 * t609 + 2880 * t170 * t600 - 640 * t170 * t603 + 1600 * t372 * t89;
      double t642 = 6784 * t200 * t379 * z + 320 * t200 * t416 * t7 - 1152 * t114 * t197 + 1600 * t139 * t409 - 4608 * t170 * t606 - 192 * t170 * t609 - 800 * t181 * t48 - 3712 * t184 * t89 + 1344 * t188 * t58 - 192 * t191 * t345 + 2560 * t194 * t96;
      double t666 = -6784 * t119 * t200 * t50 + 1152 * t114 * t213 - 3264 * t139 * t207 + 2048 * t149 * t439 - 3328 * t149 * t448 - 1792 * t153 * t442 + 1792 * t153 * t451 + 640 * t185 * t436 + 640 * t185 * t445 - 896 * t210 * t96 + 320 * t427 * t7;
      double t670 = d2chiTrunc(r,l,eps);
      double t673 = t8 * t8;
      double t686 = t462 * z;
      double t689 = d * t462;
      double t709 = t480 * z;
      double t712 = 48 * t235 * t3 * t50 + 28 * t13 * t673 + 8 * t18 * t673 - 24 * t22 * t673 - 20 * t246 * t673 - 4 * t249 * t673 + 48 * t25 * t686 - 48 * t276 * t462 + 96 * t29 * t686 + 96 * t29 * t689 + 144 * t35 * t689 - 96 * t38 * t686 - 48 * t38 * t689 - 48 * t43 * t686 - 144 * t43 * t689 + 12 * t5 * t673 + 384 * t53 * t709;
      double t713 = t235 * t50;
      double t716 = t63 * t235;
      double t747 = t513 * z;
      double t750 = t255 * t50;
      double t753 = -64 * t100 * t2 * t9 + 48 * t21 * t235 * t50 + 1152 * t103 * t747 + 192 * t103 * t750 - 192 * t235 * t85 - 96 * t483 * t88 + 384 * t519 * t95 + 48 * t53 * t713 + 480 * t62 * t709 - 96 * t62 * t713 + 288 * t62 * t716 - 288 * t709 * t71 - 480 * t709 * t78 - 96 * t71 * t713 + 192 * t71 * t716 + 48 * t713 * t78 - 288 * t716 * t78;
      double t755 = t9 * t100;
      double t758 = t118 * t9;
      double t786 = t555 * z;
      double t789 = t30 * t100;
      double t795 = 64 * t100 * t17 * t9 + 192 * t170 * t171 * t8 - 128 * t103 * t755 + 576 * t117 * t747 - 576 * t117 * t750 + 384 * t117 * t758 - 1152 * t126 * t747 - 192 * t126 * t750 + 128 * t126 * t755 - 64 * t126 * t758 - 320 * t135 * t9 - 576 * t138 * t258 + 192 * t142 * t519 + 1152 * t148 * t296 - 384 * t152 * t567 + 1536 * t156 * t786 - 384 * t156 * t789;
      double t833 = -1152 * t170 * t290 * t50 - 768 * t200 * t333 * t50 + 768 * t200 * t600 * z - 512 * t114 * t442 + 512 * t114 * t451 - 384 * t170 * t786 + 384 * t170 * t789 - 192 * t181 * t8 - 1152 * t184 * t26 + 384 * t191 * t567 + 1536 * t194 * t58 - 768 * t197 * t345 - 768 * t207 * t89 - 768 * t210 * t58 + 768 * t213 * t345 + 768 * t439 * t96 - 768 * t448 * t96;
      double t836 = d3chiTrunc(r,l,eps);
      return  (t216 + t165 + t112 + t61) * t220 * t229 * t230 / 5 + (-8 * t25 * t258 + 600 * t53 * t290 - 588 * t53 * t293 + 348 * t53 * t296 + 2040 * t62 * t290 - 952 * t62 * t293 + 504 * t62 * t296 + 952 * t71 * t290 - 56 * t71 * t293 + 504 * t71 * t296 - 600 * t78 * t290 - 640 * t103 * t345 + 4432 * t117 * t333 - 4032 * t117 * t336 + 4032 * t117 * t339 - 176 * t126 * t333 - 1184 * t126 * t336 + 320 * t126 * t339 + 640 * t126 * t345 - 976 * t135 * t48 + 1072 * t138 * t89 + 15 * t234 * t235 + 53 * t5 * t235 + 72 * t13 * t235 + 42 * t18 * t235 - t22 * t235 - 15 * t246 * t235 - 6 * t249 * t235 + 150 * t25 * t255 + t423 + t454 + 10 * t233 * t9 * z + t288 + t342 + t395) * t458 * t460 / 5 + (-800 * t126 * t561 - 320 * t126 * t567 - 1408 * t135 * t8 - 1024 * t138 * t26 + 320 * t103 * t567 + 4576 * t117 * t555 + 2784 * t117 * t558 - 2784 * t117 * t561 + 1152 * t126 * t555 - 3712 * t126 * t558 + 10 * t234 * t462 + 92 * t5 * t462 + 198 * t13 * t462 + 128 * t18 * t462 - 34 * t22 * t462 - 60 * t246 * t462 - 14 * t249 * t462 + 100 * t25 * t480 + 152 * t25 * t483 + 400 * t53 * t513 + 1640 * t53 * t516 - 224 * t53 * t519 + 2944 * t62 * t513 + 1632 * t62 * t516 - 352 * t62 * t519 + 2544 * t71 * t513 - 928 * t71 * t516 - 352 * t71 * t519 - 704 * t78 * t513 + 20 * t233 * t235 * z + t666 + t616 + t564 + t642 + t511) * t670 * t460 / 5 + (t833 + t795 + t753 + t712) * t836 * t460 / 5;
    }

    PetscScalar dzW0_tang3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double t4 = X[0]*X[0]+X[1]*X[1];
      double r  = sqrt(t4);
      double z  = X[2];

      double t1 = H * H;
      double t2 = t1 * t1;
      double t3 = t2 * M;
      double t5 = t4 * t4;
      double t6 = t5 * t4;
      double t9 = t1 * H;
      double t10 = M * M;
      double t11 = t9 * t10;
      double t14 = t10 * M;
      double t15 = t1 * t14;
      double t18 = t10 * t10;
      double t19 = H * t18;
      double t25 = t9 * M;
      double t26 = d * t5;
      double t29 = t5 * z;
      double t32 = t1 * t10;
      double t37 = H * t14;
      double t42 = t18 * d;
      double t48 = t9 * d;
      double t49 = t4 * z;
      double t53 = z * z;
      double t56 = 32 * t18 * t5 * z + 32 * t2 * t5 * z - 96 * t9 * t4 * t53 - 12 * t11 * t6 + 12 * t15 * t6 + 12 * t19 * t6 - 152 * t25 * t26 + 80 * t25 * t29 - 48 * t32 * t26 + 72 * t37 * t26 + 96 * t32 * t29 + 80 * t37 * t29 - 12 * t3 * t6 - 32 * t42 * t5 - 64 * t48 * t49;
      double t57 = t1 * M;
      double t58 = d * d;
      double t59 = t58 * t4;
      double t63 = d * t4 * z;
      double t66 = t4 * t53;
      double t69 = H * t10;
      double t76 = t14 * t58;
      double t79 = t14 * d;
      double t85 = t1 * t58;
      double t88 = t1 * d;
      double t91 = H * M;
      double t92 = t58 * d;
      double t98 = t10 * t92;
      double t100 = t10 * t58;
      double t103 = t10 * d;
      double t106 = 96 * t14 * t4 * t53 + 64 * t91 * t58 * z + 512 * t100 * z - 384 * t103 * t53 + 160 * t76 * t4 - 256 * t79 * t49 + 384 * t88 * t53 - 272 * t57 * t59 + 576 * t57 * t63 - 96 * t57 * t66 - 208 * t69 * t59 + 384 * t69 * t63 + 96 * t69 * t66 - 256 * t85 * z - 32 * t91 * t92 - 128 * t98;
      double t113 = pow(t4 * H + M * t4 + 2 * d, 2);
      double t114 = t113 * t113;
      double t115 = 0.1e1 / t114;
      double t117 = 0.1e1 / (H + M);
      double t119 = chiTrunc(r, l, eps);
      double t123 = t5 * t5;
      double t143 = d * t6;
      double t146 = t6 * z;
      double t170 = t58 * t5;
      double t173 = 5 * t2 * H * t123 - 5 * t18 * M * t123 + 40 * t2 * d * t6 - 44 * t18 * t6 * z - 44 * t2 * t6 * z + 96 * t9 * t5 * t53 + 120 * t9 * t58 * t5 + 10 * t11 * t123 - 10 * t15 * t123 - 15 * t19 * t123 + 15 * t3 * t123 + 172 * t25 * t143 + 36 * t32 * t143 - 92 * t37 * t143 - 80 * t25 * t146 - 72 * t32 * t146 - 80 * t37 * t146 + 528 * t57 * t170 - 120 * t48 * t29 + 4 * t42 * t6;
      double t174 = t26 * z;
      double t177 = t5 * t53;
      double t208 = t58 * t58;
      double t211 = H * t92;
      double t214 = H * t58;
      double t219 = M * t92;
      double t222 = M * t58;
      double t225 = 160 * t1 * t92 * t4 - 96 * t14 * t5 * t53 + 528 * t91 * t92 * t4 - 1056 * t91 * t59 * z + 80 * H * t208 + 80 * M * t208 + 48 * t100 * t49 - 72 * t69 * t170 - 552 * t57 * t174 - 360 * t69 * t174 + 96 * t57 * t177 - 96 * t69 * t177 + 224 * t211 * z - 384 * t214 * t53 - 544 * t219 * z + 384 * t222 * t53 + 72 * t79 * t29 - 208 * t98 * t4 + 48 * t85 * t49 - 96 * t76 * t5;
      double t228 = dchiTrunc(r, l, eps);
      double t232 = t123 * t4;
      double t244 = d * t123;
      double t264 = t58 * t6;
      double t267 = t143 * z;
      double t270 = t6 * t53;
      double t273 = 24 * t18 * t123 * z + 24 * t2 * t123 * z - 48 * t32 * t123 * z - 48 * t9 * t6 * t53 + 12 * t11 * t232 - 24 * t42 * t123 + 144 * t48 * t146 - 12 * t15 * t232 - 12 * t19 * t232 + 12 * t3 * t232 + 72 * t25 * t244 + 24 * t32 * t244 - 72 * t37 * t244 + 144 * t57 * t264 - 48 * t57 * t267 - 48 * t57 * t270;
      double t311 = 48 * t14 * t6 * t53 - 192 * t91 * t170 * z + 96 * t91 * t92 * t5 - 96 * t100 * t29 + 192 * t103 * t177 + 48 * t79 * t146 - 192 * t88 * t177 + 192 * t211 * t49 - 192 * t214 * t66 - 192 * t219 * t49 + 192 * t222 * t66 - 48 * t69 * t264 - 144 * t69 * t267 + 48 * t69 * t270 + 288 * t85 * t29 - 96 * t98 * t5 - 96 * t76 * t6;
      double t314 = d2chiTrunc(r, l, eps);
      return 0.2e1 / 0.5e1 * r * (t106 + t56) * t115 * t117 * t119 + 0.2e1 / 0.5e1 * r * (t225 + t173) * t228 * t115 * t117 + 0.2e1 / 0.5e1 * r * (t311 + t273) * t314 * t115 * t117;
    }

    PetscScalar ur_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double cosTheta;

      if(r2!=0)
        cosTheta = X[0]/sqrt(r2);
      else
        cosTheta = 1.;
      
      return U*U0_tang3D(X, H, M, d, U, l, eps, param)*cosTheta/2.;
    }

    PetscScalar ut_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double sinTheta;

      if(r2!=0)
        sinTheta = X[1]/sqrt(r2);
      else
        sinTheta = 0.;

      return U*V0_tang3D(X, H, M, d, U, l, eps, param)*sinTheta/2.;
    }

    PetscScalar uz_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double cosTheta;

      if(r2!=0)
        cosTheta = X[0]/sqrt(r2);
      else
        cosTheta = 1.;

      //std::cout  << " inside uz : " << W0_tang3D(X, H, M, d, U, l, eps, param)*cosTheta/2. << " r2 = " << r2 << std::endl;

      return U*W0_tang3D(X, H, M, d, U, l, eps, param)*cosTheta/2.;
    }

    PetscScalar drur_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double cosTheta;

      if(r2!=0)
        cosTheta = X[0]/sqrt(r2);
      else
        cosTheta = 1.;

      return U*drU0_tang3D(X, H, M, d, U, l, eps, param)*cosTheta/2.;
    }

    PetscScalar dtur_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double sinTheta;

      if(r2!=0)
        sinTheta = X[1]/sqrt(r2);
      else
        sinTheta = 0.;

      return -U*U0_tang3D(X, H, M, d, U, l, eps, param)*sinTheta/2.;
    }

    PetscScalar dzur_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double cosTheta;

      if(r2!=0)
        cosTheta = X[0]/sqrt(r2);
      else
        cosTheta = 1.;

      return U*dzU0_tang3D(X, H, M, d, U, l, eps, param)*cosTheta/2.;
    }

    PetscScalar drut_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double sinTheta;

      if(r2!=0)
        sinTheta = X[1]/sqrt(r2);
      else
        sinTheta = 0.;

      return U*drV0_tang3D(X, H, M, d, U, l, eps, param)*sinTheta/2.;
    }

    PetscScalar dtut_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double cosTheta;

      if(r2!=0)
        cosTheta = X[0]/sqrt(r2);
      else
        cosTheta = 1.;

      return U*V0_tang3D(X, H, M, d, U, l, eps, param)*cosTheta/2.;
    }

    PetscScalar dzut_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double sinTheta;

      if(r2!=0)
        sinTheta = X[1]/sqrt(r2);
      else
        sinTheta = 0.;

      return U*dzV0_tang3D(X, H, M, d, U, l, eps, param)*sinTheta/2.;
    }

    PetscScalar druz_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double cosTheta;

      if(r2!=0)
        cosTheta = X[0]/sqrt(r2);
      else
        cosTheta = 1.;

      return U*drW0_tang3D(X, H, M, d, U, l, eps, param)*cosTheta/2.;
    }

    PetscScalar dtuz_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double sinTheta;

      if(r2!=0)
        sinTheta = X[1]/sqrt(r2);
      else
        sinTheta = 0.;

      return - U*W0_tang3D(X, H, M, d, U, l, eps, param)*sinTheta/2.;
    }

    PetscScalar dzuz_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double cosTheta;

      if(r2!=0)
        cosTheta = X[0]/sqrt(r2);
      else
        cosTheta = 1.;

      return U*dzW0_tang3D(X, H, M, d, U, l, eps, param)*cosTheta/2.;
    }




    PetscScalar ux_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      //std::cout  << " inside ux : " << ur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param) << " r2 = " << r2 << std::endl;

      
      if(r2!=0)
        return ur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[0]/r - ut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[1]/r;
      else
        return ur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param);
    }

    PetscScalar uy_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      //std::cout  << " inside u : " << ut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param) << " r2 = " << r2 << std::endl;

      if(r2!=0)
        return ur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[1]/r + ut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[0]/r;
      else
        return ut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param);
    }















    PetscScalar dxux_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      double x2 = X[0]*X[0];
      double y2 = X[1]*X[1];
      double xy = X[0]*X[1];

      if(r2!=0)
        return drur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*x2/r2 
          +      ur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*y2/r2/r
          -    dtur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2/r
          +      ut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2/r
          -    drut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2
          +    dtut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*y2/r2/r;
      else
        return 0;
    }

    PetscScalar dyux_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      double x2 = X[0]*X[0];
      double y2 = X[1]*X[1];
      double xy = X[0]*X[1];

      if(r2!=0)
        return drur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2 
          -      ur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2/r
          +    dtur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*x2/r2/r
          -      ut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*x2/r2/r
          -    drut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*y2/r2
          -    dtut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2/r;
      else
        return 0;
    }

    PetscScalar dzux_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
        return dzur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[0]/r 
          -    dzut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[1]/r;
      else 
        return dzur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param);
    }

    PetscScalar dxuy_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      double x2 = X[0]*X[0];
      double y2 = X[1]*X[1];
      double xy = X[0]*X[1];

      if(r2!=0)
        return drur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2 
          -      ur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2/r
          -    dtur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*y2/r2/r
          +      ut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*y2/r2/r
          +    drut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*x2/r2
          -    dtut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2/r;
      else
        return 0;
    }

    PetscScalar dyuy_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      double x2 = X[0]*X[0];
      double y2 = X[1]*X[1];
      double xy = X[0]*X[1];

      if(r2!=0)
        return drur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*y2/r2 
          +      ur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*x2/r2/r
          +    dtur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2/r
          -      ut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2/r
          +    drut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*xy/r2
          +    dtut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*x2/r2/r;
      else
        return 0;
    }

    PetscScalar dzuy_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
        return dzur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[1]/r 
          +    dzut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[0]/r;
      else
	return dzut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param);
    }

    PetscScalar dxuz_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
        return druz_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[0]/r
          -    dtuz_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[1]/r2;
      else
        return 0;
    }

    PetscScalar dyuz_sing_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      if(r2!=0)
        return druz_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[1]/r
          +    dtuz_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)*X[0]/r2;
      else
        return 0;
    }

    PetscScalar DIVCYLIND_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      // return drur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)
      //   +      ur_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)/r
      //   +    dtut_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)/r
      //   +    dzuz_sing_tangMvt3D(X, H,  M, d, U, l, eps, param);

      return drU0_tang3D(X, H,  M, d, U, l, eps, param)
        +      U0_tang3D(X, H,  M, d, U, l, eps, param)/r
        +      V0_tang3D(X, H,  M, d, U, l, eps, param)/r
        +    dzW0_tang3D(X, H,  M, d, U, l, eps, param);
    }

    PetscScalar DIVCART_tangMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal d, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      return dxux_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)
        +    dyuy_sing_tangMvt3D(X, H,  M, d, U, l, eps, param)
        +    dzuz_sing_tangMvt3D(X, H,  M, d, U, l, eps, param);
    }
  }
}

#endif
