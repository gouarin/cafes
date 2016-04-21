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

    using position2d = geometry::position<2, double>;
    using position3d = geometry::position<3, double>;

    PetscScalar p_sing_withT_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double t3 = X[1]*X[1];
      double t8 = pow( H * t3 + M * t3 + 2 * a, 2);

      return -24 / (H + M) / t8 * U * mu * chiTrunc(X[1], l, eps);
    }


    PetscScalar uz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double r  = X[1];
      double z  = X[0];

      double t1 = r * U;
      double t2 = H + M;
      double t4 = 2 * z;
      double t5 = r * r;
      double t6 = M * t5;
      double t7 = 2 * a;
      double t8 = -t4 + t6 + t7;
      double t11 = H * t5;
      double t12 = t11 + t4;
      double t13 = 1. / t2;
      double t15 = t7 + t6 + t11;
      double t16 = t15 * t15;
      double t18 = 1. / t16 / t15;
      double t19 = chiTrunc(r, l, eps);
      double t23 = dchiTrunc(r, l, eps);
      return -12 * t1 * t2 * t8 * t12 * t13 * t18 * t19 - 6 * t1 * (-2 * t23 * a - t23 * M * t5 - t23 * H * t5) * t8 * t12 * t13 * t18;
    }


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

    PetscScalar dzuz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[0];
      double r   = X[1];

      double t1 = r * r;
      double t2 = H * t1;
      double t3 = a * M;
      double t4 = t3 * z;
      double t7 = H * H;
      double t8 = t7 * t1;
      double t9 = z * a;
      double t12 = t1 * t1;
      double t13 = t7 * t12;
      double t17 = H * t12;
      double t18 = M * M;
      double t22 = t18 * t12;
      double t23 = H * a;
      double t26 = t18 * t1;
      double t29 = a * a;
      double t30 = t29 * M;
      double t35 = t7 * H;
      double t39 = t12 * t1;
      double t40 = t35 * t39;
      double t43 = t7 * t39;
      double t46 = 80 * t2 * t4 + 64 * t8 * t9 - 12 * t13 * M * z + 12 * t17 * t18 * z - 8 * t22 * t23 + 16 * t26 * t9 - 24 * t2 * t30 + 4 * t13 * t3 + 12 * t35 * t12 * a + 2 * t40 * M + 4 * t43 * t18;
      double t47 = t18 * M;
      double t48 = t47 * t39;
      double t57 = H * t29 * z;
      double t59 = t30 * z;
      double t61 = z * z;
      double t65 = t23 * t61;
      double t67 = t3 * t61;
      double t71 = t35 * z;
      double t74 = t7 * t61;
      double t77 = 2 * t48 * H + 12 * t12 * t47 * z - 24 * t8 * t29 - 16 * t57 - 16 * t59 - 80 * t2 * M * t61 + 16 * t65 + 16 * t67 - 40 * t26 * t61 - 12 * t71 * t12 - 40 * t74 * t1;
      double t83 = M * t1;
      double t85 = pow((2 * a + t83 + t2), 2);
      double t86 = t85 * t85;
      double t88 = 1. / (H + M) / t86;
      double t89 = chiTrunc(r, l, eps);
      double t92 = dchiTrunc(r, l, eps);
      double t93 = t92 * t29;
      double t96 = t92 * a;
      double t99 = t12 * t12;
      double t101 = t92 * t18;
      double t105 = t92 * M;
      double t115 = t29 * a;
      double t116 = t92 * t115;
      double t119 = t92 * t12;
      double t120 = t18 * t61;
      double t131 = d2chiTrunc(r, l, eps);
      double t132 = t131 * t99;
      double t136 = t131 * t39;
      double t140 = -8 * t13 * t93 - 10 * t40 * t96 - 6 * t7 * t99 * t101 - 3 * t35 * t99 * t105 - 3 * t47 * t99 * t92 * H - 10 * t48 * t92 * z + 24 * t116 * t2 + 28 * t119 * t120 + 10 * t92 * t35 * t39 * z + 28 * t92 * t7 * t12 * t61 + 4 * t132 * t35 * a + 16 * t136 * t7 * t29;
      double t142 = t131 * t99 * t1;
      double t149 = t131 * t12;
      double t153 = t131 * t1;
      double t172 = t93 * z;
      double t175 = t96 * z;
      double t182 = 2 * t142 * t35 * M + 4 * t142 * t7 * t18 + 16 * t149 * t115 * H + 32 * t153 * t115 * z + 2 * t142 * t47 * H + 4 * t132 * t47 * z - 32 * t153 * t29 * t61 - 8 * t136 * t120 - 4 * t132 * t71 - 8 * t136 * t74 - 72 * t2 * t172 - 20 * t13 * t175 - 10 * H * t39 * t101 * z;
      double t202 = t92 * t1;
      double t207 = H * M;
      double t208 = t207 * t61;
      double t211 = a * t18;
      double t221 = 10 * t43 * t105 * z + 12 * M * t12 * t93 * H - 24 * t83 * t172 - 6 * t18 * t39 * t96 * H - 36 * t22 * t175 - 16 * t43 * t96 * M + 48 * t202 * t67 + 48 * t202 * t65 + 56 * t119 * t208 + 12 * t132 * t211 * H + 24 * t136 * t211 * z + 24 * t136 * t30 * H;
      double t258 = 48 * t149 * t59 - 8 * t136 * t7 * a * z + 16 * t149 * t57 - 4 * t132 * t7 * M * z + 4 * t132 * H * t18 * z + 16 * t132 * t3 * t7 - 32 * t149 * t67 - 32 * t149 * t65 - 16 * t136 * t208 + 16 * t116 * z - 16 * t93 * t61 - 56 * t17 * t92 * t4 + 16 * t136 * a * t207 * z;
      return 6 * U * (t46 + t77) * t88 * t89 + 6 * U * (t140 + t182 + t221 + t258) * t88;
    }


    PetscScalar dxux_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[0];
      double r   = X[1];

      double t1 = z * z;
      double t2 = r * r;
      double t3 = t1 * t2;
      double t4 = H * a;
      double t7 = a * M;
      double t10 = t2 * t2;
      double t11 = t1 * t10;
      double t12 = H * M;
      double t15 = t2 * H;
      double t16 = a * a;
      double t20 = H * H;
      double t21 = t10 * t20;
      double t22 = z * a;
      double t26 = M * z;
      double t30 = M * M;
      double t31 = t30 * z;
      double t34 = t10 * t2;
      double t35 = t34 * H;
      double t38 = t34 * t20;
      double t41 = t16 * M;
      double t42 = t10 * H;
      double t45 = a * t30;
      double t50 = 48 * t4 * t3 + 48 * t3 * t7 + 56 * t11 * t12 - 72 * t15 * t16 * z - 20 * t21 * t22 - 24 * t2 * t16 * t26 - 36 * t10 * a * t31 - 10 * t35 * t31 + 10 * t38 * t26 + 12 * t41 * t42 - 6 * t45 * t35 - 16 * t7 * t38;
      double t51 = t22 * M;
      double t54 = t1 * t16;
      double t56 = t16 * a;
      double t57 = t56 * z;
      double t59 = t56 * H;
      double t62 = t16 * t20;
      double t65 = t20 * H;
      double t66 = a * t65;
      double t69 = t30 * M;
      double t70 = t10 * t10;
      double t84 = t34 * t65;
      double t87 = t34 * t69;
      double t90 = -56 * t42 * t51 - 16 * t54 + 16 * t57 + 24 * t59 * t2 - 8 * t62 * t10 - 10 * t66 * t34 - 3 * t69 * t70 * H - 6 * t30 * t70 * t20 - 3 * M * t70 * t65 + 28 * t11 * t20 + 28 * t11 * t30 + 10 * t84 * z - 10 * t87 * z;
      double t98 = pow((2 * a + M * t2 + t15), 2);
      double t99 = t98 * t98;
      double t101 = 0.1e1 / (H + M) / t99;
      double t102 = dchiTrunc(r, l, eps);
      double t107 = t20 * t2;
      double t117 = t2 * t30;
      double t131 = 80 * t15 * t51 + 64 * t107 * t22 - 12 * t21 * t26 + 12 * t42 * t31 - 8 * t30 * t10 * t4 + 16 * t117 * t22 - 24 * t15 * t41 + 4 * t21 * t7 + 12 * t65 * t10 * a + 2 * M * t84 + 4 * t38 * t30;
      double t140 = H * t16 * z;
      double t142 = t41 * z;
      double t147 = t4 * t1;
      double t149 = t7 * t1;
      double t153 = t65 * z;
      double t156 = t20 * t1;
      double t159 = 2 * t87 * H + 12 * t10 * t69 * z - 24 * t107 * t16 - 16 * t140 - 16 * t142 - 80 * t15 * M * t1 + 16 * t147 + 16 * t149 - 40 * t117 * t1 - 12 * t153 * t10 - 40 * t156 * t2;
      double t162 = chiTrunc(r, l, eps);
      double t165 = d2chiTrunc(r, l, eps);
      double t166 = t165 * t70;
      double t169 = t165 * t34;
      double t173 = t165 * t70 * t2;
      double t180 = t165 * t10;
      double t183 = t165 * t2;
      double t201 = 4 * t166 * t66 + 16 * t169 * t62 + 2 * t173 * t65 * M + 4 * t173 * t20 * t30 + 16 * t180 * t59 + 32 * t183 * t57 + 2 * t173 * t69 * H + 4 * t166 * t69 * z - 32 * t183 * t54 - 8 * t169 * t30 * t1 - 4 * t166 * t153 - 8 * t169 * t156;
      double t241 = 12 * t166 * t45 * H + 24 * t169 * t45 * z + 24 * t169 * t41 * H + 48 * t180 * t142 - 8 * t169 * t20 * a * z + 16 * t180 * t140 - 4 * t166 * t20 * M * z + 4 * t166 * H * t30 * z + 16 * t166 * t7 * t20 - 32 * t180 * t149 - 32 * t180 * t147 - 16 * t169 * t12 * t1 + 16 * t169 * a * t12 * z;
      return -6 * U * (t50 + t90) * t101 * t102 - 6 * U * (t131 + t159) * t101 * t162 - 6 * U * (t201 + t241) * t101;


    }

    PetscScalar dxuz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double r   = X[1];
      double z   = X[0];

      double t1 = r * U;
      double t2 = 2 * a;
      double t3 = r * r;
      double t4 = M * t3;
      double t5 = H * t3;
      double t6 = -t2 - t4 - t5;
      double t10 = H + M;
      double t12 = (t5 + 4 * z - t4 - t2) / t10;
      double t13 = -t6;
      double t14 = t13 * t13;
      double t16 = 0.1e1 / t14 / t13;
      double t17 = dchiTrunc(r, l, eps);
      double t23 = chiTrunc(r, l, eps);
      return 12 * t1 * t6 * t12 * t16 * t17 + 24 * t1 * t10 * t12 * t16 * t23;
    }

    PetscScalar dzux_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double r   = X[1];
      double z   = X[0];

      double t1 = r * U;
      double t2 = r * r;
      double t3 = t2 * H;
      double t4 = a * M;
      double t5 = z * z;
      double t6 = t5 * z;
      double t7 = t4 * t6;
      double t10 = t2 * t2;
      double t11 = t10 * H;
      double t12 = M * M;
      double t13 = t12 * a;
      double t14 = t13 * t5;
      double t17 = H * H;
      double t18 = t10 * t17;
      double t19 = t4 * t5;
      double t22 = a * a;
      double t23 = H * t22;
      double t24 = M * t2;
      double t25 = t24 * t5;
      double t28 = t22 * a;
      double t29 = M * t28;
      double t30 = t29 * z;
      double t33 = t10 * t2;
      double t34 = t17 * t33;
      double t35 = t13 * z;
      double t38 = t17 * H;
      double t39 = t38 * t33;
      double t40 = t4 * z;
      double t43 = t22 * M;
      double t44 = t43 * z;
      double t47 = H * t33;
      double t48 = t12 * M;
      double t49 = t48 * a;
      double t50 = t49 * z;
      double t53 = H * t12;
      double t58 = t17 * t17;
      double t59 = t10 * t10;
      double t60 = t58 * t59;
      double t61 = M * z;
      double t64 = t17 * t2;
      double t68 = t58 * t33;
      double t69 = z * a;
      double t72 = H * t59;
      double t73 = t12 * t12;
      double t74 = t73 * z;
      double t77 = t38 * t10;
      double t78 = t22 * z;
      double t81 = -768 * t3 * t7 + 1464 * t11 * t14 + 1992 * t18 * t19 + 1152 * t23 * t25 - 864 * t3 * t30 - 264 * t34 * t35 + 264 * t39 * t40 - 1728 * t18 * t44 - 264 * t47 * t50 - 1200 * t53 * t10 * t22 * z + 24 * t60 * t61 - 1824 * t64 * t28 * z + 264 * t68 * t69 + 24 * t72 * t74 - 528 * t77 * t78;
      double t82 = t38 * t59;
      double t83 = t12 * z;
      double t86 = t17 * t59;
      double t87 = t48 * z;
      double t90 = t2 * t12;
      double t91 = a * t6;
      double t94 = t12 * t6;
      double t97 = M * t6;
      double t102 = t22 * t5;
      double t105 = t48 * t5;
      double t108 = t10 * t48;
      double t109 = a * t5;
      double t114 = M * t5;
      double t119 = t22 * t22;
      double t120 = H * t119;
      double t123 = H * t28;
      double t130 = H * a;
      double t134 = 72 * t82 * t83 + 72 * t86 * t87 - 384 * t90 * t91 - 1152 * t11 * t94 - 1152 * t18 * t97 - 384 * t64 * t91 + 1488 * t64 * t102 + 264 * t47 * t105 + 312 * t108 * t109 - 336 * t90 * t102 - 264 * t39 * t114 + 840 * t77 * t109 + 160 * t120 * t24 + 240 * t123 * t12 * t10 + 160 * t23 * t48 * t33 + 50 * t130 * t73 * t59;
      double t144 = t22 * t12;
      double t149 = t58 * t17;
      double t150 = t59 * t2;
      double t153 = t119 * a;
      double t160 = t119 * t17;
      double t163 = t28 * t38;
      double t166 = t22 * t58;
      double t169 = t58 * H;
      double t170 = a * t169;
      double t179 = -14 * t49 * t86 - 158 * t13 * t82 - 124 * t4 * t60 + 408 * t29 * t18 + 12 * t144 * t34 - 320 * t43 * t39 - 2 * t149 * t150 + 32 * H * t153 + t73 * t12 * t150 + 32 * M * t153 + 560 * t160 * t2 - 88 * t163 * t10 - 212 * t166 * t33 - 20 * t170 * t59 + 6 * t73 * t150 * t17 - 10 * t48 * t150 * t38;
      double t186 = t73 * M;
      double t206 = t23 * t6;
      double t208 = t43 * t6;
      double t210 = t123 * t5;
      double t212 = t29 * t5;
      double t221 = -21 * t12 * t150 * t58 - 12 * M * t150 * t169 + 6 * H * t186 * t150 + 80 * t119 * t12 * t2 + 80 * t28 * t48 * t10 + 40 * t22 * t73 * t33 + 10 * a * t186 * t59 - 384 * t77 * t6 - 384 * t108 * t6 + 768 * t206 + 768 * t208 - 1440 * t210 - 864 * t212 - 132 * t68 * t5 + 132 * t33 * t73 * t5 + 576 * t120 * z;
      double t226 = 0.1e1 / (H + M);
      double t228 = 2 * a + t24 + t3;
      double t229 = t228 * t228;
      double t230 = t229 * t229;
      double t232 = 0.1e1 / t230 / t228;
      double t233 = t226 * t232;
      double t234 = dchiTrunc(r, l, eps);
      double t243 = t144 * z;
      double t256 = t17 * t22 * t5;
      double t258 = t144 * t5;
      double t263 = t2 * t38;
      double t266 = t17 * a;
      double t267 = t266 * t6;
      double t272 = t13 * t6;
      double t276 = t58 * t10;
      double t283 = -1440 * t3 * t14 - 2880 * t64 * t19 + 336 * t11 * t50 + 480 * t3 * t243 + 384 * t18 * t35 - 240 * t77 * t40 + 1920 * t64 * t44 - 144 * t10 * t73 * t5 + 1152 * t256 + 576 * t258 - 576 * t17 * t28 * z + 480 * t263 * t6 - 576 * t267 + 480 * t2 * t48 * t6 - 576 * t272 - 480 * t263 * t28 + 144 * t276 * t22 + 144 * t276 * t5 - 24 * t68 * t61;
      double t300 = H * M;
      double t325 = -24 * t47 * t74 + 1440 * t263 * t78 - 288 * t276 * t69 - 72 * t39 * t83 - 72 * t34 * t87 - 576 * t123 * t61 + 1440 * t3 * t94 + 1440 * t64 * t97 - 1152 * t300 * t91 + 24 * t58 * M * a * t33 + 48 * t13 * t39 - 48 * t77 * t43 + 1728 * t300 * t102 - 288 * t11 * t105 - 192 * t144 * t18 - 1440 * t263 * t109 + 288 * t77 * t114 - 480 * t64 * t29 + 24 * t34 * t49;
      double t328 = chiTrunc(r, l, eps);
      double t331 = d3chiTrunc(r, l, eps);
      double t332 = t331 * t10;
      double t337 = t331 * t150;
      double t338 = t43 * t38;
      double t341 = t331 * t33;
      double t346 = t59 * t10;
      double t347 = t331 * t346;
      double t348 = t13 * t38;
      double t351 = t331 * t59;
      double t353 = M * t17 * t6;
      double t356 = t53 * t6;
      double t359 = t4 * t58;
      double t371 = -384 * t332 * t206 - 384 * t332 * t208 + 240 * t337 * t338 - 192 * t341 * t267 - 192 * t341 * t272 + 120 * t347 * t348 - 96 * t351 * t353 - 96 * t351 * t356 + 96 * t347 * t359 + 768 * t332 * t212 + 192 * t351 * t29 * t17 + 384 * t332 * t210 + 576 * t341 * t258;
      double t372 = t144 * t17;
      double t375 = t49 * t5;
      double t378 = t49 * t17;
      double t382 = t48 * H * t5;
      double t386 = t38 * M * t5;
      double t390 = t38 * a * t5;
      double t398 = t58 * t150;
      double t425 = 144 * t337 * t372 + 192 * t351 * t375 + 48 * t347 * t378 + 48 * t337 * t382 - 48 * t337 * t386 - 96 * t351 * t390 + 24 * t58 * t346 * t331 * M * z + 48 * t398 * t331 * a * z + 384 * t11 * t331 * t119 * z + 24 * H * t346 * t331 * t73 * z + 576 * t34 * t331 * t28 * z + 288 * t82 * t331 * t22 * z + 72 * t38 * t346 * t331 * t12 * z;
      double t432 = d2chiTrunc(r, l, eps);
      double t433 = t432 * t2;
      double t436 = t432 * t59;
      double t439 = t432 * t33;
      double t444 = t432 * t10;
      double t449 = t432 * t150;
      double t460 = t432 * t28;
      double t463 = 72 * t17 * t346 * t331 * t48 * z + 192 * t433 * t208 + 744 * t436 * t338 + 432 * t439 * t353 + 432 * t439 * t356 + 480 * t444 * t267 + 192 * t433 * t206 + 192 * t449 * t359 - 288 * t439 * t375 + 120 * t449 * t378 - 288 * t444 * t258 + 504 * t436 * t372 + 384 * t460 * t25;
      double t478 = t432 * t119;
      double t490 = H * t150;
      double t495 = t460 * z;
      double t498 = t38 * t150;
      double t507 = 864 * t460 * M * t33 * t17 - 1152 * t444 * t256 - 120 * t436 * t382 - 960 * t433 * t210 + 120 * t436 * t386 - 144 * t439 * t390 + 1344 * t3 * t478 * z - 12 * t398 * t432 * M * z - 120 * t60 * t432 * a * z - 12 * t490 * t432 * t73 * z + 864 * t18 * t495 - 36 * t498 * t432 * t12 * z - 144 * t39 * t432 * t22 * z;
      double t510 = t17 * t150;
      double t519 = t444 * H;
      double t526 = t439 * H;
      double t549 = -36 * t510 * t432 * t48 * z + 480 * t444 * t272 + 276 * t449 * t348 + 960 * t519 * t7 - 576 * t439 * M * t266 * t5 - 720 * t526 * t14 - 1440 * t519 * t43 * t5 - 144 * t82 * t432 * t40 + 72 * t86 * t432 * t35 + 720 * t34 * t432 * t44 + 96 * t72 * t432 * t50 + 864 * t526 * t243 + 1920 * t300 * t10 * t495;
      double t582 = t432 * t346;
      double t588 = t331 * t59 * t33;
      double t591 = -384 * t341 * a * t300 * t6 + 576 * t341 * t22 * t300 * t5 + 288 * t351 * t12 * t130 * t5 + 864 * t86 * t331 * t44 + 288 * t498 * t331 * t40 + 432 * t510 * t331 * t35 + 768 * t47 * t331 * t30 + 576 * t72 * t331 * t243 + 192 * t490 * t331 * t50 - 384 * t460 * t6 + 3 * t582 * t149 + 576 * t478 * t5 + 2 * t588 * t149;
      double t597 = M * t169;
      double t600 = t12 * t58;
      double t603 = t331 * t2;
      double t609 = t48 * t6;
      double t612 = t48 * t38;
      double t615 = t38 * t6;
      double t623 = t73 * t5;
      double t626 = t73 * t17;
      double t629 = 96 * t337 * t166 + 24 * t347 * t170 + 12 * t588 * t597 + 24 * t588 * t600 - 256 * t603 * t28 * t6 + 160 * t351 * t163 - 32 * t351 * t609 + 20 * t588 * t612 - 32 * t351 * t615 + 384 * t603 * t119 * t5 + 96 * t341 * t160 + 24 * t337 * t623 + 6 * t588 * t626;
      double t658 = -24 * t398 * t331 * t5 + 624 * t460 * t39 + 18 * t582 * t597 + 36 * t582 * t600 + 36 * t449 * t170 + 240 * t436 * t166 + 144 * t439 * t615 + 144 * t439 * t609 + 30 * t582 * t612 + 528 * t478 * t18 + 60 * t60 * t432 * t5 - 60 * t436 * t623 + 9 * t582 * t626;
      return -t1 * (t81 + t134 + t179 + t221) * t233 * t234 - t1 * (t283 + t325) * t233 * t328 - t1 * (t371 + t425 + t463 + t507 + t549 + t591 + t629 + t658) * t226 * t232;

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
