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

      double t1 = r * U;//std::cout << "t1 = " << t1 << std::endl;
      double t2 = r * r;//std::cout << "t2 = " << t2 << std::endl;
      double t3 = H * t2;//std::cout << "t3 = " << t3 << std::endl;
      double t4 = M * t2;//std::cout << "t4 = " << t4 << std::endl;
      double t5 = 2 * a;//std::cout << "t5 = " << t5 << std::endl;
      double t6 = -t3 - t4 - t5;//std::cout << "t6 = " << t6 << std::endl;
      double t7 = 2 * z;//std::cout << "t7 = " << t7 << std::endl;
      double t8 = t4 + t5 - t7;//std::cout << "t8 = " << t8 << std::endl;
      double t11 = t3 + t7;//std::cout << "t11= " << t11<< std::endl;
      double t12 = dchiTrunc(X[1], l, eps);//std::cout << "t12= " << t12 << std::endl;
      double t14 = H + M;//std::cout << "t14= " << t14<< std::endl;
      double t15 = 1. / t14;//std::cout << "t15= " << t15<< std::endl;
      double t16 = -t6;//std::cout << "t16= " << t16<< std::endl;
      double t17 = t16 * t16;//std::cout << "t17= " << t17<< std::endl;
      double t19 = 1. / t17 / t16;//std::cout << "t19= " << t19<< std::endl;
      double t27 = chiTrunc(X[1], l, eps);//std::cout << "t27= " << t27<< std::endl;
      //std::cout << "return " << -6 * t1 * t6 * t8 * t11 * t12 * t15 * t19 - 12 * t1 * t14 * t8 * t11 * t15 * t19 * t27 << std::endl;
      return -6 * t1 * t6 * t8 * t11 * t12 * t15 * t19 - 12 * t1 * t14 * t8 * t11 * t15 * t19 * t27;

    }


    //Attention : r cordonnee radial, donc z en 2D et donc z = x .....
    PetscScalar ux_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu = 1;
      double r  = X[1];
      double z  = X[0];

      double t1 = H * H;
      double t2 = t1 * t1;
      double t3 = t2 * H;
      double t4 = r * r;
      double t5 = t4 * t4;
      double t6 = t5 * t5;
      double t7 = t6 * t5;
      double t10 = t2 * M;
      double t13 = t1 * H;
      double t14 = M * M;
      double t15 = t13 * t14;
      double t18 = t14 * M;
      double t19 = t1 * t18;
      double t22 = t2 * a;
      double t23 = t6 * t4;
      double t26 = t13 * M;
      double t27 = a * t23;
      double t30 = t23 * z;
      double t33 = t1 * t14;
      double t38 = H * t18;
      double t41 = a * a;
      double t42 = t13 * t41;
      double t45 = t13 * a;
      double t46 = t6 * z;
      double t50 = z * z;
      double t53 = t1 * M;
      double t57 = a * t6;
      double t58 = t57 * z;
      double t61 = t6 * t50;
      double t64 = H * t14;
      double t67 = 24 * t13 * t6 * t50 - 72 * t53 * t41 * t6 - 10 * t10 * t7 - 14 * t15 * t7 - 6 * t19 * t7 - 20 * t22 * t23 - 56 * t26 * t27 - 24 * t26 * t30 - 36 * t33 * t27 - 2 * t3 * t7 - 48 * t33 * t30 - 24 * t38 * t30 - 56 * t42 * t6 - 48 * t45 * t46 - 192 * t53 * t58 + 24 * t53 * t61 - 144 * t64 * t58;
      double t73 = t41 * a;
      double t74 = t1 * t73;
      double t75 = t5 * t4;
      double t78 = t1 * t41;
      double t79 = t75 * z;
      double t82 = t1 * a;
      double t83 = t75 * t50;
      double t87 = t50 * z;
      double t90 = H * M;
      double t91 = t41 * t75;
      double t95 = a * t75;
      double t102 = t14 * a;
      double t108 = H * t73;
      double t109 = t5 * z;
      double t112 = H * t41;
      double t113 = t5 * t50;
      double t116 = H * a;
      double t117 = t5 * t87;
      double t120 = M * t41;
      double t123 = M * a;
      double t126 = t73 * t4;
      double t129 = t41 * t4;
      double t132 = 32 * t1 * t75 * t87 + 32 * t14 * t75 * t87 - 24 * t18 * t6 * t50 - 96 * t90 * t95 * t50 + 64 * t90 * t75 * t87 - 288 * t90 * t91 * z - 144 * t102 * t83 - 192 * t108 * t109 - 96 * t112 * t113 - 288 * t120 * t113 + 128 * t116 * t117 + 128 * t123 * t117 - 192 * t126 * t50 + 128 * t129 * t87 - 24 * t64 * t61 - 48 * t74 * t75 - 192 * t78 * t79 + 48 * t82 * t83;
      double t135 = d2chiTrunc(r, l, eps);
      double t137 = 1. / (H + M);
      double t143 = pow(H * t4 + M * t4 + 2 * a, 2);
      double t144 = t143 * t143;
      double t145 = 1. / t144;
      double t176 = t95 * z;
      double t183 = -60 * t13 * t75 * t50 + 5 * t10 * t23 + 7 * t15 * t23 + 192 * t53 * t176 + 72 * t64 * t176 + 3 * t19 * t23 + 6 * t22 * t6 + t3 * t23 + 36 * t26 * t46 - 12 * t26 * t57 + 72 * t33 * t46 - 18 * t33 * t57 + 36 * t38 * t46 - 52 * t42 * t75 + 120 * t45 * t79 - 60 * t53 * t83 - 108 * t53 * t91;
      double t198 = t41 * t5;
      double t202 = a * t5;
      double t213 = t4 * z;
      double t216 = t4 * t50;
      double t219 = t4 * t87;
      double t230 = -112 * t1 * t5 * t87 - 112 * t14 * t5 * t87 + 60 * t18 * t75 * t50 - 144 * t90 * t198 * z + 336 * t90 * t202 * t50 + 216 * t102 * t113 - 288 * t108 * t213 + 96 * t78 * t109 + 432 * t112 * t216 + 120 * t82 * t113 - 192 * t116 * t219 - 224 * t90 * t117 + 144 * t120 * t216 - 192 * t123 * t219 + 64 * t41 * t87 - 120 * t74 * t5 - 96 * t73 * t50 + 60 * t64 * t83;
      double t233 = dchiTrunc(r, l, eps);
      double t244 = t14 * t14;
      double t277 = t202 * z;
      double t284 = -8 * t244 * a * t75 + 72 * t13 * t5 * t50 - 144 * t45 * t109 + 72 * t53 * t113 + 48 * t53 * t198 - 72 * t64 * t198 - 48 * t53 * t277 - 48 * t33 * t79 - 24 * t38 * t79 - 32 * t38 * t95 + 96 * t42 * t5;
      double t325 = t41 * t41;
      double t338 = -480 * t90 * a * t4 * t50 - 32 * t14 * t73 * t4 + 160 * t14 * t4 * t87 - 16 * H * t325 - 16 * M * t325 - 96 * t102 * t216 + 96 * t112 * t50 - 64 * t116 * t87 + 96 * t120 * t50 - 64 * t123 * t87 + 320 * t90 * t219;
      double t342 = chiTrunc(r, l, eps);
      return  -U * (t132 + t67) * t135 * t137 * t145 / 2. + U * (t230 + t183) * t233 * t137 * t145 / 2. + U * (160 * t1 * t4 * t87 + 288 * t90 * t129 * z - 24 * t18 * t41 * t5 - 72 * t18 * t5 * t50 - 72 * t64 * t113 - 64 * t90 * t126 + 288 * t78 * t213 - 384 * t82 * t216 + 96 * t64 * t277 - 32 * t74 * t4 + t338 - 5 * H * t244 * t6 - t244 * M * t6 + 5 * t10 * t6 + 4 * t15 * t6 - 4 * t19 * t6 + 8 * t22 * t75 - 24 * t26 * t79 + 56 * t26 * t95 + t3 * t6 + 24 * t33 * t95 + t284) * t342 * t137 * t145 / 2.;
    }

    PetscScalar dzuz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[0];
      double r   = X[1];

      double t1 = H * H;
      double t2 = t1 * H;
      double t3 = t2 * M;
      double t4 = r * r;
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
      return  6 * U * (t94 + t48) * t97 * t99 * t107 + 6 * U * (t170 + t136) * t173 * t99 * t107 + 6 * U * (t228 + t201) * t99 * t107 * t232;

    }


    PetscScalar dxux_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double z   = X[0];
      double r   = X[1];

      double t1 = H * H;
      double t2 = t1 * H;
      double t3 = t2 * M;
      double t4 = r * r;
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
      return -6 * U * (t94 + t48) * t97 * t99 * t107 - 6 * U * (t170 + t136) * t173 * t99 * t107 - 6 * U * (t228 + t201) * t99 * t107 * t232;

    }

    PetscScalar dxuz_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double r   = X[1];
      double z   = X[0];

      double t1 = r * U;
      double t2 = r * r;
      double t3 = H * t2;
      double t4 = M * t2;
      double t5 = 2 * a;
      double t6 = -t3 - t4 - t5;
      double t9 = t3 - t4 - t5 + 4 * z;
      double t10 = dchiTrunc(r, l, eps);
      double t12 = H + M;
      double t13 = 0.1e1 / t12;
      double t14 = -t6;
      double t15 = t14 * t14;
      double t17 = 0.1e1 / t15 / t14;
      double t24 = chiTrunc(r, l, eps);
      return 12 * t1 * t6 * t9 * t10 * t13 * t17 + 24 * t1 * t12 * t9 * t13 * t17 * t24;

    }

    PetscScalar dzux_sing_normalMvt2D(position2d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double mu  = 1;
      double r   = X[1];
      double z   = X[0];

      double t1 = r * U;
      double t2 = H * H;
      double t3 = t2 * t2;
      double t4 = t3 * t2;
      double t5 = r * r;
      double t6 = t5 * t5;
      double t7 = t6 * t5;
      double t8 = t6 * t6;
      double t9 = t8 * t7;
      double t12 = t3 * H;
      double t13 = t12 * M;
      double t16 = M * M;
      double t17 = t3 * t16;
      double t20 = t2 * H;
      double t21 = t16 * M;
      double t22 = t20 * t21;
      double t25 = t16 * t16;
      double t26 = t2 * t25;
      double t29 = t12 * a;
      double t30 = t8 * t6;
      double t33 = t3 * M;
      double t34 = a * t30;
      double t37 = t30 * z;
      double t40 = t20 * t16;
      double t45 = t2 * t21;
      double t50 = -12 * t13 * t9 - 24 * t17 * t9 - 20 * t22 * t9 - 6 * t26 * t9 - 24 * t29 * t30 - 96 * t33 * t34 - 24 * t33 * t37 - 120 * t40 * t34 - 48 * t45 * t34 - 72 * t40 * t37 - 72 * t45 * t37 - 2 * t4 * t9;
      double t51 = H * t25;
      double t54 = a * a;
      double t55 = t3 * t54;
      double t56 = t8 * t5;
      double t59 = t3 * a;
      double t60 = t56 * z;
      double t64 = z * z;
      double t67 = t20 * M;
      double t68 = t54 * t56;
      double t71 = a * t56;
      double t72 = t71 * z;
      double t75 = t56 * t64;
      double t78 = t2 * t16;
      double t83 = H * t21;
      double t91 = t54 * a;
      double t92 = t20 * t91;
      double t95 = -24 * t25 * t56 * t64 + 24 * t3 * t56 * t64 - 24 * t51 * t37 - 96 * t55 * t56 - 48 * t59 * t60 - 240 * t67 * t68 - 288 * t67 * t72 + 48 * t67 * t75 - 144 * t78 * t68 - 432 * t78 * t72 - 192 * t83 * t72 - 48 * t83 * t75 - 160 * t92 * t8;
      double t97 = t20 * t54;
      double t98 = t8 * z;
      double t101 = t20 * a;
      double t102 = t8 * t64;
      double t106 = t64 * z;
      double t109 = t2 * M;
      double t113 = t54 * t8;
      double t114 = t113 * z;
      double t117 = t8 * t106;
      double t120 = H * t16;
      double t123 = a * t8;
      double t129 = t21 * a;
      double t135 = t54 * t54;
      double t136 = t2 * t135;
      double t139 = t2 * t91;
      double t140 = t7 * z;
      double t143 = 32 * t20 * t8 * t106 + 32 * t21 * t8 * t106 - 192 * t109 * t91 * t8 - 288 * t120 * t123 * t64 + 96 * t101 * t102 - 192 * t129 * t102 - 864 * t109 * t114 + 96 * t109 * t117 - 576 * t120 * t114 + 96 * t120 * t117 - 96 * t136 * t7 - 576 * t139 * t140 - 288 * t97 * t98;
      double t144 = t2 * a;
      double t145 = t7 * t106;
      double t148 = H * M;
      double t149 = t91 * t7;
      double t153 = t54 * t7;
      double t157 = a * t7;
      double t161 = t16 * t54;
      double t162 = t7 * t64;
      double t165 = t16 * a;
      double t168 = H * t135;
      double t169 = t6 * z;
      double t172 = H * t91;
      double t173 = t6 * t64;
      double t176 = H * t54;
      double t177 = t6 * t106;
      double t180 = M * t91;
      double t183 = M * t54;
      double t186 = t135 * t5;
      double t189 = t91 * t5;
      double t192 = 384 * t148 * t157 * t106 - 768 * t148 * t149 * z - 576 * t148 * t153 * t64 + 256 * t189 * t106 + 192 * t144 * t145 + 192 * t165 * t145 - 576 * t161 * t162 - 384 * t168 * t169 - 384 * t172 * t173 - 768 * t180 * t173 + 384 * t176 * t177 + 384 * t183 * t177 - 384 * t186 * t64;
      double t196 = d3chiTrunc(r, l, eps);
      double t198 = 0.1e1 / (H + M);
      double t203 = H * t5 + M * t5 + 2 * a;
      double t204 = t203 * t203;
      double t205 = t204 * t204;
      double t207 = 0.1e1 / t205 / t203;
      double t236 = -18 * t13 * t30 - 36 * t17 * t30 - 30 * t22 * t30 - 9 * t26 * t30 - 36 * t29 * t56 - 3 * t4 * t30 + 12 * t33 * t60 - 192 * t33 * t71 + 36 * t40 * t60 - 276 * t40 * t71 + 36 * t45 * t60 - 120 * t45 * t71 + 12 * t51 * t60;
      double t246 = t123 * z;
      double t266 = 60 * t25 * t8 * t64 - 60 * t3 * t8 * t64 - 120 * t67 * t102 + 120 * t83 * t102 - 744 * t67 * t113 - 504 * t78 * t113 + 144 * t97 * t140 + 144 * t67 * t246 - 72 * t78 * t246 - 96 * t83 * t246 - 240 * t55 * t8 + 120 * t59 * t98 - 624 * t92 * t7;
      double t275 = t153 * z;
      double t278 = t157 * t64;
      double t298 = -144 * t20 * t7 * t106 - 144 * t21 * t7 * t106 + 144 * t101 * t162 - 432 * t109 * t145 - 864 * t109 * t149 - 720 * t109 * t275 + 576 * t109 * t278 - 432 * t120 * t145 - 864 * t120 * t275 + 720 * t120 * t278 + 288 * t129 * t162 - 528 * t136 * t6 - 864 * t139 * t169;
      double t299 = t2 * t54;
      double t304 = t91 * t6;
      double t308 = t54 * t6;
      double t312 = a * t6;
      double t320 = t5 * z;
      double t323 = t5 * t64;
      double t326 = t5 * t106;
      double t337 = -960 * t148 * t312 * t106 - 1920 * t148 * t304 * z + 1440 * t148 * t308 * t64 + 384 * t91 * t106 - 576 * t135 * t64 - 480 * t144 * t177 + 288 * t161 * t173 - 480 * t165 * t177 - 1344 * t168 * t320 + 960 * t172 * t323 + 1152 * t299 * t173 - 192 * t176 * t326 - 384 * t180 * t323 - 192 * t183 * t326;
      double t341 = d2chiTrunc(r, l, eps);
      double t355 = t25 * M;
      double t377 = -6 * H * t355 * t56 - t25 * t16 * t56 + 124 * t33 * t123 + 158 * t40 * t123 + 14 * t45 * t123 - 50 * t51 * t123 + 12 * t13 * t56 + 21 * t17 * t56 + 10 * t22 * t56 - 6 * t26 * t56 + 20 * t29 * t8 - 24 * t33 * t98 + 2 * t4 * t56 - 72 * t40 * t98 - 72 * t45 * t98;
      double t392 = t157 * z;
      double t415 = -10 * t355 * a * t8 - 40 * t25 * t54 * t7 - 132 * t25 * t7 * t64 + 132 * t3 * t7 * t64 - 264 * t59 * t140 + 320 * t67 * t153 - 12 * t78 * t153 - 160 * t83 * t153 + 264 * t67 * t162 - 264 * t83 * t162 - 264 * t67 * t392 + 264 * t78 * t392 + 264 * t83 * t392 - 24 * t51 * t98 + 212 * t55 * t7 + 88 * t92 * t6;
      double t426 = t308 * z;
      double t429 = t312 * t64;
      double t454 = 384 * t20 * t6 * t106 + 384 * t21 * t6 * t106 - 80 * t21 * t91 * t6 - 840 * t101 * t173 + 1152 * t109 * t177 - 408 * t109 * t304 + 1728 * t109 * t426 - 1992 * t109 * t429 + 1152 * t120 * t177 - 240 * t120 * t304 + 1200 * t120 * t426 - 1464 * t120 * t429 - 312 * t129 * t173 - 560 * t136 * t5 + 1824 * t139 * t320 + 528 * t97 * t169;
      double t464 = t54 * t5;
      double t468 = a * t5;
      double t479 = t135 * a;
      double t494 = 768 * t148 * t468 * t106 - 80 * t16 * t135 * t5 + 864 * t148 * t189 * z - 1152 * t148 * t464 * t64 - 32 * H * t479 - 32 * M * t479 - 768 * t176 * t106 - 768 * t183 * t106 + 384 * t144 * t326 - 160 * t148 * t186 + 336 * t161 * t323 + 384 * t165 * t326 - 576 * t168 * z + 1440 * t172 * t64 + 864 * t180 * t64 - 1488 * t299 * t323;
      double t498 = dchiTrunc(r, l, eps);
      double t525 = t312 * z;
      double t543 = 144 * t25 * t6 * t64 - 144 * t3 * t6 * t64 + 24 * t33 * t140 + 72 * t40 * t140 + 72 * t45 * t140 + 24 * t51 * t140 - 24 * t33 * t157 - 48 * t40 * t157 - 24 * t45 * t157 + 288 * t59 * t169 - 288 * t67 * t173 + 288 * t83 * t173 + 48 * t67 * t308 + 192 * t78 * t308 + 480 * t92 * t5 + 240 * t67 * t525 - 384 * t78 * t525 - 336 * t83 * t525 - 144 * t55 * t6;
      double t553 = t464 * z;
      double t556 = t468 * t64;
      double t589 = 1152 * t148 * a * t106 - 480 * t20 * t5 * t106 - 480 * t21 * t5 * t106 - 1728 * t148 * t54 * t64 + 576 * t148 * t91 * z + 1440 * t101 * t323 + 576 * t144 * t106 + 576 * t165 * t106 + 480 * t109 * t189 - 1440 * t109 * t326 - 1920 * t109 * t553 + 2880 * t109 * t556 - 1440 * t120 * t326 - 480 * t120 * t553 + 1440 * t120 * t556 + 576 * t139 * z - 576 * t161 * t64 - 1152 * t299 * t64 - 1440 * t97 * t320;
      double t593 = chiTrunc(r, l, eps);
      return t1 * (t192 + t143 + t95 + t50) * t196 * t198 * t207 + t1 * (t337 + t298 + t266 + t236) * t341 * t198 * t207 + t1 * (t494 + t454 + t415 + t377) * t498 * t198 * t207 + t1 * (t589 + t543) * t198 * t207 * t593;
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
      double t11 = pow(t2 * H + t2 * M + 2 * a, 2);
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

      double t1 = r * U;
      double t3 = H * t2;
      double t4 = t2 * M;
      double t5 = 2 * a;
      double t7 = t3 - t4 - t5 + 4 * z;
      double t9 = dchiTrunc(r, l, eps);
      double t10 = t3 + t4 + t5;
      double t11 = t10 * t10;
      double t14 = H + M;
      double t15 = 0.1e1 / t14;
      double t24 = chiTrunc(r, l, eps);
      return -6 * t1 * t7 * t9 / t11 * t15 + 12 * t1 * t14 * t7 * t15 / t11 / t10 * t24;
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

      // si r=0 on a un inf .... a discuter
      // if(r2==0.0)
      //   {
      // 	r2 = 0.0000000001;
      // 	r  = 0.00001;
      //   }

      return drur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[0]/r2 + ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[1]*X[1]/r2/r;
    }

    PetscScalar dyux_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      // si r=0 on a un inf .... a discuter
      // if(r2==0)
      //   {
      // 	r2 = 0.0000000001;
      // 	r  = 0.00001;
      //   }

      return drur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[1]/r2 - ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[1]/r2/r;
    }

    PetscScalar dzux_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      // si r=0 on a un inf .... a discuter
      // if(r2==0)
      //   {
      // 	r2 = 0.0000000001;
      // 	r  = 0.00001;
      //   }

      return dzur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]/r;
    }

    PetscScalar dxuy_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      // si r=0 on a un inf .... a discuter
      // if(r2==0)
      //   {
      // 	r2 = 0.0000000001;
      // 	r  = 0.00001;
      //   }

      return drur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[1]/r2 - ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[1]/r2/r;
    }

    PetscScalar dyuy_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      // si r=0 on a un inf .... a discuter
      // if(r2==0)
      //   {
      // 	r2 = 0.0000000001;
      // 	r  = 0.00001;
      //   }

      return drur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[1]*X[1]/r2 + ur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]*X[0]/r2/r;
    }

    PetscScalar dzuy_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      // si r=0 on a un inf .... a discuter
      // if(r2==0)
      //   {
      // 	r2 = 0.0000000001;
      // 	r  = 0.00001;
      //   }

      return dzur_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[1]/r;
    }

    PetscScalar dxuz_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);

      // si r=0 on a un inf .... a discuter
      // if(r2==0)
      //   {
      // 	r2 = 0.0000000001;
      // 	r  = 0.00001;
      //   }

      return druz_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[0]/r;
    }


    PetscScalar dyuz_sing_normalMvt3D(position3d X, PetscReal H, PetscReal M, PetscReal a, double U, PetscReal l, PetscReal eps, PetscReal* param  )
    {
      double r2  = X[0]*X[0]+X[1]*X[1];
      double r   = sqrt(r2);
        
      // si r=0 on a un inf .... a discuter
      // if(r2==0)
      //   {
      // 	r2 = 0.0000000001;
      // 	r  = 0.00001;
      //   }

      return druz_sing_normalMvt3D(X, H,  M, a, U, l, eps, param)*X[1]/r;
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