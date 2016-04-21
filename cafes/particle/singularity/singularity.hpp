#ifndef CAFES_PARTICLE_SINGULARITY_SINGULARITY_HPP_INCLUDED
#define CAFES_PARTICLE_SINGULARITY_SINGULARITY_HPP_INCLUDED

#include <particle/physics/velocity.hpp>
#include <particle/geometry/position.hpp>
#include <particle/geometry/cross_product.hpp>
#include <particle/particle.hpp>

#include "particle/singularity/UandPNormal.hpp"
#include "particle/singularity/UandPTang.hpp"

#include <cmath>
#include <iostream>

namespace cafes
{
  namespace singularity
  {

    int sign(double val) {
      return (0. <= val) - (val < 0.);
    }

    template<std::size_t Dimensions, typename Shape>
    struct singularity
    {
      double alpha_;
      double threshold_;
      double H1_, H2_;
      double K_;
      double scale = 10;

      double contact_length_;
      double cutoff_dist_;
      double param_;
      double UN_, UT_;
      bool is_singularity_;

      using position_type = geometry::position<Dimensions, double>;
      std::array<position_type, Dimensions> base_;
      position_type origin_;
      
      physics::velocity<Dimensions> vector_space_;

      singularity(particle<Shape> const& p1, particle<Shape> const& p2, double h, double alpha=2, double treshold=1./5)
      : alpha_{alpha}, threshold_{treshold}
      {
        double dist = distance<Shape, Dimensions>(p1, p2);

        // Fix this if the particles are not a circle or a sphere
        double r1 = p1.shape_factors_[0];
        double r2 = p2.shape_factors_[0];
        H1_ = 1./r1;
        H2_ = 1./r2;
        contact_length_ = dist - r1 - r2;
        K_ = .5*(1./r1 + 1./r2);

        std::cout << "K = " << K_ << ", H = " << H1_ << ", M = " << H2_ << "\n"; 
        auto minr = (r1 < r2)? r1: r2;

        auto tmp = alpha*std::sqrt(contact_length_/K_);
        cutoff_dist_ = (tmp < minr)? tmp : minr;

        cutoff_dist_ = (cutoff_dist_ <= std::sqrt(2)*h)? std::sqrt(2)*h : cutoff_dist_; 

        param_ = .5*cutoff_dist_*cutoff_dist_;
        //is_singularity_ = contact_length_<threshold_*cutoff_dist_;
        is_singularity_ = contact_length_<threshold_*minr;

        if (is_singularity_)
        {
          construct_base(p1, p2);
          std::cout << "#singularity base ";
          for(std::size_t d=0; d<Dimensions; ++d)
            std::cout << base_[d] << " ";
          std::cout << "\n";

          auto origin_comp = [r1](double x, double y){return x + r1*y;};
          std::transform(p1.center_.begin(), p1.center_.end(), base_[2*Dimensions-4].begin(), origin_.begin(), origin_comp);

          std::cout << "#singularity origin " << origin_ << "\n";
        }
      }

      geometry::box<Dimensions, int> get_box(std::array<double, 2> h)
      {
        double theta1 = std::asin(cutoff_dist_*H1_);
        double theta2 = std::asin(cutoff_dist_*H2_);

        position_type p1_ref_up{ 1./H1_*(-1 + std::cos(theta1)), 1/H1_*std::sin(theta1)};
        position_type p1_ref_down{ 1./H1_*(-1 + std::cos(theta1)), -1/H1_*std::sin(theta1)};

        position_type p2_ref_up{ contact_length_ + 1./H2_*(1 - std::cos(theta2)), 1/H2_*std::sin(theta2)};
        position_type p2_ref_down{ contact_length_ + 1./H2_*(1 - std::cos(theta2)), -1/H2_*std::sin(theta2)};

        position_type p1_up{ base_[0][0]*p1_ref_up[0]+ base_[1][0]*p1_ref_up[1] + origin_[0],
                             base_[0][1]*p1_ref_up[0]+ base_[1][1]*p1_ref_up[1] + origin_[1]};
        position_type p2_up{ base_[0][0]*p2_ref_up[0]+ base_[1][0]*p2_ref_up[1] + origin_[0],
                             base_[0][1]*p2_ref_up[0]+ base_[1][1]*p2_ref_up[1] + origin_[1]};

        position_type p1_down{ base_[0][0]*p1_ref_down[0]+ base_[1][0]*p1_ref_down[1] + origin_[0],
                               base_[0][1]*p1_ref_down[0]+ base_[1][1]*p1_ref_down[1] + origin_[1]};
        position_type p2_down{ base_[0][0]*p2_ref_down[0]+ base_[1][0]*p2_ref_down[1] + origin_[0],
                               base_[0][1]*p2_ref_down[0]+ base_[1][1]*p2_ref_down[1] + origin_[1]};

        return { { std::floor(std::min({p1_up[0], p1_down[0], p2_up[0], p2_down[0]})/h[0]),
                   std::floor(std::min({p1_up[1], p1_down[1], p2_up[1], p2_down[1]})/h[1])},
                 { std::ceil(std::max({p1_up[0], p1_down[0], p2_up[0], p2_down[0]})/h[0]),
                   std::ceil(std::max({p1_up[1], p1_down[1], p2_up[1], p2_down[1]})/h[1])}
               };
      }

      geometry::box<Dimensions, int> get_box(std::array<double, 3> h)
      {
        double theta1 = std::asin(cutoff_dist_*H1_);
        double theta2 = std::asin(cutoff_dist_*H2_);

        position_type p1_ref_up_back   {  1/H1_*std::sin(theta1),  cutoff_dist_ , 1./H1_*(-1 + std::cos(theta1))};
        position_type p1_ref_up_front  {  1/H1_*std::sin(theta1), -cutoff_dist_ , 1./H1_*(-1 + std::cos(theta1))};
        position_type p1_ref_down_back { -1/H1_*std::sin(theta1),  cutoff_dist_ , 1./H1_*(-1 + std::cos(theta1))};
        position_type p1_ref_down_front{ -1/H1_*std::sin(theta1), -cutoff_dist_ , 1./H1_*(-1 + std::cos(theta1))};

        position_type p2_ref_up_back   { 1/H2_*std::sin(theta2),  cutoff_dist_, contact_length_ + 1./H2_*(1 - std::cos(theta2))};
        position_type p2_ref_up_front  { 1/H2_*std::sin(theta2), -cutoff_dist_, contact_length_ + 1./H2_*(1 - std::cos(theta2))};
        position_type p2_ref_down_back {-1/H2_*std::sin(theta2),  cutoff_dist_, contact_length_ + 1./H2_*(1 - std::cos(theta2))};
        position_type p2_ref_down_front{-1/H2_*std::sin(theta2), -cutoff_dist_, contact_length_ + 1./H2_*(1 - std::cos(theta2))};

        position_type p1_up_back{ base_[0][0]*p1_ref_up_back[0] + base_[1][0]*p1_ref_up_back[1] + base_[2][0]*p1_ref_up_back[2] + origin_[0],
                                  base_[0][1]*p1_ref_up_back[0] + base_[1][1]*p1_ref_up_back[1] + base_[2][1]*p1_ref_up_back[2] + origin_[1],
                                  base_[0][2]*p1_ref_up_back[0] + base_[1][2]*p1_ref_up_back[1] + base_[2][2]*p1_ref_up_back[2] + origin_[2]};
        position_type p1_up_front{ base_[0][0]*p1_ref_up_front[0] + base_[1][0]*p1_ref_up_front[1] + base_[2][0]*p1_ref_up_front[2] + origin_[0],
                                   base_[0][1]*p1_ref_up_front[0] + base_[1][1]*p1_ref_up_front[1] + base_[2][1]*p1_ref_up_front[2] + origin_[1],
                                   base_[0][2]*p1_ref_up_front[0] + base_[1][2]*p1_ref_up_front[1] + base_[2][2]*p1_ref_up_front[2] + origin_[2]};

        position_type p1_down_back{ base_[0][0]*p1_ref_down_back[0] + base_[1][0]*p1_ref_down_back[1] + base_[2][0]*p1_ref_down_back[2] + origin_[0],
                                    base_[0][1]*p1_ref_down_back[0] + base_[1][1]*p1_ref_down_back[1] + base_[2][1]*p1_ref_down_back[2] + origin_[1],
                                    base_[0][2]*p1_ref_down_back[0] + base_[1][2]*p1_ref_down_back[1] + base_[2][2]*p1_ref_down_back[2] + origin_[2]};
        position_type p1_down_front{ base_[0][0]*p1_ref_down_front[0] + base_[1][0]*p1_ref_down_front[1] + base_[2][0]*p1_ref_down_front[2] + origin_[0],
                                     base_[0][1]*p1_ref_down_front[0] + base_[1][1]*p1_ref_down_front[1] + base_[2][1]*p1_ref_down_front[2] + origin_[1],
                                     base_[0][2]*p1_ref_down_front[0] + base_[1][2]*p1_ref_down_front[1] + base_[2][2]*p1_ref_down_front[2] + origin_[2]};

        position_type p2_up_back{ base_[0][0]*p2_ref_up_back[0] + base_[1][0]*p2_ref_up_back[1] + base_[2][0]*p2_ref_up_back[2] + origin_[0],
                                  base_[0][1]*p2_ref_up_back[0] + base_[1][1]*p2_ref_up_back[1] + base_[2][1]*p2_ref_up_back[2] + origin_[1],
                                  base_[0][2]*p2_ref_up_back[0] + base_[1][2]*p2_ref_up_back[1] + base_[2][2]*p2_ref_up_back[2] + origin_[2]};
        position_type p2_up_front{ base_[0][0]*p2_ref_up_front[0] + base_[1][0]*p2_ref_up_front[1] + base_[2][0]*p2_ref_up_front[2] + origin_[0],
                                   base_[0][1]*p2_ref_up_front[0] + base_[1][1]*p2_ref_up_front[1] + base_[2][1]*p2_ref_up_front[2] + origin_[1],
                                   base_[0][2]*p2_ref_up_front[0] + base_[1][2]*p2_ref_up_front[1] + base_[2][2]*p2_ref_up_front[2] + origin_[2]};

        position_type p2_down_back{ base_[0][0]*p2_ref_down_back[0] + base_[1][0]*p2_ref_down_back[1] + base_[2][0]*p2_ref_down_back[2] + origin_[0],
                                    base_[0][1]*p2_ref_down_back[0] + base_[1][1]*p2_ref_down_back[1] + base_[2][1]*p2_ref_down_back[2] + origin_[1],
                                    base_[0][2]*p2_ref_down_back[0] + base_[1][2]*p2_ref_down_back[1] + base_[2][2]*p2_ref_down_back[2] + origin_[2]};
        position_type p2_down_front{ base_[0][0]*p2_ref_down_front[0] + base_[1][0]*p2_ref_down_front[1] + base_[2][0]*p2_ref_down_front[2] + origin_[0],
                                     base_[0][1]*p2_ref_down_front[0] + base_[1][1]*p2_ref_down_front[1] + base_[2][1]*p2_ref_down_front[2] + origin_[1],
                                     base_[0][2]*p2_ref_down_front[0] + base_[1][2]*p2_ref_down_front[1] + base_[2][2]*p2_ref_down_front[2] + origin_[2]};

        return { { std::floor(std::min({p1_up_back[0], p1_up_front[0], p1_down_back[0], p1_down_front[0],
                                       p2_up_back[0], p2_up_front[0], p2_down_back[0], p2_down_front[0]})/h[0]),
                   std::floor(std::min({p1_up_back[1], p1_up_front[1], p1_down_back[1], p1_down_front[1],
                                       p2_up_back[1], p2_up_front[1], p2_down_back[1], p2_down_front[1]})/h[1]),
                   std::floor(std::min({p1_up_back[2], p1_up_front[2], p1_down_back[2], p1_down_front[2],
                                       p2_up_back[2], p2_up_front[2], p2_down_back[2], p2_down_front[2]})/h[2])},
                 { std::ceil(std::max({p1_up_back[0], p1_up_front[0], p1_down_back[0], p1_down_front[0],
                                      p2_up_back[0], p2_up_front[0], p2_down_back[0], p2_down_front[0]})/h[0]),
                   std::ceil(std::max({p1_up_back[1], p1_up_front[1], p1_down_back[1], p1_down_front[1],
                                      p2_up_back[1], p2_up_front[1], p2_down_back[1], p2_down_front[1]})/h[1]),
                   std::ceil(std::max({p1_up_back[2], p1_up_front[2], p1_down_back[2], p1_down_front[2],
                                      p2_up_back[2], p2_up_front[2], p2_down_back[2], p2_down_front[2]})/h[2])}};
      }

      void construct_base(particle<Shape> const& p1, particle<Shape> const& p2,
                          std::integral_constant<int, 2>)
      {
        base_[0] = position_diff<Shape, Dimensions>(p2, p1)/distance<Shape, Dimensions>(p1, p2);
        base_[1][0] = -base_[0][1];
        base_[1][1] = base_[0][0];

        auto vel_diff = velocity_diff<Shape, Dimensions>(p2, p1);
        for(std::size_t d=0; d<Dimensions; ++d)
          vector_space_[d] = std::inner_product(vel_diff.begin(), vel_diff.end(), base_[d].begin(), 0.);
        UN_ = vector_space_[0];
        std::cout << "UN_ " << UN_ << "\n";
        UT_ = vector_space_[1];
      }

      void construct_base(particle<Shape> const& p1, particle<Shape> const& p2,
                          std::integral_constant<int, 3>)
      {
        base_[2] = position_diff<Shape, Dimensions>(p2, p1)/distance<Shape, Dimensions>(p1, p2);

        auto vel_diff = velocity_diff<Shape, Dimensions>(p1, p2);
        vector_space_[2] = std::inner_product(vel_diff.begin(), vel_diff.end(), base_[2].begin(), 0.);

        auto vel_norm = std::inner_product(vel_diff.begin(), vel_diff.end(), vel_diff.begin(), 0.);

        if (vector_space_[2]*vector_space_[2] == vel_norm)
        {
          base_[0][0]    = -sign(base_[2][2])*base_[2][2]
                           -sign(base_[2][1])*base_[2][1];
          base_[0][1]    =  sign(base_[2][1])*base_[2][0];
          base_[0][2]    =  sign(base_[2][2])*base_[0][0];
        }
        else
        {
          // Check the formula to have dU - (dU.n)n ?
          for(std::size_t d=0; d<3; ++d)
            base_[0][d] = vel_diff[d]*(1 - base_[2][d]);
        }

        auto norm = std::inner_product(base_[0].begin(), base_[0].end(), base_[0].begin(), 0.);
        base_[0] /= std::sqrt(norm);

        base_[1] = geometry::cross_product(base_[2], base_[0]);

        for(std::size_t d=0; d<2; ++d)
          vector_space_[d] = std::inner_product(vel_diff.begin(), vel_diff.end(), base_[d].begin(), 0.);

        UN_ = vector_space_[2];
        UT_ = vector_space_[0];
        //std::cout << "UN_ " << UN_ << "\n";
      }

      void construct_base(particle<Shape> const& p1, particle<Shape> const& p2)
      {
        construct_base(p1, p2, std::integral_constant<int, Dimensions>{});
      }

      auto get_pos_in_part_ref(position_type const& pos)
      {
        position_type pos_dec, pos_ref_part;

        auto diff = [](double x, double y){return x-y;};
        std::transform(pos.begin(), pos.end(), origin_.begin(), pos_dec.begin(), diff);

        for(std::size_t d=0; d<Dimensions; ++d)
          pos_ref_part[d] = std::inner_product(pos_dec.begin(), pos_dec.end(), base_[d].begin(), 0.);

        return pos_ref_part;
      }

      auto get_pos_from_part_ref(position_type const& pos_ref)
      {
        position_type pos;

        for(std::size_t d1=0; d1<Dimensions; ++d1)
        {
          pos[d1] = 0.;
          for(std::size_t d2=0; d2<Dimensions; ++d2)
            pos[d1] += pos_ref[d2]*base_[d2][d1];
        }

        auto add = [](double x, double y){return x+y;};
        std::transform(pos.begin(), pos.end(), origin_.begin(), pos.begin(), add);

        return pos;
      }

      /////////////////////////////////////////////////////////////////////////////
      //
      // U_SING
      //
      /////////////////////////////////////////////////////////////////////////////
      auto get_u_sing_ref(position_type const& pos_ref_part, std::integral_constant<int, 2>)
      {
        std::array< double, 2 > UsingRefPart{};

        UsingRefPart[0] = ux_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ ux_sing_tangMvt2D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);
        UsingRefPart[1] = uz_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ uz_sing_tangMvt2D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);

        return UsingRefPart;
      }

      auto get_u_sing_ref(position_type const& pos_ref_part, std::integral_constant<int, 3>)
      {
        std::array< double, 3 > UsingRefPart{};

        UsingRefPart[0] = ux_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ ux_sing_tangMvt3D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);
        UsingRefPart[1] = uy_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ uy_sing_tangMvt3D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);
        UsingRefPart[2] = uz_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ uz_sing_tangMvt3D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);

        return UsingRefPart;
      }

      auto get_u_sing_ref(position_type const& pos_ref_part)
      {
        return get_u_sing_ref(pos_ref_part, std::integral_constant<int, Dimensions>{});
      }

      auto get_u_sing(position_type const& pos, std::integral_constant<int, 2>)
      {
        auto pos_ref_part = get_pos_in_part_ref(pos);
        std::array< double, 2 > Using{};
        std::array< double, 2 > UsingRefPart{};

        UsingRefPart[0] = ux_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ ux_sing_tangMvt2D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);
        UsingRefPart[1] = uz_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ uz_sing_tangMvt2D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);

        for(std::size_t k=0; k<2; ++k){
          Using[0] += base_[k][0]*UsingRefPart[k];
          Using[1] += base_[k][1]*UsingRefPart[k];
        }

        return Using;
      }

      auto get_u_sing(position_type const& pos, std::integral_constant<int, 3>)
      {
        auto pos_ref_part = get_pos_in_part_ref(pos);
        std::array< double, 3 > Using{};
        std::array< double, 3 > UsingRefPart{};

        UsingRefPart[0] = ux_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ ux_sing_tangMvt3D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);
        UsingRefPart[1] = uy_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ uy_sing_tangMvt3D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);
        UsingRefPart[2] = uz_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
                        //+ uz_sing_tangMvt3D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);

        for(std::size_t k=0; k<3; ++k){
          Using[0] += base_[k][0]*UsingRefPart[k];
          Using[1] += base_[k][1]*UsingRefPart[k];
          Using[2] += base_[k][2]*UsingRefPart[k];
        }

        return Using;
      }

      auto get_u_sing(position_type const& pos)
      {
        return get_u_sing(pos, std::integral_constant<int, Dimensions>{});
      }

      // /////////////////////////////////////////////////////////////////////////////
      // //
      // // DIVU_SING
      // //
      // /////////////////////////////////////////////////////////////////////////////
      // auto get_divu_sing_ref(position_type& pos, std::integral_constant<int, 2>)
      // {
      //   return DIVCART_normalMvt2D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
      // }

      // // auto get_divu_sing_ref(position_type& pos, std::integral_constant<int, 3>)
      // // {
      // //   return p_sing_withT_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
      // // }

      // auto get_divu_sing_ref(position_type& pos)
      // {
      //   return get_divu_sing_ref(pos, std::integral_constant<int, Dimensions>{});
      // }

      // auto get_divu_sing(position_type& pos, std::integral_constant<int, 2>)
      // {
      //   auto pos_ref_part = get_pos_in_part_ref(pos);

      //   return DIVCART_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
      // }

      // // auto get_divu_sing(position_type& pos, std::integral_constant<int, 3>)
      // // {
      // //   auto pos_ref_part = get_pos_in_part_ref(pos);

      // //   return p_sing_withT_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
      // // }

      // auto get_divu_sing(position_type& pos)
      // {
      //   return get_divu_sing(pos, std::integral_constant<int, Dimensions>{});
      // }


      /////////////////////////////////////////////////////////////////////////////
      //
      // GRAD_U_SING
      //
      /////////////////////////////////////////////////////////////////////////////
      auto get_grad_u_sing_ref(position_type pos, std::integral_constant<int, 2>)
      {
        std::array< std::array<double, 2>, 2 > gradUsingRefPart{};

        gradUsingRefPart[0][0] = dxux_sing_normalMvt2D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[0][1] = dzux_sing_normalMvt2D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[1][0] = dxuz_sing_normalMvt2D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[1][1] = dzuz_sing_normalMvt2D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);

        return gradUsingRefPart;
      }

      auto get_grad_u_sing_ref(position_type& pos, std::integral_constant<int, 3>)
      {
        std::array< std::array<double, 3>, 3 > gradUsingRefPart{};

        gradUsingRefPart[0][0] = dxux_sing_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[0][1] = dyux_sing_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[0][2] = dzux_sing_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);

        gradUsingRefPart[1][0] = dxuy_sing_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[1][1] = dyuy_sing_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[1][2] = dzuy_sing_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);

        gradUsingRefPart[2][0] = dxuz_sing_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[2][1] = dyuz_sing_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[2][2] = dzuz_sing_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);

        return gradUsingRefPart;
      }

      auto get_grad_u_sing(position_type pos)
      {
        auto pos_ref_part = get_pos_in_part_ref(pos);
        std::array< std::array<double, Dimensions>, Dimensions > gradUsing{};

        auto gradUsingRefPart = get_grad_u_sing_ref(pos_ref_part, std::integral_constant<int, Dimensions>{});

        for(std::size_t k=0; k<Dimensions; ++k)
          for(std::size_t l=0; l<Dimensions; ++l)
            for(std::size_t m=0; m<Dimensions; ++m)
              for(std::size_t n=0; n<Dimensions; ++n)
                gradUsing[k][l] += base_[m][k]*base_[n][l]*gradUsingRefPart[m][n];

        return gradUsing;
      }

      auto get_grad_u_sing_ref(position_type& pos)
      {
        return get_grad_u_sing_ref(pos, std::integral_constant<int, Dimensions>{});
      }

      /////////////////////////////////////////////////////////////////////////////
      //
      // P_SING
      //
      /////////////////////////////////////////////////////////////////////////////
      auto get_p_sing_ref(position_type& pos, std::integral_constant<int, 2>)
      {
        return p_sing_withT_normalMvt2D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
      }

      auto get_p_sing_ref(position_type& pos, std::integral_constant<int, 3>)
      {
        return p_sing_withT_normalMvt3D(pos, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
      }

      auto get_p_sing_ref(position_type& pos)
      {
        return get_p_sing_ref(pos, std::integral_constant<int, Dimensions>{});
      }

      auto get_p_sing(position_type& pos, std::integral_constant<int, 2>)
      {
        auto pos_ref_part = get_pos_in_part_ref(pos);

        return p_sing_withT_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
      }

      auto get_p_sing(position_type& pos, std::integral_constant<int, 3>)
      {
        auto pos_ref_part = get_pos_in_part_ref(pos);

        return p_sing_withT_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
      }

      auto get_p_sing(position_type& pos)
      {
        return get_p_sing(pos, std::integral_constant<int, Dimensions>{});
      }

    };

  }
}
#endif
