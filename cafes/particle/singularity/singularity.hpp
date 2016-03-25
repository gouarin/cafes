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

      singularity(particle<Shape> const& p1, particle<Shape> const& p2, double h, double alpha=1, double treshold=1./5)
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
          std::cout << "singularity base " << base_[0] << " " << base_[1] << "\n";

          auto origin_comp = [r1](double x, double y){return x + r1*y;};
          std::transform(p1.center_.begin(), p1.center_.end(), base_[0].begin(), origin_.begin(), origin_comp);

          std::cout << "singularity origin " << origin_ << "\n";
        }
      }

      void construct_base(particle<Shape> const& p1, particle<Shape> const& p2,
                          std::integral_constant<int, 2>)
      {
        base_[0] = position_diff<Shape, Dimensions>(p1, p2)/distance<Shape, Dimensions>(p1, p2);
        base_[1][0] = -base_[0][1];
        base_[1][1] = base_[0][0];

        auto vel_diff = velocity_diff<Shape, Dimensions>(p1, p2);
        for(std::size_t d=0; d<Dimensions; ++d)
          vector_space_[d] = std::inner_product(vel_diff.begin(), vel_diff.end(), base_[d].begin(), 0.);
        UN_ = vector_space_[0];
        UT_ = vector_space_[1];
      }

      void construct_base(particle<Shape> const& p1, particle<Shape> const& p2,
                          std::integral_constant<int, 3>)
      {
        base_[0] = position_diff<Shape, Dimensions>(p1, p2)/distance<Shape, Dimensions>(p1, p2);

        auto vel_diff = velocity_diff<Shape, Dimensions>(p1, p2);
        vector_space_[0] = std::inner_product(vel_diff.begin(), vel_diff.end(), base_[0].begin(), 0.);

        auto vel_norm = std::inner_product(vel_diff.begin(), vel_diff.end(), vel_diff.begin(), 0.);

        if (vector_space_[0] == vel_norm)
        {
          base_[1][0]    = -sign(base_[0][2])*base_[0][2]
                           -sign(base_[0][1])*base_[0][1];
          base_[1][1]    =  sign(base_[0][1])*base_[0][0];
          base_[1][2]    =  sign(base_[0][2])*base_[0][0];
        }
        else
        {
          // Check the formula to have dU - (dU.n)n ?
          for(std::size_t d=0; d<3; ++d)
            base_[1][d] = vel_diff[d]*(1 - base_[0][d]);
        }

        auto norm = std::inner_product(base_[0].begin(), base_[0].end(), base_[0].begin(), 0.);
        base_[1] /= std::sqrt(norm);

        base_[2] = geometry::cross_product(base_[0], base_[1]);

        for(std::size_t d=0; d<2; ++d)
          vector_space_[d] = std::inner_product(vel_diff.begin(), vel_diff.end(), base_[d].begin(), 0.);

        UN_ = vector_space_[0];
        UT_ = vector_space_[1];
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

      /////////////////////////////////////////////////////////////////////////////
      //
      // U_SING
      //
      /////////////////////////////////////////////////////////////////////////////
      auto get_u_sing(position_type const& pos, std::integral_constant<int, 2>)
      {
        auto pos_ref_part = get_pos_in_part_ref(pos);
        std::array< double, 2 > Using;
        std::array< double, 2 > UsingRefPart;

        UsingRefPart[0] = ux_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL)
                        + ux_sing_tangMvt2D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);
        UsingRefPart[1] = uz_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL)
                        + uz_sing_tangMvt2D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);

        for(std::size_t k=0; k<2; ++k){
          Using[0] += base_[k][0]*UsingRefPart[k];
          Using[1] += base_[k][1]*UsingRefPart[k];
        }

        return Using;
      }

      auto get_u_sing(position_type const& pos, std::integral_constant<int, 3>)
      {
        auto pos_ref_part = get_pos_in_part_ref(pos);
        std::array< double, 3 > Using;
        std::array< double, 3 > UsingRefPart;

        UsingRefPart[0] = ux_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL)
                        + ux_sing_tangMvt3D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);
        UsingRefPart[1] = uy_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL)
                        + uy_sing_tangMvt3D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);
        UsingRefPart[2] = uz_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL)
                        + uz_sing_tangMvt3D(pos_ref_part, H1_, H2_, contact_length_, UT_, param_, param_, NULL);

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

      /////////////////////////////////////////////////////////////////////////////
      //
      // GRAD_U_SING
      //
      /////////////////////////////////////////////////////////////////////////////
      auto get_grad_u_sing(position_type pos, std::integral_constant<int, 2>)
      {
        auto pos_ref_part = get_pos_in_part_ref(pos);
        std::array< std::array<double, 2>, 2 > gradUsing;
        std::array< std::array<double, 2>, 2 > gradUsingRefPart;

        gradUsingRefPart[0][0] = dxux_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[0][1] = dzux_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[1][0] = dxuz_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[1][1] = dzuz_sing_normalMvt2D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);

        for(std::size_t k=0; k<2; ++k)
          for(std::size_t l=0; l<2; ++l)
            for(std::size_t m=0; m<2; ++m)
              for(std::size_t n=0; n<2; ++n)
                gradUsing[k][l] += base_[m][k]*base_[n][l]*gradUsingRefPart[m][n];

        return gradUsing;
      }

      auto get_grad_u_sing(position_type& pos, std::integral_constant<int, 3>)
      {
        auto pos_ref_part = get_pos_in_part_ref(pos);
        std::array< std::array<double, 3>, 3 > gradUsing;
        std::array< std::array<double, 3>, 3 > gradUsingRefPart;

        gradUsingRefPart[0][0] = dxux_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[0][1] = dyux_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[0][2] = dzux_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);

        gradUsingRefPart[1][0] = dxuy_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[1][1] = dyuy_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[1][2] = dzuy_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);

        gradUsingRefPart[2][0] = dxuy_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[2][1] = dyuy_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);
        gradUsingRefPart[2][2] = dzuy_sing_normalMvt3D(pos_ref_part, H1_, H2_, contact_length_, UN_, param_, param_, NULL);

        for(std::size_t k=0; k<3; ++k)
          for(std::size_t l=0; l<3; ++l)
            for(std::size_t m=0; m<3; ++m)
              for(std::size_t n=0; n<3; ++n)
                gradUsing[k][l] += base_[m][k]*base_[n][l]*gradUsingRefPart[m][n];

        return gradUsing;
      }

      auto get_grad_u_sing(position_type& pos)
      {
        return get_grad_u_sing(pos, std::integral_constant<int, Dimensions>{});
      }

      /////////////////////////////////////////////////////////////////////////////
      //
      // P_SING
      //
      /////////////////////////////////////////////////////////////////////////////
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