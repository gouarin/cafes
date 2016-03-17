#ifndef PARTICLE_GEOMETRY_SUPER_ELLIPSOID_HPP_INCLUDED
#define PARTICLE_GEOMETRY_SUPER_ELLIPSOID_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/geometry/box.hpp>
#include <particle/geometry/uniform_sampling.hpp>
#include <particle/geometry/quaternion.hpp>
#include <algorithm>
#include <cassert>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
namespace cafes
{
  namespace geometry
  {

      template<std::size_t Dimensions, typename T=double>
      struct super_ellipsoid
      {
        template<std::size_t N> using int_ = std::integral_constant<std::size_t,N>;
        using dimension_type = int_<Dimensions>;

        using position_type = position<Dimensions, T>;
        using shapes_type   = std::array<double, Dimensions>;

        position_type  center_;
        shapes_type    shape_factors_;

        quaternion q_;
        double perimeter = 0.;

        private:
        double         e_, n_;

        public:

        super_ellipsoid(position_type const& c, shapes_type const& s, double n, quaternion q={})
             : center_{c}, shape_factors_{s}, e_{0}, n_{2./n}, q_{q}
        {
          static_assert(Dimensions == 2, "Constructor mismatch for 3D super ellipsoid");
        }

        super_ellipsoid(position_type const& c, shapes_type const& s, double n, double e, quaternion q={})
             : center_{c}, shape_factors_{s}, e_{2./e}, n_{2./n}, q_{q}
        {
          static_assert(Dimensions == 3, "Constructor mismatch for 2D super ellipsoid");
        }

        super_ellipsoid(super_ellipsoid const&) = default;
        super_ellipsoid(super_ellipsoid&&)      = default;

        super_ellipsoid& operator=(super_ellipsoid const&) = default;
        super_ellipsoid& operator=(super_ellipsoid&&)      = default;

        box<Dimensions, double> bounding_box() const
        {
          position_type bl{center_}, ur{center_};

          for(std::size_t i=0;i<Dimensions;++i)
          {
            bl[i] -= shape_factors_[i];
            ur[i] += shape_factors_[i];
          }

          box<Dimensions, double> that{bl, ur};

          if (q_.is_rotate())
          {

            std::vector<position_type> p_list = that.get_box_points();

            for(auto& p: p_list){
              p -= center_;
              p = q_.rotate(p);
              p += center_;
            }

            bl = p_list[0];
            ur = p_list[0];
            for(std::size_t i=1; i<p_list.size(); ++i){
              for (std::size_t d=0; d<Dimensions; ++d){
                if (bl[d] > p_list[i][d]) bl[d] = p_list[i][d];
                if (ur[d] < p_list[i][d]) ur[d] = p_list[i][d];
              }
            }
            that = {bl, ur};
          }

          return that;
        }

        box<Dimensions, int> bounding_box(std::array<double, Dimensions> const& h) const
        {
            box<Dimensions, double> b = bounding_box();
            return {(b.bottom_left/h - 1.), b.upper_right/h + 1.};
        }

        std::vector<position_type> surface(std::vector<double> const& u_samples) const
        {
            std::vector<position_type> that(u_samples.size());

            std::transform( u_samples.cbegin(),u_samples.cend()
                          , that.begin()
                          , [&](double u) { return at(u); }
                          );

            return that;
        }

        double surface_area() const
        {
          return 1.;
          // if (perimeter == 0){
          //   auto pts = surface(.01);

          //   for(std::size_t i=0; i<pts.size()-1; ++i){
          //     auto tmp = pts[i]-pts[i+1];
          //     perimeter += std::sqrt(std::inner_product(tmp.begin(), tmp.end(), tmp.begin(), 0.));
          //   }
          //   auto tmp = pts[0]-pts[pts.size()-1];
          //   perimeter += std::sqrt(std::inner_product(tmp.begin(), tmp.end(), tmp.begin(), 0.));
          // }
          // return perimeter;
        }
        
        std::vector<position_type>
        surface(double const& k, double tol=1e-2)
        {
          std::vector<position_type> that;
          auto radius = shape_factors_[0];
          int n = 6*shape_factors_[0]/k;

          for(std::size_t i=0; i<n; ++i){
            auto c = radius*std::cos(2*i*M_PI/n);
            auto s = radius*std::sin(2*i*M_PI/n);
            that.push_back({center_[0] + c, center_[1] + s});
          }
          return that;
        }

        // std::vector<position_type>
        // surface(double const& k, double tol=1e-2) const
        // {
        //   auto that = uniform_sampling(shape_factors_, n_, k, tol);
        //   if (q_.is_rotate()){
        //     std::for_each(that.begin(), that.end(),[&](auto& p){ p = q_.rotate(p);});
        //   }
        //   std::for_each(that.begin(), that.end(),[&](auto& p){p += center_;});
        //   return that;
        // }

        std::vector<position_type> radial_vector(std::vector<double> const& u_samples) const
        {
            std::vector<position_type> that(u_samples.size());

            std::transform( u_samples.cbegin(),u_samples.cend()
                          , that.begin()
                          , [&](double u) { return radial(u); }
                          );

            return that;
        }

        std::vector<position_type>
        surface( std::vector<double> const& u_samples
               , std::vector<double> const& v_samples
               ) const
        {
            assert(  u_samples.size() == v_samples.size()
                  && "Incompatible sample size"
                  );

            std::vector<position_type> that(u_samples.size());

            std::transform( u_samples.cbegin(),u_samples.cend()
                          , v_samples.cbegin()
                          , that.begin()
                          , [&](double u, double v) { return at(u,v); }
                          );

            return that;
        }

        std::vector<position_type>
        surface(std::array<double, 2> const& k, double tol=1e-2)
        {
            std::vector<position_type> that;
            double omega1 = 0., omega2 = 0.;
            double eps1 = n_;
            double eps2 = e_;

            while (omega1 < M_PI/4){
              auto z1 = shape_factors_[2]*std::pow(std::sin(omega1), eps2);
              auto phi1 = std::pow(std::cos(omega1), eps2);

              auto z2 = shape_factors_[2]*std::pow(std::cos(omega2), eps2);
              auto phi2 = std::pow(std::sin(omega2), eps2);

              auto update_omega = theta(1, 1, omega1, eps2, k[1], tol);
              omega1 += update_omega;
              update_omega = theta(1, 1, omega2, eps2, k[1], tol);
              omega2 += update_omega;

              auto temp1 = uniform_sampling(shape_factors_, eps1, k[0], z1, phi1, tol);
      
              auto temp2 = uniform_sampling(shape_factors_, eps1, k[0], z2, phi2, tol);
              that.insert(that.end(), temp1.cbegin(), temp1.cend());
              that.insert(that.end(), temp2.cbegin(), temp2.cend());

            }
            std::for_each(that.begin(), that.end(),[&](auto& p){that.push_back({p[0], p[1], -p[2]});});
            std::for_each(that.begin(), that.end(),[&](auto& p){p[0] += center_[0];p[1] += center_[1];p[2] += center_[2];});
            return that;
        }

        std::vector<position_type>
        radial_vector( std::vector<double> const& u_samples
               , std::vector<double> const& v_samples
               ) const
        {
            assert(  u_samples.size() == v_samples.size()
                  && "Incompatible sample size"
                  );

            std::vector<position_type> that(u_samples.size());

            std::transform( u_samples.cbegin(),u_samples.cend()
                          , v_samples.cbegin()
                          , that.begin()
                          , [&](double u, double v) { return radial(u,v); }
                          );

            return that;
        }

        position_type at(double u) const
        {
          static_assert(Dimensions == 2, "Incompatible dimension for parametric representation");
          return  { shape_factors_[0] * c(u,n_) + center_[0]
                  , shape_factors_[1] * s(u,n_) + center_[1]
                  };
        }

        position_type radial(double u) const
        {
          static_assert(Dimensions == 2, "Incompatible dimension for parametric representation");
          return  { shape_factors_[0] * c(u,n_) 
                  , shape_factors_[1] * s(u,n_) 
                  };
        }

        position_type at(double u, double v) const
        {
          static_assert(Dimensions == 3, "Incompatible dimension for parametric representation");

          return  { shape_factors_[0] * c(v,n_) * c(u,e_) + center_[0]
                  , shape_factors_[1] * c(v,n_) * s(u,e_) + center_[1]
                  , shape_factors_[2] * s(v,n_) + center_[2]
                  };
        }

        position_type radial(double u, double v) const
        {
          static_assert(Dimensions == 3, "Incompatible dimension for parametric representation");

          return  { shape_factors_[0] * c(v,n_) * c(u,e_)
                  , shape_factors_[1] * c(v,n_) * s(u,e_)
                  , shape_factors_[2] * s(v,n_) 
                  };
        }

        double implicit(position_type const& p) const
        {
           return implicit(p, dimension_type{});
        }

        bool contains(position_type const& p) const
        {
           return implicit(p) <= 1.;
        }

        private:

        double implicit(position_type const& p, int_<2> const&) const
        {
           auto r = 2./n_;
           position_type pos = p-center_;
           if (q_.is_rotate())
            pos = q_.conj().rotate(pos);
           return std::pow( std::abs(pos[0]/shape_factors_[0]), r)
                + std::pow( std::abs(pos[1]/shape_factors_[1]), r);
        }

        double implicit(position_type const& p, int_<3> const&) const
        {
           auto r1 = 2./n_;
           auto r2 = e_/n_;
           auto r3 = 2./e_;

           return std::pow( std::pow( std::abs((p[0]-center_[0])/shape_factors_[0]), r3)
                          + std::pow( std::abs((p[1]-center_[1])/shape_factors_[1]), r3)
                          , r2
                          )

                + std::pow( std::abs((p[2]-center_[2])/shape_factors_[2]), r1);
        }

        static double c(double w, double m)
        {
          auto cw = std::cos(w);
          return std::copysign( std::pow(std::abs(cw),m), cw);
        }

        static double s(double w, double m)
        {
          auto sw = std::sin(w);
          return std::copysign( std::pow(std::abs(sw),m), sw);
        }
      };
  }
  template<std::size_t N>
  geometry::super_ellipsoid<N> make_ellipsoid(geometry::position<N,double> const& a, std::array<double, N> const& b, double d, geometry::quaternion q)
  {
    return {a,b,d, q};
  }

  template<std::size_t N>
  geometry::super_ellipsoid<N> make_ellipsoid(geometry::position<N,double> const& a, std::array<double, N> const& b, double d, double e, geometry::quaternion q)
  {
    return {a, b, d, e, q};
  }

}

#endif
