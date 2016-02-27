#ifndef PARTICLE_GEOMETRY_SUPER_ELLIPSOID_HPP_INCLUDED
#define PARTICLE_GEOMETRY_SUPER_ELLIPSOID_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/geometry/box.hpp>
#include <algorithm>
#include <cassert>
#include <vector>
#include <array>
#include <cmath>

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

        private:
        position_type  center_;
        shapes_type    shape_factors_;
        double         e_, n_;

        public:

        super_ellipsoid(position_type const& c, shapes_type const& s, double n)
             : center_{c}, shape_factors_{s}, e_{0}, n_{n}
        {
          static_assert(Dimensions == 2, "Constructor mismatch for 3D super ellipsoid");
        }

        super_ellipsoid(position_type const& c, shapes_type const& s, double n, double e)
             : center_{c}, shape_factors_{s}, e_{e}, n_{n}
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

            return {bl,ur};
        }

        box<Dimensions, int> bounding_box(std::array<double, Dimensions> const& h) const
        {
            box<Dimensions, double> b = bounding_box();
            return {b.bottom_left/h, b.upper_right/h};
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

           return std::pow( std::abs((p[0]-center_[0])/shape_factors_[0]), r)
                + std::pow( std::abs((p[1]-center_[1])/shape_factors_[1]), r);
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
}

#endif
