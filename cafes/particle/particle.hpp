#ifndef PARTICLE_PARTICLE_HPP_INCLUDED
#define PARTICLE_PARTICLE_HPP_INCLUDED

// #include <array>
// #include <algorithm>
// #include <initializer_list>
#include <particle/physics/force.hpp>
#include <particle/physics/velocity.hpp>
#include <particle/geometry/box.hpp>
#include <particle/geometry/position.hpp>

namespace cafes
{
    template<typename Shape>
    struct particle : private Shape
    {
        using shape_type     = Shape;
        using dimension_type = typename Shape::dimension_type;
        using velocity_type  = physics::velocity<dimension_type::value>;
        using force_type     = physics::force<dimension_type::value>;

        particle( Shape const& s, force_type const& f, double r )
        : Shape(s), force_(f), rho_(r)
        {}

        using Shape::surface;
        using Shape::contains;

    private:
        force_type force_;
        velocity_type velocity_, angular_velocity_;
        double rho_;
    };

    template<typename Shape>
    auto find_fluid_points_insides( particle<Shape> const& p, cafes::geometry::box<2, int> const& b, std::array<double,2> const& h)
    {
        std::vector<geometry::position<2, int>> that;
        that.reserve(b.length());

        double x, y;
        std::size_t ix, iy;
        for(y=b.bottom_left[1]*h[1], iy=b.bottom_left[1]; iy<b.upper_right[1]; y+=h[1], ++iy)
            for(x=b.bottom_left[0]*h[0], ix=b.bottom_left[0]; ix<b.upper_right[0]; x+=h[0], ++ix)
            {
                cafes::geometry::position<2, double> pt {x,y};
                if(p.contains(pt)) that.push_back({ix, iy});
            }

        return that;
    }

    template<typename Shape>
    auto find_fluid_points_insides( particle<Shape> const& p, cafes::geometry::box<3, int> const& b, std::array<double,3> const& h)
    {
        std::vector<geometry::position<3, int>> that;
        that.reserve(b.length());

        double x, y, z;
        std::size_t ix, iy, iz;
        for(z=b.bottom_left[2]*h[2], iz=b.bottom_left[2]; iz<b.upper_right[2]; z+=h[2], ++iz)
            for(y=b.bottom_left[1]*h[1], iy=b.bottom_left[1]; iy<b.upper_right[1]; y+=h[1], ++iy)
                for(x=b.bottom_left[0]*h[0], ix=b.bottom_left[0]; ix<b.upper_right[0]; x+=h[0], ++ix)
                {
                   cafes::geometry::position<3, double> pt {x, y, z};
                   if(p.contains(pt)) that.push_back({ix, iy, iz});
                }   

        return that;
    }

    template<std::size_t Dimensions>
    auto find_surf_points_insides( std::vector<geometry::position<Dimensions, double>> const& surf_p, cafes::geometry::box<Dimensions, int> const& b, std::array<double,Dimensions> const& h)
    {
        std::vector<std::pair<geometry::position<Dimensions, int>, geometry::position<Dimensions, double>>> that;
        that.reserve(surf_p.size());

        for(std::size_t i=0; i<surf_p.size(); ++i){
            auto surf_pi = static_cast<geometry::position<Dimensions, int>>(surf_p[i]/h);
            if (cafes::geometry::point_inside(b, surf_pi))
                that.push_back(std::make_pair(surf_pi, surf_p[i]));
        }
        return that;
    }
}



#endif
