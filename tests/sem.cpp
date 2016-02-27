#include <particle/particle.hpp>
#include <problem/sem.hpp>
#include <particle/geometry/super_ellipsoid.hpp>
#include <particle/geometry/linspace.hpp>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <petsc.h>

void zeros(const PetscReal x[], PetscScalar *u){
  *u = 0.;
}

void haut(const PetscReal x[], PetscScalar *u){
  *u = 0.;
}

void bas(const PetscReal x[], PetscScalar *u){
  *u = 0.;
}

void left(const PetscReal x[], PetscScalar *u){
  *u = 0.;
}

void right(const PetscReal x[], PetscScalar *u){
  *u = 0.;
}

template<std::size_t Dimensions>
struct dirichlet_conditions
{
  using conditon_fn = void(*)(const PetscReal*, PetscScalar*);

  dirichlet_conditions( std::initializer_list<std::array<conditon_fn, Dimensions>> il )
  {
    std::copy(il.begin(), il.end(), conditions_.begin());
  }

  std::array<std::array<conditon_fn, Dimensions>, (1<<Dimensions)> conditions_;
};

template<typename ST>
cafes::particle<ST> make_particle(ST const& se, cafes::physics::force<ST::dimension_type::value> const& a, double d)
{
  return {se,a,d};
}

template<std::size_t N>
cafes::geometry::super_ellipsoid<N> make_ellipsoid(cafes::geometry::position<N,double> const& a, std::array<double, N> const& b, double d)
{
  return {a,b,d};
}

template<typename PL> auto make_SEM(PL const& pt)
{
  using s_t = typename PL::value_type::shape_type;
  return cafes::problem::SEM<s_t>{pt};
}

/*
template<typename P> auto make_SEM(P* const pt)
{
  using s_t = typename P::shape_type;
  return cafes::problem::SEM<s_t>{pt};
}
*/

int main(int argc, char **argv)
{
  {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);


    dirichlet_conditions<2> dc = {  {{haut , zeros}} 
                                 ,  {{bas  , zeros}}
                                 ,  {{zeros, left }}
                                 ,  {{zeros, right}}
                                 };
    auto se = make_ellipsoid<2>( {1.,1.}, {.5,.5}, 1. );
    std::vector<cafes::particle<decltype(se)>> pt{ make_particle(se, {0.,150.}, 0.25)
                                                 , make_particle(se, {0.,150.}, 0.25)
                                                 };

    auto s = make_SEM(pt);

    return 0;
  }
}
