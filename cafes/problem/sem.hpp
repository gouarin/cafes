#ifndef PARTICLE_PROBLEM_SEM_HPP_INCLUDED
#define PARTICLE_PROBLEM_SEM_HPP_INCLUDED

namespace cafes
{
    namespace problem
    {
        template<typename Shape>
        struct SEM{
          std::vector<particle<Shape>> parts_;

          SEM(std::vector<particle<Shape>> parts):
          parts_{parts}
          {}
           
        private:
            std::size_t scale=4;
        };
    }
}
#endif