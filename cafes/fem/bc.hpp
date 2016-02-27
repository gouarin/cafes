#ifndef CAFES_FEM_BC_HPP_INCLUDED
#define CAFES_FEM_BC_HPP_INCLUDED

#include <petsc.h>
#include <array>

namespace cafes
{
  namespace fem
  {
    using condition_fn = void(*)(const PetscReal*, PetscScalar*);

    template<std::size_t Dimensions>
    struct dirichlet_conditions
    {
      dirichlet_conditions( std::initializer_list<std::array<condition_fn, Dimensions>> il )
      {
        std::copy(il.begin(), il.end(), conditions_.begin());
      }

      std::array<std::array<condition_fn, Dimensions>, 2*Dimensions> conditions_;

      dirichlet_conditions(dirichlet_conditions const&) = default;
      dirichlet_conditions(dirichlet_conditions&&)      = default;

      dirichlet_conditions& operator=(dirichlet_conditions const&) = default;
      dirichlet_conditions& operator=(dirichlet_conditions&&)      = default;

    };

    void set_vec_with_vec_on_bc(std::array<std::array<PetscInt, 2>, 2> const& index, 
                                PetscScalar ***px, PetscScalar ***py, 
                                std::size_t const& ndof, std::array<condition_fn, 2> const& bc_conditions)
    {
      for (std::size_t j=index[1][0]; j<index[1][1]; ++j)
        for (std::size_t i=index[0][0]; i<index[0][1]; ++i)
          for(std::size_t dof=0; dof<ndof; ++dof)
            if (bc_conditions[dof])
              py[j][i][dof] = px[j][i][dof];
    }

    void set_vec_with_vec_on_bc(std::array<std::array<PetscInt, 2>, 3> const& index, 
                                PetscScalar ****px, PetscScalar ****py, 
                                std::size_t const& ndof, std::array<condition_fn, 3> const& bc_conditions)
    {
      for (std::size_t k=index[2][0]; k<index[2][1]; ++k)
        for (std::size_t j=index[1][0]; j<index[1][1]; ++j)
          for (std::size_t i=index[0][0]; i<index[0][1]; ++i)
            for(std::size_t dof=0; dof<ndof; ++dof)
              if (bc_conditions[dof])
                py[k][j][i][dof] = px[k][j][i][dof];
    }

    void set_vec_with_func_on_bc(std::array<std::array<PetscInt, 2>, 2> const& index, 
                                 PetscScalar ***px, std::array<double, 2> const& h,
                                 std::size_t const& ndof, std::array<condition_fn, 2> const& bc_conditions)
    {
      PetscReal coord[2];
      for (std::size_t j=index[1][0]; j<index[1][1]; ++j){
        coord[1] = j*h[1];
        for (std::size_t i=index[0][0]; i<index[0][1]; ++i){
          coord[0] = i*h[0];
          for(std::size_t dof=0; dof<ndof; ++dof)
            if (bc_conditions[dof])
              bc_conditions[dof](coord, &px[j][i][dof]);
        }
      }
    }

    void set_vec_with_func_on_bc(std::array<std::array<PetscInt, 2>, 3> const& index, 
                                 PetscScalar ****px, std::array<double, 3> const& h,
                                 std::size_t const& ndof, std::array<condition_fn, 3> const& bc_conditions)
    {
      PetscReal coord[3];
      for (std::size_t k=index[2][0]; k<index[2][1]; ++k){
        coord[2] = k*h[2];
        for (std::size_t j=index[1][0]; j<index[1][1]; ++j){
          coord[1] = j*h[1];
          for (std::size_t i=index[0][0]; i<index[0][1]; ++i){
            coord[0] = i*h[0];
            for(std::size_t dof=0; dof<ndof; ++dof)
              if (bc_conditions[dof])
                bc_conditions[dof](coord, &px[k][j][i][dof]);
          }
        }
      }
    }

    #undef __FUNCT__
    #define __FUNCT__ "SetDirichletOnVec"
    PetscErrorCode SetDirichletOnVec(DM dm, dirichlet_conditions<2> bc, Vec x, Vec y){
      PetscErrorCode ierr;
      PetscScalar ***px, ***py;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      // warning: we assume that both u and v on a border have a Dirichlet condition

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);

      ierr = DMDAVecGetArrayDOFRead(dm, x, &px);
      ierr = DMDAVecGetArrayDOF(dm, y, &py);

      using array2d = std::array<std::array<PetscInt, 2>, 2>;

      if (info.bx != DM_BOUNDARY_PERIODIC){
        // left border
        if (info.xs == 0){
          array2d index{{ {{       0,               1 }}, 
                          {{ info.ys, info.ys+info.ym }} 
                       }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[0]);
        }

        // right border
        if (info.xs + info.xm == info.mx){
          array2d index{{ {{ info.mx - 1,         info.mx }}, 
                          {{     info.ys, info.ys+info.ym }} 
                       }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[1]);
        }
      }

      if (info.by != DM_BOUNDARY_PERIODIC){
        // bottom border
        if (info.ys == 0){
          array2d index{{ {{ info.xs, info.xs+info.xm }}, 
                          {{       0,               1 }} 
                        }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[2]);
        }

        // top border
        if (info.ys + info.ym == info.my){
          array2d index{{ {{     info.xs, info.xs+info.xm }}, 
                          {{ info.my - 1,         info.my }} 
                        }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[3]);
        }
      }

      // Restore x and y
      ierr = DMDAVecRestoreArrayDOF(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(dm, y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "SetDirichletOnVec"
    PetscErrorCode SetDirichletOnVec(DM dm, dirichlet_conditions<3> bc, Vec x, Vec y){
      PetscErrorCode ierr;
      PetscScalar ****px, ****py;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);

      ierr = DMDAVecGetArrayDOFRead(dm, x, &px);
      ierr = DMDAVecGetArrayDOF(dm, y, &py);

      using array3d = std::array<std::array<PetscInt, 2>, 3>;

      if (info.bx != DM_BOUNDARY_PERIODIC){
        // left border
        if (info.xs == 0){
          array3d index{{ {{       0,               1 }}, 
                          {{ info.ys, info.ys+info.ym }},
                          {{ info.zs, info.zs+info.zm }}
                        }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[0]);
        }

        // right border
        if (info.xs + info.xm == info.mx){
          array3d index{{ {{ info.mx - 1,         info.mx }}, 
                          {{     info.ys, info.ys+info.ym }},
                          {{     info.zs, info.zs+info.zm }}
                        }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[1]);
        }
      }

      if (info.by != DM_BOUNDARY_PERIODIC){
        // bottom border
        if (info.ys == 0){
          array3d index{{ {{ info.xs, info.xs+info.xm }},
                          {{       0,               1 }}, 
                          {{ info.zs, info.zs+info.zm }}
                        }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[2]);
        }

        // top border
        if (info.ys + info.ym == info.my){
          array3d index{{ {{     info.xs, info.xs+info.xm }},
                          {{ info.my - 1,         info.my }}, 
                          {{     info.zs, info.zs+info.zm }}
                        }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[3]);
        }
      }

      if (info.bz != DM_BOUNDARY_PERIODIC){
        // front border
        if (info.zs == 0){
          array3d index{{ {{ info.xs, info.xs+info.xm }},
                          {{ info.ys, info.ys+info.ym }},
                          {{       0,               1 }}
                        }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[4]);
        }

        // back border
        if (info.zs + info.zm == info.mz){
          array3d index{{ {{     info.xs, info.xs+info.xm }},
                          {{     info.ys, info.ys+info.ym }},
                          {{ info.mz - 1,         info.mz }} 
                        }};
          set_vec_with_vec_on_bc(index, px, py, info.dof, bc.conditions_[5]);
        }
      }

      // Restore x and y
      ierr = DMDAVecRestoreArrayDOF(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(dm, y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "SetDirichletOnRHS"
    PetscErrorCode SetDirichletOnRHS(DM dm, dirichlet_conditions<2> const& bc, Vec x, std::array<double, 2> const& h){
      PetscErrorCode ierr;
      PetscScalar ***px;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      // warning: we assume that both u and v on a border have a Dirichlet condition

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);

      ierr = DMDAVecGetArrayDOF(dm, x, &px);

      using array2d = std::array<std::array<PetscInt, 2>, 2>;

      if (info.bx != DM_BOUNDARY_PERIODIC){
        // left border
        if (info.xs == 0){
          array2d index{{ {{       0,               1 }}, 
                          {{ info.ys, info.ys+info.ym }} 
                        }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[0]);
        }

        // right border
        if (info.xs + info.xm == info.mx){
          array2d index{{ {{ info.mx - 1,         info.mx }}, 
                          {{     info.ys, info.ys+info.ym }} 
                        }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[1]);
        }
      }

      if (info.by != DM_BOUNDARY_PERIODIC){
        // bottom border
        if (info.ys == 0){
          array2d index{{ {{ info.xs, info.xs+info.xm }}, 
                          {{       0,               1 }} 
                        }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[2]);
        }

        // top border
        if (info.ys + info.ym == info.my){
          array2d index{{ {{     info.xs, info.xs+info.xm }}, 
                          {{ info.my - 1,         info.my }} 
                       }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[3]);
        }
      }

      // Restore x
      ierr = DMDAVecRestoreArrayDOF(dm, x, &px);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "SetDirichletOnRHS"
    PetscErrorCode SetDirichletOnRHS(DM dm, dirichlet_conditions<3> const& bc, Vec x, std::array<double, 3> const& h){
      PetscErrorCode ierr;
      PetscScalar ****px;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      // warning: we assume that both u and v on a border have a Dirichlet condition

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);

      ierr = DMDAVecGetArrayDOF(dm, x, &px);

      using array3d = std::array<std::array<PetscInt, 2>, 3>;

      if (info.bx != DM_BOUNDARY_PERIODIC){
        // left border
        if (info.xs == 0){
          array3d index{{ {{       0,               1 }}, 
                          {{ info.ys, info.ys+info.ym }},
                          {{ info.zs, info.zs+info.zm }}
                        }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[0]);
        }

        // right border
        if (info.xs + info.xm == info.mx){
          array3d index{{ {{ info.mx - 1,         info.mx }}, 
                          {{     info.ys, info.ys+info.ym }},
                          {{     info.zs, info.zs+info.zm }}
                        }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[1]);
        }
      }

      if (info.by != DM_BOUNDARY_PERIODIC){
        // bottom border
        if (info.ys == 0){
          array3d index{{ {{ info.xs, info.xs+info.xm }},
                          {{       0,               1 }}, 
                          {{ info.zs, info.zs+info.zm }}
                        }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[2]);
        }

        // top border
        if (info.ys + info.ym == info.my){
          array3d index{{ {{     info.xs, info.xs+info.xm }},
                          {{ info.my - 1,         info.my }}, 
                          {{     info.zs, info.zs+info.zm }}
                        }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[3]);
        }
      }

      if (info.bz != DM_BOUNDARY_PERIODIC){
        // front border
        if (info.zs == 0){
          array3d index{{ {{ info.xs, info.xs+info.xm }},
                          {{ info.ys, info.ys+info.ym }},
                          {{       0,               1 }}
                        }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[4]);
        }

        // back border
        if (info.zs + info.zm == info.mz){
          array3d index{{ {{     info.xs, info.xs+info.xm }},
                          {{     info.ys, info.ys+info.ym }},
                          {{ info.mz - 1,         info.mz }} 
                        }};
          set_vec_with_func_on_bc(index, px, h, info.dof, bc.conditions_[5]);
        }
      }

      // Restore x
      ierr = DMDAVecRestoreArrayDOF(dm, x, &px);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  }
}
#endif