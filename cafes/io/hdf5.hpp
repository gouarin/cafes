// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef IO_HDF5_HPP_INCLUDED
#define IO_HDF5_HPP_INCLUDED

#include <fstream>
#include <hdf5.h>
#include <iostream>
#include <petsc.h>
#include <petscviewerhdf5.h>
// #include <highfive/H5File.hpp>

namespace cafes
{
    namespace io
    {
        void create_xdmf(const char *path, const char *filename,
                         DMDALocalInfo const &info,
                         std::array<double, 2> const &h)
        {
            std::stringstream output;
            std::stringstream output_h5;

            output << path << "/" << filename << ".xdmf";
            output_h5 << filename << ".h5";

            std::ofstream myfile;
            myfile.open(output.str());

            myfile << "<Xdmf>\n";
            myfile << "<Domain>\n";
            myfile << "  <Grid CollectionType=\"Spatial\" "
                      "GridType=\"Collection\" Name=\"Collection\">\n";
            myfile
                << "    <Grid Name=\"pressure_grid\" GridType=\"Uniform\">\n";
            myfile << "      <Topology TopologyType=\"2DCoRectMesh\" "
                      "NumberOfElements=\""
                   << info.mx / 2 + 1 << " " << info.my / 2 + 1 << "\"/>\n";
            myfile << "      <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
            myfile << "        <DataItem DataType=\"Float\" Dimensions= \"2\" "
                      "Format=\"XML\">0.0 0.0 </DataItem>\n";
            myfile << "        <DataItem DataType=\"Float\" Dimensions= \"2\" "
                      "Format=\"XML\">"
                   << 2 * h[0] << " " << 2 * h[1] << "</DataItem>\n";
            myfile << "      </Geometry>\n";
            myfile << "      <Attribute Name=\"pressure\" "
                      "AttributeType=\"Scalar\" Center=\"Node\">\n";
            myfile << "        <DataItem DataType=\"Float\" Precision=\"8\" "
                      "Format=\"HDF\" Dimensions=\""
                   << info.mx / 2 + 1 << " " << info.my / 2 + 1 << "\">\n";
            myfile << "          " << output_h5.str() << ":/pressure\n";
            myfile << "        </DataItem>\n";
            myfile << "      </Attribute>\n";
            myfile << "    </Grid>\n";
            myfile
                << "    <Grid Name=\"velocity_grid\" GridType=\"Uniform\">\n";
            myfile << "      <Topology TopologyType=\"2DCoRectMesh\" "
                      "NumberOfElements=\""
                   << info.mx << " " << info.my << "\"/>\n";
            myfile << "      <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
            myfile << "        <DataItem DataType=\"Float\" Dimensions= \"2\" "
                      "Format=\"XML\">0.0 0.0 </DataItem>\n";
            myfile << "        <DataItem DataType=\"Float\" Dimensions= \"2\" "
                      "Format=\"XML\">"
                   << h[0] << " " << h[1] << "</DataItem>\n";
            myfile << "      </Geometry>\n";
            myfile << "      <Attribute Name=\"velocity\" "
                      "AttributeType=\"Vector\" Center=\"Node\">\n";
            myfile << "        <DataItem DataType=\"Float\" Precision=\"8\" "
                      "Format=\"HDF\" Dimensions=\""
                   << info.mx << " " << info.my << " 2\">\n";
            myfile << "          " << output_h5.str() << ":/velocity\n";
            myfile << "        </DataItem>\n";
            myfile << "      </Attribute>\n";
            myfile << "    </Grid>\n";
            myfile << "  </Grid>\n";
            myfile << "</Domain>\n";
            myfile << "</Xdmf>\n";

            myfile.close();
        }

#undef __FUNCT__
#define __FUNCT__ "save_hdf5"
        template<std::size_t dim>
        PetscErrorCode save_hdf5(const char *path, const char *filename,
                                 Vec sol, DM dm, std::array<double, dim> h)
        {
            PetscErrorCode ierr;
            int rank;
            DM dau, dap;
            DMDALocalInfo infou, infop;
            Vec solu, solp;
            PetscViewer H5viewer;

            std::stringstream output;
            output << path << "/" << filename << ".h5";

            ierr = DMCompositeGetEntries(dm, &dau, &dap);
            CHKERRQ(ierr);
            ierr = DMCompositeGetAccess(dm, sol, &solu, &solp);
            CHKERRQ(ierr);

            ierr = DMDAGetLocalInfo(dau, &infou);
            CHKERRQ(ierr);
            ierr = DMDAGetLocalInfo(dap, &infop);
            CHKERRQ(ierr);

            ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, output.str().data(),
                                       FILE_MODE_WRITE, &H5viewer);
            CHKERRQ(ierr);
            ierr = PetscViewerSetFromOptions(H5viewer);
            CHKERRQ(ierr);

            ierr = PetscObjectSetName((PetscObject)solu, "velocity");
            CHKERRQ(ierr);
            ierr = VecView(solu, H5viewer);
            CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)solp, "pressure");
            CHKERRQ(ierr);
            ierr = VecView(solp, H5viewer);
            CHKERRQ(ierr);

            ierr = DMCompositeRestoreAccess(dm, sol, &solu, &solp);
            CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&H5viewer);
            CHKERRQ(ierr);

            MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
            if (rank == 0)
            {
                create_xdmf(path, filename, infou, h);
            }

            PetscFunctionReturn(0);
        }


#undef __FUNCT__
#define __FUNCT__ "save_particles"
        template<typename Shape, std::size_t dim, typename Surf>
        PetscErrorCode save_particles(const char *path, const char *filename,
                                      std::array<double, dim> h,
                                      const std::vector<particle<Shape>>& particles,
                                      const Surf& surf_points)
        {
            PetscErrorCode ierr;
            std::stringstream output;
            output << path << "/" << filename << ".h5";

            // HighFive::File file(output.str().data(), HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

            // for(std::size_t ipart=0; ipart< particles.size(); ++ipart)
            // {
            //     std::vector<double> surf(surf_points[ipart].size()*dim);
            //     for(std::size_t isurf=0; isurf < surf_points[ipart].size(); ++isurf)
            //     {
            //         auto index = problem::get_index(surf_points[ipart][isurf]);
            //         auto pos = problem::get_position(surf_points[ipart][isurf]);
            //         auto new_pos = index*h + pos;
            //         // std::cout << "[" << new_pos[0] << ", " << new_pos[1] << "],\n";
            //         // surf[isurf] = new_pos[0];
            //         // surf[surf_points[ipart].size() + isurf] = new_pos[1];
            //     }

            //     // HighFive::DataSet dataset = file.createDataSet<double>("/particle_" + ipart,  HighFive::DataSpace{{dim, surf_points[ipart].size()}});
            //     // dataset.write(surf);
            // }

            PetscFunctionReturn(0);
        }
    } // namespace io
} // namespace cafes

#endif