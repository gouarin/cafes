#ifndef IO_VTK_HPP_INCLUDED
#define IO_VTK_HPP_INCLUDED

#include <iostream>
#include <sstream>

#include <petsc.h>

#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkCell.h"
#include "vtkPointData.h"
#include "vtkStructuredGrid.h"
#include "vtkImageData.h"
#include "vtkXMLStructuredGridWriter.h"
#include "vtkXMLImageDataWriter.h"
#include "vtkXMLPStructuredGridWriter.h"
#include "vtkStructuredGridOutlineFilter.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkAxes.h"
#include "vtkContourFilter.h"
#include "vtkArrowSource.h"
#include "vtkGlyph3D.h"
#include "vtkScalarBarActor.h"
#include "vtkTubeFilter.h"
#include "vtkProperty.h"
#include "vtkHedgeHog.h"
#include "vtkShrinkFilter.h"
#include "vtkAssignAttribute.h"
#include "vtkMaskPoints.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPolyData.h"

#include <array>

namespace cafes
{
  namespace io
  {
    #undef __FUNCT__
    #define __FUNCT__ "save_VTK"
    PetscErrorCode save_VTK(const char* path, const char* filename, Vec sol, DM dm, std::array<double, 2> h)
    {
      PetscErrorCode ierr;
      DM             dau, dap;
      DMDALocalInfo  infou, infop;
      Vec            solu, solp;
      Vec            locsolu, locsolp;
      PetscScalar    ***psolu, **psolp;

      ierr = DMCompositeGetEntries(dm, &dau, &dap);CHKERRQ(ierr);

      ierr = DMDAGetLocalInfo(dau, &infou);CHKERRQ(ierr);
      ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);

      ierr = DMCompositeGetAccess(dm,sol,&solu,&solp);CHKERRQ(ierr);

      ierr = DMGetLocalVector(dau, &locsolu);CHKERRQ(ierr);
      ierr = DMGetLocalVector(dap, &locsolp);CHKERRQ(ierr);

      ierr = DMGlobalToLocalBegin(dau, solu, INSERT_VALUES, locsolu);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dau, solu, INSERT_VALUES, locsolu);CHKERRQ(ierr);

      ierr = DMGlobalToLocalBegin(dap, solp, INSERT_VALUES, locsolp);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dap, solp, INSERT_VALUES, locsolp);CHKERRQ(ierr);

      ierr = DMDAVecGetArrayDOF(dau, locsolu, &psolu);CHKERRQ(ierr);
      ierr = DMDAVecGetArray(dap, locsolp, &psolp);CHKERRQ(ierr);
      
      vtkStructuredGrid* pressureDataSet = vtkStructuredGrid::New();
      pressureDataSet->SetExtent(infop.gxs, infop.gxs+infop.gxm-1, infop.gys, infop.gys+infop.gym-1, 0, 0);

      vtkStructuredGrid* velocityDataSet = vtkStructuredGrid::New();
      velocityDataSet->SetExtent(infou.gxs, infou.gxs+infou.gxm-1, infou.gys, infou.gys+infou.gym-1, 0, 0);
      
      vtkPoints* pressurePoints = vtkPoints::New();
      vtkPoints* velocityPoints = vtkPoints::New();
      
      vtkDoubleArray* velocity = vtkDoubleArray::New();
      velocity->SetNumberOfComponents(3);
      velocity->SetName("velocity");
      
      vtkDoubleArray* pressure = vtkDoubleArray::New();
      pressure->SetName("pressure");
      
      for (int j=infop.gys; j<infop.gys+infop.gym; j++)
        for (int i=infop.gxs; i<infop.gxs+infop.gxm; i++)
          {
            pressurePoints->InsertNextPoint(2*i*h[0], 2*j*h[1], 0.);
            pressure->InsertNextValue(psolp[j][i]);
          }

      int rank;
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      for (int j=infou.gys; j<infou.gys+infou.gym; j++)
        for (int i=infou.gxs; i<infou.gxs+infou.gxm; i++)
          {
            velocityPoints->InsertNextPoint(i*h[0], j*h[1], 0.);
            velocity->InsertNextTuple3(psolu[j][i][0], psolu[j][i][1], 0.);
          }

      pressureDataSet->SetPoints(pressurePoints);
      pressureDataSet->GetPointData()->SetScalars(pressure);

      velocityDataSet->SetPoints(velocityPoints);
      velocityDataSet->GetPointData()->SetScalars(velocity);
      
      vtkXMLStructuredGridWriter* pressureDataWriter = vtkXMLStructuredGridWriter::New();
      vtkXMLStructuredGridWriter* velocityDataWriter = vtkXMLStructuredGridWriter::New();
      std::stringstream op, ov;
      int size;//, rank;
      MPI_Comm_size(PETSC_COMM_WORLD, &size);
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      int x[size], y[size], z[size];
      int sx[size], sy[size], sz[size];

      MPI_Gather(&infop.gxs, 1, MPI_INT, x, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infop.gxm, 1, MPI_INT, sx, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      MPI_Gather(&infop.gys, 1, MPI_INT, y, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infop.gym, 1, MPI_INT, sy, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      MPI_Gather(&infop.gzs, 1, MPI_INT, z, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infop.gzm, 1, MPI_INT, sz, 1, MPI_INT, 0, PETSC_COMM_WORLD);


      op << path << "/" << filename << "_pressure_" << rank << ".vts";
      pressureDataWriter->SetFileName(op.str().data());
      //dataWriter->SetDataModeToAscii();
    #if VTK_MAJOR_VERSION <= 5
      pressureDataWriter->SetInput(pressureDataSet);
    #else
      pressureDataWriter->SetInputData(pressureDataSet);
    #endif
      pressureDataWriter->Write();

      ov << path << "/" << filename << "_velocity_" << rank << ".vts";
      velocityDataWriter->SetFileName(ov.str().data());
      //dataWriter->SetDataModeToAscii();
    #if VTK_MAJOR_VERSION <= 5
      velocityDataWriter->SetInput(velocityDataSet);
    #else
      velocityDataWriter->SetInputData(velocityDataSet);
    #endif
      velocityDataWriter->Write();

      if (rank == 0){
        std::stringstream oall;
        oall << path << "/" << filename << "_pressure.pvts";
        ofstream pvts;
        pvts.open(oall.str().data(), ios::out | ios::trunc);
        pvts << "<?xml version=\"1.0\"?>" << endl;
        pvts << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
        pvts << "<PStructuredGrid WholeExtent=\"0 " << infop.mx-1 << " 0 " << infop.my-1 << " 0 0\" GhostLevel=\"0\">";
        pvts << "<PPoints>" << endl;
        pvts << "<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << endl;
        pvts << "</PPoints>" << endl;
        pvts << "<PPointData Scalars=\"pressure\">" << endl;
        pvts << "<PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
        pvts << "</PPointData>" << endl;
        pvts << "<PCells>" << endl;
        pvts << "</PCells>" << endl;
        // Attention!! ajout d'un point apres source pour aller dans le repertoire ..
        // changer les arguments d'entree en dir et filename
        for(int i=0; i<size; i++){
          pvts << "<Piece Extent=\"" << x[i] << " " << x[i]+sx[i]-1 ;
          pvts << " " << y[i] << " " << y[i]+sy[i]-1 << " ";
          pvts << " 0 0\" Source=\"./" << filename << "_pressure_" << i << ".vts\"/>" << endl;
        }
            pvts << "</PStructuredGrid>" << endl;
        pvts << "</VTKFile>" << endl;
        pvts.close();
      }

      MPI_Gather(&infou.gxs, 1, MPI_INT, x, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infou.gxm, 1, MPI_INT, sx, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      MPI_Gather(&infou.gys, 1, MPI_INT, y, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infou.gym, 1, MPI_INT, sy, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      MPI_Gather(&infou.gzs, 1, MPI_INT, z, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infou.gzm, 1, MPI_INT, sz, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      if (rank == 0){
        std::stringstream oall;
        oall << path << "/" << filename << "_velocity.pvts";
        ofstream pvts;
        pvts.open(oall.str().data(), ios::out | ios::trunc);
        pvts << "<?xml version=\"1.0\"?>" << endl;
        pvts << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
        pvts << "<PStructuredGrid WholeExtent=\"0 " << infou.mx-1 << " 0 " << infou.my-1 << " 0 0\" GhostLevel=\"0\">";
        pvts << "<PPoints>" << endl;
        pvts << "<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << endl;
        pvts << "</PPoints>" << endl;
        pvts << "<PPointData Vectors=\"velocity\">" << endl;
        pvts << "<PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
        pvts << "</PPointData>" << endl;
        pvts << "<PCells>" << endl;
        pvts << "</PCells>" << endl;
        // Attention!! ajout d'un point apres source pour aller dans le repertoire ..
        // changer les arguments d'entree en dir et filename
        for(int i=0; i<size; i++){
          pvts << "<Piece Extent=\"" << x[i] << " " << x[i]+sx[i]-1;
          pvts << " " << y[i] << " " << y[i]+sy[i]-1 << " ";
          pvts << " 0 0\" Source=\"./" << filename << "_velocity_" << i << ".vts\"/>" << endl;
        }
        pvts << "</PStructuredGrid>" << endl;
        pvts << "</VTKFile>" << endl;
        pvts.close();
      }

      pressurePoints->Delete();
      velocityPoints->Delete();
      velocity->Delete();
      pressure->Delete();
      pressureDataSet->Delete();
      velocityDataSet->Delete();
      pressureDataWriter->Delete();
      velocityDataWriter->Delete();

      ierr = DMDAVecRestoreArray(dau, locsolu, &psolu);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(dap, locsolp, &psolp);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(dau, &locsolu);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(dap, &locsolp);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(dm,sol,&solu,&solp);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

// PetscErrorCode saveInVTK3dVTI(const char* path, const char* filename, Vec sol, StokesCtx *user){
//   PetscErrorCode ierr;
//   DM             dau, dap;
//   DMDALocalInfo  infou, infop;
//   Vec            solu, solp;
//   Vec            locsolu, locsolp;
//   PetscScalar    ****psolu, ***psolp;

//   ierr = DMCompositeGetEntries(user->dm, &dau, &dap);CHKERRQ(ierr);
//   ierr = DMDAGetLocalInfo(dau, &infou);CHKERRQ(ierr);
//   ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);

//   ierr = DMCompositeGetAccess(user->dm,sol,&solu,&solp);CHKERRQ(ierr);

//   ierr = DMGetLocalVector(dau, &locsolu);CHKERRQ(ierr);
//   ierr = DMGetLocalVector(dap, &locsolp);CHKERRQ(ierr);

//   ierr = DMGlobalToLocalBegin(dau, solu, INSERT_VALUES, locsolu);CHKERRQ(ierr);
//   ierr = DMGlobalToLocalEnd(dau, solu, INSERT_VALUES, locsolu);CHKERRQ(ierr);

//   ierr = DMGlobalToLocalBegin(dap, solp, INSERT_VALUES, locsolp);CHKERRQ(ierr);
//   ierr = DMGlobalToLocalEnd(dap, solp, INSERT_VALUES, locsolp);CHKERRQ(ierr);

//   ierr = DMDAVecGetArrayDOF(dau, locsolu, &psolu);CHKERRQ(ierr);
//   ierr = DMDAVecGetArray(dap, locsolp, &psolp);CHKERRQ(ierr);
  
//   vtkImageData* pressureDataSet = vtkImageData::New();
//   pressureDataSet->SetExtent(infop.gxs, infop.gxs+infop.gxm-1, infop.gys, infop.gys+infop.gym-1, infop.gzs, infop.gzs+infop.gzm-1);
//   pressureDataSet->SetOrigin(0., 0., 0.);
//   pressureDataSet->SetSpacing(2*user->hx, 2*user->hy, 2*user->hz);
//   pressureDataSet->SetExtent(infop.gxs, infop.gxs+infop.gxm-1, infop.gys, infop.gys+infop.gym-1, infop.gzs, infop.gzs+infop.gzm-1);

//   vtkImageData* velocityDataSet = vtkImageData::New();
//   velocityDataSet->SetExtent(infou.gxs, infou.gxs+infou.gxm-1, infou.gys, infou.gys+infou.gym-1, infou.gzs, infou.gzs+infou.gzm-1);
//   velocityDataSet->SetOrigin(0., 0., 0.);
//   velocityDataSet->SetSpacing(user->hx, user->hy, user->hz);
//   velocityDataSet->SetExtent(infou.gxs, infou.gxs+infou.gxm-1, infou.gys, infou.gys+infou.gym-1, infou.gzs, infou.gzs+infou.gzm-1);
    
//   vtkDoubleArray* velocity = vtkDoubleArray::New();
//   velocity->SetNumberOfComponents(3);
//   velocity->SetName("velocity");
  
//   vtkDoubleArray* pressure = vtkDoubleArray::New();
//   pressure->SetName("pressure");
  
//   for (int k=infop.gzs; k<infop.gzs+infop.gzm; k++)
//     for (int j=infop.gys; j<infop.gys+infop.gym; j++)
//       for (int i=infop.gxs; i<infop.gxs+infop.gxm; i++)
//   pressure->InsertNextValue(psolp[k][j][i]);

//   int rank;
//   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//   for (int k=infou.gzs; k<infou.gzs+infou.gzm; k++)
//     for (int j=infou.gys; j<infou.gys+infou.gym; j++)
//       for (int i=infou.gxs; i<infou.gxs+infou.gxm; i++){
//   velocity->InsertNextTuple3(psolu[k][j][i][0], psolu[k][j][i][1], psolu[k][j][i][2]);
//       }
//   pressureDataSet->GetPointData()->SetScalars(pressure);

//   velocityDataSet->GetPointData()->SetVectors(velocity);
  
//   vtkXMLImageDataWriter* pressureDataWriter = vtkXMLImageDataWriter::New();
//   vtkXMLImageDataWriter* velocityDataWriter = vtkXMLImageDataWriter::New();
//   stringstream op, ov;
//   //int rank;
//   int size;
//   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//   MPI_Comm_size(PETSC_COMM_WORLD, &size);

//   op << path << "/" << filename << "_pressure_" << rank << ".vti";
//   ov << path << "/" << filename << "_velocity_" << rank << ".vti";

//   pressureDataWriter->SetFileName(op.str().data());

// #if VTK_MAJOR_VERSION <= 5
//   pressureDataWriter->SetInput(pressureDataSet);
// #else
//   pressureDataWriter->SetInputData(pressureDataSet);
// #endif

//   pressureDataWriter->Write();

//   velocityDataWriter->SetFileName(ov.str().data());

// #if VTK_MAJOR_VERSION <= 5
//   velocityDataWriter->SetInput(velocityDataSet);
// #else
//   velocityDataWriter->SetInputData(velocityDataSet);
// #endif

//   velocityDataWriter->Write();

//   int x[size], y[size], z[size];
//   int sx[size], sy[size], sz[size];
//   MPI_Gather(&infop.gxs, 1, MPI_INT, x, 1, MPI_INT, 0, PETSC_COMM_WORLD);
//   MPI_Gather(&infop.gxm, 1, MPI_INT, sx, 1, MPI_INT, 0, PETSC_COMM_WORLD);

//   MPI_Gather(&infop.gys, 1, MPI_INT, y, 1, MPI_INT, 0, PETSC_COMM_WORLD);
//   MPI_Gather(&infop.gym, 1, MPI_INT, sy, 1, MPI_INT, 0, PETSC_COMM_WORLD);

//   MPI_Gather(&infop.gzs, 1, MPI_INT, z, 1, MPI_INT, 0, PETSC_COMM_WORLD);
//   MPI_Gather(&infop.gzm, 1, MPI_INT, sz, 1, MPI_INT, 0, PETSC_COMM_WORLD);

//   if (rank == 0){
//     stringstream oall;
//     oall << path << "/" << filename << "_pressure.pvti";
//     ofstream pvts;
//     pvts.open(oall.str().data(), ios::out | ios::trunc);
//     pvts << "<?xml version=\"1.0\"?>" << endl;
//     pvts << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
//     pvts << "<PImageData WholeExtent=\"0 " << infop.mx-1 << " 0 " << infop.my-1 << " 0 " << infop.mz-1 << "\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"" << 2*user->hx << "," << 2*user->hy << "," << 2*user->hz << "\">";
//     pvts << "<PPoints>" << endl;
//     pvts << "<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << endl;
//     pvts << "</PPoints>" << endl;
//     pvts << "<PPointData Scalars=\"pressure\">" << endl;
//     pvts << "<PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
//     pvts << "</PPointData>" << endl;
//     pvts << "<PCells>" << endl;
//     pvts << "</PCells>" << endl;
//     for(int i=0; i<size; i++)
//       pvts << "<Piece Extent=\"" << x[i] << " " << x[i]+sx[i]-1 << " " << y[i] << " " << y[i]+sy[i]-1 << " " << z[i] << " " << z[i]+sz[i]-1 << "\" Source=\"" << filename << "_pressure_" << i << ".vti\"/>" << endl;
//         pvts << "</PImageData>" << endl;
//     pvts << "</VTKFile>" << endl;
//     pvts.close();
//   }

//   MPI_Gather(&infou.gxs, 1, MPI_INT, x, 1, MPI_INT, 0, PETSC_COMM_WORLD);
//   MPI_Gather(&infou.gxm, 1, MPI_INT, sx, 1, MPI_INT, 0, PETSC_COMM_WORLD);

//   MPI_Gather(&infou.gys, 1, MPI_INT, y, 1, MPI_INT, 0, PETSC_COMM_WORLD);
//   MPI_Gather(&infou.gym, 1, MPI_INT, sy, 1, MPI_INT, 0, PETSC_COMM_WORLD);

//   MPI_Gather(&infou.gzs, 1, MPI_INT, z, 1, MPI_INT, 0, PETSC_COMM_WORLD);
//   MPI_Gather(&infou.gzm, 1, MPI_INT, sz, 1, MPI_INT, 0, PETSC_COMM_WORLD);

//   if (rank == 0){
//     stringstream oall;
//     oall << path << "/" << filename << "_velocity.pvti";
//     ofstream pvts;
//     pvts.open(oall.str().data(), ios::out | ios::trunc);
//     pvts << "<?xml version=\"1.0\"?>" << endl;
//     pvts << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
//     pvts << "<PImageData WholeExtent=\"0 " << infou.mx-1 << " 0 " << infou.my-1 << " 0 " << infou.mz-1 << "\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"" << user->hx << "," << user->hy << "," << user->hz << "\">";
//     pvts << "<PPoints>" << endl;
//     pvts << "<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << endl;
//     pvts << "</PPoints>" << endl;
//     pvts << "<PPointData Vectors=\"velocity\">" << endl;
//     pvts << "<PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
//     pvts << "</PPointData>" << endl;
//     pvts << "<PCells>" << endl;
//     pvts << "</PCells>" << endl;
//     for(int i=0; i<size; i++)
//       pvts << "<Piece Extent=\"" << x[i] << " " << x[i]+sx[i]-1 << " " << y[i] << " " << y[i]+sy[i]-1 << " " << z[i] << " " << z[i]+sz[i]-1 << "\" Source=\"" << filename << "_velocity_" << i << ".vti\"/>" << endl;
//         pvts << "</PImageData>" << endl;
//     pvts << "</VTKFile>" << endl;
//     pvts.close();
//   }

//   velocity->Delete();
//   pressure->Delete();
//   pressureDataSet->Delete();
//   pressureDataWriter->Delete();
//   velocityDataSet->Delete();
//   velocityDataWriter->Delete();

//   ierr = DMDAVecRestoreArray(dau, locsolu, &psolu);CHKERRQ(ierr);
//   ierr = DMDAVecRestoreArray(dap, locsolp, &psolp);CHKERRQ(ierr);
//   ierr = DMRestoreLocalVector(dau, &locsolu);CHKERRQ(ierr);
//   ierr = DMRestoreLocalVector(dap, &locsolp);CHKERRQ(ierr);
//   ierr = DMCompositeRestoreAccess(user->dm,sol,&solu,&solp);CHKERRQ(ierr);
//   PetscFunctionReturn(0);
// }

  #undef __FUNCT__
  #define __FUNCT__ "save_VTK"
  PetscErrorCode save_VTK(const char* path, const char* filename, Vec sol, DM dm, std::array<double, 3> h)
  {
      PetscErrorCode ierr;
      DM             dau, dap;
      DMDALocalInfo  infou, infop;
      Vec            solu, solp;
      Vec            locsolu, locsolp;
      PetscScalar    ****psolu, ***psolp;

      ierr = DMCompositeGetEntries(dm, &dau, &dap);CHKERRQ(ierr);
      ierr = DMDAGetLocalInfo(dau, &infou);CHKERRQ(ierr);
      ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);

      ierr = DMCompositeGetAccess(dm,sol,&solu,&solp);CHKERRQ(ierr);

      ierr = DMGetLocalVector(dau, &locsolu);CHKERRQ(ierr);
      ierr = DMGetLocalVector(dap, &locsolp);CHKERRQ(ierr);

      ierr = DMGlobalToLocalBegin(dau, solu, INSERT_VALUES, locsolu);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dau, solu, INSERT_VALUES, locsolu);CHKERRQ(ierr);

      ierr = DMGlobalToLocalBegin(dap, solp, INSERT_VALUES, locsolp);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dap, solp, INSERT_VALUES, locsolp);CHKERRQ(ierr);

      ierr = DMDAVecGetArrayDOF(dau, locsolu, &psolu);CHKERRQ(ierr);
      ierr = DMDAVecGetArray(dap, locsolp, &psolp);CHKERRQ(ierr);
      
      vtkStructuredGrid* pressureDataSet = vtkStructuredGrid::New();
      pressureDataSet->SetExtent(infop.gxs, infop.gxs+infop.gxm-1, infop.gys, infop.gys+infop.gym-1, infop.gzs, infop.gzs+infop.gzm-1);

      vtkStructuredGrid* velocityDataSet = vtkStructuredGrid::New();
      velocityDataSet->SetExtent(infou.gxs, infou.gxs+infou.gxm-1, infou.gys, infou.gys+infou.gym-1, infou.gzs, infou.gzs+infou.gzm-1);
      
      vtkPoints* pressurePoints = vtkPoints::New();
      vtkPoints* velocityPoints = vtkPoints::New();
      
      vtkDoubleArray* velocity = vtkDoubleArray::New();
      velocity->SetNumberOfComponents(3);
      velocity->SetName("velocity");
      
      vtkDoubleArray* pressure = vtkDoubleArray::New();
      pressure->SetName("pressure");
      
      for (int k=infop.gzs; k<infop.gzs+infop.gzm; k++)
        for (int j=infop.gys; j<infop.gys+infop.gym; j++)
          for (int i=infop.gxs; i<infop.gxs+infop.gxm; i++)
            {
              pressurePoints->InsertNextPoint(2*i*h[0], 2*j*h[1], 2*k*h[2]);
              pressure->InsertNextValue(psolp[k][j][i]);
          }

      for (int k=infou.gzs; k<infou.gzs+infou.gzm; k++)
        for (int j=infou.gys; j<infou.gys+infou.gym; j++)
          for (int i=infou.gxs; i<infou.gxs+infou.gxm; i++)
            {
              velocityPoints->InsertNextPoint(i*h[0], j*h[1], k*h[2]);
              velocity->InsertNextTuple3(psolu[k][j][i][0], psolu[k][j][i][1], psolu[k][j][i][2]);
          }

      pressureDataSet->SetPoints(pressurePoints);
      pressureDataSet->GetPointData()->SetScalars(pressure);

      velocityDataSet->SetPoints(velocityPoints);
      velocityDataSet->GetPointData()->SetVectors(velocity);
      
      vtkXMLStructuredGridWriter* pressureDataWriter = vtkXMLStructuredGridWriter::New();
      vtkXMLStructuredGridWriter* velocityDataWriter = vtkXMLStructuredGridWriter::New();
      std::stringstream op, ov;
      int rank;
      int size;
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      MPI_Comm_size(PETSC_COMM_WORLD, &size);

      op << path << "/" << filename << "_pressure_" << rank << ".vts";
      ov << path << "/" << filename << "_velocity_" << rank << ".vts";

      pressureDataWriter->SetFileName(op.str().data());
    #if VTK_MAJOR_VERSION <= 5
      pressureDataWriter->SetInput(pressureDataSet);
    #else
      pressureDataWriter->SetInputData(pressureDataSet);
    #endif
      pressureDataWriter->Write();

      velocityDataWriter->SetFileName(ov.str().data());
    #if VTK_MAJOR_VERSION <= 5
      velocityDataWriter->SetInput(velocityDataSet);
    #else
      velocityDataWriter->SetInputData(velocityDataSet);
    #endif
      velocityDataWriter->Write();

      int x[size], y[size], z[size];
      int sx[size], sy[size], sz[size];
      MPI_Gather(&infop.gxs, 1, MPI_INT, x, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infop.gxm, 1, MPI_INT, sx, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      MPI_Gather(&infop.gys, 1, MPI_INT, y, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infop.gym, 1, MPI_INT, sy, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      MPI_Gather(&infop.gzs, 1, MPI_INT, z, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infop.gzm, 1, MPI_INT, sz, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      if (rank == 0){
        std::stringstream oall;
        oall << path << "/" << filename << "_pressure.pvts";
        ofstream pvts;
        pvts.open(oall.str().data(), ios::out | ios::trunc);
        pvts << "<?xml version=\"1.0\"?>" << endl;
        pvts << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
        pvts << "<PStructuredGrid WholeExtent=\"0 " << infop.mx-1 << " 0 " << infop.my-1 << " 0 " << infop.mz-1 << "\" GhostLevel=\"0\">";
        pvts << "<PPoints>" << endl;
        pvts << "<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << endl;
        pvts << "</PPoints>" << endl;
        pvts << "<PPointData Scalars=\"pressure\">" << endl;
        pvts << "<PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
        pvts << "</PPointData>" << endl;
        pvts << "<PCells>" << endl;
        pvts << "</PCells>" << endl;
        for(int i=0; i<size; i++)
          pvts << "<Piece Extent=\"" << x[i] << " " << x[i]+sx[i]-1 << " " << y[i] << " " << y[i]+sy[i]-1 << " " << z[i] << " " << z[i]+sz[i]-1 << "\" Source=\"" << filename << "_pressure_" << i << ".vts\"/>" << endl;
            pvts << "</PStructuredGrid>" << endl;
        pvts << "</VTKFile>" << endl;
        pvts.close();
      }

      MPI_Gather(&infou.gxs, 1, MPI_INT, x, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infou.gxm, 1, MPI_INT, sx, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      MPI_Gather(&infou.gys, 1, MPI_INT, y, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infou.gym, 1, MPI_INT, sy, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      MPI_Gather(&infou.gzs, 1, MPI_INT, z, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Gather(&infou.gzm, 1, MPI_INT, sz, 1, MPI_INT, 0, PETSC_COMM_WORLD);

      if (rank == 0){
        std::stringstream oall;
        oall << path << "/" << filename << "_velocity.pvts";
        ofstream pvts;
        pvts.open(oall.str().data(), ios::out | ios::trunc);
        pvts << "<?xml version=\"1.0\"?>" << endl;
        pvts << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
        pvts << "<PStructuredGrid WholeExtent=\"0 " << infou.mx-1 << " 0 " << infou.my-1 << " 0 " << infou.mz-1 << "\" GhostLevel=\"0\">";
        pvts << "<PPoints>" << endl;
        pvts << "<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << endl;
        pvts << "</PPoints>" << endl;
        pvts << "<PPointData Vectors=\"velocity\">" << endl;
        pvts << "<PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
        pvts << "</PPointData>" << endl;
        pvts << "<PCells>" << endl;
        pvts << "</PCells>" << endl;
        for(int i=0; i<size; i++)
          pvts << "<Piece Extent=\"" << x[i] << " " << x[i]+sx[i]-1 << " " << y[i] << " " << y[i]+sy[i]-1 << " " << z[i] << " " << z[i]+sz[i]-1 << "\" Source=\"" << filename << "_velocity_" << i << ".vts\"/>" << endl;
            pvts << "</PStructuredGrid>" << endl;
        pvts << "</VTKFile>" << endl;
        pvts.close();
      }

      pressurePoints->Delete();
      velocityPoints->Delete();
      velocity->Delete();
      pressure->Delete();
      pressureDataSet->Delete();
      pressureDataWriter->Delete();
      velocityDataSet->Delete();
      velocityDataWriter->Delete();

      ierr = DMDAVecRestoreArray(dau, locsolu, &psolu);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(dap, locsolp, &psolp);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(dau, &locsolu);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(dap, &locsolp);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(dm,sol,&solu,&solp);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "saveParticles"
    template<typename Shape>
    PetscErrorCode saveParticles(const char* path, const char* filename,
                                 std::vector<particle<Shape>>const& particles)
    {
      PetscErrorCode ierr;
      int rank;
      ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);

      if (rank == 0){

        vtkPoints* spheresPoints = vtkPoints::New();
        vtkDoubleArray* shperesRadius = vtkDoubleArray::New();
        vtkPolyData* data = vtkPolyData::New();
        PetscFunctionBeginUser;

        data->Allocate(1, 1);

        spheresPoints->SetDataTypeToDouble();
        spheresPoints->SetNumberOfPoints(particles.size());
        data->SetPoints(spheresPoints);
        shperesRadius->SetName("radius");
        shperesRadius->SetNumberOfComponents(1);
        shperesRadius->SetNumberOfTuples(particles.size());
        data->GetPointData()->AddArray(shperesRadius);


        std::size_t i = 0;
        for(auto& p: particles){
          spheresPoints->SetPoint(i, p.center_[0], p.center_[1], p.center_[2]);
          shperesRadius->SetValue(i, p.shape_factors_[0]);
          i++;
        }

        //data->Update();

        std::stringstream output;
        output << path << "/" << filename << ".vtp";
        vtkXMLPolyDataWriter* writer = vtkXMLPolyDataWriter::New();
        writer->SetFileName(output.str().data());
      #if VTK_MAJOR_VERSION <= 5
        writer->SetInput(data);
      #else
        writer->SetInputData(data);
      #endif

        writer->Write();
      }
      PetscFunctionReturn(0);
    }

  }
}
#endif