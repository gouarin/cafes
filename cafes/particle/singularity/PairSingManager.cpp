#ifndef PAIRSINGMANAGER_CPP
#define PAIRSINGMANAGER_CPP

#include <petsc.h>
#include <DtoN/DtoNProblem.h>

#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <cmath>

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
#include "stokes/stokes.h"
#include "DtoN/UandPNormal.h"
#include "DtoN/UandPTang.h"
#include "DtoN/Truncation.h"
#include "../src/DtoN/PairSingInfo.cpp"
#include "../src/DtoN/forcesManager.cpp"


namespace DtoN{

  template<typename PTYPE, std::size_t Dimensions>
  struct Display_impl{};

  template<typename PTYPE, std::size_t Dimensions>
  struct  ComputesingularST_impl{};

  template<typename PTYPE, std::size_t Dimensions>
  struct  ComputesingularBC_impl{};

  template<typename PTYPE, std::size_t Dimensions>
  class PairSingManager : public PairSingInfo<PTYPE, Dimensions>{

    static double mu;

    std::vector<PetscInt>  _nMin;
    std::vector<PetscInt>  _nMax;

    std::vector<double>    _h;
    
    bool                   _testbox;
    bool                   _isSingularity;
    double                 _cutoffdist;
    double                 _a;
    double                 _UN;
    double                 _UT;
    double                 _param;
    int                    _NLOC;

    const char* _path;
    const char* _filename;

  public:
    PairSingManager(){std::cout << "DEFAULT" << std::endl;}
    PairSingManager(PTYPE& part1, PTYPE& part2, std::vector<double> h, const char* path, const char* filename);

    //Simple accessor
    double geta(){return _a;}
    

    void initializebox(DMDALocalInfo& infoU);
    bool isTestbox(){return _testbox;}
    bool isSingularity(){return _isSingularity;}

    friend struct Display_impl<PTYPE, Dimensions>;
    PetscErrorCode display();

    friend struct ComputesingularST_impl<PTYPE, Dimensions>;
    PetscErrorCode computesingularST(DM& dau,Vec& ulocal, Vec& u){return ComputesingularST_impl<PTYPE,Dimensions>::computesingularST(*this,dau,ulocal,u);}

    friend struct ComputesingularBC_impl<PTYPE, Dimensions>;
    PetscErrorCode computesingularBC(std::vector<PTYPE> pair){return ComputesingularBC_impl<PTYPE,Dimensions>::computesingularForces(*this,pair);}


    PetscErrorCode computesingularForces();

    void getForces(std::vector<double>& F);	 
    void printForces();
    void clearForces(ForcesManager<Dimensions>& Forces);

    void printForcesNormal();
    void printForcesTang();

    // void printNormalDIV();
    // void printTangDIV();
    // void printDIV();
  };

  template<class PTYPE, std::size_t Dimensions> double PairSingManager<PTYPE, Dimensions>::mu = 1.0;

#undef __FUNCT__
#define __FUNCT__ "display_"
  template<class PTYPE, std::size_t Dimensions>
  PetscErrorCode PairSingManager<PTYPE, Dimensions>::display(){

    PetscErrorCode ierr;
    PetscFunctionBegin;
    vtkStructuredGrid* velocity_singDataSet;
    vtkStructuredGrid* pressure_singDataSet;

    velocity_singDataSet = vtkStructuredGrid::New();
    pressure_singDataSet = vtkStructuredGrid::New();
    if(Dimensions==2){
      velocity_singDataSet->SetExtent(_nMin[0]*_NLOC, _nMax[0]*_NLOC-1, _nMin[1]*_NLOC, _nMax[1]*_NLOC-1, 0, 0);
      pressure_singDataSet->SetExtent(_nMin[0]*_NLOC, _nMax[0]*_NLOC-1, _nMin[1]*_NLOC, _nMax[1]*_NLOC-1, 0, 0);
    }
    else{
      velocity_singDataSet->SetExtent(_nMin[0]*_NLOC, _nMax[0]*_NLOC-1, _nMin[1]*_NLOC, _nMax[1]*_NLOC-1,_nMin[2]*_NLOC, _nMax[2]*_NLOC-1);
      pressure_singDataSet->SetExtent(_nMin[0]*_NLOC, _nMax[0]*_NLOC-1, _nMin[1]*_NLOC, _nMax[1]*_NLOC-1, _nMin[2]*_NLOC, _nMax[2]*_NLOC-1);
    }

    vtkPoints* pressure_singPoints = vtkPoints::New();
    vtkPoints* velocity_singPoints = vtkPoints::New();

    vtkDoubleArray* velocity_sing = vtkDoubleArray::New();
    velocity_sing->SetNumberOfComponents(3);
    velocity_sing->SetName("velocity_sing");

    vtkDoubleArray* pressure_sing = vtkDoubleArray::New();
    pressure_sing->SetName("pressure_sing");

    
    Display_impl<PTYPE, Dimensions>::display(*this,pressure_singPoints,velocity_singPoints,pressure_sing,velocity_sing);

    velocity_singDataSet->SetPoints(velocity_singPoints);
    velocity_singDataSet->GetPointData()->SetScalars(velocity_sing);

    pressure_singDataSet->SetPoints(pressure_singPoints);
    pressure_singDataSet->GetPointData()->SetScalars(pressure_sing);

    vtkXMLStructuredGridWriter* velocity_singDataWriter = vtkXMLStructuredGridWriter::New();
    stringstream oc;

    oc << _path << "/" << _filename << "_velocity_sing_" << 0 << "_" << 1 << ".vts";//a changer
    velocity_singDataWriter->SetFileName(oc.str().data());
    //dataWriter->SetDataModeToAscii();
#if VTK_MAJOR_VERSION <= 5
    velocity_singDataWriter->SetInput(velocity_singDataSet);
#else
    velocity_singDataWriter->SetInputData(velocity_singDataSet);
#endif
    velocity_singDataWriter->Write();

    vtkXMLStructuredGridWriter* pressure_singDataWriter = vtkXMLStructuredGridWriter::New();
    stringstream occ;

    occ << _path << "/" << _filename << "_pressure_sing_" << 0 << "_" << 1 << ".vts";//a changer

    pressure_singDataWriter->SetFileName(occ.str().data());
    //dataWriter->SetDataModeToAscii();
#if VTK_MAJOR_VERSION <= 5
    pressure_singDataWriter->SetInput(pressure_singDataSet);
#else
    pressure_singDataWriter->SetInputData(pressure_singDataSet);
#endif
    pressure_singDataWriter->Write();

    velocity_singPoints->Delete();
    velocity_sing->Delete();
    velocity_singDataSet->Delete();
    velocity_singDataWriter->Delete();

    pressure_singPoints->Delete();
    pressure_sing->Delete();
    pressure_singDataSet->Delete();
    pressure_singDataWriter->Delete();

    PetscFunctionReturn(0);
  }

  template<typename PTYPE>
  struct Display_impl<PTYPE, 2>{
    static void display(PairSingManager<PTYPE, 2>& p, vtkPoints*& pressure_singPoints, vtkPoints*& velocity_singPoints, vtkDoubleArray*& pressure_sing, vtkDoubleArray*& velocity_sing){

      std::vector<double> posNetwork(2);
      std::vector<double> posNetworkDec(2); // On place l'origine
      std::vector<double> posNetworkRefPart(2); // On utilise la base dans ref particules
      std::vector<double> posC1Diff(2);
      std::vector<double> posC2Diff(2);

      double divU = 0;

      p._param = p._cutoffdist*p._cutoffdist/2.;
      std::vector<double> U;
      p.getU(U);
      p._UN = U[0];
      p._UT = U[1];

      for(int j=p._nMin[1]; j<p._nMax[1]; j++){
        for (int jloc=0; jloc<p._NLOC; jloc++){
          for(int i=p._nMin[0]; i<p._nMax[0]; i++){
            for (int iloc=0; iloc<p._NLOC; iloc++){

              posNetwork.clear();
              posNetwork.push_back( (i+iloc/double(p._NLOC))*p._h[0]);
              posNetwork.push_back( (j+jloc/double(p._NLOC))*p._h[1]);


              //distance between network dot and particles' centers
              double distC1, distC2;
              auto diff = [](double x, double y){return x-y;};
              std::transform(p._center_p1.cbegin(), p._center_p1.cend(), posNetwork.begin(), posC1Diff.begin(), diff); 
              distC1 = std::inner_product(posC1Diff.begin(), posC1Diff.end(),posC1Diff.begin(), 0.);

              std::transform(p._center_p2.cbegin(), p._center_p2.cend(), posNetwork.begin(), posC2Diff.begin(), diff); 
              distC2 = std::inner_product(posC2Diff.begin(), posC2Diff.end(),posC2Diff.begin(), 0.);

              if((distC1>=1./p._H1/p._H1)&&(distC2>=1./p._H2/p._H2))
              {
                std::transform(posNetwork.begin(), posNetwork.end(), p._posOrigin.begin(), posNetworkDec.begin(), diff); 
                for(int i_dim=0; i_dim<2;i_dim++)
                  posNetworkRefPart[i_dim] = std::inner_product(posNetworkDec.begin(), posNetworkDec.end(), p._base[i_dim].begin(), 0.);

                //attention : ur = uy et uz= ux
                auto ux = ux_sing_normalMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL) + ux_sing_tangMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL);
                auto uz = uz_sing_normalMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL) + uz_sing_tangMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL);


                // Add points to vtk + singular value to vtk
                velocity_singPoints->InsertNextPoint(posNetwork[0], posNetwork[1], 0.);
                velocity_sing->InsertNextTuple3(ux*p._base[0][0]+uz*p._base[1][0],
                                                ux*p._base[0][1]+uz*p._base[1][1], 0.);


                // Add points to vtk + singular value to vtk
                pressure_singPoints->InsertNextPoint(posNetwork[0], posNetwork[1], 0.);
                pressure_sing->InsertNextValue(p_sing_withT_normalMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL  ) + p_sing_withT_tangMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL  ));


                // A CHANGER : FAIRE FONCTION POUR RECUPERER LE RESEAU A PART
                divU = max(divU, DIVCART_normalMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL  ) + DIVCART_tangMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL  ));



              }
              else{
                // Add points to vtk + singular value to vtk
                velocity_singPoints->InsertNextPoint(posNetwork[0], posNetwork[1], 0.);
                velocity_sing->InsertNextTuple3(0., 0., 0.);

                // Add points to vtk + singular value to vtk
                pressure_singPoints->InsertNextPoint(posNetwork[0], posNetwork[1], 0.);
                pressure_sing->InsertNextValue(0.);          
              }
            }
          }
        }
      }
      std::cout << "Divergence value : " << divU << std::endl;
    }
  };

  template<typename PTYPE>
  struct Display_impl<PTYPE, 3>{
    static void display(PairSingManager<PTYPE, 3>& p, vtkPoints*& pressure_singPoints, vtkPoints*& velocity_singPoints, vtkDoubleArray*& pressure_sing, vtkDoubleArray*& velocity_sing){

      std::vector<double> posNetwork(3);
      std::vector<double> posNetworkDec(3); // On place l'origine
      std::vector<double> posNetworkRefPart(3); // On utilise la base dans ref particules
      std::vector<double> posC1Diff(3);
      std::vector<double> posC2Diff(3);

      p._param = p._cutoffdist*p._cutoffdist/2.;
      std::vector<double> U;
      p.getU(U);
      p._UN = U[2];std::cout << "UN = " << p._UN << std::endl;
      p._UT = U[0];std::cout << "UT = " << p._UT << std::endl;

      for(int k=p._nMin[2]; k<p._nMax[2]; k++){
       for (int kloc=0; kloc<p._NLOC; kloc++){
         for(int j=p._nMin[1]; j<p._nMax[1]; j++){
           for (int jloc=0; jloc<p._NLOC; jloc++){
             for(int i=p._nMin[0]; i<p._nMax[0]; i++){
              for (int iloc=0; iloc<p._NLOC; iloc++){

                posNetwork.clear();
                posNetwork.push_back( (i+iloc/double(p._NLOC))*p._h[0]);
                posNetwork.push_back( (j+jloc/double(p._NLOC))*p._h[1]);
                posNetwork.push_back( (k+kloc/double(p._NLOC))*p._h[2]);

                //distance between network dot and particles' centers
                double distC1, distC2;
                auto diff = [](double x, double y){return x-y;};
                std::transform(p._center_p1.cbegin(), p._center_p1.cend(), posNetwork.begin(), posC1Diff.begin(), diff); 
                distC1 = std::inner_product(posC1Diff.begin(), posC1Diff.end(),posC1Diff.begin(), 0.);

                std::transform(p._center_p2.cbegin(), p._center_p2.cend(), posNetwork.begin(), posC2Diff.begin(), diff); 
                distC2 = std::inner_product(posC2Diff.begin(), posC2Diff.end(),posC2Diff.begin(), 0.);

                if((distC1>=1./p._H1/p._H1)&&(distC2>=1./p._H2/p._H2))
                {
                  std::transform(posNetwork.begin(), posNetwork.end(), p._posOrigin.begin(), posNetworkDec.begin(), diff); 
                  for(int i_dim=0; i_dim<3;i_dim++)
                   posNetworkRefPart[i_dim] = std::inner_product(posNetworkDec.begin(), posNetworkDec.end(), p._base[i_dim].begin(), 0.);

                 auto ux =  ux_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL) + ux_sing_tangMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL);
                 auto uy =  uy_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL) + uy_sing_tangMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL);
                 auto uz =  uz_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL) + uz_sing_tangMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL);

                 //std::cout  << "Display : " << ux << " " << uy << " " << uz << std::endl;
                 // Add points to vtk + singular value to vtk 
                 //velocity_singPoints->InsertNextPoint(posNetwork[0], posNetwork[1], posNetwork[2]);
                 velocity_sing->InsertNextTuple3(ux*p._base[0][0]+uy*p._base[1][0]+uz*p._base[2][0],
                                                 ux*p._base[0][1]+uy*p._base[1][1]+uz*p._base[2][1],
                                                 ux*p._base[0][2]+uy*p._base[1][2]+uz*p._base[2][2]);
                 //velocity_sing->InsertNextTuple3(drW0_tang3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL),0,0);
                 // Add points to vtk + singular value to vtk
                 pressure_singPoints->InsertNextPoint(posNetwork[0], posNetwork[1], posNetwork[2]);
                 pressure_sing->InsertNextValue(p_sing_withT_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL  ) + p_sing_withT_tangMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL  ));
               }
               else{
                // Add points to vtk + singular value to vtk
                velocity_singPoints->InsertNextPoint(posNetwork[0], posNetwork[1], posNetwork[2]);		    
                velocity_sing->InsertNextTuple3(0., 0., 0.);

                // Add points to vtk + singular value to vtk
                pressure_singPoints->InsertNextPoint(posNetwork[0], posNetwork[1], posNetwork[2]);
                pressure_sing->InsertNextValue(0.);          
              }
            }
          }
        }
      }
    }
  }
}
};


  template<typename PTYPE>
struct ComputesingularST_impl<PTYPE, 2>{

#undef __FUNCT__
#define __FUNCT__ "computesingularST_"
  static PetscErrorCode computesingularST(PairSingManager<PTYPE, 2>& p, DM& dau,Vec& ulocal, Vec& u){

    PetscErrorCode ierr;
    PetscFunctionBegin;

    PetscScalar ***pu;
    ierr = DMCreateLocalVector(dau, &ulocal);CHKERRQ(ierr);
    ierr = VecSet(ulocal, 0.);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(dau, ulocal, &pu);CHKERRQ(ierr);

    std::vector<double> posNetwork       (2);	
      std::vector<double> posNetworkDec    (2); // On place l'origine
      std::vector<double> posNetworkRefPart(2); // On utilise la base dans ref particules

      std::vector<double> posLoc           (2);

      std::vector<double> posC1Diff        (2);
      std::vector<double> posC2Diff        (2);


      for(int j=p._nMin[1]; j<p._nMax[1]; j++){
       for (int jloc=0; jloc<p._NLOC; jloc++){
         for(int i=p._nMin[0]; i<p._nMax[0]; i++){
           for (int iloc=0; iloc<p._NLOC; iloc++){

              posNetwork.clear();
              posNetwork.push_back( (i+iloc/double(p._NLOC))*p._h[0]);
              posNetwork.push_back( (j+jloc/double(p._NLOC))*p._h[1]);

              posLoc.clear();
              posLoc.push_back(p._h[0]*iloc/double(p._NLOC));
              posLoc.push_back(p._h[1]*jloc/double(p._NLOC));


              //distance between network dot and particles' centers
              double distC1, distC2;
              auto diff = [](double x, double y){return x-y;};
              std::transform(p._center_p1.cbegin(), p._center_p1.cend(), posNetwork.begin(), posC1Diff.begin(), diff); 
              distC1 = std::inner_product(posC1Diff.begin(), posC1Diff.end(),posC1Diff.begin(), 0.);

              std::transform(p._center_p2.cbegin(), p._center_p2.cend(), posNetwork.begin(), posC2Diff.begin(), diff); 
              distC2 = std::inner_product(posC2Diff.begin(), posC2Diff.end(),posC2Diff.begin(), 0.);

              if((distC1>=1./p._H1/p._H1)&&(distC2>=1./p._H2/p._H2))
              {
                std::transform(posNetwork.begin(), posNetwork.end(), p._posOrigin.begin(), posNetworkDec.begin(), diff); 
                for(int i_dim=0; i_dim<2;i_dim++)
                  posNetworkRefPart[i_dim] = std::inner_product(posNetworkDec.begin(), posNetworkDec.end(), p._base[i_dim].begin(), 0.);

                std::vector< std::vector<double> > gradUsing       (2);
                std::vector< std::vector<double> > gradUsingRefPart(2);
                double                             divUsing;

                gradUsingRefPart[0].resize(2);
                gradUsingRefPart[0][0]     = dxux_sing_normalMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);
                gradUsingRefPart[0][1]     = dzux_sing_normalMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);

                gradUsingRefPart[1].resize(2);
                gradUsingRefPart[1][0]     = dxuz_sing_normalMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);
                gradUsingRefPart[1][1]     = dzuz_sing_normalMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);

                //changement de base GradUsing 
                gradUsing[0].resize(2);
                gradUsing[1].resize(2);
                for(int k=0; k<2;k++)
                  for(int l=0; l<2;l++)
                    for(int m=0; m<2;m++)
                     for(int n=0; n<2;n++)
                       gradUsing[k][l] += p._base[m][k]*p._base[n][l]*gradUsingRefPart[m][n];



                //fonction test (x-hx)*(y-hy) /(hx*hy)--------------------------------
                getGrad(p,posLoc[1]-p._h[1], posLoc[0]-p._h[0], pu, posNetworkRefPart, gradUsing,i,j);
                //--------------------------------------------------------------------


                //fonction test -x*(y-hy) /(hx*hy)------------------------------------
                getGrad(p,-(posLoc[1]-p._h[1]), -posLoc[0], pu, posNetworkRefPart, gradUsing,i+1,j);
                //--------------------------------------------------------------------


                //fonction test -y*(x-hx) /(hx*hy)------------------------------------
                getGrad(p,-posLoc[1], -(posLoc[0]-p._h[0]), pu, posNetworkRefPart, gradUsing,i,j+1);
                //--------------------------------------------------------------------


                //fonction test x*y/(hx*hy)-------------------------------------------
                getGrad(p,posLoc[1], posLoc[0], pu, posNetworkRefPart, gradUsing,i+1,j+1);
                //--------------------------------------------------------------------

               }
               }
             }
           }
         }
         ierr = DMDAVecRestoreArrayDOF(dau, ulocal, &pu);CHKERRQ(ierr);
         ierr = DMLocalToGlobalBegin(dau, ulocal, ADD_VALUES, u);CHKERRQ(ierr);
         ierr = DMLocalToGlobalEnd(dau, ulocal, ADD_VALUES, u);CHKERRQ(ierr);
         ierr = VecDestroy(&ulocal);CHKERRQ(ierr);

         PetscFunctionReturn(0);
       }

       static void getGrad(PairSingManager<PTYPE, 2>& p,double gradvX,double gradvY, PetscScalar*** pu, std::vector<double>& posNetworkRefPart, std::vector< std::vector<double> >& gradUsing, int i, int j){ 

        std::vector< std::vector<double> > gradv(2);
        gradv[0].resize(2);
        gradv[1].resize(2);


        gradv[0][0]     = gradvX/(p._NLOC*p._NLOC);
        gradv[0][1]     = gradvY/(p._NLOC*p._NLOC);
        gradv[1][0]     = gradv[0][0];
        gradv[1][1]     = gradv[1][0];


        for(int l=0; l<2;l++){
         for(int k=0; k<2;k++)
           pu[j][i][l] += p.mu*gradUsing[l][k]*gradv[l][k];
         pu[j][i][l] -= p_sing_withT_normalMvt2D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL)*gradv[l][l];
       }
     }
   };

  template<typename PTYPE>
   struct ComputesingularST_impl<PTYPE, 3>{
#undef __FUNCT__
#define __FUNCT__ "computesingularST_"
    static PetscErrorCode computesingularST(PairSingManager<PTYPE, 3>& p, DM& dau,Vec& ulocal, Vec u){

      PetscErrorCode ierr;
      PetscFunctionBegin;
      PetscScalar ****pu;
      ierr = DMCreateLocalVector(dau, &ulocal);CHKERRQ(ierr);
      ierr = VecSet(ulocal, 0.);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(dau, ulocal, &pu);CHKERRQ(ierr);

      std::vector<double> posNetwork       (3);
      std::vector<double> posNetworkDec    (3); // On place l'origine
      std::vector<double> posNetworkRefPart(3); // On utilise la base dans ref particules

      std::vector<double> posLoc           (3);

      std::vector<double> posC1Diff        (3);
      std::vector<double> posC2Diff        (3);

      for(int k=p._nMin[2]; k<p._nMax[2]; k++){
       for (int kloc=0; kloc<p._NLOC; kloc++){
         for(int j=p._nMin[1]; j<p._nMax[1]; j++){
           for (int jloc=0; jloc<p._NLOC; jloc++){
             for(int i=p._nMin[0]; i<p._nMax[0]; i++){
              for (int iloc=0; iloc<p._NLOC; iloc++){

                posNetwork.clear();
                posNetwork.push_back( (i+iloc/double(p._NLOC))*p._h[0]);
                posNetwork.push_back( (j+jloc/double(p._NLOC))*p._h[1]);
                posNetwork.push_back( (k+kloc/double(p._NLOC))*p._h[2]);

                posLoc.clear();
                posLoc.push_back(p._h[0]*iloc/double(p._NLOC));
                posLoc.push_back(p._h[1]*jloc/double(p._NLOC));
                posLoc.push_back(p._h[1]*kloc/double(p._NLOC));

		  //distance between network dot and particles' centers
                double distC1, distC2;
                auto diff = [](double x, double y){return x-y;};
                std::transform(p._center_p1.cbegin(), p._center_p1.cend(), posNetwork.begin(), posC1Diff.begin(), diff); 
                distC1 = std::inner_product(posC1Diff.begin(), posC1Diff.end(),posC1Diff.begin(), 0.);

                std::transform(p._center_p2.cbegin(), p._center_p2.cend(), posNetwork.begin(), posC2Diff.begin(), diff); 
                distC2 = std::inner_product(posC2Diff.begin(), posC2Diff.end(),posC2Diff.begin(), 0.);

                if((distC1>=1./p._H1/p._H1)&&(distC2>=1./p._H2/p._H2))
                {

                  std::vector< std::vector<double> > gradUsing       (3);
                  std::vector< std::vector<double> > gradUsingRefPart(3);
                  double                             divUsing;

                  gradUsingRefPart[0].resize(3);
                  gradUsingRefPart[0][0]     = dxux_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);
                  gradUsingRefPart[0][1]     = dyux_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);
                  gradUsingRefPart[0][2]     = dzux_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);

                  gradUsingRefPart[1].resize(3);
                  gradUsingRefPart[1][0]     = dxuy_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);
                  gradUsingRefPart[1][1]     = dyuy_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);
                  gradUsingRefPart[1][2]     = dzuy_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);

                  gradUsingRefPart[2].resize(3);
                  gradUsingRefPart[2][0]     = dxuy_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);
                  gradUsingRefPart[2][1]     = dyuy_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);
                  gradUsingRefPart[2][2]     = dzuy_sing_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL);

		      //changement de base GradUsing 
                  gradUsing[0].resize(3);
                  gradUsing[1].resize(3);
                  gradUsing[2].resize(3);

                  for(int l=0; l<3;l++)
                   for(int m=0; m<3;m++)
                     for(int n=0; n<3;n++)
                       for(int o=0; o<3;o++)
                         gradUsing[l][m] += p._base[n][l]*p._base[o][m]*gradUsingRefPart[n][o];



		      //fonction test -(x-hx)*(y-hy)*(z-hz) /(hx*hy*hz)-------------------------
                       getGrad(p,-(posLoc[1]-p._h[1])*(posLoc[2]-p._h[2]), -(posLoc[0]-p._h[0])*(posLoc[2]-p._h[2]), -(posLoc[0]-p._h[0])*(posLoc[1]-p._h[1]),pu, posNetworkRefPart, gradUsing,i,j,k);
		      //-------------------------------------------------------------------------

		      //fonction test x*(y-hy)*(z-hz) /(hx*hy*hz)--------------------------------
                       getGrad(p,(posLoc[1]-p._h[1])*(posLoc[2]-p._h[2]), posLoc[0]*(posLoc[2]-p._h[2]), posLoc[0]*(posLoc[1]-p._h[1]),pu, posNetworkRefPart, gradUsing,i+1,j,k);
		      //-------------------------------------------------------------------------

		      //fonction test (x-hx)*y*(z-hz) /(hx*hy*hz)--------------------------------
                       getGrad(p,posLoc[1]*(posLoc[2]-p._h[2]), (posLoc[0]-p._h[0])*(posLoc[2]-p._h[2]), (posLoc[0]-p._h[0])*posLoc[1],pu, posNetworkRefPart, gradUsing,i,j+1,k);
		      //--------------------------------------------------------------------

		      //fonction test (x-hx)*(y-hy)*z /(hx*hy*hz)--------------------------------
                       getGrad(p,(posLoc[1]-p._h[1])*posLoc[2], (posLoc[0]-p._h[0])*posLoc[2], (posLoc[0]-p._h[0])*(posLoc[1]-p._h[1]),pu, posNetworkRefPart, gradUsing,i,j,k+1);
		      //--------------------------------------------------------------------

		      //fonction test - (x)*(y)*(z-hz) /(hx*hy*hz)--------------------------------
                       getGrad(p,-posLoc[1]*(posLoc[2]-p._h[2]), -posLoc[0]*(posLoc[2]-p._h[2]), -posLoc[0]*posLoc[1],pu, posNetworkRefPart, gradUsing,i+1,j+1,k);
		      //--------------------------------------------------------------------

		      //fonction test -(x)*(y-hy)*(z) /(hx*hy*hz)--------------------------------
                       getGrad(p,-(posLoc[1]-p._h[1])*posLoc[2], -posLoc[0]*posLoc[2], -posLoc[0]*(posLoc[1]-p._h[1]),pu, posNetworkRefPart, gradUsing,i+1,j,k+1);
		      //--------------------------------------------------------------------

		      //fonction test -(x-hx)*(y)*(z) /(hx*hy*hz)--------------------------------
                       getGrad(p,-posLoc[1]*posLoc[2], -(posLoc[0]-p._h[0])*posLoc[2], -(posLoc[0]-p._h[0])*posLoc[1],pu, posNetworkRefPart, gradUsing,i,j+1,k+1);
		      //--------------------------------------------------------------------

		      //fonction test (x)*(y)*(z) /(hx*hy*hz)--------------------------------
                       getGrad(p,posLoc[1]*posLoc[2], posLoc[0]*posLoc[2], posLoc[0]*posLoc[1],pu, posNetworkRefPart, gradUsing,i+1,j+1,k+1);
		      //--------------------------------------------------------------------
                     }
                   }
                 }
               }
             }
           }
         }


         ierr = DMDAVecRestoreArrayDOF(dau, ulocal, &pu);CHKERRQ(ierr);
         ierr = DMLocalToGlobalBegin(dau, ulocal, ADD_VALUES, u);CHKERRQ(ierr);
         ierr = DMLocalToGlobalEnd(dau, ulocal, ADD_VALUES, u);CHKERRQ(ierr);
         ierr = VecDestroy(&ulocal);CHKERRQ(ierr);

         PetscFunctionReturn(0);
       }

       void getGrad(PairSingManager<PTYPE, 3>& p, double gradvX,double gradvY, double gradvZ, PetscScalar**** pu, std::vector<double>& posNetworkRefPart, std::vector< std::vector<double> >& gradUsing, int i, int j, int k){ 

        std::vector< std::vector<double> > gradv(3);
        gradv[0].resize(3);
        gradv[1].resize(3);
        gradv[2].resize(3);

        gradv[0][0]     = gradvX/(p._NLOC*p._NLOC*p._NLOC);
        gradv[0][1]     = gradvY/(p._NLOC*p._NLOC*p._NLOC);
        gradv[0][2]     = gradvZ/(p._NLOC*p._NLOC*p._NLOC);

        gradv[1][0]     = gradv[0][0];
        gradv[1][1]     = gradv[0][1];
        gradv[1][2]     = gradv[0][2];

        gradv[2][0]     = gradv[0][0];
        gradv[2][1]     = gradv[0][1];
        gradv[2][2]     = gradv[0][2];

        for(int m=0; m<3;m++){
         for(int l=0; l<3;l++)
           pu[k][j][i][m] += p.mu*gradUsing[m][l]*gradv[m][l];
         pu[k][j][i][m] -= p_sing_withT_normalMvt3D(posNetworkRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL)*gradv[m][m];
       }
     }
   };


  template<class PTYPE, std::size_t Dimensions>
   PairSingManager<PTYPE, Dimensions>::PairSingManager(PTYPE& part1, PTYPE& part2, std::vector<double> h, const char* path, const char* filename) : PairSingInfo<PTYPE, Dimensions>(part1,part2), _path(path), _h(h),_filename(filename), _NLOC(10) {
    double R1 = part1->_radius;
    double R2 = part2->_radius;
    double minR = (R1 < R2)? R1: R2;
    double K  = .5*(1./R1 + 1./R2);
    double alpha = 4.;//coef for radial cutoff
    double threshold = 1./5;
    _a    = PairSingInfo<PTYPE, Dimensions>::_distC1C2 - R1 - R2;
    _cutoffdist = minR;
    //_cutoffdist = (alpha*sqrt(_a/K)<minR)? alpha*sqrt(_a/K) : minR;
    _isSingularity = _a<threshold*minR;
  }

  template<typename PTYPE>
  struct ComputesingularBC_impl<PTYPE, 2>{
#undef __FUNCT__
#define __FUNCT__ "computesingularBC_"
    static PetscErrorCode computesingularBC(PairSingManager<PTYPE, 2>& p, std::vector<PTYPE> pair){
      PetscErrorCode ierr;
      PetscFunctionBegin;
      
      for (auto it = pair.begin(); it != pair.end(); ++it){
       if(it->hasPart){
         double *pg;
         double *pmsphere;
         std::vector<double> pos(2);
         std::vector<double> posDec(2);
         std::vector<double> posRefPart(2);
         std::vector<double> u(2);	
         std::vector<double> uRefPart(2);
         int size;

         ierr  = VecGetArray(it->g, &pg);CHKERRQ(ierr);
         ierr  = VecGetArray(it->_msphere, &pmsphere);CHKERRQ(ierr);

         ierr  = VecGetSize(it->_msphere, &size);CHKERRQ(ierr);
         size /= 2;

	  // parcourt les points du bord et on interpole

         for(int k=0; k<size; k++){
	    // Coordonees du point
           pos[0] = pmsphere[2*k];
           pos[1] = pmsphere[2*k + 1];

           if(p._testbox){
             auto diff = [](double x, double y){return x-y;};
             std::transform(pos.begin(), pos.end(), p._posOrigin.begin(), posDec.begin(), diff); 
             for(int i_dim=0; i_dim<2;i_dim++)
              posRefPart[i_dim] = std::inner_product(posDec.begin(), posDec.end(), p._base[i_dim].begin(), 0.);

            uRefPart[0] = ux_sing_normalMvt2D(posRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL  ) + ux_sing_TanglMvt2D(posRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL  );
            uRefPart[1] = uz_sing_normalMvt2D(posRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL  ) + uz_sing_TangMvt2D(posRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL  );

            std:cout << "computeBC "  << posRefPart[0] << " " << uRefPart[0] << " " << posRefPart[1] << " " << uRefPart[1] << std::endl;


            u.clear();
            u.resize(2);
            for(int l=0; l<2;l++){
              u[0] += p._base[l][0]*uRefPart[l];
              u[1] += p._base[l][1]*uRefPart[l];
            }

	      //add rigid motion - singularity to simple layer
	      pg[k]      += u[0]; // -ctx->_particles[i]->_velocity[0] + ctx->_particles[i]->_angular_velocity[2]*r[1];
	      pg[k+size] += u[1]; // -ctx->_particles[i]->_velocity[1] - ctx->_particles[i]->_angular_velocity[2]*r[0]   ;

	    }
	  }
	  ierr = VecRestoreArray(it->g, &pg);CHKERRQ(ierr);
	  ierr = VecRestoreArray(it->_msphere, &pmsphere);CHKERRQ(ierr);
	}
}
PetscFunctionReturn(0);
}
};

  template<typename PTYPE>
struct ComputesingularBC_impl<PTYPE, 3>{
#undef __FUNCT__
#define __FUNCT__ "computesingularBC_"
  static PetscErrorCode computesingularBC(PairSingManager<PTYPE, 3>& p, std::vector<PTYPE> pair){
    PetscErrorCode ierr;
    PetscFunctionBegin;

    for (auto it = pair.begin(); it != pair.end(); ++it){
     if(it->hasPart){
       double *pg;
       double *pmsphere;
       std::vector<double> pos(3);
       std::vector<double> posDec(3);
       std::vector<double> posRefPart(3);
       std::vector<double> u(3);	
       std::vector<double> uRefPart(3);
       int size;

       ierr  = VecGetArray(it->g, &pg);CHKERRQ(ierr);
       ierr  = VecGetArray(it->_msphere, &pmsphere);CHKERRQ(ierr);

       ierr  = VecGetSize(it->_msphere, &size);CHKERRQ(ierr);
       size /= 3;

	  // parcourt les points du bord et on interpole

       for(int k=0; k<size; k++){
	    // Coordonees du point
         pos[0] = pmsphere[2*k];
         pos[1] = pmsphere[2*k + 1];
         pos[2] = pmsphere[2*k + 1];

         if(p._testbox){
           auto diff = [](double x, double y){return x-y;};
           std::transform(pos.begin(), pos.end(), p._posOrigin.begin(), posDec.begin(), diff); 
           for(int i_dim=0; i_dim<3;i_dim++)
            posRefPart[i_dim] = std::inner_product(posDec.begin(), posDec.end(), p._base[i_dim].begin(), 0.);

          uRefPart[0] = ux_sing_normalMvt2D(posRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL  ) + ux_sing_TangMvt2D(posRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL  );
          uRefPart[1] = uy_sing_normalMvt2D(posRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL  ) + uy_sing_TangMvt2D(posRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL  );
          uRefPart[2] = uz_sing_normalMvt3D(posRefPart, p._H1, p._H2, p._a, p._UN, p._param, p._param, NULL  ) + uz_sing_TangMvt3D(posRefPart, p._H1, p._H2, p._a, p._UT, p._param, p._param, NULL  );

          u.clear();
          u.resize(2);
          for(int l=0; l<3;l++){
            u[0] += p._base[l][0]*uRefPart[l];
            u[1] += p._base[l][1]*uRefPart[l];
            u[2] += p._base[l][2]*uRefPart[l];
          }

	      //add rigid motion - singularity to simple layer
	      pg[k]        += u[0]; // -ctx->_particles[i]->_velocity[0] + ctx->_particles[i]->_angular_velocity[2]*r[1];
	      pg[k+ size]  += u[1]; // -ctx->_particles[i]->_velocity[1] - ctx->_particles[i]->_angular_velocity[2]*r[0]   ;
	      pg[k+2*size] += u[2];		
	    }
	  }
	  ierr = VecRestoreArray(it->g, &pg);CHKERRQ(ierr);
	  ierr = VecRestoreArray(it->_msphere, &pmsphere);CHKERRQ(ierr);
	}
}
PetscFunctionReturn(0);
}
};

#undef __FUNCT__
#define __FUNCT__ "computesingularForces_"
 template<typename PTYPE, std::size_t Dimensions>
PetscErrorCode computesingularForces(PairSingManager<PTYPE, 2>& p, std::vector<PTYPE> pair, DtoNContext* ctx){
  PetscErrorCode ierr;
  PetscFunctionBegin;

  std::vector<double> F;
  std::vector<double> forces_sing;
  p.getForces(F);


  forces_sing.resize(Dimensions);
  for(int l=0; l<Dimensions;l++){
   for(int m=0; l<Dimensions;l++)
     forces_sing[m] = p._base[l][m]*F[l];

   ctx->forces_sing[l]     += forces_sing[l];
   ctx->forces_sing[Dimensions+l] -= forces_sing[l];
 }

 PetscFunctionReturn(0);
}

  template<class PTYPE, std::size_t Dimensions>
void PairSingManager<PTYPE, Dimensions>::initializebox(DMDALocalInfo& infoU)
{
  std::vector<PetscInt>    indBeginVec;
  std::vector<PetscInt>    indEndVec;

  indBeginVec.push_back(infoU.xs);
  indBeginVec.push_back(infoU.ys);
  if(Dimensions==3) indBeginVec.push_back(infoU.zs);

  indEndVec.push_back(infoU.xs + infoU.xm);
  indEndVec.push_back(infoU.ys + infoU.ym);
  if(Dimensions==3) indEndVec.push_back(infoU.zs + infoU.zm);

    // For the not periodic boundary condition.
  if(indEndVec[0] == infoU.mx && infoU.bx != DM_BOUNDARY_PERIODIC)             indEndVec[0] --;
  if(indEndVec[1] == infoU.my && infoU.by != DM_BOUNDARY_PERIODIC)             indEndVec[1] --;
  if(indEndVec[2] == infoU.mz && infoU.bz != DM_BOUNDARY_PERIODIC && Dimensions == 3) indEndVec[2] --;

  _testbox = true;

  for(int i=0; i<Dimensions; i++)
  {
   _nMin.push_back( min(floor((PairSingInfo<PTYPE, Dimensions>::_center_p1[i]-_cutoffdist)/_h[i]),floor((PairSingInfo<PTYPE, Dimensions>::_center_p2[i]-_cutoffdist)/_h[i])) );
   _nMax.push_back( max(ceil ((PairSingInfo<PTYPE, Dimensions>::_center_p1[i]+_cutoffdist)/_h[i]),ceil ((PairSingInfo<PTYPE, Dimensions>::_center_p2[i]+_cutoffdist)/_h[i])) );
   if(_nMin[i]>=indEndVec[i]||_nMax[i]<indBeginVec[i])
     _testbox = false;
   else{
     _nMin[i] = (_nMin[i]>indBeginVec[i])?_nMin[i]:indBeginVec[i];
     _nMax[i] = (_nMax[i]<indEndVec[i])  ?_nMax[i]:indEndVec[i];
   }	  
 }
}

  template<class PTYPE,std::size_t Dimensions>
void PairSingManager<PTYPE, Dimensions>::getForces(std::vector<double>& F){
  ForcesManager<Dimensions> Forces(PairSingInfo<PTYPE, Dimensions>::_H1, PairSingInfo<PTYPE, Dimensions>::_H2, _a, _UN, _UT, _param, _param, _cutoffdist);
  Forces.getForces(F);
}	 

  template<class PTYPE,std::size_t Dimensions>
void PairSingManager<PTYPE, Dimensions>::printForces(){
  ForcesManager<Dimensions> Forces(PairSingInfo<PTYPE, Dimensions>::_H1, PairSingInfo<PTYPE, Dimensions>::_H2, _a, _UN, _UT, _param, _param, _cutoffdist);
  Forces.printForces();
}

  template<class PTYPE,std::size_t Dimensions>
void PairSingManager<PTYPE, Dimensions>::clearForces(ForcesManager<Dimensions>& Forces){
  Forces.clearForces();
}

  template<class PTYPE,std::size_t Dimensions>
void PairSingManager<PTYPE, Dimensions>::printForcesNormal(){
  ForcesManager<Dimensions> Forces(PairSingInfo<PTYPE, Dimensions>::_H1, PairSingInfo<PTYPE, Dimensions>::_H2, _a, _UN, _UT, _param, _param, _cutoffdist);
  Forces.printForcesNormal();
}

  template<class PTYPE,std::size_t Dimensions>
void PairSingManager<PTYPE, Dimensions>::printForcesTang(){
  ForcesManager<Dimensions> Forces(PairSingInfo<PTYPE, Dimensions>::_H1, PairSingInfo<PTYPE, Dimensions>::_H2, _a, _UN, _UT, _param, _param, _cutoffdist);
  Forces.printForcesTang();
}

  // template<class PTYPE,std::size_t Dimensions>
  // void PairSingManager<PTYPE, Dimensions>::printNormalDIV(){
  //   std::cout << "Normal divergence value : " << DIVCART_normalMvt3D(PairSingInfo<PTYPE, Dimensions>::_H1, PairSingInfo<PTYPE, Dimensions>::_H2, _a, _UN, _UT, _param, _param, _cutoffdist) << std::endl;
  // }


  // template<class PTYPE,std::size_t Dimensions>
  // void PairSingManager<PTYPE, Dimensions>::printTangDIV(){
  //   std::cout << "Tangential divergence value : " << DIVCART_tangMvt3D(PairSingInfo<PTYPE, Dimensions>::_H1, PairSingInfo<PTYPE, Dimensions>::_H2, _a, _UN, _UT, _param, _param, _cutoffdist) << std::endl;
  //     }

  // template<class PTYPE,std::size_t Dimensions>
  // void PairSingManager<PTYPE, Dimensions>::printDIV(){
  //   std::cout << "Divergence value : " << DIVCART_normalMvt3D(PairSingInfo<PTYPE, Dimensions>::_H1, PairSingInfo<PTYPE, Dimensions>::_H2, _a, _UN, _UT, _param, _param, _cutoffdist) + DIVCART_tangMvt3D(PairSingInfo<PTYPE, Dimensions>::_H1, PairSingInfo<PTYPE, Dimensions>::_H2, _a, _UN, _UT, _param, _param, _cutoffdist) << std::endl;
  //     }
};
#endif
