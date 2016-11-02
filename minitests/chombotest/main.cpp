/**
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2016  Michael Welter and Thierry Fredrich

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <iostream>
using std::endl;
#include <fenv.h>


#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "BoxIterator.H"
#include "FABView.H"
#include "AMRIO.H"

#include "VCAMRPoissonOp2.H"
#include "AMRPoissonOp.H"
#include "BCFunc.H"
#include "BiCGStabSolver.H"
#include "CH_Timer.H"

#include "math_ext.h"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testMultiGrid" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

template<class T>
inline RefCountedPtr<T> rcp(T *p) { return RefCountedPtr<T>(p); }



static Box domain = Box(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(127,127,127)));
static Real dx = 0.0125;
static Real xshift = 0.0;


int
testMultiGrid();


int
main(int argc ,char* argv[])
{
  int except =  FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW ;
  feenableexcept(except);

  parseTestOptions( argc ,argv ) ;
  if( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int overallStatus = 0;
  int status = testMultiGrid();

  if( status == 0 )
  {
    pout() << indent << pgmname << " passed." << endl ;
  }
  else
  {
    overallStatus = 1;
    pout() << indent << pgmname << " failed with return code " << status << endl ;
  }

  xshift = 0.2;

  return overallStatus;
}


extern "C" {
  void Parabola_neum(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    switch (*dir) {
    case 0:
      a_values[0]=2*(pos[0]-xshift);
      return;
    case 1:
      a_values[0]=2*pos[1];
    return;
    case 2:
      a_values[0]=2*pos[2];
      return;
    default:
      MayDay::Error("no such dimension");
    };
  }

  void Parabola_diri(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    //a_values[0] = D_TERM((pos[0]-xshift)*(pos[0]-xshift),+pos[1]*pos[1],+pos[2]*pos[2]);
    a_values[0] = 0.;
  }

  void DirParabolaBC(FArrayBox& a_state,
                     const Box& valid,
                     const ProblemDomain& a_domain,
                     Real a_dx,
                     bool a_homogeneous)
  {

    for(int i=0; i<CH_SPACEDIM; ++i)
      {
        DiriBC(a_state,
               valid,
               dx,
               a_homogeneous,
               Parabola_diri,
               i,
               Side::Lo);
        DiriBC(a_state,
               valid,
               dx,
               a_homogeneous,
               Parabola_diri,
               i,
               Side::Hi);
      }
  }

  void NeumParabolaBC(FArrayBox& a_state,
                      const Box& valid,
                      const ProblemDomain& a_domain,
                      Real a_dx,
                      bool a_homogeneous)
  {

    for(int i=0; i<CH_SPACEDIM; ++i)
      {
        NeumBC(a_state,
               valid,
               dx,
               a_homogeneous,
               Parabola_neum,
               i,
               Side::Lo);
        NeumBC(a_state,
               valid,
               dx,
               a_homogeneous,
               Parabola_neum,
               i,
               Side::Hi);
      }
  }
}

static BCValueFunc pointFunc = Parabola_diri;

void parabola(const Box& box, int comps, FArrayBox& t)
{
  RealVect pos;
  Side::LoHiSide side;
  int dir;
  int num = 1;
  ForAllXBNN(Real,t, box, 0, comps)
    {
      num=nR;
      D_TERM(pos[0]=dx*(iR+0.5);, pos[1]=dx*(jR+0.5);, pos[2]=dx*(kR+0.5));
      pointFunc(&(pos[0]), &dir, &side, &tR);
    }EndFor;
}

void fillFluxCoeff(FArrayBox &a, int dir, Real dx, const ProblemDomain &domain)
{
  Real diam = domain.domainBox().longside() * dx;
  RealVect pos;
  const Box bb = a.box();
  ForAllXBNN(Real, a, bb, 0, 1) // type, name of the array, box, component start index, number of components
  {
    pos[0] = dx * (iR + 0.5) - diam*0.5;
    pos[1] = dx * (jR + 0.5) - diam*0.5;
    pos[2] = dx * (kR + 0.5) - diam*0.5;
    pos[dir] -= dx * 0.5;
    Real l = pos.vectorLength();
    Real f = my::smooth_heaviside(l - diam * 0.25, 3 * dx); // 0 in center, 1 outside
    aR = my::lerp(f, 10., 1.);
  }
  EndFor;
}



int
testMultiGrid()
{
  ProblemDomain regularDomain(domain);
  Vector<DisjointBoxLayout> vectGrids;
  Vector<int> refRatio;
  {
    // just one box
    Vector<Box> boxVec; boxVec.push_back(domain);
    Vector<int> procVec; procVec.push_back(0);
    DisjointBoxLayout dbl(boxVec, procVec, regularDomain);
    vectGrids.push_back(dbl);
  }
  
  pout()<<"\n single grid MultiGrid solver \n";
  // single grid solver test
  {
    Vector<LevelData<FArrayBox>*> phi, rhs;
    for (int i=0; i<vectGrids.size(); ++i)
    {
      // 1 component, 2 sites ghost boundary
      phi.push_back(new LevelData<FArrayBox>(vectGrids[i], 1, 2*IntVect::Unit));
      // 1 component, 0 ghost boudnary
      rhs.push_back(new LevelData<FArrayBox>(vectGrids[i], 1, IntVect::Zero));
      for (DataIterator dit = vectGrids[i].dataIterator(); dit.ok(); ++dit)
      {
        (*phi[i])[dit()].setVal(0.);
        (*rhs[i])[dit()].setVal(2*CH_SPACEDIM);
      }
    }

    #if 0
    AMRPoissonOpFactory opFactory;
    opFactory.define(regularDomain, vectGrids, refRatio, dx, DirParabolaBC);
    #else
    Vector<RefCountedPtr<LevelData<FArrayBox> > > vecIdentCoeff;
    Vector<RefCountedPtr<LevelData<FluxBox> > > vecFluxCoeff;
    for (int i=0; i<vectGrids.size(); ++i)
    {
      vecIdentCoeff.push_back(rcp(new LevelData<FArrayBox>(vectGrids[i], 1, IntVect::Zero)));
      vecFluxCoeff.push_back(rcp(new LevelData<FluxBox>(vectGrids[i], 1, IntVect::Zero)));
      for (DataIterator dit = vectGrids[i].dataIterator(); dit.ok(); ++dit)
      {
        (*vecIdentCoeff[i])[dit()].setVal(1.);
        for (int dir = 0; dir < SpaceDim; ++dir)
        {
          FArrayBox& arr = (*vecFluxCoeff[i])[dit()][dir];
          fillFluxCoeff(arr, dir, dx, regularDomain);
        }
      }
    }
    VCAMRPoissonOp2Factory opFactory;
    opFactory.define(regularDomain, vectGrids, refRatio, dx, DirParabolaBC, 0, vecIdentCoeff, 1., vecFluxCoeff);
    #endif

    BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
    bottomSolver.m_verbosity=0;
    bottomSolver.m_eps = 1.e-12;

    AMRMultiGrid<LevelData<FArrayBox> > solver;
    solver.define(regularDomain, opFactory, &bottomSolver, 1); // one amr level
    solver.m_numMG = 1;
    solver.m_post = 3;
    solver.m_pre = 3;
    solver.m_eps = 1.e-9;
    solver.m_verbosity = 10;
    
    solver.solve(phi, rhs, 0, 0);
    {
      Real time = 1.;
      int numLevels = 1;
      Vector<string> name; name.push_back("phi");
      WriteAMRHierarchyHDF5(string("testresult.h5"),
                            vectGrids,
                            phi,
                            name,
                            regularDomain.domainBox(),
                            dx, time, time,
                            refRatio,
                            numLevels);
    }
  }
  return 0;
}
