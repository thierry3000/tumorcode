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
#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "common.h"
#include <boost/scoped_ptr.hpp>
#include <boost/smart_ptr.hpp>

//#define MULTIGRID_DEBUG

#ifdef MULTIGRID_DEBUG
#define MGDB(a) a
#else
#define MGDB(a)
#endif



template<class Operator_, class Vector_, class MgOperations_>
class Multigrid
{
  typedef Operator_ Operator;
  typedef Vector_ Vector;
  typedef MgOperations_ MgOperations;
  typedef typename Operator::value_type value_type;
  typedef boost::shared_ptr<Operator> OperatorPtr;
  typedef boost::scoped_ptr<Vector> VectorPtr;
  /* MgOperations provides
   * - prolongation
   * - restriction
   * - bottom solve
   * - smoothing
   * Prolongation turns coarse grid vectors into fine grid vectors.
   * Restriction turns fine grid operators and vectors into coarse grid version.
   * The bottom solve takes a coarse grid system and solves it "exactly". The bottom
   * operator never changes, so it can store factorized matrices, etc. for reuse.
   */
  
public:
  int numSmootherIterations;
  int numVCycleIterations;
  T minResidual;
  T lastSolutionResidual;
  MgOperations mgops;
  DynArray<boost::shared_ptr<Operator> > operator_hierarchy;

  Multigrid(OperatorPtr initial_operator)
  {
    MGDB(printf("*generating grids*\n"));

    this->numSmootherIterations = 5;
    this->numVCycleIterations = 1;
    this->minResidual = 1.e-18;
    this->lastSolutionResidual = 0.;

    GenerateGrids(initial_operator);

    MGDB(printf("*done generating grids*\n"));
  }

  void GenerateGrids(const OperatorPtr initial_operator)
  {    
    operator_hierarchy.push_back(initial_operator);
    while(true)
    {
      const int level = operator_hierarchy.size();
      const OperatorPtr op = operator_hierarchy.back();
      const OperatorPtr coarse = mgops.restrict(*op, level);
      if (!coarse.get()) break;
      operator_hierarchy.push_back(coarse);
    }
  }

  VectorPtr AllocateVector(int level)
  {
    return mgops.allocateVector(*operator_hierarchy[level], level);
  }

  void BottomSolve(Vector &lhs, const Vector &rhs, int level)
  {
    MGDB(printf("*direct solve level %i*\n",level));
    myAssert(level == operator_hierarchy.size()-1);

    mgops.bottomSolve(*operator_hierarchy[level], lhs, rhs, level);

    MGDB(
      double r = operator_hierarchy[level]->residualNorm(lhs, rhs);
      printf("direct solve level %i residual norm = %e\n",level,r);
    )
  }

  void Smooth(Vector &lhs, const Vector &rhs, int level)
  {
    for(int iter=0; iter<numSmootherIterations; ++iter)
      mgops.smooth(*operator_hierarchy[level], lhs, rhs, level);
  }

  void VCycle(V &lhs, const V &rhs, int level=0)
  {    
    if(level == operator_hierarchy.size()-1)
    {
      BottomSolve(lhs,rhs,level);
    }
    else
    {
      MGDB(printf("*vcycle level %i*\n",level));
      Smooth(lhs,rhs,level);

      VectorPtr res = AllocateVector(level);
      VectorPtr resc = AllocateVector(level+1);
      
      operator_hierarchy[level]->residual(*res, lhs, rhs);
      mgops.restrict(*resc, res, level); // coarse rhs computed

      res.reset();
      
      VectorPtr lhsc = AllocateVector(level+1);
      VCycle(*lhsc, *resc, level+1);   

      resc.reset();
      
      mgops.addProlongation(lhs, lhsc);

      lhsc.reset();
      
      Smooth(lhs,rhs,level);

      MGDB(
        double r = operator_hierarchy[level]->residualNorm(lhs, rhs);
        printf("vcycle level %i residual norm = %e\n",level,r);
      )
    }
  }

  void VCycles(V &lhs, const V &rhs, int level=0)
  {
    if(level < operator_hierarchy.size()-1)
    {
      for(int i=0; i<numVCycleIterations; ++i)
      {
        VCycle(lhs,rhs,level);
        if(minResidual > 0.)
        {
          double r = lastSolutionResidual = operator_hierarchy[level]->residualNorm(lhs, rhs);
          if(r < minResidual) break;
        }
      }
    }
    else
    {
      BottomSolve(lhs,rhs,level);
    }
  }

  void FMG(V &lhs, const V &rhs, int level=0)
  {
    if(level == operator_hierarchy.size()-1)
    {
      BottomSolve(lhs,rhs,level);
    }
    else
    {
      MGDB(printf("*fmg level %i*\n",level));
      VectorPtr rhsc = AllocateVector(level+1);
      mgops.restrict(*rhsc, rhs, level);
      VectorPtr lhsc = AllocateVector(level+1);

      FMG(*lhsc, *rhsc, level+1);

      rhsc.reset();

      mgops.prolongate(lhs, *lhsc);
      
      lhsc.reset();
      
      VCycles(lhs, rhs, level);

      MGDB(
        double r = operator_hierarchy[level]->residualNorm(lhs, rhs);
        printf("fmg level %i residual norm = %e\n",level,r);
      )
    }
  }
};

#undef MGDB

#endif
