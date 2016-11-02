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
#ifndef _LEVELSET_H_
#define _LEVELSET_H_

#include "mwlib/field.h"
#include "mwlib/refcounted.h"
#include "continuum-utils.h"

#include <boost/property_tree/ptree.hpp>

#define USE_HIGH_ORDER_LEVELSET_METHODS
//#define USE_CONSTRAINT_LEVELSET_METHOD


template<class T>
void AddScaled(double afactor, Array3d<T> a, double bfactor, ConstArray3d<T> b); /// multithreaded


/**
@brief Represents the distancemap function used in the levelset method. 

At each point in space it gives approximately the distance to its zero-level. The Zero-level defines a d-1 dimensional surface in d dimensional space. This class also comprises various mathematical operations and helpers.
*/
struct Levelset : public RefCounted
{
  Array3d<float> phi; /// a signed distance function
  double max_val;     /// maximal value found in phi

  Levelset() : max_val(-1) {}
  
  void init(const LatticeDataQuad3d &grid, int ndims, int border = 0)
  {
    phi.initFromBox(ExtendForDim(grid.Box(), ndims, border));
    phi.setBox(grid.Box());
    max_val = -1.; // invalid
  }
  
  void addScaled(double s, const Levelset &u, double s_self = 1.)
  {
    //phi.addScaled<>(s, u.phi, s_self);
    AddScaled(s_self, phi, s, u.phi);
  }

  void initCloned(const Levelset &other)
  {
    phi.initDeepCopy<>(other.phi);
    max_val = other.max_val;
  }

  void Draw(Image &img, const LatticeDataQuad3d &ld, const DrawArrayOpts &opts, const Array3d<char>* active = NULL);

  /**
   @brief  Returns the application of the Heaviside function on the distancemap phi. 
   
    The result is positive inside the volume described by the levelset and zero outside. The transition region between the two has the width given by the width parameter.
  */
  void ToHeaviside(Array3d<float> theta, const BBox3 &bbox, double width) const;
  static float ToHeaviside(float phi, float width);
  
  void initLike(const Levelset &other)
  {
    phi.initFromBox(other.phi.getBox());
    max_val = -1;
  }

  void clear()
  {
    phi.clear();
  }
  
  void setId(int id_) {  }
  int getId() const { return 0; }
};

void AddScaled(double f0, Levelset& l0, double f1, const Levelset &l1);

/*------------------------------------------------------
------------------------------------------------------*/

/**
 * @brief Mathematical operations used by my time stepper classes
 */
class LevelsetOps
{
  Array3dOps<float> ops;

public:
  typedef Levelset state_type;
  
  LevelsetOps() : ops() {}
  LevelsetOps(const DomainDecomposition &mtboxes_, const BBox3 &bbox_, int dim_, int border_) : ops(mtboxes_, bbox_, dim_, border_) {}
  void init(const DomainDecomposition &mtboxes_, const BBox3 &bbox_, int dim_, int border_) { ops.init(mtboxes_, bbox_, dim_, border_); }

  void init(Levelset& u, bool clear = true) const
  {
   ops.init(u.phi, clear);
   u.max_val = -1.;
  }

  void addScaled(double fa, Levelset &u, double fb, const Levelset &v) const
  {
    ops.addScaled(fa, u.phi, fb, v.phi);
  }
  void initFrom(Levelset &u, const Levelset &v, ConsMode mode) const
  {
    ops.initFrom(u.phi, v.phi, mode);
    u.max_val = v.max_val;
  }
};

/*------------------------------------------------------
------------------------------------------------------*/

/**
  @brief Does the magic to turn a distorted distancemap back into a proper one where each point in space gives 
 */
boost::property_tree::ptree reinit_levelset(Levelset& levelset, const LatticeDataQuad3d& ld, int dim, const DomainDecomposition &mtboxes, int width, const boost::property_tree::ptree& params = boost::property_tree::ptree());
boost::property_tree::ptree reinit_levelset(Levelset& levelset, const LatticeDataQuad3d& ld, int dim, int width, const boost::property_tree::ptree& params = boost::property_tree::ptree());

/**
 @brief How badly distorted is the levelset function?
 */
double calc_levelset_badness(const Levelset &levelset, const LatticeDataQuad3d &ld, int dim, const DomainDecomposition &mtboxes);
double calc_levelset_badness(const Levelset &levelset, const LatticeDataQuad3d &ld, int dim);
/**
 @brief The badness beyond which the levelset function should be "redistanced" i.e. reinit_levelset be called.
*/
double levelset_ok_badness();

#endif