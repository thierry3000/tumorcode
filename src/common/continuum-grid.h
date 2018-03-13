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
#ifndef CONTINUUM_GRID_H
#define CONTINUUM_GRID_H

#include "common.h"
#include "mwlib/lattice-data.h"

#include <boost/array.hpp>
#include <boost/range/iterator_range.hpp>

//-----------------------------------
//-----------------------------------

struct LatticeIndexRanges
{
  static LatticeIndexRanges FromCellRange(const BBox3 &b, int ndim = 3);
  BBox3 cells, faces[3], verts;
};

/**
@brief For a regular grid we can consider the faces of and vertices of grid cells and
recognize that these faces and vertices are themselves arranged on a lattice. 

These functions are used to derive corresponding Lattices from the base computational grid.
*/
BBox3 CellRangeToFaceRange(const BBox3 &b, int ndim, int axis);
BBox3 CellRangeToVertexRange(const BBox3 &b, int ndim);
LatticeDataQuad3d CellLdToFaceLd(const LatticeDataQuad3d &ld, int ndim, int axis);
LatticeDataQuad3d CellLdToVertexLd(const LatticeDataQuad3d &ld, int ndim);
Bool3 CellCenteringForDim(int ndim);

//-----------------------------------
//-----------------------------------
/*
 * makes a list of blocks which form a grid and their union is bb
 * bs -> hint for block size
 */
DynArray<BBox3> MakeMtBoxGrid(const BBox3 &bb, const Int3 &bs);
DynArray<BBox3> MakeMtBoxGrid(const BBox3 &bb); // smallish boxes 32^dim
DynArray<BBox3> MakeMtBoxGridLarge(const BBox3 &bb, int max_size = 128); // fairly large boxes

class DomainDecomposition
{
public:
  struct ThreadBox : public BBox3
  {
    ThreadBox(const BBox3 &bbox, int global_index_) : BBox3(bbox), global_index(global_index_) {}
    int global_index;
  };
  DynArray<ThreadBox> boxes;
  boost::array<DynArray<ThreadBox>, 32> by_thread;
  ThreadBox largestBox;
public:
  DomainDecomposition(const DynArray<BBox3> &boxes);
  DomainDecomposition();
  void init(const DynArray<BBox3> &boxes); //, boost::optional<DynArray<double>&> weights);

  typedef boost::iterator_range<DynArray<ThreadBox>::const_iterator> range_type;
  range_type getCurrentThreadRange() const;
  const ThreadBox getBoxWithLargestVolume() const;

  std::size_t size() const { return boxes.size(); }
  const ThreadBox& operator[](int i) const { return boxes[i]; }

  void insert(int thread_num, const BBox3 &bbox); // manually insert a box
  std::ostream& print(std::ostream &os) const;
};



//-----------------------------------
//-----------------------------------

/*
 * All the commonly needed info to work on staggered grids with variable
 * dimension. Also includes a subdivision into smaller boxes (mtboxes)
 * for multithreading.
 */
struct ContinuumGrid
{
  ContinuumGrid(const LatticeDataQuad3d &ld_, int dim_ = -1);
  ContinuumGrid();
  void init(const LatticeDataQuad3d &ld_, int dim_ = -1);
  
  LatticeDataQuad3d ld;
  LatticeIndexRanges ir; // index ranges for cell, vertex and face variables
  int dim;
  BBox3 Box() const { return ld.Box(); }
  Int3 Size() const { return ld.Size(); }
  int Dim() const { return dim; }
  float Spacing() const { return ld.Scale(); }
};


//-----------------------------------
//-----------------------------------
/**
 * @brief Collect three arrays that represent function values 
 * on the faces of grid cells into a single structure.
 * 
 * Note: 
 * Face index 0 corresponds to the left/bottom face of the first grid cell.
 * Face index 1 corresponds to the face between the first and the second cell. 
 * And so on. See info on staggered grids.
 */
class FaceVarArrays : public boost::array<Array3d<float>, 3>
{
public:
  typedef boost::array<Array3d<float>, 3> Super;
  FaceVarArrays() {}
  // the Box bb is cell centered! Space for the outer left,right,top,bottom,etc. boundary faces is allocated
  FaceVarArrays(const BBox3 &bb, int dim) { init(bb, dim); }
  FaceVarArrays(const LatticeIndexRanges &ir, int dim) { init(ir, dim); }
  void init(const BBox3 &bb, int dim);
  void init(const LatticeIndexRanges &ir, int dim);

  float operator()(int axis, const Int3 &p) const { return Super::operator[](axis)(p); }
  float& operator()(int axis, const Int3 &p) { return Super::operator[](axis)(p); }
};

//typedef boost::array<Array3d<float>, 3> FaceVarArrays;
//void InitFaceVar(const BBox3 &bb, int dim, FaceVarArrays &arr);

//-----------------------------------
//-----------------------------------

enum MakeArraySetBoxStyle {
  MAKE_ARRAY3D_BOX_INNER,
  MAKE_ARRAY3D_BOX_OUTER
};

/** @brief
 * Allocate a Array3d for a bounding box but actually allocate more space which
 * comprises a ghost boundary around bb of width border. If set_box_style is
 * MAKE_ARRAY3D_BOX_OUTER then the arrays box is set to include the boundary.
 */
template<class T>
inline Array3d<T> MakeArray3dWithBorder(const BBox3 &bb, int dim, int border, MakeArraySetBoxStyle set_box_style = MAKE_ARRAY3D_BOX_INNER)
{
  Array3d<T> a(ExtendForDim(bb, dim, border));
  if (set_box_style == MAKE_ARRAY3D_BOX_INNER)
    a.setBox(bb);
  return a;
}


//-----------------------------------
//-----------------------------------


/** @brief
 * multithreaded linear combination function
 * a = afactor * a + bfactor * b
 */
template<class T>
void AddScaled(double afactor, Array3d<T> a, double bfactor, ConstArray3d<T> b);


/*------------------------------------------------------
------------------------------------------------------*/

/** @brief 
 * Rules:
 *  + boxindex may indicate the threads sub-box, otherwise it is -1 and indicates the full domain
 *  + bbox corresponds to the boxindex and the indexrange of the thread box is a subset of box,
 *    but bbox may be larger to incorporate boundary items
 *  + it is allowed to return mixed vertex/face/cell variables when the caller and callee
 *    agree on a common policy
 *  + The Array3d instances in values should be pre-allocated and will filled by the callee
 */
template<class T, int COUNT_>
class VirtualGridFunctions
{
public:
  enum { COUNT = COUNT_ };
  typedef boost::array<Array3d<T>, COUNT> List;
  virtual ~VirtualGridFunctions() {}
  virtual void GetValues(int boxindex, const BBox3 &bbox, List &values) const = 0;
};

/*------------------------------------------------------
------------------------------------------------------*/

template<class T, int COUNT_>
class VirtualArray3dGridFunctions_ : public VirtualGridFunctions<T, COUNT_>
{
public:
  typedef VirtualGridFunctions<T, COUNT_> Super;
  typename Super::List fields;
  VirtualArray3dGridFunctions_() {}
  Array3d<T>& operator[](int i) { return fields[i]; }
  const ConstArray3d<T>& operator[](int i) const { return fields[i]; }
  void GetValues(int boxindex, const BBox3 &bbox, typename Super::List &values) const
  {
    for (int i=0; i<fields.size(); ++i)
      values[i][bbox].fill(fields[i][bbox]);
  }
};

template<class T, int COUNT_>
class VirtualArray3dGridFunctions : public VirtualArray3dGridFunctions_<T, COUNT_>
{
};

template<>
class VirtualArray3dGridFunctions<float, 1> : public VirtualArray3dGridFunctions_<float, 1>
{
public:
  VirtualArray3dGridFunctions<float, 1>(Array3d<float> a1) { fields[0] = a1; }
};

/*------------------------------------------------------
------------------------------------------------------*/
#endif
