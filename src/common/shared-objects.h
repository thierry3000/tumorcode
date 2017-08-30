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
#ifndef SHARED_OBJECTS_H
#define SHARED_OBJECTS_H

#include "common.h"
#include "hdfio.h"

#include "mwlib/lattice-data.h"
#include "mwlib/listgraph.h"
#include "mwlib/random.h"
#include <memory>

#include <boost/function.hpp>
#include <boost/array.hpp>
#include <boost/filesystem/convenience.hpp>

class  Vessel;
class  VesselNode;
class  VesselList3d;
struct HemodynamicBounds;
class  ContinuumTumor;
namespace polymorphic_latticedata {
  class LatticeData;
}
class ContinuumGrid;
class DomainDecomposition;
class Levelset;

enum TissueId {
  TCS = 0,
  DEAD = 1,
  TISSUE = 2,
};

struct TissuePhases
{
  TissuePhases() : count(0) {}
  TissuePhases(int count_, const BBox3 &bbox_) : count(count_)
  {
    for (int i=0; i<count_; ++i)
      phase_arrays[i].initFromBox(bbox_);
  }

  Float3 operator()(const Int3 &p) const
  {
    Float3 result(count>0 ? phase_arrays[0](p) : 0,
                  count>1 ? phase_arrays[1](p) : 0,
                  count>2 ? phase_arrays[2](p) : 0);
    return result;
  }
  
  int count;
  Array3df phase_arrays[3];
};
void SetupTissuePhases(TissuePhases &phases, const ContinuumGrid &grid, DomainDecomposition &mtboxes, boost::optional<h5cpp::Group> tumorgroup);
enum TumorTypes
{
  FAKE = 0,
  BULKTISSUE = 1,
  NONE = 3,
};
TumorTypes determineTumorType(h5cpp::Group *tumorgroup);
/*-----------------------------------------------------------
-----------------------------------------------------------*/

#if 0
// not currently used
void InitCoarseLD( LatticeDataQuad3d &ld_coarse, // result size includes border
                   const LatticeDataQuad3d &ld_fine,
                   const int border_fine,
                   const int border_coarse,
                   const int subdiv );
#endif


//-----------------------------------
double EstimateTumorRadius(const ContinuumGrid& grid, const DomainDecomposition& mtboxes, const Levelset& ls);

//-----------------------------------
enum
{
  TUM_NUM_REGIONS = 4,
  TUMBIT_REGION_CENTER = 1,
  TUMBIT_REGION_PERIPHERY = 2,
  TUMBIT_REGION_BORDER = 4,
  TUMBIT_REGION_OUTSIDE = 8
};

template<class CheckInside>
void ClassifyVesselsByRegion( const VesselList3d* list, const CheckInside &idfield, DynArray<my::Bitfield<uchar> > &vessid );


class CheckInside_ContTumor3d
{
  LatticeDataQuad3d ld;
  BBox3 bbox;
  Array3d<uchar> f;
public:
  CheckInside_ContTumor3d();
  void Init(const ContinuumTumor &tl, float boundary_radius, bool bPeripherySeperate);
  my::Bitfield<uchar>  operator()(const Float3 &p) const;
};

//-----------------------------------
//-----------------------------------

float GetInitialThickness( float rad );
bool FindHemodynamicBounds( const VesselList3d *vl, HemodynamicBounds &hb, bool bAutoFill );

/* ReadVesselList3d - parameters
 * "scale override"
 * "hematocrit"
 * "calcflow"
 * "filter"
 */
std::auto_ptr<VesselList3d> ReadVesselList3d(h5cpp::Group vesselgroup, const ptree &params);
/* WriteVesselList3dVolume - parameters
 * w_all
 * w_pressure
 * w_flow
 * -----------------
 * this function writes nodes/, and edges/ subgroups and adds datasets into them
 */
void WriteVesselList3d(const VesselList3d& vl, h5cpp::Group vesselgroup, const ptree& params = ptree());

void CenterVesselListAndSetupFieldLattice(VesselList3d &vl, int dim, float spacing, LatticeDataQuad3d &field_ld);
void SetupFieldLattice(const FloatBBox3 &wbbox, int dim, float spacing, float safety_spacing, LatticeDataQuad3d &ld);
void printPtree(ptree const &pt);

int FindPositionOnVessel(const polymorphic_latticedata::LatticeData &ld, const Vessel* v, const Int3 &p );
int FindDistanceToJunction( const Vessel* vstart, int vpos, const VesselNode* vcstart, int maxdist );


//-----------------------------------
//-----------------------------------

//picks a random
inline int RandomPick( Random& random, int cnt, double *prob )
{
  double r = random.Get01();
  double pp=0.0;
  for( int k=0; k<cnt; ++k )
  {
    pp += prob[k];
    if( r<pp ) return k;
  }
  myAssert( std::abs(pp-1.0)<1.0e-6 );
  return cnt-1;
}


inline double NormalizeSumNorm( double* x, int cnt )
{
  double n = 0.0;
  for( int i=0; i<cnt; ++i ) n+=std::abs(x[i]);
  if( n>0.0 ) {
    const double nn = 1.0/n;
    for( int i=0; i<cnt; ++i ) x[i]*=nn;
  }
  return n;
}


template<class T>
inline void RandomPermute( Random &random, T* x, int cnt )
{
  int i,j;
  for( i=0; i<cnt; ++i )
  {
    j = random.Get( cnt );
    std::swap( x[i], x[j] );
  }
}


namespace FieldInterpolate
{

template<class T>
struct OutOfDomainConstT
{
  T m_value;
  OutOfDomainConstT(const T &value_) : m_value(value_) {}
};

struct OutOfDomainExtrapolateT
{
};


template<class T>
inline OutOfDomainConstT<T> Const(const T &value) { return OutOfDomainConstT<T>(value); }
inline OutOfDomainExtrapolateT Extrapolate() { return OutOfDomainExtrapolateT(); }

template<class T, class OOD>
T ValueAveraged(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Float3& pos);
template<class T, class OOD>
T Value(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Float3& pos);

template<class T, class OOD>
Vec<T,3> Gradient(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Float3 &pos);

template<class T, class OOD>
T Value(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Int3& ip);

}

//-----------------------------------
//-----------------------------------

/** @brief
 * As the name suggests, this guy will create samples.
 * You need to provide:
 * 1) world coodinates of the vesselsegment under consideration
 * 2) radius 
 * 3) how many sample points you want
 * 
 * There are different types. You can sample on
 *    SURFACE 
 *    LINE  
 *    VOLUME
 */
class CylinderNetworkSampler
{
public:
  typedef my::eqpair<Float3> WorldEdge;
  struct Sample {
    /*
     * fpos is a fractional coordinate,
     * component 0 and 1 are cartesian coordinates perpendicular to the cylinder axis,
     * and comp. 2 is the axial coordinate
     * all components \in [0,1]
     *
     * wpos is the world coordinate
     */
    Float3 wpos, fpos;
  };

  /*
   * this works as a stream. once the arrays are set, 
   * you pull samples with these functions until it returns 0
   * each iteration you access the current sample buffer with GetSample()
   * Reset() restarts the stream from the beginning
   */
  int GenerateLineSamples(); // clears the buffer, generates new samples, returns number of samples in the buffer
  int GenerateSurfaceSamples();
  int GenerateVolumeSamples();
  const Sample& GetSample(int i) const { return buffer[i]; } // access buffer
  void Restart(); // start completely from the beginning

  float weight; // length, area or volume, depending on the generating function
  float weight_per_volume;
  int index; // index of the edge from which it was generated

  CylinderNetworkSampler() {}

  void Set(const WorldEdge *edges_, const float *radii_, int count_);
  void Set(const Float3 &posa, const Float3 &posb, float rad);

  void Init(double spacing, const ptree &params);

  int samples_per_cell;
  uint seed;

private:
  DynArray<Sample> buffer;

  const WorldEdge *edges;
  const float *radii;
  int edge_count;
  my::eqpair<Float3> single_edge;
  float single_rad;

  double spacing;
  Random rnd;
};

void Set(CylinderNetworkSampler &sampler, const Vessel* v, const polymorphic_latticedata::LatticeData &ld);

/*------------------------------------------------------
------------------------------------------------------*/

typedef DynArray< DynArray<const Vessel*> > VesselsInBoxes;
void SortVesselsIntoMtBoxGrid(const LatticeDataQuad3d &ld, const VesselList3d &vl, int boxborder, const DynArray<BBox3> &boxes, VesselsInBoxes &grid);
void SortVesselsIntoMtBoxGrid(const LatticeDataQuad3d &ld, const VesselList3d &vl, int boxborder, const DomainDecomposition &mtboxes, VesselsInBoxes &grid);

/*------------------------------------------------------
------------------------------------------------------*/

/* Fills a box with the volume fraction occupied by vessels. It is NOT thread safe!
 */
class VesselVolumeGenerator
{
  const VesselList3d &vl;
  const polymorphic_latticedata::LatticeData& ld;
  CylinderNetworkSampler sampler;
  bool use_smooth_delta;
  float upper_volume_bound;
  const LatticeDataQuad3d &fieldld;
  int dim;

public:
  VesselVolumeGenerator(const VesselList3d &vl_, const LatticeDataQuad3d &fieldld_, int dim_, const ptree &params);
  void Fill(const BBox3 &bbox, Array3df vessel_volume, const DynArray<const Vessel*> vessels, int &box_sample_cnt);
};

/*
 * parameters:
 *  smooth : bool (true) - use smothended delta function with 4 cell support width
 *  cut_at : float (1.) - limit resulting values to this upper bound. The limiting is applied to the values in the supplied vessel_volume array, use < 0 to disable
 *  samples_per_cell : int (2) - the number of samples (approximately) taken within one cell volume of the array lattice
 * return:
 *  num_samples : int - total number of samples taken
 */
ptree AddVesselVolumeFraction(const VesselList3d &vl, const LatticeDataQuad3d &field_ld, int dim, Array3d<float> vessel_volume, const ptree &params = ptree());



/*------------------------------------------------------
------------------------------------------------------*/

/*
 * T.F. 20.10.16. this was needed for the legacy cpp side.
 * My standard program arguments handling routine.
 * Arguments are -p and -pa
 *  -p is either a filename of a .info file or "show"
 *  -pa is as .info formated data as string
 * -returns empty ptree if program should be stopped
 * -throws if there is a error parsing, or parameter which is not in the defaults
 */
boost::optional<ptree> HandleSimulationProgramArguments(const ptree &default_params, int argc, char **argv, const ptree &vararg_params = ptree());
/*------------------------------------------------------
------------------------------------------------------*/

/** @brief tests if ReadVesselList3d was successful
 */
void isVesselListGood(VesselList3d &vl);

#endif
