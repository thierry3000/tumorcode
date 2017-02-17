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
#ifndef LATTICE_DATA_POLYMORPHIC
#define LATTICE_DATA_POLYMORPHIC

#include "mwlib/lattice-data.h"
#include "hdf_wrapper.h"
#include <boost/noncopyable.hpp>

// forward declaration, their location may change
void WriteHdfLd( h5cpp::Group f, const LatticeDataQuad3d &ld );
void ReadHdfLd( h5cpp::Group f, LatticeDataQuad3d &ld );
void WriteHdfLd( h5cpp::Group f, const LatticeDataFCC &ld );
void ReadHdfLd( h5cpp::Group f, LatticeDataFCC &ld );


namespace polymorphic_latticedata
{
/*
 * polymorphic variant, for performance reasons this is a seperate class
 */

class LatticeData : boost::noncopyable
{
  protected:
    LatticeData() {}
  public:
    typedef Int3 LatticeIndexType;
    typedef int64 SiteType;
    
    virtual ~LatticeData() {}
    virtual boost::shared_ptr<LatticeData> Clone() const = 0;
    virtual void Init(const BBox3 &bb, float scale) = 0;

    virtual float Scale() const = 0;
    virtual void SetScale(float s) = 0;
    
    virtual BBox3 Box() const = 0;

    virtual void SetOriginPosition(const Float3 &pos) = 0;
    virtual Float3 GetOriginPosition() const = 0;
    virtual FloatBBox3 GetWorldBox() const = 0;
    virtual void SetCellCentering(const Vec<bool, 3> &cc) = 0;

    virtual int NbCount() const = 0;
    virtual SiteType NbSite(SiteType site, int dir) const  = 0;
    virtual Int3 NbLattice(const Int3 &p, int dir) const = 0;
    virtual Float3 LatticeToWorld( const Int3 &p) const = 0;
    virtual Int3 WorldToLattice(const Float3 &p) const =  0;
    
    virtual SiteType LatticeToSite(const Int3 &p) const = 0;
    virtual const Int3 SiteToLattice(SiteType site_) const = 0;
    virtual bool IsInsideLattice(const Int3 &p) const = 0;

    virtual AxisDirLen GetAxisDirLen(const Int3 &p1, const Int3 &p2 ) const = 0;
    virtual Int3 GetLatticeIndexOnRefinedGrid(const Int3 &pos, int refinement_subdivision) const = 0;

    virtual void print(std::ostream &os) const = 0;

    /*
     * ldtype = quad or fcc
     */
    static boost::shared_ptr<LatticeData> Make(const char* ldtype, const BBox3 &bb, float scale);

    // hdf 5 support
    static boost::shared_ptr<LatticeData> ReadHdf(h5cpp::Group g);
    virtual void WriteHdf(h5cpp::Group g) const = 0;
};

template<class Ld>
struct fwd_cell_centering {
  static void set(Ld& ld, const Bool3 &b) { throw std::runtime_error("not impl."); } 
};

template<>
struct fwd_cell_centering<LatticeDataQuad3d> {
  static void set(LatticeDataQuad3d& ld, const Bool3 &b) { ld.SetCellCentering(b); }
};


template<class LD>
Int3 WorldToLatticeWrapper(const LD &ld, const Float3 &p);


template<class Ld>
class Derived : public LatticeData
{
  Ld ld;
public:
  // trivially forward all calls
  Derived(const Ld &ld) : ld(ld) {}
  Derived(const BBox3 &bb, float scale) { ld.Init(bb, scale); }
  Derived(const Derived &other) : ld(other.ld) {}
  Ld& get() { return ld; }
  Ld get() const { return ld; }

  virtual boost::shared_ptr<LatticeData> Clone() const { return boost::shared_ptr<LatticeData>(new Derived(*this)); }
  virtual void Init(const BBox3 &bb, float scale) { ld.Init(bb, scale); }

  virtual float Scale() const { return ld.Scale(); }
  virtual void SetScale(float s) { ld.Scale(s); }
  virtual BBox3 Box() const { return ld.Box(); }

  virtual void SetOriginPosition(const Float3 &pos) { ld.SetOriginPosition(pos); }
  virtual Float3 GetOriginPosition() const { return ld.GetOriginPosition(); }
  virtual FloatBBox3 GetWorldBox() const { return ld.GetWorldBox(); }
  virtual void SetCellCentering(const Vec<bool, 3> &cc) { fwd_cell_centering<Ld>::set(ld, cc); }

  virtual int NbCount() const { return ld.NbCount(); }
  virtual SiteType NbSite(SiteType site, int dir) const { return ld.NbSite(site, dir); }
  virtual Int3 NbLattice(const Int3 &p, int dir) const { return ld.NbLattice(p, dir); }
  virtual Float3 LatticeToWorld( const Int3 &p) const { return ld.LatticeToWorld(p); }
  virtual Int3 WorldToLattice(const Float3 &p) const { return WorldToLatticeWrapper(ld, p); }

  virtual SiteType LatticeToSite(const Int3 &p) const { return ld.LatticeToSite(p); }
  virtual const Int3 SiteToLattice(SiteType site) const { return ld.SiteToLattice(site); }
  virtual bool IsInsideLattice(const Int3 &p) const { return ld.IsInsideLattice(p); }

  virtual AxisDirLen GetAxisDirLen(const Int3 &p1, const Int3 &p2 ) const { return ld.GetAxisDirLen(p1, p2); }
  virtual Int3 GetLatticeIndexOnRefinedGrid(const Int3 &pos, int refinement_subdivision) const { return ld.GetLatticeIndexOnRefinedGrid(pos, refinement_subdivision); }

  virtual void print(std::ostream &os) const { ld.print(os); }
  // hdf5 support
  virtual void WriteHdf(h5cpp::Group g) const { WriteHdfLd(g, ld); }
};


inline std::ostream& operator<<(std::ostream &os, const LatticeData &ld)
{
  ld.print(os);
  return os;
}


}

#endif