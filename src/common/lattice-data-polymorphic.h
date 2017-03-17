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
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>

// forward declaration, their location may change
void WriteHdfLd( h5cpp::Group f, const LatticeDataQuad3d &ld );
void ReadHdfLd( h5cpp::Group f, LatticeDataQuad3d &ld );
void WriteHdfLd( h5cpp::Group f, const LatticeDataFCC &ld );
void ReadHdfLd( h5cpp::Group f, LatticeDataFCC &ld );
namespace polymorphic_latticedata{
template<class Ld>
class Derived;
}
//declaration in order to friend it
// namespace boost{ namespace serialization{
// template<class Archive >
// inline void save_construct_data(
//   Archive & ar, const polymorphic_latticedata::Derived<LatticeDataFCC> *t, const unsigned int file_version);
// template<class Archive>
// inline void load_construct_data(
//   Archive & ar, polymorphic_latticedata::Derived<LatticeDataFCC> *t, const unsigned int file_version);
// template<class Archive >
// inline void save_construct_data(
//   Archive & ar, const polymorphic_latticedata::Derived<LatticeDataQuad3d> *t, const unsigned int file_version);
// template<class Archive>
// inline void load_construct_data(
//   Archive & ar, polymorphic_latticedata::Derived<LatticeDataQuad3d> *t, const unsigned int file_version);
// }}//namespace boost{ namespace serialization{

namespace polymorphic_latticedata
{
/*
 * polymorphic variant, for performance reasons this is a seperate class
 */

class LatticeData : boost::noncopyable
{
  protected:
    LatticeData() = default;
  public:
    typedef Int3 LatticeIndexType;
    typedef int64 SiteType;
    
    virtual ~LatticeData() = default;
    virtual boost::shared_ptr<LatticeData> Clone() const {};
    virtual void Init(const BBox3 &bb, float scale) {};

    virtual float Scale() const {};
    virtual void SetScale(float s) {};
    
    virtual BBox3 Box() const {};

    virtual void SetOriginPosition(const Float3 &pos) {};
    virtual Float3 GetOriginPosition() const {};
    virtual FloatBBox3 GetWorldBox() const {};
    virtual void SetCellCentering(const Vec<bool, 3> &cc){};

    virtual int NbCount() const{};
    virtual SiteType NbSite(SiteType site, int dir) const{};
    virtual Int3 NbLattice(const Int3 &p, int dir) const{};
    virtual Float3 LatticeToWorld( const Int3 &p) const{};
    virtual Int3 WorldToLattice(const Float3 &p) const{};
    
    virtual SiteType LatticeToSite(const Int3 &p) const{};
    virtual const Int3 SiteToLattice(SiteType site_) const{};
    virtual bool IsInsideLattice(const Int3 &p) const{};

    virtual AxisDirLen GetAxisDirLen(const Int3 &p1, const Int3 &p2 ) const{};
    virtual Int3 GetLatticeIndexOnRefinedGrid(const Int3 &pos, int refinement_subdivision) const{};

    virtual void print(std::ostream &os) const{};

    /*
     * ldtype = quad or fcc
     */
    static boost::shared_ptr<LatticeData> Make(const char* ldtype, const BBox3 &bb, float scale);

    // hdf 5 support
    static boost::shared_ptr<LatticeData> ReadHdf(h5cpp::Group g);
    virtual void WriteHdf(h5cpp::Group g) const {};
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar.template register_type<polymorphic_latticedata::Derived<LatticeDataFCC>>();
      //ar & boost::serialization::base_object<boost::noncopyable>(*this);
    };
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
  Derived() = default;
  Derived(const Ld &ld) : ld(ld) {
    std::printf("create one of this\n");
  }
  Derived(const BBox3 &bb, float scale) { ld.Init(bb, scale); }
  Derived(const Derived &other) : ld(other.ld) {}
  ~Derived() { }
  //Ld& get() { return ld; }
  Ld get() const { return ld; }

  boost::shared_ptr<LatticeData> Clone() const { return boost::shared_ptr<LatticeData>(new Derived(*this)); }
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
private:
  friend class boost::serialization::access;
//   template<class Archive> 
//   inline friend void boost::serialization::save_construct_data(
// 	Archive & ar, const polymorphic_latticedata::Derived<Ld> *t, const unsigned int file_version);
//   template<class Archive> 
//   inline friend void boost::serialization::load_construct_data(
// 	Archive & ar, polymorphic_latticedata::Derived<Ld> *t, const unsigned int file_version);
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    //boost::serialization::void_cast_register<LatticeData,Derived<Ld>>();
    //ar.template register_type<Derived<Ld>>();
    //ar.template register_type<Ld>();
    ar & boost::serialization::base_object<LatticeData>(*this);
    //ar & boost::serialization::base_object<LatticeData>(Derived<Ld>);
    ar & ld;
  };
};

inline std::ostream& operator<<(std::ostream &os, const LatticeData &ld)
{
  ld.print(os);
  return os;
}

// template<class Ld>
// template<class Archive>
// void Derived<Ld>::serialize(Archive & ar, const unsigned int version)
// {
//    //ar & boost::serialization::base_object<LatticeData>(*this);  //serialize base class
//    //ar & BOOST_SERIALIZATION_NVP(T);
//    ar & ld;
// }

}//namespace polymorphic_latticedata

// namespace boost { namespace serialization {
//   template<class Archive>
//   inline void save_construct_data(
//     Archive &ar, const polymorphic_latticedata::Derived<LatticeDataQuad3d> *t, const unsigned int file_version)
//   {
//     ar & t->ld;
//   }
//   template <class Archive>
//   inline void load_construct_data(
//     Archive &ar, polymorphic_latticedata::Derived<LatticeDataQuad3d> *t, const unsigned int file_version)
//   {
// #ifdef DEBUG
//     std::printf("in load_construct_data at polymorphic_latticedata::Derived<LatticeDataQuad3d");
// #endif
//     LatticeDataQuad3d attribute;
//     ar & attribute;
//     ::new(t)polymorphic_latticedata::Derived<LatticeDataQuad3d>(attribute);
//   }
//   template<class Archive>
//   inline void save_construct_data(
//     Archive &ar, const polymorphic_latticedata::Derived<LatticeDataFCC> *t, const unsigned int file_version)
//   {
//     ar & t->ld;
//   }
//   template <class Archive>
//   inline void load_construct_data(
//     Archive &ar, polymorphic_latticedata::Derived<LatticeDataFCC> *t, const unsigned int file_version)
//   {
// #ifdef DEBUG
//     std::printf("in load_construct_data at polymorphic_latticedata::Derived<LatticeDataFCC");
// #endif
//     LatticeDataFCC attribute;
//     ar & attribute;
//     ::new(t)polymorphic_latticedata::Derived<LatticeDataFCC>(attribute);
//   }
// }}//namespace boost { namespace serialization {

BOOST_SERIALIZATION_ASSUME_ABSTRACT(polymorphic_latticedata::LatticeData)
BOOST_CLASS_EXPORT_KEY(polymorphic_latticedata::LatticeData)
BOOST_CLASS_EXPORT_KEY(polymorphic_latticedata::Derived<LatticeDataQuad3d>)
BOOST_CLASS_EXPORT_KEY(polymorphic_latticedata::Derived<LatticeDataFCC>)


//typedef foo<Ld>polymorphic_latticedata::Derived<Ld>
//BOOST_CLASS_EXPORT_KEY(foo)
//BOOST_CLASS_EXPORT_KEY(polymorphic_latticedata::Derived)
//BOOST_CLASS_EXPORT_KEY(polymorphic_latticedata::Derived<LatticeDataQuad3d>)
//BOOST_CLASS_EXPORT_KEY(polymorphic_latticedata::Derived<LatticeDataFCC>)
//BOOST_CLASS_EXPORT(polymorphic_latticedata::LatticeData)
//BOOST_CLASS_EXPORT(polymorphic_latticedata::Derived<LatticeDataQuad3d>)
//BOOST_CLASS_EXPORT_KEY(polymorphic_latticedata::Derived<LatticeDataFCC>)
//BOOST_CLASS_EXPORT_GUID(polymorphic_latticedata::Derived<LatticeDataQuad3d>, "Quad3d")
//BOOST_CLASS_EXPORT_GUID(polymorphic_latticedata::Derived<LatticeDataFCC>, "FCC")
#endif