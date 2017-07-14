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

#include "calcflow_common.h"
#include "mwlib/histogram.h"
#include "vessels3d.h"
#include "mwlib/ptree_ext.h"
#include <boost/property_tree/info_parser.hpp>

//Clin Hemorheol Microcirc. 2008;39(1-4):243-6.
//Plasma viscosity: a forgotten variable.
// they give the viscosity as 1.1 - 1.3 mPa at 37 deg C as normal value!!!!

//Clin Chem. 1996 Aug;42(8 Pt 1):1189-95.
//Distribution of blood viscosity values and biochemical correlates in healthy adults.
// gives 1.39 +/- 0.08

BloodFlowParameters::BloodFlowParameters()
{
  viscosityPlasma = 4.0e-6; // [kPa s] // i think this is actually the blood viscosity measured in-vitro perhaps?!
  rheology = RheologyForHuman;
  inletHematocrit = 0.45;
  includePhaseSeparationEffect = false;
}

void BloodFlowParameters::assign(const ptree &pt)
{
  viscosityPlasma = pt.get<double>("viscosityPlasma", 4.0e-6); // this is the WRONG!!! value, for backward compatibility!!!!!
  rheology        = pt.get<Rheology>("rheology", RheologyForRats);
  inletHematocrit = pt.get<double>("inletHematocrit", 0.45);
  includePhaseSeparationEffect = pt.get<bool>("includePhaseSeparationEffect", false);
}

ptree BloodFlowParameters::as_ptree() const
{
  return make_ptree("viscosityPlasma", viscosityPlasma)
                   ("rheology", rheology)
                   ("inletHematocrit", inletHematocrit)
                   ("includePhaseSeparationEffect", includePhaseSeparationEffect);
}

ostream& operator<<(ostream &os, Rheology rheology)
{
  switch (rheology)
  {
    case RheologyForHuman: os << "RheologyForHuman"; break;
    case RheologyForRats: os << "RheologyForRats"; break;
    case RheologySecomb2005: os << "RheologySecomb2005"; break;
    default:
      throw std::runtime_error("not implemented");
  }
  return os;
}

istream& operator>>(istream &is, Rheology &rheology)
{
  string s;
  is >> s;
  if (s == "RheologyForHuman")     rheology = RheologyForHuman;
  else if (s == "RheologyForRats") rheology = RheologyForRats;
  else if (s == "RheologySecomb2005") rheology = RheologySecomb2005;
  else throw std::runtime_error("reading enum Rheology error: unknown "+s);
  return is;
}


ostream& operator<<(ostream &os, const BloodFlowParameters &bfparams)
{
  const ptree pt = bfparams.as_ptree();
  boost::property_tree::write_info(os, pt);
  return os;
}

/*------------------------------------------------------------------
--------------------------------------------------------------------*/

double CalcFahraeusEffect(double h, double r, Rheology rheology)
{
  double d = r * 2.;
  if (rheology == RheologyForRats) d /= 0.84;
  const double e1 = 0.415;
  const double e2 = 0.011;
  double t1 = 1 + 1.7 * std::exp(-e1 * d) - 0.6 * std::exp(-e2 * d);
  return h*(h + (1.-h)*t1);
}

namespace ImprovedRelativeViscosityInternal 
{
// Code by Thierry; implements relative viscosity formula from Secomb et al. (2005)
// The main change is reduced viscosity for thinnest possible capillaries (r < 4 mu)
// where the previous formula produced a sharp rise in viscosity.
typedef  double myacc;
inline myacc etha_45(myacc d)
{
  return 220.*std::exp(-1.3*d)+3.2-2.44*std::exp(-0.06*std::pow(d,0.645));
}
inline myacc ca(myacc d)
{
  return (0.8+std::exp(-0.075*d))*(-1+1/(1+std::pow(10,-11)*std::pow(d,12)))+1/(1+std::pow(10,-11)*std::pow(d,12));
}

inline myacc etha_vitro(myacc d, myacc h)
{
  return 1+(etha_45(d)-1)*(std::pow(1-h,ca(d))-1)/(std::pow(1-0.45,ca(d))-1);
}
inline myacc was(myacc d)
{
  if( d<2.4)
  {
    return 0.;
  }
  else
  {
    return (d-2.4)/(d+100.-2.*2.4)*2.6;
  }
}
inline myacc wpeak(myacc d)
{
  if(d<2.4)
  {
    return 0.;
  }
  else if ((2.4<d) and (d<10.5))
  {
    return 1.1*(d-2.4)/(10.5-2.4);
  }
  else
  {
    return 1.1*std::exp(-0.03*(d-10.5));
  }
}
inline myacc weff(myacc d, myacc h)
{
  return was(d)+wpeak(d)*(1+h*1.18);
}
inline myacc Deff(myacc d, myacc h)
{
  return d-2*weff(d,h);
}
inline myacc EtaRel(myacc r, myacc h)
{
  return etha_vitro(2*r,h)*std::pow(2*r/Deff(2*r,h),4);
}
} //ImprovedRelativeViscosityInternal

// "In Vivo" Law from Secomb et al. (1994) Resistance to blood flow in microvessels in vivo
double CalcRelViscosity( double r, double h, Rheology rheology)
{
  h = my::cut(h, 0., 0.999); // hematocrit could be slightly beyond allowed limits [0,1]
  r = my::cut(r, 0.1, 2000.); // protect against division by zero due to r=0; flow resistance of smaller vessels (r<0.1) hence depends on hematocrit and r^4, however eta ist constant in this case.
  if (rheology == RheologyForHuman || rheology == RheologyForRats)
  {
    if (rheology == RheologyForRats) r /= 0.84;  // increases the radius that the viscosity formula "sees" because RBC of rats are smaller
    const double u = 6.0*std::exp( -0.17*r ) + 3.2 - 2.44*std::exp(-0.06*std::pow(2.0*r,0.645) );
    const double rr1 = r/5.0;
    const double rr2 = rr1*rr1;
    const double rr4 = rr2*rr2;
    const double rr8 = rr4*rr4;
    const double rr12 = rr4*rr8;
    const double rc = 1.0/(1.0+10.0*rr12);
    const double C = (0.8+std::exp(-0.15*r))*(-1.0+rc)+rc;
    const double f = (2.0*r)/(2.0*r-1.1);
    const double ff = f*f;
    return (1.0 + (u-1.0)*ff*(std::pow(1.0-h,C)-1.0)/(std::pow(1.0-0.45,C)-1.0))*ff;
  }
  else
  {
    return ImprovedRelativeViscosityInternal::EtaRel(r, h);
  }
}

/*------------------------------------------------------------------
  lookup table for viscosity
--------------------------------------------------------------------*/
#if 0
class ViscosityLUT : private LookupTable2D<FlReal,FlReal>
{
  typedef LookupTable2D<FlReal,FlReal> SUPER;
public:
  ViscosityLUT()
  {
    SUPER::Init( 1.0, 400.0, 200, true,
                 0.01, 0.99, 100, false );

    for( uint y=0; y<ny; ++y )
    {
      const float h = yAxis.GetBucketBegin(y);
      for( uint x=0; x<nx; ++x )
      {
        const float r = xAxis.GetBucketBegin(x);
        SUPER::operator()(x,y) = CalcRelViscosity( r, h, RheologyForRats);
      }
    }
  }
  FlReal DoLookup(FlReal r, FlReal h) const { return SUPER::DoLookup(r,h); }
};


FlReal CalcRelViscosityByTable( FlReal r, FlReal h, Rheology rheology)
{
  myAssert(r > 1.e-3 && h > -1.e-3);
  static ViscosityLUT viscosLut;
  if (rheology == RheologyForRats) r /= 0.84;
  return viscosLut.DoLookup( r, h );
}

void CalcViscosities( const FlArray &rad, const FlArray &hema, FlArray &visc )
{
  // init static table
  static ViscosityLUT viscosLut;

  int ecnt = (int)rad.size();
  if(visc.size()!=ecnt) visc.resize(ecnt);
  
  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    float x = bloodFlowParameters.viscosityPlasma;
    //x *= viscosLut.DoLookup(rad[i], hema[i]);
    
    myAssert(x > 0. && isFinite(x));
    visc[i] = x;
  }
}
#endif

void CalcConductivities(const FlArray &rad, const FlArray &len, const FlArray &visc, FlArray &cond)
{
  int ecnt = (int)rad.size();
  if(cond.size()!=ecnt) cond.resize(ecnt);

  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    FlReal coeff = CalcFlowCoeff(visc[i],rad[i],len[i]);
    myAssert(coeff > 0.);
    myAssert(isFinite(coeff));
    cond[i] = coeff;
  }
}


void CalcViscosities( const FlArray &rad, const FlArray &hema, const BloodFlowParameters &bloodFlowParameters, FlArray &visc)
{
  int ecnt = (int)rad.size();
  if(visc.size()!=ecnt) visc.resize(ecnt);
  
  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    float x = bloodFlowParameters.viscosityPlasma;
    x *= CalcRelViscosity(rad[i], hema[i], bloodFlowParameters.rheology);
    myAssert(x > 0.);
    myAssert(isFinite(x));
    visc[i] = x;
  }
}


/*------------------------------------------------------------------
--------------------------------------------------------------------*/
/*
 * pries_secomb_1995: design principles of vascular beds; fig5.
 * fit to the plotted data
 * fits quite well and goes indeed from 2.5kPa to 10kPa
 * from 100 micron diameter veins to arteries
 */
#if 1
FlReal PressureRadiusRelation( FlReal rad, bool bArtery )
{
  if( bArtery ) rad = -rad;
  const FlReal A2 = 18.;
  const FlReal A1 = 89.;
  const FlReal x0 = -21.;
  const FlReal dx = 16;
  FlReal p = A2 + (A1-A2)/(1.+std::exp((rad-x0)/dx));
  p *= 0.1f*1.33f;
  return p;
}
#else


/* from Marieb Human Anatomy & Physiology,
 * blood pressure [mmHg]
 * arterioles: 80 - 35 
 * capillaries: 35 - 15
 * venules:     15 - 10
 * veins:       10 - 0
 */

FlReal PressureRadiusRelation( FlReal rad, bool bArtery )
{
  // see testpressureradiusrelation.py
  if( bArtery ) rad = -rad;
  FlReal A2 = 0.;
  FlReal A1 = 89.;
  FlReal x0 = -17.;
  FlReal dx = 10.;
  FlReal p1 = A2 + (A1-A2)/(1.+std::exp((rad-x0)/dx));
  x0 = -21.;
  dx = 2000;
  FlReal p2 = A2 + (A1-A2)/(1.+std::exp((rad-x0)/dx));
  FlReal p = 0.2*p2+0.8*p1;
  p *= 0.1*1.33;
  return p;
}
#endif

/*------------------------------------------------------------------
--------------------------------------------------------------------*/
uint GetFlowNetwork(CompressedFlowNetwork &fl,
                    const VesselList3d* vl,
                    const BloodFlowParameters &bfparams,
                    bool keepTheVesselHematocrit
 		  )
{
  int ecnt = vl->GetECount();
  int ncnt = vl->GetNCount();

  fl.edges.reserve(ecnt);
  fl.org2new_vertex.resize(ncnt);

  // copy only perfused edges
  for(int i=0; i<ecnt; ++i)
  {
    const Vessel* v = vl->GetEdge(i);
    if (!v->IsCirculated()) 
      continue;
    if (!vl->HasLattice())
    {
      if (v->NodeA()->worldpos == v->NodeB()->worldpos) 
	continue; //well obviously this happens :-(
    }
    int a = fl.org2new_vertex.add(v->NodeA()->Index());
    int b = fl.org2new_vertex.add(v->NodeB()->Index());
    myAssert(a!=b);
    fl.edges.push_back(my::make_eqpair(a, b));
  }
  /*sometimes this happens for to small system sizes, 
   *when there are independant arteries and veins
   * ComputeCirculatedComponents fails an all vessels
   * are labeled as uncirculated. For the next hierachical 
   * iteration this fail here!
   */
  if(fl.org2new_vertex.num_added()== 0 or fl.edges.size() == 0)
  {
    //that does not make sence, but could happen during optimization
    return 1;
  }
  myAssert(fl.org2new_vertex.num_added()>0);
  myAssert(fl.edges.size() > 0);

#if 0//this was micheal way-->freak
  // copy bcs and map index numbers
  auto mapKey = [&](const VesselNode* nd) -> int 
  { 
    return fl.org2new_vertex[nd->Index()];
  };
  remap_keys(vl->GetBCMap(), fl.bcs, mapKey);
#endif
  fl.bcs.clear();
  for(int i=0; i<ncnt; ++i)
  {
    const VesselNode* vc= vl->GetNode(i);
    // insert boundary nodes into the bcs array
    if (vc->IsBoundary())
    {
      int id = fl.org2new_vertex[vc->Index()];
      //myAssert(id<=fl.org2new_vertex.num_added());
      /*
       * fancy new stuff
       * use the BCMap if node is in
       */
      if (id != IntegerMap::unused())//and fl.bcs.find(vc->Index()) == fl.bcs.end())
      {
	if( vl->GetBCMap().find(vc) == vl->GetBCMap().end())//not present
	{
	  fl.bcs[id] = FlowBC(FlowBC::PIN, vc->press);
	}
	else //is present
	{
	  fl.bcs[id] = vl->GetBCMap().at(vc);
	}
      }//if not unused
    }//if boundary
  }//for all nodes
#ifdef DEBUG
  for(auto bc: fl.bcs)
  {
    printf("first: %i, second: %f\n", bc.first, bc.second.val);
  }
#endif
  //if (!(flags & FLOWNET_NOPROPERTIES))
  {
    // copy edge properties
    
    fl.len.resize(fl.edges.size());
    fl.rad.resize(fl.edges.size());
    fl.hema.resize(fl.edges.size());
    const VesselList3d::LatticeData &ld = vl->Ld();
    for(int i=0, k=0; i<ecnt; ++i)
    {
      const Vessel* v = vl->GetEdge(i);
      if (!v->IsCirculated()) continue;
      fl.hema[k] = keepTheVesselHematocrit ? v->hematocrit : bfparams.inletHematocrit;
      fl.rad[k] = v->r;
      double wl=0;
      if (!vl->HasLattice())
      {
        wl = (v->NodeA()->worldpos.transpose()-v->NodeB()->worldpos.transpose()).norm();
	myAssert(wl>0);
#ifdef DEBUG
cout<<"vessel no.: "<<v->Index()<<" length: "<<wl<<endl;
cout<<"posA : "<<v->NodeA()->worldpos<<endl;
cout<<"posB : "<<v->NodeB()->worldpos<<endl;
#endif
      }
      else
      {
	wl = v->WorldLength(ld);
        myAssert(wl>0);
      }
      fl.len[k] = wl;
      ++k;
    }
  }
  return 0;
}


void SetFlowValues( VesselList3d* vl,
                    const CompressedFlowNetwork &fl,
                    const FlArray &cond,
                    const FlArray &press,
                    const FlArray &hema
		    )
{
  int ecnt = vl->GetECount();
  int ncnt = vl->GetNCount();
  
  //set node stuff
  for(int i=0; i<ncnt; ++i)
  {
    VesselNode *nd= vl->GetNode(i);
    if (fl.org2new_vertex[nd->Index()] != IntegerMap::unused())
    {
      nd->press = press[fl.org2new_vertex[nd->Index()]];
    }
  }
  //set edge stuff
  for(int i=0, k=0; i<ecnt; ++i)
  {
    Vessel* v = vl->GetEdge(i);
    if(v->IsCirculated())
    {
      int a = fl.org2new_vertex[v->NodeA()->Index()];
      int b = fl.org2new_vertex[v->NodeB()->Index()];
      myAssert(a != IntegerMap::unused() && b != IntegerMap::unused());
      if(hema[k]>0)
        v->hematocrit = hema[k];
#if DEBUG
      else
        cout<< "bad hema at " << k << "... " << hema[k] <<endl;
#endif
      double wl;
      if( !vl->HasLattice() )
      {
	wl = (v->NodeA()->worldpos.transpose()-v->NodeB()->worldpos.transpose()).norm();
      }
      else
      {
	const VesselList3d::LatticeData &ld = vl->Ld();
	wl = v->WorldLength(ld);
      }
      myAssert(wl>0);
      v->q = cond[k] * std::abs(press[a]-press[b]);
      v->f = CalcShearFromFlowAndCond(v->q, v->r, wl, cond[k]);
      //cout<<v->f*10000.<<endl;
      ++k;
    }
    else
    {
      v->q = 0;
      v->f = 0;
      v->hematocrit = 0.;
    }
  }
}

