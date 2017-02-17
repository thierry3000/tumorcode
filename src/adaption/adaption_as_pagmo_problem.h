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

#ifndef PAGMO_PROBLEM_ADAPTION_H
#define PAGMO_PROBLEM_ADAPTION_H

#include <list>
#include <memory>
#include <fstream>
#include <string>

#include <cstdio> // remove, std::autoptr inteface wrong in dinkumware
#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::remove;
}
#endif

#include <boost/archive/tmpdir.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/split_free.hpp>



//#define PAGMO_ENABLE_MPI
#include <pagmo/src/pagmo.h>
#include <pagmo/src/config.h>
#include <pagmo/src/serialization.h>
#include <pagmo/src/types.h>
#include <pagmo/src/problem/base.h>

//#include "../serialization.h"
//#include "../types.h"
//#include "base.h"
#include <adaption/adaption_model2.h>

namespace pagmo{ namespace problem {
//const std::auto_ptr<VesselList3d> vl;
// const BloodFlowParameters bfparams;
// const Adaption::Parameters params;
// void set_vl(std::auto_ptr<VesselList3d> vl)
// {
//   pagmo::problem::vl = vl;
// }
class __PAGMO_VISIBLE adaption_problem : public base
{
	//static boost::shared_ptr<VesselList3d> s_vl;
	//static Adaption::Parameters s_params;
	//static BloodFlowParameters s_bfparams;
	public:
		adaption_problem(int n = 3);
		adaption_problem(VesselList3d vl,Adaption::Parameters params_, BloodFlowParameters bfparams, int n = 3);
		//adaption_problem();
		//adaption_problem(int n=3);
		base_ptr clone() const;
		std::string get_name() const;
		//void set_static_members(Adaption::Parameters params, BloodFlowParameters bfparams, std::auto_ptr<VesselList3d> vl);
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, unsigned int)
		{
		  ar & boost::serialization::base_object<base>(*this);
		  ar & params;
		  ar & bfparams;
		  ar & vl;	
		}
		Adaption::Parameters params;
		BloodFlowParameters bfparams;
		mutable VesselList3d vl;
		//VesselList3d vl;
		//mutable boost::shared_ptr<BloodFlowParameters> bfparams;
		//mutable boost::shared_ptr<Adaption::Parameters> params;
		//mutable boost::shared_ptr<VesselList3d> vl;
// 		mutable BloodFlowParameters bfparams;
// 		mutable Adaption::Parameters params;
// 		mutable std::auto_ptr<VesselList3d> vl;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::adaption_problem)


//from http://www.boost.org/doc/libs/1_54_0/libs/serialization/example/demo_auto_ptr.cpp
namespace boost { 
namespace serialization {

/////////////////////////////////////////////////////////////
// implement serialization for auto_ptr< T >
// note: this must be added to the boost namespace in order to
// be called by the library
template<class Archive, class T>
inline void save(
    Archive & ar,
    const std::auto_ptr< T > &t,
    const unsigned int file_version
){
    // only the raw pointer has to be saved
    // the ref count is rebuilt automatically on load
    const T * const tx = t.get();
    ar << tx;
}

template<class Archive, class T>
inline void load(
    Archive & ar,
    std::auto_ptr< T > &t,
    const unsigned int file_version
){
    T *pTarget;
    ar >> pTarget;
    // note that the reset automagically maintains the reference count
    #if BOOST_WORKAROUND(BOOST_DINKUMWARE_STDLIB, == 1)
        t.release();
        t = std::auto_ptr< T >(pTarget);
    #else
        t.reset(pTarget);
    #endif
}

// split non-intrusive serialization function member into separate
// non intrusive save/load member functions
template<class Archive, class T>
inline void serialize(
    Archive & ar,
    std::auto_ptr< T > &t,
    const unsigned int file_version
){
    boost::serialization::split_free(ar, t, file_version);
}

} // namespace serialization
} // namespace boost

#endif
