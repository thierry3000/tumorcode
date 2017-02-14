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

#include <string>

#include <pagmo/src/config.h>
#include <pagmo/src/serialization.h>
#include <pagmo/src/types.h>
#include <pagmo/src/problem/base.h>

//#include "../serialization.h"
//#include "../types.h"
//#include "base.h"
#include <adaption/adaption_model2.h>

namespace pagmo{ namespace problem {


class __PAGMO_VISIBLE adaption_problem : public base
{
	public:
		//adaption_problem(int n = 3);
		adaption_problem(
		  Adaption::Parameters , 
		  BloodFlowParameters,
		  std::auto_ptr<VesselList3d>,
		   int n);
		//adaption_problem(int n=3);
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		//Adaption::Parameters *params;
		//BloodFlowParameters *bfparams;
		//VesselList3d *vl;
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, unsigned int)
// 		void serialize(Archive &ar, unsigned int)
		{
		  ar & boost::serialization::base_object<base>(*this);
		  ar & params;
		  ar & bfparams;
		  ar & vl;
// 		  ar & boost::serialization::base_object<base>(*this) & bfparams & params & vl;
			
		}
		//mutable boost::shared_ptr<BloodFlowParameters> bfparams;
		//mutable boost::shared_ptr<Adaption::Parameters> params;
		//mutable boost::shared_ptr<VesselList3d> vl;
		mutable BloodFlowParameters bfparams;
		mutable Adaption::Parameters params;
		mutable std::auto_ptr<VesselList3d> vl;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::adaption_problem)

#endif
