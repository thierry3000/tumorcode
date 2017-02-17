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
#include "adaption_as_pagmo_problem.h"
#include "adaption_model2.h"

#include <boost/foreach.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional De Jong's problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
/* hm, since Adaption::Parameters are a real class,
 * we we to treat it differently than the BloodFlowParams
 */
//Adaption::Parameters adaption_problem::s_params;
//BloodFlowParameters adaption_problem::s_bfparams;
//boost::shared_ptr<VesselList3d> adaption_problem::s_vl;
// void adaption_problem::set_static_members(Adaption::Parameters params,BloodFlowParameters bfparams, std::auto_ptr< VesselList3d > vl)
// {
//   this->s_params = params;
//   this->s_bfparams = bfparams;
//   this->s_vl = vl;
//   //this->sr_vl = *vl.get();
// }

adaption_problem::adaption_problem(int n): 
base(n)
{
	set_lb(0.5);
	set_ub(4.0);
// 	VesselList3d tmp_vl = VesselList3d();
// 	tmp_vl = *s_vl;
//  	this->p_vl = s_vl;
// 	this->params = params;
// 	this->bfparams = bfparams;
}
adaption_problem::adaption_problem(VesselList3d vl_, Adaption::Parameters params_, BloodFlowParameters bfparams_, int n_):
vl(vl_),params(params_),bfparams(bfparams_), base(n_)
{
  set_lb(0.5);
  set_ub(4.0);
  //this->vl(vl);
  //this->params = params_;
  //this->bfparams = bfparams;
}

/// Clone method.
base_ptr adaption_problem::clone() const {
//	return base_ptr(new adaption_problem(*this));
	VesselList3d aNewVl(this->vl);
	return base_ptr(new adaption_problem(aNewVl,params,bfparams,3));
//	return base_ptr(new adaption_problem(3));
  // 	VesselList3d vl_;
// 	vl_.Init(vl.Ld());
// 	return base_ptr(new adaption_problem());
}

/// Implementation of the objective function.
void adaption_problem::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	decision_vector::size_type n = x.size();
	double retval = 0.0;

	for (decision_vector::size_type i=0; i<n; i++){
		retval += x[i]*x[i];
	}
	f[0]=retval;
// 	this->params->k_m = x[0];
// 	this->params->k_c = x[1];
// 	this->params->k_s = x[2];
	std::tuple<uint,FlReal> adaption_loop_return;
	#pragma omp single
  {
	adaption_loop_return = Adaption::runAdaption_Loop(this->params, this->bfparams,&this->vl, false);
  }
	
	f[0] = std::get<1>(adaption_loop_return);
}

std::string adaption_problem::get_name() const
{
	return "secomb adaption";
}
//adaption_problem::
// template<class Archive>
// void adaption_problem::serialize(Archive &ar, adaption_problem &apt_prob, unsigned int)
// {
//   ar & apt_prob.params;
//   ar & apt_prob.bfparams;
//   ar & apt_prob.vl;
// }


}} //namespaces
// namespace boost {namespace serialization {
// template<class Archive>
// inline void save_construct_data(
// Archive &ar, const pagmo::problem::adaption_problem *t, const unsigned int file_version)
// {
//   ar << t->params;
//   ar << t->bfparams;
//   ar << t->vl;
// }
// }}//namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::adaption_problem)