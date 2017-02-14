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
adaption_problem::adaption_problem( Adaption::Parameters param_, BloodFlowParameters bfparams_, std::auto_ptr<VesselList3d> vl_, int n): 
base(n)
{
	set_lb(0.5);
	set_ub(4.0);
	this->params = param_;
	this->bfparams = bfparams_;
	this->vl = vl_;
}
// adaption_problem::adaption_problem( int n): 
// base(n)
// {
// 	set_lb(0.5);
// 	set_ub(4.0);
// 	//this->params = params;
// 	//this->bfparams = bfparams;
// 	//vl = VesselList3d();
// }

/// Clone method.
base_ptr adaption_problem::clone() const {
// 	return base_ptr(new adaption_problem(*this));
// 	return base_ptr(new adaption_problem(3));
// 	VesselList3d vl_;
// 	vl_.Init(vl.Ld());
	return base_ptr(new adaption_problem(params,bfparams,vl,3));
}

/// Implementation of the objective function.
void adaption_problem::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 1);
	decision_vector::size_type n = x.size();
// 	double retval = 0.0;
// 
// 	for (decision_vector::size_type i=0; i<n; i++){
// 		retval += x[i]*x[i];
// 	}
	//this->params->k_m = x[0];
	//this->params->k_c = x[1];
	//this->params->k_s = x[2];
	using namespace boost::accumulators;
	uint return_state = Adaption::runAdaption_Loop(this->params, this->bfparams, this->vl, false);
	accumulator_set<double, features<tag::mean, tag::variance>> acc;
#pragma omp parallel for
  for(int i =0;i<vl->GetECount();++i)
  {
    Vessel* v = vl->GetEdge(i);
    acc(v->q);
  }
  double mean_value = mean(acc);
  double mean_std = sqrt(variance(acc));
#pragma omp barrier
	
	f[0] = mean_std;
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