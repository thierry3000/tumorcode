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



namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional De Jong's problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
adaption_problem::adaption_problem(std::auto_ptr<VesselList3d> vl_, Adaption::Parameters params_, BloodFlowParameters bfparams_):params(params_),bfparams(bfparams_), base(3)
{
  set_lb(0.5);
  set_ub(4.0);
  this->vl=vl_->Clone();
}
/// Clone method.
base_ptr adaption_problem::clone() const {
	//use the copy constructor of a boost::noncopyable
	//VesselList3d aNewVl(*this->vl);
	//std::auto_ptr<VesselList3d> my_vl(&aNewVl);
	//return base_ptr(new adaption_problem(*this));
	//VesselList3d my_vl(new VesselList3d);
	std::auto_ptr<VesselList3d> my_vl;// = this->vl->Clone();
	//my_vl->Init(this->vl->Ld());
// 	my_vl->init_from_other_vl(*this->vl);
	my_vl = this->vl->Clone();
	return base_ptr(new adaption_problem(my_vl,params,bfparams));
}
/// Implementation of the objective function.
void adaption_problem::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
  pagmo_assert(f.size() == 1);
  this->params.k_m = x[0];
  this->params.k_c = x[1];
  this->params.k_s = x[2];
  std::tuple<uint,FlReal> adaption_loop_return;
  {
    adaption_loop_return = Adaption::runAdaption_Loop(this->params, this->bfparams,this->vl, false);
  }
  double average_cap_flow = std::get<1>(adaption_loop_return);
  double suggested_cap_flow = 200000.;
  f[0] = (average_cap_flow-suggested_cap_flow)*(average_cap_flow-suggested_cap_flow);
#ifdef DEBUG
	cout<<format("f[0]: %f\n") % f[0];
#endif
}

std::string adaption_problem::get_name() const
{
	return "secomb adaption";
}
//getter and setter
Adaption::Parameters adaption_problem::get_params() const
{
  return this->params;
}
BloodFlowParameters adaption_problem::get_bfparams() const
{
  return this->bfparams;
}
std::auto_ptr<VesselList3d> adaption_problem::get_vl() const
{
  return this->vl;
}

//serialization
template <class Archive>
void adaption_problem::serialize(Archive &ar, unsigned int)
{
  ar & boost::serialization::base_object<pagmo::problem::base>(*this);
  ar & params;
  ar & bfparams;
  ar & vl;
}
}} //namespaces namespace pagmo { namespace problem {


BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::adaption_problem)