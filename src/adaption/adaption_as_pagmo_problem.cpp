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


namespace boost{ namespace serialization{
template<class Archive>
inline void save_construct_data(
  Archive &ar, const pagmo::problem::adaption_problem *t, const unsigned int file_version)
{
  //save data required to construct instances
  ar & t->params;
  ar & t->bfparams;
  ar & t->vl;
}
template<class Archive>
inline void load_construct_data(
  Archive &ar, pagmo::problem::adaption_problem *t, const unsigned int file_version)
{
  //retrieve data from archive required to construct new instance
  Adaption::Parameters params;
  ar & params;
  BloodFlowParameters bfparams;
  ar & bfparams;
  boost::shared_ptr<VesselList3d> vl;
  cout<<"I am here "<<endl;
  //vl->init_from_other_vl();
  ar & vl;
  // invoke inplace constructor to initialize instance of adaption_problem
  ::new(t)pagmo::problem::adaption_problem(vl,params,bfparams);
  //t(vl,params,bfparams);
}
}}//namespace boost{ namespace serialization{

namespace pagmo { namespace problem {
/// Constructor from dimension.
/**
 * Will construct an n dimensional De Jong's problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
 

adaption_problem::adaption_problem(boost::shared_ptr<VesselList3d> vl_, Adaption::Parameters params_, BloodFlowParameters bfparams_):
params(params_),bfparams(bfparams_), base(3)
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
	//boost::shared_ptr<VesselList3d> my_vl(new VesselList3d(vl->getLD()));// = this->vl->Clone();
	boost::shared_ptr<VesselList3d> my_vl = this->vl->Clone();
  //my_vl->Init(this->vl->Ld());
// 	my_vl->init_from_other_vl(*this->vl);
	//my_vl = this->vl->Clone();
	return base_ptr(new adaption_problem(my_vl,params,bfparams));
	//return base_ptr(new ackley(*this));
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
template<class Archive> 
void adaption_problem::serialize(Archive& ar, unsigned int)
{
  //boost::serialization::void_cast_register<pagmo::problem::adaption_problem, pagmo::problem::base>();
  ar & boost::serialization::base_object<base>(*this);
  ar & params;
  ar & bfparams;
  ar & vl;
}
}} //namespaces namespace pagmo { namespace problem {
BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::adaption_problem)



