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
  //std::printf("save_construct_data: adaption_as_pagmo_problem\n");
  ar & t->params;
  ar & t->bfparams;
}
template<class Archive>
inline void load_construct_data(
  Archive &ar, pagmo::problem::adaption_problem *t, const unsigned int file_version)
{
  //retrieve data from archive required to construct new instance
#ifdef DEBUG
  std::printf("load_construct_data: adaption_as_pagmo_problem\n");
#endif
  Adaption::Parameters params;
  ar & params;
  BloodFlowParameters bfparams;
  ar & bfparams;
  //knowing the filename will allow each instance to load the vesselList 
  ::new(t)pagmo::problem::adaption_problem(params,bfparams);
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
 

adaption_problem::adaption_problem(Adaption::Parameters params_, BloodFlowParameters bfparams_):
params(params_),bfparams(bfparams_), base(3)
{
#ifdef DEBUG
  printf("invoking construtor: adaption_as_pagmo_problem\n");
#endif
  //std::mt19937 mt(rd());
  //std::uniform_real_distribution<float> adapt_dist(0.5,4.0);
//   h5cpp::File *readInFile = new h5cpp::File(this->vessel_fn,"r");
//   h5cpp::Group vl_grp = h5cpp::Group(readInFile->root().open_group("adaption/vessels_after_adaption"));
//   //h5cpp::Group vl_grp = h5cpp::Group(readInFile->root().open_group("adaption/recomputed"));
//   this->vl = ReadVesselList3d(vl_grp, make_ptree("filter", false));
  set_lb(0.5);
  set_ub(4.0);
  // is it really necessary to initialize that?
  //params.k_m = adapt_dist(mt);
  //params.k_c = adapt_dist(mt);
  //params.k_s = adapt_dist(mt);
//   for(int i=0;i<3;i++)
//   {
//     printf("x[0]: %f, ", get_best_x()[i]);
//   }
//   printf("\n");
}

/// Clone method.
base_ptr adaption_problem::clone() const {
#ifdef DEBUG
  printf("clone adaption_as_pagmo_problem\n");
#endif
//   std::mt19937 mt(rd());
//   params.k_m = adapt_dist(mt);
//   params.k_c = adapt_dist(mt);
//   params.k_s = adapt_dist(mt);
	return base_ptr(new adaption_problem(params,bfparams));
}

/// Implementation of the objective function.
void adaption_problem::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
  pagmo_assert(f.size() == 1);
  
#ifdef DEBUG
  printf("before: this->params.k_m: %f\n", this->params.k_m);
#endif
  this->params.k_m = x[0];
  this->params.k_c = x[1];
  this->params.k_s = x[2];
#ifdef DEBUG
  printf("after: this->params.k_m: %f\n", this->params.k_m);
#endif
  std::tuple<uint,FlReal> adaption_loop_return;
  {
    //adaption_loop_return = Adaption::runAdaption_Loop(this->params, this->bfparams, *vl, false);
    adaption_loop_return = Adaption::runAdaption_Loop(this->params, this->bfparams, false);
  }
  int returnValue = std::get<0>(adaption_loop_return);
  bool isConvergent = false;
  if(returnValue == 0)
  {
    isConvergent = true;
  }
  if(returnValue == 1)
  {
    isConvergent = false;
  }
  if( isConvergent )
  {
    double average_cap_flow = std::get<1>(adaption_loop_return);
    double suggested_cap_flow = 5000.;
    f[0] = (average_cap_flow-suggested_cap_flow)*(average_cap_flow-suggested_cap_flow);
#ifdef DEBUG
	cout<<format("f[0]: %f\n") % f[0];
#endif
  }
  else
  {
    f[0] = std::numeric_limits< double >::max();
  }
}

std::string adaption_problem::get_name() const
{
	return "secomb adaption";
}
template<class Archive> 
void adaption_problem::serialize(Archive& ar, unsigned int)
{
  std::printf("serialize: adaption_as_pagmo_problem\n");
  ar & boost::serialization::base_object<base>(*this);
  ar & params;
  ar & bfparams;
}
}} //namespaces namespace pagmo { namespace problem {

int archi_best_idx(pagmo::archipelago archi) {
	double min = archi.get_island(0)->get_population().champion().f[0];
	int idx=0;
	for (size_t i=1;i<archi.get_size();++i) {
		double cur = archi.get_island(i)->get_population().champion().f[0];
		if (cur < min) {
			min=cur;
			idx=i;
		}
	}
	return idx;
}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::adaption_problem)



