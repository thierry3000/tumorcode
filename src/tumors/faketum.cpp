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
#include "faketum.h"

#include "../common/calcflow.h"
#include "../common/shared-objects.h"
#include "../python_krebsutils/python_helpers.h"

//Parameters::Parameters()
FakeTum::Parameters::Parameters()
{
  rGf = 200;
  out_intervall = 100;
  tend = 1000;
  dt = 1;
  //num_threads = 1;
  tumor_radius = 200.;
  tumor_speed = 2.;  // \mu m / hour
  vesselfile_ensemble_index = 0;
  tissuePressureDistribution = TISSUE_PRESSURE_SPHERE;
  tissuePressureWidth = 500.;
  tissuePressureCenterFraction = 0.;
  stopping_radius_fraction = 0.6;
  paramset_name = "aname";
  message = "";
  //isRerun = false;
  vessel_path = "vesesls";
#ifdef USE_ADAPTION
  apply_adaption_intervall = 1;// earlier adaption was done in each step, so for backward compatibility, default in 1
#endif
}

void FakeTum::Parameters::assign(const ptree &pt)
{
  #define DOPT(name) boost::property_tree::get(name, #name, pt);
  DOPT(message);
  //DOPT(num_threads);
  DOPT(out_intervall);
  DOPT(tend);
  DOPT(dt);
  DOPT(fn_out);
  DOPT(fn_vessel);
  DOPT(vessel_path);
  DOPT(rGf);
  DOPT(tumor_radius);
  DOPT(tumor_speed);
  DOPT(stopping_radius_fraction);
  DOPT(tissuePressureWidth);
  DOPT(tissuePressureCenterFraction);
  //DOPT(isRerun);
//   std::vector<string> myOptions = {"message", 
//     "paramset_name", 
//     "num_threads",
//     "out_intervall",
//     "tend",
//     "dt",
//     "fn_out",
//     "fn_vessel",
//     "rGf",
//     "tumor_radius",
//     "tumor_speed",
//     "stopping_radius_fraction",
//     "tissuePressureWidth",
//     "tissuePressureCenterFraction",
//   };
//   for( auto aOption: myOptions)
//   {
//     ptree::const_assoc_iterator it = pt.find(aOption);
//     if( it == pt.not_found())
//       std::cout<< format("warning: %s not found") % aOption <<std::endl;
//     else
//       boost::property_tree::get(aOption, #aOption, pt);
//   }

  //DOPT(tissuePressureDistribution);
  string s = pt.get<string>("tissuePressureDistribution");
  if (s == "sphere") tissuePressureDistribution = TISSUE_PRESSURE_SPHERE;
  else if (s == "shell") tissuePressureDistribution = TISSUE_PRESSURE_SHELL;
  else throw std::runtime_error("unknown tissuePressureDistribution "+s);
  const auto bfparamsPtree = pt.get_child_optional("calcflow");
  if (bfparamsPtree) bfparams.assign(*bfparamsPtree);
#ifdef USE_ADAPTION
  const auto adapt_paramsPtree = pt.get_child_optional("adaption");
  if (adapt_paramsPtree) adap_params.assign(*adapt_paramsPtree);
  DOPT(apply_adaption_intervall);
#endif
  #undef DOPT
}

ptree FakeTum::Parameters::as_ptree() const
{
  boost::property_tree::ptree pt;
  #define DOPT(name) pt.put(#name, name)
  DOPT(paramset_name);
  DOPT(message);
  //DOPT(num_threads);
  DOPT(out_intervall);
  DOPT(tend);
  DOPT(dt);
//   DOPT(message);
  DOPT(fn_out);
  DOPT(fn_vessel);
  DOPT(vessel_path);
  DOPT(rGf);
  DOPT(tumor_radius);
  DOPT(tumor_speed);
  //DOPT(isRerun);
  DOPT(stopping_radius_fraction);
  if (tissuePressureDistribution == TISSUE_PRESSURE_SPHERE) pt.put("tissuePressureDistribution", "sphere");
  else if (tissuePressureDistribution == TISSUE_PRESSURE_SHELL) pt.put("tissuePressureDistribution", "shell");
  //DOPT(tissuePressureDistribution);
  DOPT(tissuePressureWidth);
  DOPT(tissuePressureCenterFraction);
  pt.put_child("calcflow", bfparams.as_ptree());
#if USE_ADAPTION
  pt.put_child("adaption", adap_params.as_ptree());
  DOPT(apply_adaption_intervall);
#endif
  #undef DOPT
  return pt;
}

void FakeTum::Parameters::update_ptree(ptree &dst, const ptree &src)
{
  boost::property_tree::update(dst, src);
}


float FakeTum::FakeTumorSim::getGf(const Float3 &pos) const
{
  float r = pos.norm();
  return std::max(tumor_radius + params.rGf - r, 0.);
}

float FakeTum::FakeTumorSim::getPress(const Float3 &pos) const
{
  float r = pos.norm();
  if (params.tissuePressureDistribution == TISSUE_PRESSURE_SPHERE)
  {
    return my::smooth_heaviside<float>(tumor_radius - r, params.tissuePressureWidth);
  }
  else
  {
    return                                       2.*params.tissuePressureWidth*my::smooth_delta_cos<float>(tumor_radius - r, params.tissuePressureWidth) +
           params.tissuePressureCenterFraction * my::smooth_heaviside_sin<float>(tumor_radius - r, params.tissuePressureWidth);
  }
}

float FakeTum::FakeTumorSim::getTumorDens(const Float3 &pos) const
{
  float r = pos.norm();
  return my::smooth_heaviside<float>(tumor_radius - r, 30.);
}

Float3 FakeTum::FakeTumorSim::getGfGrad(const Float3 &pos) const
{
  float r = pos.norm();
  if (r > 0.)
    return (-1./r) * pos;
  else
    return Float3(0.);
}

int FakeTum::FakeTumorSim::run()
{
#ifndef NDEBUG
  std::cout << "starting FakeTumor run in c++" << std::endl;
#endif

  // direct cout through log
  cout.rdbuf(my::log().rdbuf());
  
// #ifdef USE_ADAPTION
//     ptree adaptionSettings = params.adap_params.as_ptree();
//     boost::property_tree::update(adaptionSettings,(pt_params.get_child("adaption")).get_child("adaption"));
//     params.adap_params.assign(adaptionSettings);
//     // BAD HACK
//     // do be done
//     //this->params.adap_params.radMin_for_kill = this->model.params.radMin;
// #endif
  
  H5::H5File file;
  H5::Group h5_vessels;
  //H5::Group h5params;
  try{
    file = H5::H5File(params.fn_vessel, H5F_ACC_RDONLY);
    h5_vessels = file.openGroup(params.vessel_path);
  }
  catch(H5::Exception &e)
  {
    e.printErrorStack();
  }
  
  ptree pt;
  if( ! mySystemParameters.isRerun)
  {
    pt.put("scale subdivide", 10.);
  }
  vl = ReadVesselList3d(h5_vessels, pt);
  
  // adjust vessel list ld to have the center inside the
  // simulation domain, an not at the edge of th cube
  const Float3 c = 0.5 * (vl->Ld().GetWorldBox().max + vl->Ld().GetWorldBox().min);
  vl->SetDomainOrigin(vl->Ld().LatticeToWorld(Int3(0))-c);

  cout << "--------------------"<< endl;
  cout << "Vessel Lattice is: " << endl;
  vl->Ld().print(cout); cout  << endl;
  cout << "--------------------"<< endl;

  h5_vessels.close();
  //h5params.close();
  file.close();
  
  //h5params = file.openGroup("/parameters");
  
  /* this part is for multiple files
    * -> spared for later
    */
  /*
  try{
    string message;
    readAttrFromH5(h5params, string("MESSAGE"),message);
    params.vesselfile_message = message;
    int index;
    readAttrFromH5(h5params, string("ENSEMBLE_INDEX"),index);
    params.vesselfile_ensemble_index = index;
  }
  catch()
  {
    e.printErrorStack();
  
  */
  
  
  
  
//     params.vesselfile_message = file.root().open_group("parameters").attrs().get<string>("MESSAGE");
//     params.vesselfile_ensemble_index = file.root().open_group("parameters").attrs().get<int>("ENSEMBLE_INDEX");
  
  VesselModel1::Callbacks callbacks;
  callbacks.getGf = boost::bind(&FakeTumorSim::getGf, boost::ref(*this), _1);
  callbacks.getPress = boost::bind(&FakeTumorSim::getPress, boost::ref(*this), _1);
  callbacks.getTumorDens = boost::bind(&FakeTumorSim::getTumorDens, boost::ref(*this), _1);
  callbacks.getGfGrad = boost::bind(&FakeTumorSim::getGfGrad, boost::ref(*this), _1);

  /* need to compute flow because shearforce must be
    * known and be consistent with current parameters.
    * Shear force is used e.g. in model.Init to initialize
    * f_initial. */
  CalcFlow(*vl, params.bfparams); 
  vessel_model.Init(vl.get(), this->vessel_model.params, callbacks);
  

  tumor_radius = params.tumor_radius;
  if( !mySystemParameters.isRerun)
  {
    time = 0.;
    num_iteration = 0.;
    output_num = 0;
    next_output_time = 0;
    next_adaption_time = 0;
  }
  else
  {
    //this is store before the increment in the last run
    num_iteration++;
    output_num++;
    time = time+params.dt;
  }
  
  //for the adaption it could be usefull to have the
  //vessel network after the adaption in the beginning   ---> done automatically since adaption is in tum-only-vessels
//   bool writeVesselsafter_initial_adaption = false;
//   if(model.params.badaption_on_off)
//   {
//     writeVesselsafter_initial_adaption = true;
//   }

  int iteration_in_this_rerun=0;
  
  while (iteration_in_this_rerun <= max_iteration_per_rerun and not PyCheckAbort() )
  {
    if (time > params.tend)
    {
      cout << "stopped because of time" << endl;
      break;
    }
    
    double size_limit = 0.5*maxCoeff(Size(vl->Ld().GetWorldBox())) * params.stopping_radius_fraction; 
    //cout << format("size_limit = %f vs tumor_radius = %f\n") % size_limit % tumor_radius;
    
    if (tumor_radius >  size_limit)
    {
      cout << "stopped because of size" << endl;
      break;
    }
#ifdef USE_ADAPTION
#if 1
    if (time >= next_adaption_time - params.dt * 0.1)
    {
      //do adaption if wanted
      //GenerateSprouts();
      //if (IS_DEBUG) vl->IntegrityCheck();
      //VesselModel1::myprint(params.adap_params.as_ptree());
      //note: not yet adaption ready
      std::tuple<uint,FlReal,FlReal, FlReal> return_state;
  
      return_state = Adaption::runAdaption_Loop(*vl, params.adap_params, params.bfparams, false);
      if(std::get<0>(return_state)>0)
      {
        cout<<"adaption not successfull --> bad"<<endl;
      }
      else
      {
        cout<< "adaption good" << endl;
      }
      
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<vl->GetECount();++i)
        {
          Vessel *v = vl->GetEdge(i);
          v->reference_r = v->r;
        }
      }
      next_adaption_time += params.apply_adaption_intervall;
    }
#endif
#endif
    
    if (time >= next_output_time - params.dt * 1.0)
    {
      //this happens only for fixed instances of time
      writeOutput(true);
      next_output_time += params.out_intervall;
    }
    //for a rerun we need to access the latest instant of time
    params.latest_executed_timepoint = time;
    writeOutput(false);

    doStep(params.dt);
    
    //depricated since adaption is now in tum-only-vessls
//     if(writeVesselsafter_initial_adaption)
//     {
//       writeOutput();
//       writeVesselsafter_initial_adaption = false;
//     }
    time += params.dt;
    ++output_num;
    ++num_iteration;
    ++iteration_in_this_rerun;
  }

  return 0;
}



void FakeTum::FakeTumorSim::doStep(double dt)
{
  cout << format("step %i, t=%f") % num_iteration % time << endl;
  CalcFlow(*vl, params.bfparams);

  // adaption switch inside this function   
  vessel_model.DoStep(dt, &params.bfparams);

  tumor_radius += dt * params.tumor_speed;
}


void FakeTum::FakeTumorSim::writeOutput(bool doPermanentSafe)
{
  if( doPermanentSafe )
  {
    cout << format("permanent output %i -> %s") % output_num % params.fn_out << endl;
  }
  else
  {
    cout << format("buffer output %i -> %s") % output_num % params.fn_out << endl;
  }
  H5::H5File f_out;
  H5::Group root, gout, h5_tum, h5_parameters, h5_vessel_parameters,h5_system_parameters;
  
  try{
    if( !mySystemParameters.isRerun)
    {
      f_out = H5::H5File(params.fn_out, output_num==0 ? H5F_ACC_TRUNC : H5F_ACC_RDWR);
    }
    else
    {
      f_out = H5::H5File(params.fn_out, H5F_ACC_RDWR );
    }
    root = f_out.openGroup("/");
  }
  catch(H5::Exception &e)
  {
    e.printErrorStack();
  }

  if (output_num == 0)
  {
    root.createGroup("last_state");
    
    try 
    {
      h5_parameters = root.createGroup("parameters");
      h5_vessel_parameters = h5_parameters.createGroup("vessels");
      h5_system_parameters = h5_parameters.createGroup("system");
      writeAttrToH5(root, string("MESSAGE"), params.message);
      writeAttrToH5(root, string("VESSELTREEFILE"), params.fn_vessel);
      writeAttrToH5(root, string("VESSELTREEPATH"), params.vessel_path);
      writeAttrToH5(root, string("OUTPUT_NAME"), params.fn_out);
      writeAttrToH5(root, string("VESSELFILE_MESSAGE"), params.vesselfile_message);
      writeAttrToH5(root, string("VESSELFILE_ENSEMBLE_INDEX"), params.vesselfile_ensemble_index);
      WriteHdfPtree(h5_vessel_parameters,vessel_model.params.as_ptree());
      WriteHdfPtree(h5_system_parameters, mySystemParameters.as_ptree());
      WriteHdfPtree(h5_parameters, params.as_ptree());
    }
    catch(H5::Exception &e)
    {
      e.printErrorStack();
    }
  }
  
  try
  {
    if(!doPermanentSafe)
    {
      root.unlink("last_state");
      gout=root.createGroup("last_state");
      writeAttrToH5(gout, "CURRENT_RERUN_NUMBER", mySystemParameters.reRunNumber);
    }
    else
    {
      gout = root.createGroup(str(format("out%04i") % output_num));
    }
    writeAttrToH5(gout, string("time"), time);
    writeAttrToH5(gout, string("OUTPUT_NUM"), output_num);
    writeAttrToH5(gout, string("NUM_ITERATION"), num_iteration);
    writeAttrToH5(gout, string("NEXT_OUTPUT_TIME"), next_output_time);
    writeAttrToH5(gout, string("NEXT_ADAPTION_TIME"), next_adaption_time);
    
    H5::Group h5_current_vessels = gout.createGroup("vessels");
    WriteVesselList3d(*vl, h5_current_vessels);
    LatticeDataQuad3d ld;
    SetupFieldLattice(vl->Ld().GetWorldBox(), 3, 100., 0.1 * vl->Ld().Scale(), ld);
    Array3d<float> tum_field(ld.Box());
    FOR_BBOX3(p, ld.Box())
    {
      float t = getTumorDens(ld.LatticeToWorld(p));
      tum_field(p) = t;
    }
    h5_tum = gout.createGroup("tumor");
    writeAttrToH5(h5_tum, string("TYPE"), string("faketumor"));
    writeAttrToH5(h5_tum, string("TUMOR_RADIUS"), tumor_radius);
    // could be done, but since it is a sphere, you can easily calculate the tc_density from the radius
    //WriteScalarField(h5_tum, string("tc_density"), tum_field, ld, field_ld_group);
  }
  catch(H5::Exception &e)
  {
    e.printErrorStack();
  }
  f_out.close();
  /* close all open groups!!! otherwise hdf5 library will not close the file
   * see: https://support.hdfgroup.org/HDF5/doc/RM/RM_H5F.html#File-Close
   */
  root.close();
  gout.close();
  h5_tum.close();
  h5_parameters.close();
  h5_vessel_parameters.close();
  h5_system_parameters.close();
}
