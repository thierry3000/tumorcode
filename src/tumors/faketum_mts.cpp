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
#include "faketum_mts.h"

#include "../common/calcflow.h"

#include "../common/shared-objects.h"

//Parameters::Parameters()
FakeTumMTS::Parameters::Parameters()
{
  rGf = 200;
  out_intervall = 100;
  apply_adaption_intervall = 1;// earlier adaption was done in each step, so for backward compatibility, default in 1
  tend = 1000;
  dt = 1;
  num_threads = 1;
  tumor_radius = 200.;
  tumor_speed = 2.;  // \mu m / hour
  vesselfile_ensemble_index = 0;
  tissuePressureDistribution = TISSUE_PRESSURE_SPHERE;
  tissuePressureWidth = 500.;
  tissuePressureCenterFraction = 0.;
  stopping_radius_fraction = 0.6;
  paramset_name = "aname";
}

void FakeTumMTS::Parameters::assign(const ptree &pt)
{
  #define DOPT(name) boost::property_tree::get(name, #name, pt)
  DOPT(paramset_name);
  DOPT(num_threads);
  DOPT(out_intervall);
  DOPT(apply_adaption_intervall);
  DOPT(tend);
  DOPT(dt);
  DOPT(message);
  DOPT(fn_out);
  DOPT(fn_vessel);
  DOPT(rGf);
  DOPT(tumor_radius);
  DOPT(tumor_speed);
  DOPT(stopping_radius_fraction);
  //DOPT(tissuePressureDistribution);
  string s = pt.get<string>("tissuePressureDistribution");
  if (s == "sphere") tissuePressureDistribution = TISSUE_PRESSURE_SPHERE;
  else if (s == "shell") tissuePressureDistribution = TISSUE_PRESSURE_SHELL;
  else throw std::runtime_error("unknown tissuePressureDistribution "+s);
  DOPT(tissuePressureWidth);
  DOPT(tissuePressureCenterFraction);
  #undef DOPT
  const auto bfparamsPtree = pt.get_child_optional("calcflow");
  if (bfparamsPtree) bfparams.assign(*bfparamsPtree);
#ifdef USE_ADAPTION
  const auto adapt_paramsPtree = pt.get_child_optional("adaption");
  if (adapt_paramsPtree) adap_params.assign(*adapt_paramsPtree);
#endif
}

ptree FakeTumMTS::Parameters::as_ptree() const
{
  boost::property_tree::ptree pt;
  #define DOPT(name) pt.put(#name, name)
  DOPT(paramset_name);
  DOPT(num_threads);
  DOPT(out_intervall);
  DOPT(apply_adaption_intervall);
  DOPT(tend);
  DOPT(dt);
  DOPT(message);
  DOPT(fn_out);
  DOPT(fn_vessel);
  DOPT(rGf);
  DOPT(tumor_radius);
  DOPT(tumor_speed);
  DOPT(stopping_radius_fraction);
  if (tissuePressureDistribution == TISSUE_PRESSURE_SPHERE) pt.put("tissuePressureDistribution", "sphere");
  else if (tissuePressureDistribution == TISSUE_PRESSURE_SHELL) pt.put("tissuePressureDistribution", "shell");
  //DOPT(tissuePressureDistribution);
  DOPT(tissuePressureWidth);
  DOPT(tissuePressureCenterFraction);
  #undef DOPT
  pt.put_child("calcflow", bfparams.as_ptree());
#if USE_ADAPTION
  pt.put_child("adaption", adap_params.as_ptree());
#endif
  return pt;
}

void FakeTumMTS::Parameters::update_ptree(ptree &dst, const ptree &src)
{
  boost::property_tree::update(dst, src);
}


float FakeTumMTS::FakeTumorSimMTS::getGf(const Float3 &pos) const
{
  float r = pos.norm();
  return std::max(tumor_radius + params.rGf - r, 0.);
}

float FakeTumMTS::FakeTumorSimMTS::getPress(const Float3 &pos) const
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

float FakeTumMTS::FakeTumorSimMTS::getTumorDens(const Float3 &pos) const
{
  float r = pos.norm();
  return my::smooth_heaviside<float>(tumor_radius - r, 30.);
}

Float3 FakeTumMTS::FakeTumorSimMTS::getGfGrad(const Float3 &pos) const
{
  float r = pos.norm();
  if (r > 0.)
    return (-1./r) * pos;
  else
    return Float3(0.);
}

int FakeTumMTS::FakeTumorSimMTS::run(const ptree &pt_params)
{
  {
    FakeTumMTS::Parameters::update_ptree(all_pt_params, pt_params);
    this->params.assign(all_pt_params);
  }
  // direct cout through log
  cout.rdbuf(my::log().rdbuf());
  {
#ifdef USE_ADAPTION
    // BAD HACK
    this->params.adap_params.radMin_for_kill = this->model.params.radMin;
#endif
    my::SetNumThreads(params.num_threads);
    
    h5cpp::File file(params.fn_vessel, "r");
    ptree pt;
    pt.put("scale subdivide", 10.);
    vl = ReadVesselList3d(file.root().open_group("vessels"), pt);
    
    // adjust vessel list ld
    const Float3 c = 0.5 * (vl->Ld().GetWorldBox().max + vl->Ld().GetWorldBox().min);
    vl->SetDomainOrigin(vl->Ld().LatticeToWorld(Int3(0))-c);

    cout << "--------------------"<< endl;
    cout << "Vessel Lattice is: " << endl;
    vl->Ld().print(cout); cout  << endl;
    cout << "--------------------"<< endl;

    params.vesselfile_message = file.root().open_group("parameters").attrs().get<string>("MESSAGE");
    params.vesselfile_ensemble_index = file.root().open_group("parameters").attrs().get<int>("ENSEMBLE_INDEX");
    
    VesselModel1::Callbacks callbacks;
    callbacks.getGf = boost::bind(&FakeTumorSimMTS::getGf, boost::ref(*this), _1);
    callbacks.getPress = boost::bind(&FakeTumorSimMTS::getPress, boost::ref(*this), _1);
    callbacks.getTumorDens = boost::bind(&FakeTumorSimMTS::getTumorDens, boost::ref(*this), _1);
    callbacks.getGfGrad = boost::bind(&FakeTumorSimMTS::getGfGrad, boost::ref(*this), _1);

    /* need to compute flow because shearforce must be
     * known and be consistent with current parameters.
     * Shear force is used e.g. in model.Init to initialize
     * f_initial. */
    CalcFlow(*vl, params.bfparams); 
    model.Init(vl.get(), this->model.params, callbacks);
  }

  tumor_radius = params.tumor_radius;
  time = 0.;
  num_iteration = 0.;
  output_num = 0;
  
  double next_output_time = 0.;
  double next_adaption_time = 0.;
  
  /**   INIT Milotti   */
  int run_type = 1; //command file
  bool terminal = false;
  string run_name;
  //CellsSystem CellsSystem;	// Standard allocation of the CellsSystem (in this case, the initial dynamic reserve is 2000000)
  currentCellsSystem.Set_Commands( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/commands.txt" );
  currentCellsSystem.Set_CellTypeFile( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/CellType.txt" );
  currentCellsSystem.Set_CellTypeFileAlt( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/CellType.txt" );
  currentCellsSystem.Set_EnvironmentFile( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/Environment.txt" );
  currentCellsSystem.InitializeCellsSystem( terminal );
  cout << "Initialization milotti completed" << endl;
  currentCellsSystem.RunDefinition( );// Run number and output directory output directory & output file opening for metabolism
  currentCellsSystem.Set_nconfiguration( 0 ); // The configuration number is initialized to 0
  currentCellsSystem.Geometry( );// Initial calculation of cluster geometry
  currentCellsSystem.Set_time_from_CGAL(0.);	// Timer reset from last call to CGAL
  if(run_type == 0 || run_type == 1)
    currentCellsSystem.Print2logfile("Cell status at the end of initialization");
  else if (run_type == 2)
    currentCellsSystem.Print2logfile("Cell status at restart of simulation");
  

  currentCellsSystem.CPU_timer(Start_timer);		// start del CPU timer (e reset del timer degli intertempi)
  currentCellsSystem.Timing( true );				// reset del timer
  currentCellsSystem.StepStat( true );			// reset delle statistiche (azzera anche il vettore convergence_fail)

  cout << "\nStartup milotti completed" << endl;
  
  //for the adaption it could be usefull to have the
  //vessel network after the adaption in the beginning   ---> done automatically since adaption is in tum-only-vessels
//   bool writeVesselsafter_initial_adaption = false;
//   if(model.params.badaption_on_off)
//   {
//     writeVesselsafter_initial_adaption = true;
//   }
  while (true)
  {
    if (time >= next_adaption_time - params.dt * 0.1)
    {
      //do adaption if wanted
      if(model.params.badaption_on_off)
      {
	//GenerateSprouts();
	//if (IS_DEBUG) vl->IntegrityCheck();
	//VesselModel1::myprint(params.adap_params.as_ptree());
	//note: not yet adaption ready
	//Adaption::runAdaption_Loop(params.adap_params, params.bfparams, vl, false);
        #pragma omp parallel
        {
          #pragma omp for
          for(int i=0;i<vl->GetECount();++i)
          {
            Vessel *v = vl->GetEdge(i);
            v->reference_r = v->r;
          }
        }
      }
      next_adaption_time += params.apply_adaption_intervall;
    }
    if (time >= next_output_time - params.dt * 0.1)
    {
      
      writeOutput();
      next_output_time += params.out_intervall;
    }

    if (time > params.tend) break;
    
    double size_limit = 0.5*maxCoeff(Size(vl->Ld().GetWorldBox())) * params.stopping_radius_fraction; 
    //cout << format("size_limit = %f vs tumor_radius = %f\n") % size_limit % tumor_radius;
    
    if (tumor_radius >  size_limit) break;

    doStep(params.dt);
    time += params.dt;
    cout << boost::format("advance milotti until: %f") % time;
    //currentCellsSystem.Set_tmax(time);
    doMilottiStep();
    ++num_iteration;
  }
  
  currentCellsSystem.CloseOutputFiles();						// Closing output files
  
  return 0;
}



void FakeTumMTS::FakeTumorSimMTS::doStep(double dt)
{
  cout << format("step %i, t=%f") % num_iteration % time << endl;
  CalcFlow(*vl, params.bfparams);
#ifdef USE_ADAPTION
  model.DoStep(dt, &params.adap_params,&params.bfparams);
#else
  //do be implemented
#endif
  tumor_radius += dt * params.tumor_speed;
}

void FakeTumMTS::FakeTumorSimMTS::doMilottiStep()
{
  cout << format("start mts at tumor time: %f\n" ) % time;
  //   //from milotti


  uint returnValue = currentCellsSystem.runMainLoop( time * 3600 );
  if ( currentCellsSystem.Get_ready2start() )
  {
    currentCellsSystem.Print2logfile("Cells at the end of the run");
    currentCellsSystem.WriteCellsSystem( );					// dump of the final configuration
  }
  
  cout << format(" mts at tumor time: %f\n" ) % time;
//   //end milotti
}

namespace h5 = h5cpp;

void FakeTumMTS::FakeTumorSimMTS::writeOutput()
{
  cout << format("output %i -> %s") % output_num % params.fn_out << endl;
  h5::File f(params.fn_out, output_num==0 ? "w" : "a");
  h5::Group g, root = f.root();

  h5::Attributes a = root.attrs();
  
  if (output_num == 0)
  {
    a.set("MESSAGE",params.message);
    a.set("VESSELTREEFILE",params.fn_vessel);
    a.set("OUTPUT_NAME", params.fn_out);
    a.set("VESSELFILE_MESSAGE", params.vesselfile_message);
    a.set("VESSELFILE_ENSEMBLE_INDEX", params.vesselfile_ensemble_index);
    g = root.create_group("parameters");
    WriteHdfPtree(g.create_group("vessels"), model.params.as_ptree());
    WriteHdfPtree(g, params.as_ptree());
  }

  h5::Group gout = root.create_group(str(format("out%04i") % output_num));
  a = gout.attrs();
  a.set("time", time);
  a.set("OUTPUT_NUM",output_num);
  
  WriteVesselList3d(*vl, gout.create_group("vessels"));
  {
    LatticeDataQuad3d ld;
    SetupFieldLattice(vl->Ld().GetWorldBox(), 3, 100., 0.1 * vl->Ld().Scale(), ld);
    Array3d<float> tum_field(ld.Box());
    FOR_BBOX3(p, ld.Box())
    {
      float t = getTumorDens(ld.LatticeToWorld(p));
      tum_field(p) = t;
    }

    h5::Group field_ld_group = root.require_group("field_ld");
    if (output_num==0) WriteHdfLd(field_ld_group, ld);

    h5::Group gtum = gout.create_group("tumor");
    gtum.attrs().set("TYPE", "faketumor");
    gtum.attrs().set("TUMOR_RADIUS", tumor_radius);
    WriteScalarField(gtum, "tc_density", tum_field, ld, field_ld_group);
  }
  ++output_num;
}
