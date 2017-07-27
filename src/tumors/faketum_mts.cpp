/**
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2016 Thierry Fredrich

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
  gf_production_threshold = 0.1;
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
  int lattice_size_per_single_dim;
  DOPT(paramset_name);
  DOPT(num_threads);
  DOPT(out_intervall);
  DOPT(tend);
  DOPT(dt);
  DOPT(message);
  DOPT(fn_out);
  DOPT(fn_vessel);
  DOPT(rGf);
  DOPT(gf_production_threshold);
  DOPT(tumor_radius);
  DOPT(tumor_speed);
  DOPT(stopping_radius_fraction);
  DOPT(lattice_scale);
  boost::property_tree::get(lattice_size_per_single_dim, "lattice_size", pt);
  lattice_size = {lattice_size_per_single_dim,lattice_size_per_single_dim,lattice_size_per_single_dim};
  //DOPT(tissuePressureDistribution);
  string s = pt.get<string>("tissuePressureDistribution");
  if (s == "sphere") tissuePressureDistribution = TISSUE_PRESSURE_SPHERE;
  else if (s == "shell") tissuePressureDistribution = TISSUE_PRESSURE_SHELL;
  else throw std::runtime_error("unknown tissuePressureDistribution "+s);
  DOPT(tissuePressureWidth);
  DOPT(tissuePressureCenterFraction);
  
  const auto bfparamsPtree = pt.get_child_optional("calcflow");
  if (bfparamsPtree) bfparams.assign(*bfparamsPtree);
#ifdef USE_ADAPTION
  DOPT(apply_adaption_intervall);
  const auto adapt_paramsPtree = pt.get_child_optional("adaption");
  if (adapt_paramsPtree) adap_params.assign(*adapt_paramsPtree);
#endif
  #undef DOPT
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
  DOPT(gf_production_threshold);
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
    cout<<"print inside"<<endl;
    printPtree(pt_params);
    this->params.assign(pt_params);
  }
  // direct cout through log
  cout.rdbuf(my::log().rdbuf());
  {

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
#ifndef undo
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
#endif  

  
  //for the adaption it could be usefull to have the
  //vessel network after the adaption in the beginning   ---> done automatically since adaption is in tum-only-vessels
//   bool writeVesselsafter_initial_adaption = false;
//   if(model.params.badaption_on_off)
//   {
//     writeVesselsafter_initial_adaption = true;
//   }
  
  boost::optional<h5cpp::Group> lastTumorGroupWrittenByFakeTum;
  std::string lastTumorGroupWrittenByFakeTumName;
//   h5cpp::Group lastTumorGroupWrittenByFakeTum;
  
#ifdef USE_DETAILED_O2
//   h5cpp::Group lastTumorGroupWrittenByFakeTum;
  //set up oxygen calculation
  
  // this should be read from file later on
  //DetailedPO2::Parameters oxy_params;
  double grid_lattice_const = 40;
  double safety_layer_size = 120;
  boost::optional<Int3> grid_lattice_size;
#endif
  //init continum lattice stuff
//     Int3 s = params.lattice_size;
    Int3 s = params.lattice_size;
    int dim = s[2]<=1 ? (s[1]<=1 ? 1 : 2) : 3;
    LatticeDataQuad3d field_ld;
    Bool3 centering = Bool3::mapIndex([=](int i) { return i<dim; });
    field_ld.Init(params.lattice_size, params.lattice_scale);
    field_ld.SetCellCentering(centering);
    field_ld.SetOriginPosition(-field_ld.GetWorldBox().max.cwiseProduct(centering.cast<float>()) * 0.5); // set origin = lower left side
    grid.init(field_ld, dim);
    mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
//     SetupTissuePhases(phases, grid, mtboxes, lastTumorGroupWrittenByFakeTum);//filling
  
//   gf_model.init(grid, mtboxes, all_pt_params);
//   gf_model.initField(state.gffield);
    
  oxyops.init(mtboxes, grid.Box(), grid.dim, 2);
  oxyops.init(state.o2field);
  glucoseOps.init(mtboxes, grid.Box(), grid.dim, 2);
  glucoseOps.init(state.glucoseField);
// //   ptree aNewPtree = glucoseParams.as_ptree();
//   boost::property_tree::ptree ptr2;
//   glucoseParams.assigen(ptr2);
//   ptr2.put("name", "d6");
  all_pt_params.add_child("glucose" , glucoseParams.as_ptree());
  last_chem_update = -1;
  last_vessels_checksum = -1;
  UpdateVesselVolumeFraction();
  
  while (true)
  {
    if (time >= next_adaption_time - params.dt * 0.1)
    {
#ifdef USE_ADAPTION
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
#endif
    }
#pragma omp barrier
    if (time >= next_output_time - params.dt * 0.1)
    {
      lastTumorGroupWrittenByFakeTumName = writeOutput();
      h5::File f(params.fn_out + ".h5", "a");
      lastTumorGroupWrittenByFakeTum = f.root().open_group(lastTumorGroupWrittenByFakeTumName+"/tumor");
#ifdef USE_DETAILED_O2
      {
        o2_sim.init(o2_params, params.bfparams,*vl,grid_lattice_const, safety_layer_size, grid_lattice_size, lastTumorGroupWrittenByFakeTum);
        o2_sim.run(*vl);
      }
#else
      
      
#endif 
      SetupTissuePhases(phases, grid, mtboxes, lastTumorGroupWrittenByFakeTum);//filling
      next_output_time += params.out_intervall;
    }

    if (time > params.tend) break;
    
    double size_limit = 0.5*maxCoeff(Size(vl->Ld().GetWorldBox())) * params.stopping_radius_fraction; 
    //cout << format("size_limit = %f vs tumor_radius = %f\n") % size_limit % tumor_radius;
    
    if (tumor_radius >  size_limit) break;

    doStep(params.dt);
    time += params.dt;
    cout << boost::format("advance milotti until: %f") % time;
#ifndef undo
    //currentCellsSystem.Set_tmax(time);
    doMilottiStep();
#endif
    ++num_iteration;
  }
#ifndef undo
  currentCellsSystem.CloseOutputFiles();						// Closing output files
#endif 
  return 0;
}



void FakeTumMTS::FakeTumorSimMTS::doStep(double dt)
{
  cout << format("step %i, t=%f") % num_iteration % time << endl;
  CalcFlow(*vl, params.bfparams);
  /* do not use adaptation stuff for mts*/
// #ifdef USE_ADAPTION
//   model.DoStep(dt, &params.adap_params,&params.bfparams);
// #else
//   //do be implemented
// #endif
  model.DoStep(dt, &params.bfparams);
  calcChemFields();
  tumor_radius += dt * params.tumor_speed;
}

#ifndef undo
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
#endif

std::string FakeTumMTS::FakeTumorSimMTS::writeOutput()
{
  cout << format("output %i -> %s") % output_num % params.fn_out << endl;
  //if this is the first output, we have to create the file, otherwise we append
  h5::File f(params.fn_out + ".h5", output_num==0 ? "w" : "a");
  h5::Group g, root = f.root();
  h5::Group g_o2;

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
    g_o2= g.create_group("o2_params");
    WriteHdfPtree(g_o2, o2_params.as_ptree());
  }
  std::string tumOutName = str(format("out%04i") % output_num);
  h5::Group gout = root.create_group(tumOutName);
  a = gout.attrs();
  a.set("time", time);
  a.set("OUTPUT_NUM",output_num);
  
  WriteVesselList3d(*vl, gout.create_group("vessels"));
  {
//     LatticeDataQuad3d ld;
//     SetupFieldLattice(vl->Ld().GetWorldBox(), 3, 100., 0.1 * vl->Ld().Scale(), ld);
//     Array3d<float> tum_field(ld.Box());
    Array3d<float> tum_field(grid.ld.Box());
    FOR_BBOX3(p, grid.ld.Box())
    {
      float t = getTumorDens(grid.ld.LatticeToWorld(p));
      tum_field(p) = t;
    }

    h5::Group field_ld_group = root.require_group("field_ld");
    if (output_num==0) WriteHdfLd(field_ld_group, grid.ld);

    h5::Group gtum = gout.create_group("tumor");
    gtum.attrs().set("TYPE", "faketumor");
    gtum.attrs().set("TUMOR_RADIUS", tumor_radius);
    WriteScalarField(gtum, "tc_density", tum_field, grid.ld, field_ld_group);
//     WriteScalarField(gtum, "fieldGf", state.gffield, grid.ld, field_ld_group);
#ifdef USE_DETAILED_O2
    //something is different here, I do not know how to handle this 
//      WriteScalarField(gtum, "fieldDetailedOxy", o2_sim.po2field, o2_sim.grid.ld, field_ld_group);
#else
    WriteScalarField(gtum, "fieldOxy", state.o2field, grid.ld, field_ld_group);
#endif
    WriteScalarField(gtum, "fieldGlucose", state.glucoseField, grid.ld, field_ld_group);
//     WriteScalarField(gtum, "fieldOxy", state.o2field, ld, field_ld_group);
    // vessel continuum
    UpdateVesselVolumeFraction();
    WriteScalarField(gtum, "vessel_volume_fraction", vessel_volume_fraction, grid.ld, field_ld_group);
    WriteScalarField(gtum, "oxy_source_lin", vessel_o2src_clin, grid.ld, field_ld_group);
  }
  ++output_num;
  return tumOutName;
}
void FakeTumMTS::FakeTumorSimMTS::calcChemFields()
{
  {
    my::log().push("o2:");

    ptree pt_params = make_ptree("preconditioner","multigrid")("output", 1)("max_iter", 100);
#ifndef USE_DETAILED_O2
    StationaryDiffusionSolve(grid, mtboxes, state.o2field, boost::bind(&FakeTumorSimMTS::insertO2Coefficients, this, _1, _2, boost::cref(state), _3) , pt_params);
#endif
    StationaryDiffusionSolve(grid, mtboxes, state.glucoseField, boost::bind(&FakeTumorSimMTS::insertGlucoseCoefficients, this, _1, _2, boost::cref(state), _3) , pt_params);
    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
      {
        CopyBorder(state.o2field, bbox, grid.Box(), grid.dim, 2);
      }
    }
#ifdef DEBUG
    //cout << "stats: " << state.o2field.valueStatistics() << endl;
#endif
    my::log().pop();
  }

  {
//     my::log().push("gf:");
//     cout << "update" << endl;

    Array3df src;
    Array3dOps<float>(mtboxes, grid.Box(), grid.dim, 0).init(src, false);

    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bb, mtboxes.getCurrentThreadRange())
      {
        //Array3df phases[3];
        float phases_tcs;
        //tie(phases[TISSUE], phases[TCS], phases[DEAD]) = tumor_model.getTissuePhases(bb, state.tumor);
        //phases[TCS] = getTumorDens();
        FOR_BBOX3(p, bb)
        {
          float o2 = state.o2field(p);
//           phases_tcs = phases.phase_arrays[TCS](p);
          //phases_tcs = getTumorDens(grid.ld.LatticeToWorld(p));
//           phases_tcs = getTumorDens(p);
          src(p) = (o2 < params.gf_production_threshold) ? phases.phase_arrays[TCS](p) : 0;
//           src(p) = (o2 < params.gf_production_threshold) ? phases_tcs : 0;
        }
      }
    }
//     gf_model.update(state.gffield, src);
// #ifdef DEBUG
//     cout << "stats: " << state.gffield.valueStatistics() << endl;
// #endif
//     my::log().pop();
  }
  state.chem_checksum++;
}

#ifndef USE_DETAILED_O2
void FakeTumMTS::FakeTumorSimMTS::insertO2Coefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder &mb)
{
  Array3d<float> l_coeff(bbox),
                 rhs(bbox);
  l_coeff[bbox] += vessel_o2src_clin[bbox];
  rhs[bbox] += vessel_o2src_crhs[bbox];

//   Array3df phases[3];
//   tie(phases[TISSUE], phases[TCS], phases[DEAD]) = tumor_model.getTissuePhases(bbox, state.tumor);
// 
  float hemostatic_cell_frac_norm[3];
//   hemostatic_cell_frac_norm[TISSUE] = 1./tumor_model.params.ncells_norm;
//   hemostatic_cell_frac_norm[TCS] = 1./tumor_model.params.ncells_tumor;
//   hemostatic_cell_frac_norm[DEAD] = 1./tumor_model.params.ncells_tumor;
  hemostatic_cell_frac_norm[TISSUE] = 1./0.4;
  hemostatic_cell_frac_norm[TCS] = 1./0.6;
  hemostatic_cell_frac_norm[DEAD] = 1./0.6;
//   
  mb.AddDiffusion(bbox, ConstValueFunctor<float>(1.), -1.);
  
  FOR_BBOX3(p, bbox)
  {
    const float loc_phases[3] = { phases.phase_arrays[0](p), phases.phase_arrays[1](p), phases.phase_arrays[2](p) };
//     const float loc_phases = phases.phase_arrays[0](p);
    //const float loc_phases[3] = { phases[0](p), phases[1](p), phases[2](p) };
    // we must ensure that the resulting diffusion distance agrees with the parameters which define it
    // by l = sqrt(c/D). So the effect that the cell volume fraction is about 0.5 must be canceled out.
    float loc_l_coeff = 0.;
    for (int i=0; i<3; ++i)
    {
      loc_l_coeff -= hemostatic_cell_frac_norm[i] * loc_phases[i] * o2_params.o2_cons_coeff[i];
    }
//     loc_l_coeff = hemostatic_cell_frac_norm[0] * loc_phases * o2_params.o2_cons_coeff[0];
    mb.AddLocally(p, -loc_l_coeff - l_coeff(p), -rhs(p)); /// this was it before
//     mb.AddLocally(p, - l_coeff(p), -rhs(p));
  }
}
#endif

void FakeTumMTS::FakeTumorSimMTS::UpdateVesselVolumeFraction()
{
  my::LogScope log_push_(my::log(), "ves:");
  
  // update only when the vessel system state has changed since last update 
  if (!vessel_volume_fraction.empty() && last_vessels_checksum==state.vessels_checksum)
    return;

  cout << "volume fraction and o2 sources update!" << endl;
  
  VesselsInBoxes vessboxes;
  SortVesselsIntoMtBoxGrid(grid.ld, *vl, 2, mtboxes, vessboxes);
  
  // reinit and fill the array
  Array3dOps<float> ops(mtboxes, grid.Box(), grid.dim, 0);
  ops.init(vessel_volume_fraction, true);
  ops.init(vessel_o2src_clin);
  ops.init(vessel_o2src_crhs);
  ops.init(vessel_glucosesrc_clin);
  ops.init(vessel_glucosesrc_crhs);

#ifndef USE_DETAILED_O2
  #pragma omp parallel
  {
    VesselVolumeGenerator volumegen(*vl, grid.ld, grid.dim, make_ptree("samples_per_cell", 100));
    
    BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, mtboxes.getCurrentThreadRange())
    {
      int dummy_;
      volumegen.Fill(bbox, vessel_volume_fraction, vessboxes[bbox.global_index], dummy_);
      O2Model::AddSourceDistribution( bbox, grid.ld, 
                                      grid.dim, 
                                      vessel_o2src_clin, 
                                      vessel_o2src_crhs, 
                                      vl->Ld(), 
                                      vessboxes[bbox.global_index], 
                                      all_pt_params.get_child("simple_o2"));
    }
  }
#endif
  #pragma omp parallel
  {
    VesselVolumeGenerator volumegen(*vl, grid.ld, grid.dim, make_ptree("samples_per_cell", 100));
    
    BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, mtboxes.getCurrentThreadRange())
    {
      int dummy_;
      volumegen.Fill(bbox, vessel_volume_fraction, vessboxes[bbox.global_index], dummy_);
      //ptree pt = make_ptree("seed", 123456);
      GlucoseModel::AddSourceDistribution( bbox, grid.ld, 
                                      grid.dim, 
                                      vessel_glucosesrc_clin, 
                                      vessel_glucosesrc_crhs, 
                                      vl->Ld(), 
                                      vessboxes[bbox.global_index], 
                                      all_pt_params.get_child("glucose"));
    }
  }

  last_vessels_checksum = state.vessels_checksum;
}
void FakeTumMTS::FakeTumorSimMTS::insertGlucoseCoefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder &mb)
{
  Array3d<float> l_coeff(bbox),
                 rhs(bbox);
  l_coeff[bbox] += vessel_glucosesrc_clin[bbox];
  rhs[bbox] += vessel_glucosesrc_crhs[bbox];

//   Array3df phases[3];
//   tie(phases[TISSUE], phases[TCS], phases[DEAD]) = tumor_model.getTissuePhases(bbox, state.tumor);
// 
//   float hemostatic_cell_frac_norm[3];
// //   hemostatic_cell_frac_norm[TISSUE] = 1./tumor_model.params.ncells_norm;
// //   hemostatic_cell_frac_norm[TCS] = 1./tumor_model.params.ncells_tumor;
// //   hemostatic_cell_frac_norm[DEAD] = 1./tumor_model.params.ncells_tumor;
//   hemostatic_cell_frac_norm[TISSUE] = 1./0.4;
//   hemostatic_cell_frac_norm[TCS] = 1./0.6;
//   hemostatic_cell_frac_norm[DEAD] = 1./0.6;
//   
  //mb.AddDiffusion(bbox, ConstValueFunctor<float>(1.), -1.);
  //mb.AddDiffusion<> (bbox, ConstValueFunctor<float>(1.), -params.po2_kdiff);
  mb.AddDiffusion<> (bbox, ConstValueFunctor<float>(1.), -0.001);
  
  FOR_BBOX3(p, bbox)
  {
    const float loc_phases[3] = { phases.phase_arrays[0](p), phases.phase_arrays[1](p), phases.phase_arrays[2](p) };
//     const float loc_phases = phases.phase_arrays[0](p);
    //const float loc_phases[3] = { phases[0](p), phases[1](p), phases[2](p) };
    // we must ensure that the resulting diffusion distance agrees with the parameters which define it
    // by l = sqrt(c/D). So the effect that the cell volume fraction is about 0.5 must be canceled out.
    float loc_l_coeff = 0.;
    //assume similar behaviour for all phases
    for (int i=0; i<3; ++i)
    {
      //loc_l_coeff -= hemostatic_cell_frac_norm[i] * loc_phases[i] * o2_params.o2_cons_coeff[i];
      loc_l_coeff -= 1 * loc_phases[i] * 1;
    }
//     loc_l_coeff = hemostatic_cell_frac_norm[0] * loc_phases * o2_params.o2_cons_coeff[0];
    mb.AddLocally(p, -loc_l_coeff - l_coeff(p), -rhs(p)); /// this was it before
//     mb.AddLocally(p, - l_coeff(p), -rhs(p));
  }
}

