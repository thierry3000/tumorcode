/**
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2017 Thierry Fredrich

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have receivqed a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "faketum_mts.h"

void FakeTumMTS::FakeTumorSimMTS::initMilotti()
{
  /* let the CellsSystem know the continuous grid of tumorcode */
#if VBL_USE_TUMORCODE
  std::cout<< "use tumorcod" << std::endl;
  tumorcode_pointer_to_currentCellsSystem->Set_Tumorcode_O2_uptake_model(o2_uptake_model);
#endif
    /**   INIT Milotti   */
  int run_type = 1; //command file
  bool terminal = false;
  string run_name;
  bool all_parameter_files_there= true;
  std::vector<std::string> possibleFiles = {"commands.txt", "CellType.txt", "Environment.txt"};
  for( auto afileName:possibleFiles)
  {
    all_parameter_files_there = boost::filesystem::exists(afileName) && all_parameter_files_there;
  }
  if( not all_parameter_files_there )
  {
    std::cout<<"not all parameter files found"<< std::endl;
    exit(EXIT_FAILURE);
  }
//   currentCellsSystem->Set_Commands( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/commands.txt" );
//   currentCellsSystem->Set_CellTypeFile( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/CellType.txt" );
//   currentCellsSystem->Set_CellTypeFileAlt( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/CellType.txt" );
//   currentCellsSystem->Set_EnvironmentFile( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/Environment.txt" );
  tumorcode_pointer_to_currentCellsSystem->Set_Commands( "commands.txt" );
  tumorcode_pointer_to_currentCellsSystem->Set_CellTypeFile( "CellType.txt" );
  tumorcode_pointer_to_currentCellsSystem->Set_CellTypeFileAlt( "CellType.txt" );
  tumorcode_pointer_to_currentCellsSystem->Set_EnvironmentFile( "Environment.txt" );
  
  if( !mySystemParameters.isRerun )
  {
    tumorcode_pointer_to_currentCellsSystem->InitializeCellsSystem( terminal );
    cout << "Initialization milotti completed" << endl;
    tumorcode_pointer_to_currentCellsSystem->RunDefinition( );// Run number and output directory output directory & output file opening for metabolism
    tumorcode_pointer_to_currentCellsSystem->Set_nconfiguration( 0 ); // The configuration number is initialized to 0
    if(tumorcode_pointer_to_currentCellsSystem->Get_ncells() > 1)
    {
      tumorcode_pointer_to_currentCellsSystem->CleanCellsSystem();
      tumorcode_pointer_to_currentCellsSystem->Geometry();
      
      tumorcode_pointer_to_currentCellsSystem->Set_time_from_CGAL(0.);		// Timer reset from last call to CGAL
    }
    else 
    {
      tumorcode_pointer_to_currentCellsSystem->NoGeometry( );
    }
    
    tumorcode_pointer_to_currentCellsSystem->Set_time_from_CGAL(0.);	// Timer reset from last call to CGAL
    if(run_type == 0 || run_type == 1)
      tumorcode_pointer_to_currentCellsSystem->Print2logfile(std::string("Cell status at the end of initialization"));
    else if (run_type == 2)
      tumorcode_pointer_to_currentCellsSystem->Print2logfile(std::string("Cell status at restart of simulation"));
    

    tumorcode_pointer_to_currentCellsSystem->CPU_timer(vbl::timer_button::Start_timer);		// start del CPU timer (e reset del timer degli intertempi)
    tumorcode_pointer_to_currentCellsSystem->Timing( true );				// reset del timer
    tumorcode_pointer_to_currentCellsSystem->StepStat( true );			// reset delle statistiche (azzera anche il vettore convergence_fail)
    if (tumorcode_pointer_to_currentCellsSystem->Get_alive()==0)
    {
      throw std::runtime_error(" no alive cell present");
    }
  }
  else
  {
    //tumorcode_pointer_to_currentCellsSystem->InitializeCellsSystem( terminal );
    cout << "Initialization milotti completed" << endl;
    //tumorcode_pointer_to_currentCellsSystem->RunDefinition( );// Run number and output directory output directory & output file opening for metabolism
    //tumorcode_pointer_to_currentCellsSystem->Set_nconfiguration( 0 ); // The configuration number is initialized to 0
//     if(tumorcode_pointer_to_currentCellsSystem->Get_ncells() > 1)
//     {
//       tumorcode_pointer_to_currentCellsSystem->CleanCellsSystem();
//       tumorcode_pointer_to_currentCellsSystem->Geometry();
//       
//       tumorcode_pointer_to_currentCellsSystem->Set_time_from_CGAL(0.);		// Timer reset from last call to CGAL
//     }
//     else 
//     {
//       tumorcode_pointer_to_currentCellsSystem->NoGeometry( );
//     }
    
    tumorcode_pointer_to_currentCellsSystem->Set_time_from_CGAL(0.);	// Timer reset from last call to CGAL
//     if(run_type == 0 || run_type == 1)
//       tumorcode_pointer_to_currentCellsSystem->Print2logfile(std::string("Cell status at the end of initialization"));
//     else if (run_type == 2)
//       tumorcode_pointer_to_currentCellsSystem->Print2logfile(std::string("Cell status at restart of simulation"));
    

    tumorcode_pointer_to_currentCellsSystem->CPU_timer(vbl::timer_button::Start_timer);		// start del CPU timer (e reset del timer degli intertempi)
    tumorcode_pointer_to_currentCellsSystem->Timing( false );				// reset del timer
    tumorcode_pointer_to_currentCellsSystem->StepStat( true );			// reset delle statistiche (azzera anche il vettore convergence_fail)
    if (tumorcode_pointer_to_currentCellsSystem->Get_alive()==0)
    {
      throw std::runtime_error(" no alive cell present");
    }
  }
  
  cout << "\nStartup milotti completed" << endl; 
}

// FakeTumMTS::SystemParameters::SystemParameters()
// {
//   num_threads = 1;
//   cluster = "local";
//   computing_node = "local";
//   num_threads_queuing = 1;
//   mem_in_GB = 1;
// }
// void FakeTumMTS::SystemParameters::assign(const boost::property_tree::ptree& pt)
// {
//   #define DOPT(name) boost::property_tree::get(name, #name, pt)
//   DOPT(num_threads);
//   DOPT(cluster);
//   DOPT(computing_node);
//   DOPT(num_threads_queuing);
//   DOPT(mem_in_GB);
//   #undef DOPT
// }
// boost::property_tree::ptree FakeTumMTS::SystemParameters::as_ptree() const
// {
//   boost::property_tree::ptree pt;
//   #define DOPT(name) pt.put(#name, name)
//   DOPT(cluster);
//   DOPT(computing_node);
//   DOPT(num_threads);
//   DOPT(num_threads_queuing);
//   DOPT(mem_in_GB);
//   #undef DOPT
//   return pt;
// }

FakeTumMTS::Parameters::Parameters()
{
  out_intervall = 100;
  tend = 1000;
  dt = 1;
  latest_executed_timepoint = 0;
  
  apply_adaption_intervall = 1;// earlier adaption was done in each step, so for backward compatibility, default in 1
  message = "";
  max_iteration_per_rerun = -1;
  fn_out = "";
  fn_vessel = "";
  vessel_path = "vessels";
  vesselfile_message = "";
  paramset_name = "aname";
  vesselfile_ensemble_index = 0;
  rGf = 200;
  gf_production_threshold = 0.1;
  rO2Consumtion = 10;
  tumor_radius = 200.;
  tumor_speed = 2.;  // \mu m / hour
  stopping_radius_fraction = 0.6;
  tissuePressureDistribution = TISSUE_PRESSURE_SPHERE;
  tissuePressureWidth = 500.;
  tissuePressureCenterFraction = 0.;
  lattice_size = {42,42,42};
  lattice_scale = 42;
  useConstO2 = true;
  useTumorcodeVessels = true;
  output_num = 0;
  time = 0;
}


void FakeTumMTS::Parameters::assign(const ptree &pt)
{
  #define DOPT(name) boost::property_tree::get(name, #name, pt)
  DOPT(max_iteration_per_rerun);
  //int lattice_size_per_single_dim;
  DOPT(time);
  DOPT(useConstO2);
  DOPT(useTumorcodeVessels);
  DOPT(paramset_name);
  DOPT(out_intervall);
  DOPT(latest_executed_timepoint);
  DOPT(tend);
  DOPT(dt);
  DOPT(message);
  DOPT(fn_out);
  DOPT(fn_vessel);
  DOPT(vessel_path);
  DOPT(rGf);
  DOPT(rO2Consumtion);
  DOPT(gf_production_threshold);
  DOPT(tumor_radius);
  DOPT(tumor_speed);
  DOPT(stopping_radius_fraction);
  DOPT(lattice_scale);
  //boost::property_tree::get(lattice_size_per_single_dim, "lattice_size", pt);
  //lattice_size = {lattice_size_per_single_dim,lattice_size_per_single_dim,lattice_size_per_single_dim};
  //DOPT(tissuePressureDistribution);
  string s = pt.get<string>("tissuePressureDistribution");
  if (s == "sphere") tissuePressureDistribution = TISSUE_PRESSURE_SPHERE;
  else if (s == "shell") tissuePressureDistribution = TISSUE_PRESSURE_SHELL;
  else throw std::runtime_error("unknown tissuePressureDistribution "+s);
  DOPT(tissuePressureWidth);
  DOPT(tissuePressureCenterFraction);
  
//   const auto bfparamsPtree = pt.get_child_optional("calcflow");
//   if (bfparamsPtree) bfparams.assign(*bfparamsPtree);
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
  DOPT(useConstO2);
  DOPT(useTumorcodeVessels);
  DOPT(paramset_name);
  DOPT(out_intervall);
  DOPT(apply_adaption_intervall);
  DOPT(tend);
  DOPT(time);
  DOPT(latest_executed_timepoint);
  DOPT(lattice_scale);
  DOPT(dt);
  DOPT(message);
  DOPT(max_iteration_per_rerun);
  DOPT(fn_out);
  DOPT(fn_vessel);
  DOPT(vessel_path);
  DOPT(rGf);
  DOPT(gf_production_threshold);
  DOPT(rO2Consumtion);
  DOPT(tumor_radius);
  DOPT(tumor_speed);
  DOPT(stopping_radius_fraction);
  if (tissuePressureDistribution == TISSUE_PRESSURE_SPHERE) pt.put("tissuePressureDistribution", "sphere");
  else if (tissuePressureDistribution == TISSUE_PRESSURE_SHELL) pt.put("tissuePressureDistribution", "shell");
  //DOPT(tissuePressureDistribution);
  DOPT(tissuePressureWidth);
  DOPT(tissuePressureCenterFraction);
  #undef DOPT
//   pt.put_child("calcflow", bfparams.as_ptree());
#if USE_ADAPTION
  pt.put_child("adaption", adap_params.as_ptree());
#endif
  return pt;
}

void FakeTumMTS::Parameters::update_ptree(ptree &dst, const ptree &src)
{
  boost::property_tree::update(dst, src);
}

float FakeTumMTS::FakeTumorSimMTS::estimateTumorRadiusFromCells()
{
  /** for now, we keep it simple, simplifications implied:
   * 1) neglect volume in between cells
   * 2) assume spherical shape
   */ 
  std::vector<double> allCellVolumina = tumorcode_pointer_to_currentCellsSystem->Get_volume();
  double sum = 0;
#pragma omp parallel for reduction(+:sum)
  for(std::size_t i = 0; i<allCellVolumina.size();i++)
  {
    sum+=allCellVolumina[i];
  }
  // estimate radius from Volume ((3*V)/(4*pi))^(1/3): 3/(4*pi) = 0.238732414637843
  float estimated_tum_radius = pow(0.238732414637843 * sum, 0.3333333333333333333333);
  return estimated_tum_radius;
}

//float FakeTumMTS::FakeTumorSimMTS::getGf(const Float3 &pos) const
float FakeTumMTS::FakeTumorSimMTS::getGf(const Float3 &pos)
{
// #define someInteraction
// #ifdef someInteraction // do the super cool new interaction stuff
//   float estimated_tum_radius = estimateTumorRadiusFromCells();
//   float currentShell = pos.norm();
//   return std::max(estimated_tum_radius + params.rGf - currentShell, 0.);
// #else // do it like fake tum has done it before
//   float r = pos.norm();
//   return std::max(tumor_radius + params.rGf - r, 0.);
// #endif
  return FieldInterpolate::ValueAveraged(state.gffield, grid.ld, FieldInterpolate::Const(0.f), pos);
}


float FakeTumMTS::FakeTumorSimMTS::getPress(const Float3 &pos)
{
#ifdef someInteraction
  float estimated_tum_radius = estimateTumorRadiusFromCells();
  float currentShell = pos.norm();
  if (params.tissuePressureDistribution == TISSUE_PRESSURE_SPHERE)
  {
    return my::smooth_heaviside<float>(estimated_tum_radius - currentShell, params.tissuePressureWidth);
  }
  else
  {
    return 2.*params.tissuePressureWidth*my::smooth_delta_cos<float>(estimated_tum_radius - currentShell, params.tissuePressureWidth) +
           params.tissuePressureCenterFraction * my::smooth_heaviside_sin<float>(estimated_tum_radius - currentShell, params.tissuePressureWidth);
  }
#else
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
#endif
}

//float FakeTumMTS::FakeTumorSimMTS::getTumorDens(const Float3 &pos) const
float FakeTumMTS::FakeTumorSimMTS::getTumorDens(const Float3 &pos)
{
#ifdef someInteraction
  float estimated_tum_radius = estimateTumorRadiusFromCells();
  float currentShell = pos.norm();
  return my::smooth_heaviside<float>(estimated_tum_radius - currentShell, 30.);
#else
  float r = pos.norm();
  return my::smooth_heaviside<float>(tumor_radius - r, 30.);
#endif
}

Float3 FakeTumMTS::FakeTumorSimMTS::getGfGrad(const Float3 &pos) const
{
  float r = pos.norm();
  if (r > 0.)
    return (-1./r) * pos;
  else
    return Float3(0.);
}

int FakeTumMTS::FakeTumorSimMTS::run()
{
  // direct cout through log
  cout.rdbuf(my::log().rdbuf());
  
  //read vessels from file
  H5::H5File file;
  H5::Group h5_vessels;
  //H5::Group h5params;
  try{
    file = H5::H5File(params.fn_vessel, H5F_ACC_RDONLY);
    h5_vessels = file.openGroup(params.vessel_path);
  }
  catch(H5::Exception &e)
  {
    cout << "Error reading vessels from hdf" << endl;
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

  // I/O params
  //h5params = file.openGroup("/parameters");
  /*
  try{
    string message;
    readAttrFromH5(h5params, string("MESSAGE"),message);
    params.vesselfile_message = message;
    int index;
    readAttrFromH5(h5params, string("ENSEMBLE_INDEX"),index);
    params.vesselfile_ensemble_index = index;
  }
  catch(H5::Exception &e)
  {
    e.printErrorStack();
  }
  */
  
  /* READ IN DONE */
  
  /* define callback providing information about the simulation */
  VesselModel1::Callbacks callbacks;
  callbacks.getGf = boost::bind(&FakeTumorSimMTS::getGf, boost::ref(*this), _1);
  callbacks.getPress = boost::bind(&FakeTumorSimMTS::getPress, boost::ref(*this), _1);
  callbacks.getTumorDens = boost::bind(&FakeTumorSimMTS::getTumorDens, boost::ref(*this), _1);
  callbacks.getGfGrad = boost::bind(&FakeTumorSimMTS::getGfGrad, boost::ref(*this), _1);

  /* need to compute flow because shearforce must be
    * known and be consistent with current parameters.
    * Shear force is used e.g. in model.Init to initialize
    * f_initial. */
  CalcFlow(*vl, bfparams);
  // initialize vessel model
  vessel_model.Init(vl.get(), this->vessel_model.params, callbacks);
  
  
  
  

  
  double grid_lattice_const = 30;
  double safety_layer_size = 2*grid_lattice_const;
  //boost::optional<Int3> grid_lattice_size;
  /* continum lattice stuff
   * set up grid for calculating diffusion equations
   * needed for solving diffusion equations, here Growthfactors
   */
  //Int3 s = params.lattice_size;
  //int dim = s[2]<=1 ? (s[1]<=1 ? 1 : 2) : 3;
  LatticeDataQuad3d field_ld;
  int dim = 3;
  Bool3 centering = Bool3::mapIndex([=](int i) { return i<dim; });
  
#ifdef DEBUG
  printf("params.lattice_size: %i %i %i\n" , params.lattice_size[0],params.lattice_size[1],params.lattice_size[2]);
  printf("params.lattice_scale: %f\n" , params.lattice_scale);
#endif
  
  SetupFieldLattice(vl->Ld().GetWorldBox(), dim, grid_lattice_const, safety_layer_size, field_ld);
  grid.init(field_ld, dim);
  mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(8, 8, 8)));
  //mtboxes.init(MakeMtBoxGridLarge(grid.Box(), 128));
  
  gf_model.init(field_ld, mtboxes, params.as_ptree());
  gf_model.initField(state.gffield);
  
  o2_uptake_model.init(field_ld,mtboxes, params.as_ptree());
  o2_uptake_model.initField(state.cell_O2_consumption);
  
  last_chem_update = -1;
  last_vessels_checksum = -1;
  
  
  // this is needed to pass tumor information to the DetailedPO2 simulation
  boost::optional<H5::Group> lastTumorGroupWrittenByFakeTum;
  std::string lastTumorGroupWrittenByFakeTumName;
  /* set initial conditions of the FakeTumMTS simulation */
  //tumor_radius = params.tumor_radius; not needed for cell based simulation
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
    
    if( !params.useConstO2 )
    {
      //read detailed o2 from previous run
      H5::Group h5_previous_po2vessels;
      H5::Group h5_previous_tumor;
      H5::DataSet h5_po2field;
      H5::DataSet h5_fieldGf;
      H5::DataSet h5_fieldO2Consumption;
      try
      {
        //file.reOpen();
        file = H5::H5File(params.fn_vessel, H5F_ACC_RDONLY);
        h5_previous_po2vessels = file.openGroup("last_state/po2/vessels");
        h5_previous_tumor = file.openGroup("last_state/tumor");
        h5_po2field = file.openDataSet("last_state/po2/vessels/po2field");
        h5_fieldGf = h5_previous_tumor.openDataSet("fieldGf");
        h5_fieldO2Consumption = h5_previous_tumor.openDataSet("fieldO2Consumption");
      }
      catch(H5::Exception &e)
      {
        cout << "Error reading o2 data from previous simulation" << endl;
        e.printErrorStack();
      }
      
      ReadArray3D<float>(h5_fieldGf, state.gffield);
      ReadArray3D<float>(h5_fieldO2Consumption, state.cell_O2_consumption);
      ReadArray3D<float>(h5_po2field, state.previous_po2field);
      readDataSetFromGroup<Float2>(h5_previous_po2vessels,string("po2vessels"), state.previous_po2vessels);
  #ifndef NDEBUG
      cout<<(*state.previous_po2field)(42,42,42) << endl;
      cout<<(*state.previous_po2vessels)[42] << endl;
  #endif
      
      h5_previous_po2vessels.close();
      h5_previous_tumor.close();
      h5_po2field.close();
      h5_fieldGf.close();
      h5_fieldO2Consumption.close();
      file.close();
      cout<<"reloaded o2 from previous simulation" << endl;
    }
  }
  
  /* init the cell system */
  try{
    initMilotti();
  }
  catch( const std::runtime_error& e)
  {
    std::cerr << "exception in initMilotti: " << e.what() << std::endl;
  }
  
  
  /** 
   * MAIN LOOP
   *  1) calculated vessel po2 and p2 field
   *  2) hand vessels and their po2 to cells simulation
   *  3) OUTPUT
   *  4) run cells simulation
   *  5) remodel vessel tree
   */
  int iteration_in_this_rerun=0;
  while (not PyCheckAbort())
  {
    if (iteration_in_this_rerun >= params.max_iteration_per_rerun)
    {
      std::cout << ">>>>>>> max_iteration_per_rerun reached" << std::endl;
      break;
    }
    if (time > params.tend) 
    {
      std::cout << ">>>>>>> stopped because time limit" << std::endl;
      break;
    }
  
    /* stop if tumor reaches certain fraction of volume */
    double size_limit = 0.5*maxCoeff(Size(vl->Ld().GetWorldBox())) * params.stopping_radius_fraction; 
    //cout << format("size_limit = %f vs tumor_radius = %f\n") % size_limit % tumor_radius;
    if (tumor_radius >  size_limit) 
    {
      std::cout << ">>>>>>> stopped because of size limit" << std::endl;
      break;
    }
    // time =0, next_output_time = 0 at begining
    // 0>=-1 if dt=1 --> true
    // first one is negative in time, so we are preparing
    //maybe we need the output intervalls later again?
    //if (time >= next_output_time - params.dt )
    
    if( ! params.useConstO2)
    {
      /* right now, the o2 simulation is called at every output time */
#ifdef W_timing
      currentTiming.begin_o2 = std::chrono::steady_clock::now();
#endif
      /* O2 stuff is initialized
      * NOTE: size of po2Store is set here
      */
      o2_sim.vl = vl;
      //o2_sim.init(bfparams,grid_lattice_const, safety_layer_size, grid_lattice_size, lastTumorGroupWrittenByFakeTum, state.previous_po2field,state.previous_po2vessels,state.cell_O2_consumption);
      o2_sim.init(bfparams,grid_lattice_const, safety_layer_size, lastTumorGroupWrittenByFakeTum, state.previous_po2field,state.previous_po2vessels,state.cell_O2_consumption);
      cout << "\nInit O2 completed" << endl;
      o2_sim.run();
      cout << "\n mts run finished" << endl;
      //we store the results, to achive quicker convergence in consecutive runs
      state.previous_po2field = o2_sim.getPo2field();
      cout << "\nAccquired po2Field  completed" << endl;
      state.previous_po2vessels = o2_sim.getVesselPO2Storrage();
      cout << "\nDetailed O2 completed" << endl;
#ifdef W_timing
      currentTiming.end_o2 = std::chrono::steady_clock::now();
      currentTiming.time_diff = currentTiming.end_o2 - currentTiming.begin_o2;
      currentTiming.run_o2 = currentTiming.run_o2 + currentTiming.time_diff.count();
#endif
    }
    else 
    {
      state.previous_po2field = boost::none;
      state.previous_po2vessels = boost::none;
    }
    
#ifdef W_timing
    currentTiming.begin_ann = std::chrono::steady_clock::now();
#endif
    findNearestVessel(state.previous_po2vessels);// to have the information for the first output
#ifdef W_timing
    currentTiming.end_ann = std::chrono::steady_clock::now();
    currentTiming.time_diff = currentTiming.end_ann-currentTiming.begin_ann;
    currentTiming.run_ann = currentTiming.run_ann + currentTiming.time_diff.count();
#endif


    
    
    if (time >= next_output_time - params.dt*1.0)
    {
      //this happens only for fixed instances of time
      lastTumorGroupWrittenByFakeTumName = writeOutput(true);//detailedO2 should be calculated prior to this call
      next_output_time += params.out_intervall;
    }
    //for a rerun we need to access the latest instant of time
    params.latest_executed_timepoint = time;
    writeOutput(false);
  
    //lastTumorGroupWrittenByFakeTumName = writeOutput(true);//detailedO2 should be calculated prior to this call
    currentTiming.reset();


    /**
      * milotti vessel structure is initialized with tumorcodes vessels
      */
    /* increment tumor time */
    time += params.dt;
    
    /* propergate cells in time until current fake tumor time */
    cout << boost::format("advance milotti until: %f\n") % time;
    cout.flush();
#ifdef W_timing
    currentTiming.begin_doMilottiStep = std::chrono::steady_clock::now();
#endif
    doMilottiStep();
#ifdef W_timing
    currentTiming.end_doMilottiStep = std::chrono::steady_clock::now();
    currentTiming.time_diff = currentTiming.end_doMilottiStep-currentTiming.begin_doMilottiStep;
    currentTiming.run_doMilottiStep = currentTiming.run_doMilottiStep + currentTiming.time_diff.count();
#endif

    
  
    

    /**
    * do a vessel model remodeling step 
    */
#ifndef NDEBUG
    cout << boost::format("start vessel remodel step! \n");
#endif

#ifdef W_timing
    currentTiming.begin_doStep = std::chrono::steady_clock::now();
#endif
    doStep(params.dt);
#ifdef W_timing
    currentTiming.end_doStep = std::chrono::steady_clock::now();
    currentTiming.time_diff = currentTiming.end_doStep-currentTiming.begin_doStep;
    currentTiming.run_doStep = currentTiming.run_doStep + currentTiming.time_diff.count();
#endif
  
    cout << boost::format("finished vessel remodel step! \n");
    ++output_num;
    ++num_iteration;
    ++iteration_in_this_rerun;
  }

  tumorcode_pointer_to_currentCellsSystem->CloseOutputFiles();						// Closing output files
  return 0;
}

void FakeTumMTS::FakeTumorSimMTS::doStep(double dt)
{
  auto start_here = std::chrono::steady_clock::now();
  cout << format("step %i, t=%f") % num_iteration % time << endl;
  CalcFlow(*vl, bfparams);
  /* do not use adaptation stuff for mts*/
  vessel_model.DoStep(dt, &bfparams);
  auto end_here = std::chrono::steady_clock::now();
  std::cout << "timing 1: " << (end_here-start_here).count() << endl;
  /* this calculates the simple diffusion of substances */
  
  if(isFirstIterationInRerun or num_iteration%21 == 0)
  {
    isFirstIterationInRerun = false;
    calcChemFields();
  }
  auto end_here_2 = std::chrono::steady_clock::now();
  std::cout << "timing 2: " << (end_here_2-end_here).count() << endl;
  
  /* maybe not the best estimate, but first approach */
  tumor_radius = estimateTumorRadiusFromCells();
}


void FakeTumMTS::FakeTumorSimMTS::doMilottiStep()
{
  cout << format("start mts at tumor time: %f\n" ) % time;
  cout.flush();
  //std::cout << "o2_uptake before milotti" << std::endl;
  //std::cout << state.cell_O2_consumption(0,0,0) << std::endl;
  /** 
   * the tumor time is given in hours 
   * evolve cells until that
   */
  uint returnValue = tumorcode_pointer_to_currentCellsSystem->runMainLoop( time * 3600 );
  
  //std::cout << "o2_uptake after milotti" << std::endl;
  //std::cout << state.cell_O2_consumption(0,0,0) << std::endl;
  
  /** 
   * for safety reasons we use both output structures (hdf and vbl)
   * this one the the output of milotti, see WriteCellsSystemHDF for the hdf output
   */
  if ( tumorcode_pointer_to_currentCellsSystem->Get_ready2start() )
  {
    tumorcode_pointer_to_currentCellsSystem->Print2logfile(std::string("Cells at the end of the run"));
    tumorcode_pointer_to_currentCellsSystem->WriteCellsSystem( );					// dump of the final configuration
  }
  
  cout << format("finished FakeTumMTS::FakeTumorSimMTS::doMilottiStep(),  mts at tumor time: %f\n" ) % time;
  //end milotti
}

void FakeTumMTS::FakeTumorSimMTS::writeVBLDataToHDF(H5::Group &h5_vbl)
{
  writeDataSetToGroup(h5_vbl, string("name"), tumorcode_pointer_to_currentCellsSystem->Get_name());
  writeDataSetToGroup(h5_vbl, string("mark"), tumorcode_pointer_to_currentCellsSystem->Get_mark());
  
  writeDataSetToGroup(h5_vbl, string("cell_type"), tumorcode_pointer_to_currentCellsSystem->get_CellTypeIndexVector());
  
  writeDataSetToGroup(h5_vbl, string("phase"), tumorcode_pointer_to_currentCellsSystem->Get_phase_int());
  writeDataSetToGroup(h5_vbl, string("death_condition"), tumorcode_pointer_to_currentCellsSystem->Get_death_condition());
  writeDataSetToGroup(h5_vbl, string("age"), tumorcode_pointer_to_currentCellsSystem->Get_age());
  writeDataSetToGroup(h5_vbl, string("phase_age"), tumorcode_pointer_to_currentCellsSystem->Get_phase_age());
  writeDataSetToGroup(h5_vbl, string("age_mother"), tumorcode_pointer_to_currentCellsSystem->Get_age_mother());
  writeDataSetToGroup(h5_vbl, string("n_mitosis"), tumorcode_pointer_to_currentCellsSystem->Get_n_mitosis());
  
  writeDataSetToGroup(h5_vbl, string("Temperature"), tumorcode_pointer_to_currentCellsSystem->Get_Temperature());
  writeDataSetToGroup(h5_vbl, string("x"), tumorcode_pointer_to_currentCellsSystem->Get_x());
  writeDataSetToGroup(h5_vbl, string("y"), tumorcode_pointer_to_currentCellsSystem->Get_y());
  writeDataSetToGroup(h5_vbl, string("z"), tumorcode_pointer_to_currentCellsSystem->Get_z());
  writeDataSetToGroup(h5_vbl, string("vx"), tumorcode_pointer_to_currentCellsSystem->Get_vx());
  writeDataSetToGroup(h5_vbl, string("vy"), tumorcode_pointer_to_currentCellsSystem->Get_vy());
  writeDataSetToGroup(h5_vbl, string("vz"), tumorcode_pointer_to_currentCellsSystem->Get_vz());
  writeDataSetToGroup(h5_vbl, string("r"), tumorcode_pointer_to_currentCellsSystem->Get_r());
  writeDataSetToGroup(h5_vbl, string("surface"), tumorcode_pointer_to_currentCellsSystem->Get_surface());
  writeDataSetToGroup(h5_vbl, string("volume"), tumorcode_pointer_to_currentCellsSystem->Get_volume());
  writeDataSetToGroup(h5_vbl, string("mass"), tumorcode_pointer_to_currentCellsSystem->Get_mass());
  
  writeDataSetToGroup(h5_vbl, string("volume_extra"), tumorcode_pointer_to_currentCellsSystem->Get_volume_extra());
  
  writeDataSetToGroup(h5_vbl, string("M"), tumorcode_pointer_to_currentCellsSystem->Get_M());
  
  writeDataSetToGroup(h5_vbl, string("G"), tumorcode_pointer_to_currentCellsSystem->Get_G());
  writeDataSetToGroup(h5_vbl, string("G6P"), tumorcode_pointer_to_currentCellsSystem->Get_G6P());
  writeDataSetToGroup(h5_vbl, string("O2"), tumorcode_pointer_to_currentCellsSystem->Get_O2());
  writeDataSetToGroup(h5_vbl, string("store"), tumorcode_pointer_to_currentCellsSystem->Get_store());
  writeDataSetToGroup(h5_vbl, string("A"), tumorcode_pointer_to_currentCellsSystem->Get_A());
  writeDataSetToGroup(h5_vbl, string("AcL"), tumorcode_pointer_to_currentCellsSystem->Get_AcL());
  writeDataSetToGroup(h5_vbl, string("h"), tumorcode_pointer_to_currentCellsSystem->Get_h());
  writeDataSetToGroup(h5_vbl, string("pHi"), tumorcode_pointer_to_currentCellsSystem->Get_pHi());
  writeDataSetToGroup(h5_vbl, string("protein"), tumorcode_pointer_to_currentCellsSystem->Get_protein());
  writeDataSetToGroup(h5_vbl, string("prot_rate"), tumorcode_pointer_to_currentCellsSystem->Get_prot_rate());
  writeDataSetToGroup(h5_vbl, string("DNA"), tumorcode_pointer_to_currentCellsSystem->Get_DNA());
  writeDataSetToGroup(h5_vbl, string("DNA_rate"), tumorcode_pointer_to_currentCellsSystem->Get_DNA_rate());
  
  writeDataSetToGroup(h5_vbl, string("ATP_St"), tumorcode_pointer_to_currentCellsSystem->Get_ATP_St());
  writeDataSetToGroup(h5_vbl, string("ATP_Ox"), tumorcode_pointer_to_currentCellsSystem->Get_ATP_Ox());
  writeDataSetToGroup(h5_vbl, string("ATP_NOx"), tumorcode_pointer_to_currentCellsSystem->Get_ATP_NOx());
  writeDataSetToGroup(h5_vbl, string("ATP2"), tumorcode_pointer_to_currentCellsSystem->Get_ATP2());
  writeDataSetToGroup(h5_vbl, string("ATP3"), tumorcode_pointer_to_currentCellsSystem->Get_ATP3());
  writeDataSetToGroup(h5_vbl, string("ConsATP"), tumorcode_pointer_to_currentCellsSystem->Get_ConsATP());
  writeDataSetToGroup(h5_vbl, string("ConsATP_1"), tumorcode_pointer_to_currentCellsSystem->Get_ConsATP_1());
  writeDataSetToGroup(h5_vbl, string("ConsATP_2"), tumorcode_pointer_to_currentCellsSystem->Get_ConsATP_2());
  writeDataSetToGroup(h5_vbl, string("ConsATP_3"), tumorcode_pointer_to_currentCellsSystem->Get_ConsATP_3());
  writeDataSetToGroup(h5_vbl, string("ConsATP_4"), tumorcode_pointer_to_currentCellsSystem->Get_ConsATP_4());
  writeDataSetToGroup(h5_vbl, string("ConsATP_5"), tumorcode_pointer_to_currentCellsSystem->Get_ConsATP_5());
  writeDataSetToGroup(h5_vbl, string("ATPtot"), tumorcode_pointer_to_currentCellsSystem->Get_ATPtot());
  writeDataSetToGroup(h5_vbl, string("ATPp"), tumorcode_pointer_to_currentCellsSystem->Get_ATPp());
  writeDataSetToGroup(h5_vbl, string("ATPmin"), tumorcode_pointer_to_currentCellsSystem->Get_ATPmin());
  writeDataSetToGroup(h5_vbl, string("ATPstart"), tumorcode_pointer_to_currentCellsSystem->Get_ATPstart());
  writeDataSetToGroup(h5_vbl, string("ATPprod"), tumorcode_pointer_to_currentCellsSystem->Get_ATPprod());
  writeDataSetToGroup(h5_vbl, string("ATPcons"), tumorcode_pointer_to_currentCellsSystem->Get_ATPcons());
  
  writeDataSetToGroup(h5_vbl, string("G_extra"), tumorcode_pointer_to_currentCellsSystem->Get_G_extra());
  writeDataSetToGroup(h5_vbl, string("A_extra"), tumorcode_pointer_to_currentCellsSystem->Get_A_extra());
  writeDataSetToGroup(h5_vbl, string("AcL_extra"), tumorcode_pointer_to_currentCellsSystem->Get_AcL_extra());
  
  writeDataSetToGroup(h5_vbl, string("SensO2"), tumorcode_pointer_to_currentCellsSystem->Get_SensO2());
  writeDataSetToGroup(h5_vbl, string("ConsO"), tumorcode_pointer_to_currentCellsSystem->Get_ConsO());
  
  writeDataSetToGroup(h5_vbl, string("DNA_spread"), tumorcode_pointer_to_currentCellsSystem->Get_DNA_spread());
  writeDataSetToGroup(h5_vbl, string("M_T"), tumorcode_pointer_to_currentCellsSystem->Get_M_T());
  writeDataSetToGroup(h5_vbl, string("pRb"), tumorcode_pointer_to_currentCellsSystem->Get_pRb());
  
  writeDataSetToGroup(h5_vbl, string("ConcS"), tumorcode_pointer_to_currentCellsSystem->Get_ConcS());
  
  writeDataSetToGroup(h5_vbl, string("cyclinD"), tumorcode_pointer_to_currentCellsSystem->Get_cyclinD());
  writeDataSetToGroup(h5_vbl, string("cyclinE"), tumorcode_pointer_to_currentCellsSystem->Get_cyclinE());
  writeDataSetToGroup(h5_vbl, string("cyclinX"), tumorcode_pointer_to_currentCellsSystem->Get_cyclinX());
  writeDataSetToGroup(h5_vbl, string("NpRbk"), tumorcode_pointer_to_currentCellsSystem->Get_NpRbk());
  
  writeDataSetToGroup(h5_vbl, string("isonAS"),tumorcode_pointer_to_currentCellsSystem->Get_isonAS());
  writeDataSetToGroup(h5_vbl, string("isonCH"),tumorcode_pointer_to_currentCellsSystem->Get_isonCH());
  //writeDataSetToGroup(h5_vbl, string("DNA_spread"), tumorcode_pointer_to_currentCellsSystem->Get_DNA_spread());
  std::cout << "finished writeVBLDataToHDF" << std::endl;
}

void FakeTumMTS::FakeTumorSimMTS::readVBLDataFromHDF(H5::Group &h5_vbl)
{
  std::vector<unsigned long> ul_buffer;
  std::vector<int> int_buffer;
  std::vector<double> double_buffer;
  std::vector<float> float_buffer;
  readDataSetFromGroup(h5_vbl, string("cell_type"), ul_buffer);
  tumorcode_pointer_to_currentCellsSystem->set_CellTypeFromIndexVector(ul_buffer);
  readDataSetFromGroup(h5_vbl, string("name"), ul_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_name(ul_buffer);
  readDataSetFromGroup(h5_vbl, string("mark"), int_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_mark(int_buffer);
  readDataSetFromGroup(h5_vbl, string("Temperature"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_Temperature(double_buffer);
  readDataSetFromGroup(h5_vbl, string("phase"), int_buffer);
  tumorcode_pointer_to_currentCellsSystem->set_CellPhaseFromIntVector(int_buffer);
  readDataSetFromGroup(h5_vbl, string("death_condition"), int_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_death_condition(int_buffer);
  readDataSetFromGroup(h5_vbl, string("age"), float_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_age(float_buffer);
  readDataSetFromGroup(h5_vbl, string("phase_age"), float_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_phase_age(float_buffer);
  readDataSetFromGroup(h5_vbl, string("n_mitosis"), int_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_n_mitosis(int_buffer);
  
  readDataSetFromGroup(h5_vbl, string("x"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_x(double_buffer);
  readDataSetFromGroup(h5_vbl, string("y"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_y(double_buffer);
  readDataSetFromGroup(h5_vbl, string("z"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_z(double_buffer);
  readDataSetFromGroup(h5_vbl, string("vx"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_vx(double_buffer);
  readDataSetFromGroup(h5_vbl, string("vy"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_vy(double_buffer);
  readDataSetFromGroup(h5_vbl, string("vz"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_vz(double_buffer);
  readDataSetFromGroup(h5_vbl, string("r"), double_buffer);
//   for(int i=0;i<double_buffer.size();i++)
//   {
//     tumorcode_pointer_to_currentCellsSystem->Set_r(i, double_buffer[i]);
//   }
  tumorcode_pointer_to_currentCellsSystem->Set_r(double_buffer);
  readDataSetFromGroup(h5_vbl, string("surface"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_surface(double_buffer);
  readDataSetFromGroup(h5_vbl, string("volume"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_volume(double_buffer);
//   for(int i=0;i<double_buffer.size();i++)
//   {
//     tumorcode_pointer_to_currentCellsSystem->Set_volume(i, double_buffer[i]);
//   }
  readDataSetFromGroup(h5_vbl, string("mass"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_mass(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("volume_extra"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_volume_extra(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("M"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_M(double_buffer);
//   for(int i=0;i<double_buffer.size();i++)
//   {
//     tumorcode_pointer_to_currentCellsSystem->Set_M(i, double_buffer[i]);
//   }
  readDataSetFromGroup(h5_vbl, string("G"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_G(double_buffer);
  readDataSetFromGroup(h5_vbl, string("G6P"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_G6P(double_buffer);
  readDataSetFromGroup(h5_vbl, string("O2"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_O2(double_buffer);
  readDataSetFromGroup(h5_vbl, string("store"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_store(double_buffer);
  readDataSetFromGroup(h5_vbl, string("A"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_A(double_buffer);
  readDataSetFromGroup(h5_vbl, string("AcL"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_AcL(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("h"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_h(double_buffer);
  readDataSetFromGroup(h5_vbl, string("pHi"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_pHi(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("protein"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_protein(double_buffer);
  readDataSetFromGroup(h5_vbl, string("prot_rate"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_prot_rate(double_buffer);
  readDataSetFromGroup(h5_vbl, string("DNA"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_DNA(double_buffer);
  readDataSetFromGroup(h5_vbl, string("DNA_rate"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_DNA_rate(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("ATP_St"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATP_St(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATP_Ox"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATP_Ox(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATP_NOx"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATP_NOx(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATP2"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATP2(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATP3"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATP3(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ConsATP"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ConsATP(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ConsATP_1"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ConsATP_1(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ConsATP_2"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ConsATP_2(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ConsATP_3"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ConsATP_3(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ConsATP_4"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ConsATP_4(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ConsATP_5"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ConsATP_5(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATPtot"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATPtot(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATPp"), double_buffer);
//   for(int i=0;i<double_buffer.size();i++)
//   {
//     tumorcode_pointer_to_currentCellsSystem->Set_ATPp(i, double_buffer[i]);
//   }
  tumorcode_pointer_to_currentCellsSystem->Set_ATPp(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATPmin"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATPmin(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATPstart"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATPstart(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATPprod"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATPprod(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ATPcons"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ATPcons(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("G_extra"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_G_extra(double_buffer);
  readDataSetFromGroup(h5_vbl, string("A_extra"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_A_extra(double_buffer);
  readDataSetFromGroup(h5_vbl, string("AcL_extra"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_AcL_extra(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("SensO2"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_SensO2(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ConsO"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ConsO(double_buffer);
  readDataSetFromGroup(h5_vbl, string("DNA_spread"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_DNA_spread(double_buffer);
  readDataSetFromGroup(h5_vbl, string("M_T"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_M_T(double_buffer);
  readDataSetFromGroup(h5_vbl, string("pRb"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_pRb(double_buffer);
  readDataSetFromGroup(h5_vbl, string("ConcS"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_ConcS(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("cyclinD"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_cyclinD(double_buffer);
  readDataSetFromGroup(h5_vbl, string("cyclinE"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_cyclinE(double_buffer);
  readDataSetFromGroup(h5_vbl, string("cyclinX"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_cyclinX(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("NpRbk"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_NpRbk(double_buffer);
  
  readDataSetFromGroup(h5_vbl, string("DNA_spread"), double_buffer);
  tumorcode_pointer_to_currentCellsSystem->Set_DNA_spread(double_buffer);
  
  
  
#if 0
  writeDataSetToGroup(h5_vbl, string("mark"), tumorcode_pointer_to_currentCellsSystem->Get_mark());
  
  writeDataSetToGroup(h5_vbl, string("cell_type"), tumorcode_pointer_to_currentCellsSystem->get_CellTypeIndexVector());
  
  writeDataSetToGroup(h5_vbl, string("phase"), tumorcode_pointer_to_currentCellsSystem->Get_phase_int());
  writeDataSetToGroup(h5_vbl, string("death_condition"), tumorcode_pointer_to_currentCellsSystem->Get_death_condition());
  writeDataSetToGroup(h5_vbl, string("age"), tumorcode_pointer_to_currentCellsSystem->Get_age());
  writeDataSetToGroup(h5_vbl, string("phase_age"), tumorcode_pointer_to_currentCellsSystem->Get_phase_age());
  writeDataSetToGroup(h5_vbl, string("age_mother"), tumorcode_pointer_to_currentCellsSystem->Get_age_mother());
  writeDataSetToGroup(h5_vbl, string("n_mitosis"), tumorcode_pointer_to_currentCellsSystem->Get_n_mitosis());
  
  writeDataSetToGroup(h5_vbl, string("Temperature"), tumorcode_pointer_to_currentCellsSystem->Get_Temperature());
  writeDataSetToGroup(h5_vbl, string("x"), tumorcode_pointer_to_currentCellsSystem->Get_x());
  writeDataSetToGroup(h5_vbl, string("y"), tumorcode_pointer_to_currentCellsSystem->Get_y());
  writeDataSetToGroup(h5_vbl, string("z"), tumorcode_pointer_to_currentCellsSystem->Get_z());
  writeDataSetToGroup(h5_vbl, string("vx"), tumorcode_pointer_to_currentCellsSystem->Get_vx());
  writeDataSetToGroup(h5_vbl, string("vy"), tumorcode_pointer_to_currentCellsSystem->Get_vy());
  writeDataSetToGroup(h5_vbl, string("vz"), tumorcode_pointer_to_currentCellsSystem->Get_vz());
  writeDataSetToGroup(h5_vbl, string("r"), tumorcode_pointer_to_currentCellsSystem->Get_r());
  writeDataSetToGroup(h5_vbl, string("surface"), tumorcode_pointer_to_currentCellsSystem->Get_surface());
  writeDataSetToGroup(h5_vbl, string("volume"), tumorcode_pointer_to_currentCellsSystem->Get_volume());
  writeDataSetToGroup(h5_vbl, string("mass"), tumorcode_pointer_to_currentCellsSystem->Get_mass());
  
  writeDataSetToGroup(h5_vbl, string("volume_extra"), tumorcode_pointer_to_currentCellsSystem->Get_volume_extra());
  
  writeDataSetToGroup(h5_vbl, string("M"), tumorcode_pointer_to_currentCellsSystem->Get_M());
  
  writeDataSetToGroup(h5_vbl, string("G"), tumorcode_pointer_to_currentCellsSystem->Get_G());
  writeDataSetToGroup(h5_vbl, string("G6P"), tumorcode_pointer_to_currentCellsSystem->Get_G6P());
  
  writeDataSetToGroup(h5_vbl, string("NpRbk"), tumorcode_pointer_to_currentCellsSystem->Get_NpRbk());
  writeDataSetToGroup(h5_vbl, string("DNA_spread"), tumorcode_pointer_to_currentCellsSystem->Get_DNA_spread());
#endif 
  std::cout << "finished readVBLDataFromHDF" << std::endl;
}

std::string FakeTumMTS::FakeTumorSimMTS::writeOutput(bool doPermanentSafe)
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
  H5::Group root, gout, h5_tum, h5_cells_out, h5_ld_last_state, h5_parameters, h5_vessel_parameters, h5_system_parameters, h5_o2_parameters, h5_calcflow_parameters, h5_field_ld_group, h5_timing, h5_current_vessels, po2outputGroup, ldgroup,
  h5_vbl,h5_memory;
  
  std::string tumOutName = "nothing";
  try
  {
    if( !mySystemParameters.isRerun )
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
      h5_o2_parameters = h5_parameters.createGroup("o2_sim");
      h5_calcflow_parameters = h5_parameters.createGroup("calcflow");
      
      writeAttrToH5(root, string("MESSAGE"), params.message);
      writeAttrToH5(root, string("VESSELTREEFILE"), params.fn_vessel);
      writeAttrToH5(root, string("OUTPUT_NAME"), params.fn_out);
      writeAttrToH5(root, string("VESSELFILE_MESSAGE"), params.vesselfile_message);
      writeAttrToH5(root, string("VESSELFILE_ENSEMBLE_INDEX"), params.vesselfile_ensemble_index);
      WriteHdfPtree(h5_vessel_parameters,vessel_model.params.as_ptree());
      WriteHdfPtree(h5_parameters, params.as_ptree());
      WriteHdfPtree(h5_system_parameters, mySystemParameters.as_ptree());
      WriteHdfPtree(h5_calcflow_parameters, o2_sim.bfparams.as_ptree());
      o2_sim.WriteParametersToHDF(h5_o2_parameters);
      //WriteHdfPtree(h5_o2_parameters, o2_sim.params.as_ptree());
      /* on first occasion, we write field_ld to the root folder */
      h5_field_ld_group = root.createGroup("field_ld");
      grid.ld.WriteHdfLd(h5_field_ld_group);
    }
    catch(H5::Exception &e)
    {
      e.printErrorStack();
    }
  }
  
  try
  {
    if( !doPermanentSafe)
    {
      root.unlink("last_state");
      //the flush deletes the unlinked object immediatelly
      //somehow messes up other stuff
      //root.flush(H5F_SCOPE_LOCAL);
      gout = root.createGroup("last_state");
      writeAttrToH5(gout, "CURRENT_RERUN_NUMBER", mySystemParameters.reRunNumber);
      //data needed for rerun:
      //   lattice, detailed o2, fieldGf, fieldO2Consumption
      
    }
    else
    {
      tumOutName = str(format("out%04i") % output_num);
      gout = root.createGroup(tumOutName);
    }
    
    h5_vbl = gout.createGroup("vbl");
    h5_current_vessels = gout.createGroup("vessels");
    h5_tum = gout.createGroup("tumor");
    
    if( doPermanentSafe)
    {
      /* write timing */
      h5_timing = gout.createGroup("timing");
      ///// tumorcode
      writeAttrToH5(h5_timing, string("run_init_o2"), currentTiming.run_init_o2);
      writeAttrToH5(h5_timing, string("run_o2"), currentTiming.run_o2);
      writeAttrToH5(h5_timing, string("run_ann"), currentTiming.run_ann);
      writeAttrToH5(h5_timing, string("run_doStep"), currentTiming.run_doStep);
      writeAttrToH5(h5_timing, string("run_doMilottiStep"), currentTiming.run_doMilottiStep);
      ///// vbl
      writeAttrToH5(h5_timing, string("run_vbl_diff"), tumorcode_pointer_to_currentCellsSystem->myTiming.diff);
      writeAttrToH5(h5_timing, string("run_vbl_diff_loop_1"), tumorcode_pointer_to_currentCellsSystem->myTiming.diff_loop_1);
      writeAttrToH5(h5_timing, string("run_vbl_diff_loop_2"), tumorcode_pointer_to_currentCellsSystem->myTiming.diff_loop_2);
      writeAttrToH5(h5_timing, string("run_vbl_diff_loop_3"), tumorcode_pointer_to_currentCellsSystem->myTiming.diff_loop_3);
      writeAttrToH5(h5_timing, string("run_vbl_dynamics"), tumorcode_pointer_to_currentCellsSystem->myTiming.dynamics);
      writeAttrToH5(h5_timing, string("run_vbl_geometry"), tumorcode_pointer_to_currentCellsSystem->myTiming.geometry);
      writeAttrToH5(h5_timing, string("run_vbl_cellEvents"), tumorcode_pointer_to_currentCellsSystem->myTiming.cellEvents);
      writeAttrToH5(h5_timing, string("run_vbl_writeToFile"), tumorcode_pointer_to_currentCellsSystem->myTiming.writeToFile);
      writeAttrToH5(h5_timing, string("run_vbl_bico_call"), tumorcode_pointer_to_currentCellsSystem->myTiming.bico_call);
      writeAttrToH5(h5_timing, string("geometry_neighborhood"), tumorcode_pointer_to_currentCellsSystem->myTiming.geometry_neighborhood);
      ///// global timing
      const auto now = std::chrono::system_clock::now();
      const auto epoch   = now.time_since_epoch();
      const auto seconds = std::chrono::duration_cast<std::chrono::seconds>(epoch);
      writeAttrToH5(h5_timing, string("secondsSinceEpoch"), (int)seconds.count());
      /* write memory */
      h5_memory = gout.createGroup("memory");
      MemUsage m = GetMemoryUsage_();
      writeAttrToH5(h5_memory, string("vmem"), m.vmem);
      writeAttrToH5(h5_memory, string("vmem_peak"), m.vmem_peak);
      writeAttrToH5(h5_memory, string("rss"), m.rss);
      writeAttrToH5(h5_memory, string("rss_peak"), m.rss_peak);
    }
    
    
    
    /** runtime data
     */
    writeAttrToH5(gout, string("time"), time);
    writeAttrToH5(gout, string("OUTPUT_NUM"), output_num);
    writeAttrToH5(gout, string("NUM_ITERATION"), num_iteration);
    writeAttrToH5(gout, string("NEXT_OUTPUT_TIME"), next_output_time);
    writeAttrToH5(gout, string("NEXT_ADAPTION_TIME"), next_adaption_time);
    
    /* writes the vessel list */
    WriteVesselList3d(*vl, h5_current_vessels);
    
    /* tumor information, might change in future? */
    writeAttrToH5(h5_tum, string("TYPE"), string("faketumor"));
    writeAttrToH5(h5_tum, string("TUMOR_RADIUS"), tumor_radius);
    // could be done, but since it is a sphere, you can easily calculate the tc_density from the radius
    /** BIG DATA HERE 
     * so we output this only every 5.th time step 
     */
    if( ( output_num > 0 and (output_num % 5 == 0) ) or !doPermanentSafe)
    {
      /** VBL 
      */
      WriteHdfPtree(h5_vbl, tumorcode_pointer_to_currentCellsSystem->get_params_pointer()->as_ptree());
      
      boost::property_tree::ptree return_from_vbl = tumorcode_pointer_to_currentCellsSystem->as_ptree();
      for( auto it:return_from_vbl)
      {
  #ifndef NDEBUG
	cout << it.first << endl;
  #endif
	H5::Group h5_vbl_type = h5_vbl.createGroup(it.first);
	WriteHdfPtree(h5_vbl_type, return_from_vbl.get_child(it.first));
      }
      writeVBLDataToHDF(h5_vbl);
      /** END VBL 
      */
      WriteScalarField(h5_tum, string("fieldGf"), state.gffield, grid.ld, root.openGroup("field_ld"));
      WriteScalarField(h5_tum, string("fieldO2Consumption"), state.cell_O2_consumption, grid.ld, root.openGroup("field_ld"));
    }
    /* needs to calculate, before output! */
    //UpdateVesselVolumeFraction();
    //WriteScalarField(h5_tum, "vessel_volume_fraction", vessel_volume_fraction, grid.ld, root.openGroup("field_ld"));
    if (tumorcode_pointer_to_currentCellsSystem->Get_alive()>0)
    {
      h5_cells_out = gout.createGroup("cells");
      WriteCellsSystemHDF_with_nearest_vessel_index(h5_cells_out);
      h5_cells_out.close();
    }
    /* write oxygen stuff */
    if(!params.useConstO2)
    {
      // copied from python-oxygen2.cpp
      // write the detailed o2 stuff necessary for the cells
      po2outputGroup = gout.createGroup("po2");
      o2_sim.WriteDataOutput(po2outputGroup);
      //ldgroup = po2outputGroup.createGroup("field_ld");
      //o2_sim.grid.ld.WriteHdfLd(ldgroup);
      //WriteScalarField<float>(po2outputGroup, string("po2field"), o2_sim.po2field, o2_sim.grid.ld, ldgroup);
      //writeDataSetToGroup(po2outputGroup, string("po2vessels"),o2_sim.po2vessels);
      writeAttrToH5(po2outputGroup, string("simType"), string("MTS"));
      //WriteHdfPtree(po2outputGroup, o2_sim.metadata, HDF_WRITE_PTREE_AS_ATTRIBUTE);
      po2outputGroup.close();
    }
  }
  catch(H5::Exception &e)
  {
    e.printErrorStack();
  }
  
  root.close();
  gout.close();
  h5_tum.close();
  
  h5_ld_last_state.close();
  h5_parameters.close();
  h5_vessel_parameters.close();
  h5_system_parameters.close();
  h5_o2_parameters.close();
  h5_field_ld_group.close();
  h5_timing.close();
  h5_current_vessels.close();
  
  ldgroup.close();
  h5_vbl.close();
  h5_memory.close();
  
  //f_out.flush(H5F_SCOPE_LOCAL);
  f_out.close();
  //++output_num;
  
  cout << format("files %s flushed and closed")  % params.fn_out << endl;
  return tumOutName;
}

//This checks if value is in the range [low, high)
template <typename T>
    bool IsInBounds(const T& value, const T& low, const T& high) {
    return !(value < low) && (value < high);
}

//This checks if value is in the range [low, high]
// template <typename T>
//     bool IsInBounds(const T& value, const T& low, const T& high) {
//     return !(value < low) && !(high < value);
// }


/** 
 *          - origin
 *      VBL  ----> tumorcode
 *  tumorcode <---- VBL 
 *              + origin
 */
void FakeTumMTS::FakeTumorSimMTS::calcChemFields()
{
  /** needs to calculate, before output-> any field operation!
   *  initializes cell_GFsrc and cell_O2src 
   */
  UpdateVesselVolumeFraction();
  //*********** gf ****************
  
  my::log().push("gf:");
  
  /** identify boxes containing cells
   */
  Float3 origin_of_tum = grid.ld.GetOriginPosition();
  //cout << "origin_of_tum:\t" << origin_of_tum << endl;
  
  const long ncells = tumorcode_pointer_to_currentCellsSystem->Get_ncells();

  float max_cell_distance_to_seeding = 0.0;
  std::array<float, 3> seeding_pos_in_vbl_coordinates = tumorcode_pointer_to_currentCellsSystem->get_seeding_position();
  std::array<float, 3> seeding_pos_in_vessel_coordinates= tumorcode_pointer_to_currentCellsSystem->get_seeding_position();
  for(int i=0;i<3;i++)
  {
    seeding_pos_in_vessel_coordinates[i] = seeding_pos_in_vbl_coordinates[i]-origin_of_tum[i];
#ifndef NDEBUG
    cout << "seeding_pos_in_vbl_coordinates " << i << ":\t" << seeding_pos_in_vbl_coordinates[i] << endl;
    cout << "seeding_pos_in_vessel_coordinates " << i << ":\t" << seeding_pos_in_vessel_coordinates[i] << endl;
#endif
  }
  
  //seeding_pos[0]=seeding_pos[1]=seeding_pos[2] = 0.0;
#pragma omp parallel for reduction( max: max_cell_distance_to_seeding)
  for(uint k=0;k<ncells; ++k)
  {
    float x = tumorcode_pointer_to_currentCellsSystem->Get_x(k);
    x= x-seeding_pos_in_vbl_coordinates[0];
    x=x*x;
    float y = tumorcode_pointer_to_currentCellsSystem->Get_y(k);
    y= y-seeding_pos_in_vbl_coordinates[1];
    y=y*y;
    float z = tumorcode_pointer_to_currentCellsSystem->Get_z(k);
    z= z-seeding_pos_in_vbl_coordinates[2];
    z=z*z;
    float distance_this_cell = x+y+z;
    
    if( distance_this_cell > max_cell_distance_to_seeding)
    {
      max_cell_distance_to_seeding = distance_this_cell;
    }
  }
  max_cell_distance_to_seeding = std::sqrt(max_cell_distance_to_seeding);
  //DynArray<DomainDecomposition::ThreadBox> thread_boxes_with_cells;
  std::vector<DomainDecomposition::ThreadBox> thread_boxes_with_cells;
  
#if 0  // this is the maximum norm case!
  float max_cell_x = tumorcode_pointer_to_currentCellsSystem->Get_max_x()-origin_of_tum[0];
  float max_cell_y = tumorcode_pointer_to_currentCellsSystem->Get_max_y()-origin_of_tum[1];
  float max_cell_z = tumorcode_pointer_to_currentCellsSystem->Get_max_z()-origin_of_tum[2];
  float min_cell_x = tumorcode_pointer_to_currentCellsSystem->Get_min_x()-origin_of_tum[0];
  float min_cell_y = tumorcode_pointer_to_currentCellsSystem->Get_min_y()-origin_of_tum[1];
  float min_cell_z = tumorcode_pointer_to_currentCellsSystem->Get_min_z()-origin_of_tum[2];
  
#if 0
  //an example
  float extention=30.;
  max_cell_x = extention-origin[0];
  min_cell_x =  -1.*extention -origin[0];
  max_cell_y = extention-origin[1];
  min_cell_y =  -1.*extention-origin[1];
  max_cell_z = extention-origin[2];
  min_cell_z =  -1.*extention -origin[2];
#endif
  
#ifndef NDEBUG
  cout  << "x: [ " << min_cell_x << "," << max_cell_x << " ] \t"
        << "y: [ " << min_cell_y << "," << max_cell_y << " ] \t"
        << "z: [ " << min_cell_z << "," << max_cell_z << " ]" << endl;
#endif
  
  //push back is not thread safe!
  //#pragma omp parallel
  //{
    BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, mtboxes.boxes)
    {
      float multi_box_x_min = bbox[0][0]*grid.ld.Scale();
      float multi_box_x_max = bbox[1][0]*grid.ld.Scale();
      float multi_box_y_min = bbox[0][1]*grid.ld.Scale();
      float multi_box_y_max = bbox[1][1]*grid.ld.Scale();
      float multi_box_z_min = bbox[0][2]*grid.ld.Scale();
      float multi_box_z_max = bbox[1][2]*grid.ld.Scale();
      
      /** in principal all those cases can overlap in 
       * each dimension, but this would be a mess 
       */
      if
	(
	  (
          IsInBounds(min_cell_x, multi_box_x_min, multi_box_x_max)
          &&
          IsInBounds(min_cell_y, multi_box_y_min,multi_box_y_max)
          &&
          IsInBounds(min_cell_z,multi_box_z_min,multi_box_z_max)
	  )
	or
	  (
          IsInBounds(max_cell_x, multi_box_x_min, multi_box_x_max)
          &&
          IsInBounds(max_cell_y, multi_box_y_min,multi_box_y_max)
          &&
          IsInBounds(max_cell_z,multi_box_z_min,multi_box_z_max)
	  )
	or
	  (
          IsInBounds(multi_box_x_min,min_cell_x, max_cell_x)
          &&
          IsInBounds(multi_box_y_min, min_cell_y, max_cell_y)
          &&
          IsInBounds(multi_box_z_min,min_cell_z, max_cell_z)
	  )
	or
	  (
          IsInBounds(multi_box_x_max,min_cell_x, max_cell_x)
          &&
          IsInBounds(multi_box_y_max, min_cell_y, max_cell_y)
          &&
          IsInBounds(multi_box_z_max,min_cell_z, max_cell_z)
	  )
	)
        {
#ifndef NDEBUG
          cout << "bbox: " << bbox << endl;
          cout << "bbox[0]: " << bbox[0] << endl;
          cout << "bbox[1]: " << bbox[1] << endl;
          cout << "multi_box_x_min" << multi_box_x_min<< endl;
          cout << "multi_box_x_max" << multi_box_x_max<< endl;
          cout << "multi_box_y_min" << multi_box_y_min<< endl;
          cout << "multi_box_y_max" << multi_box_y_max<< endl;
          cout << "multi_box_z_min" << multi_box_z_min<< endl;
          cout << "multi_box_z_max" << multi_box_z_max<< endl;
#endif
          thread_boxes_with_cells.push_back(bbox);
        }
    }
  //}

  cout << "found:\t " << thread_boxes_with_cells.size() << " boxes." << endl;
  thread_boxes_with_cells.erase(std::unique( thread_boxes_with_cells.begin(), thread_boxes_with_cells.end() ), thread_boxes_with_cells.end() );
  cout << "found:\t " << thread_boxes_with_cells.size() << " unique boxes." << endl;
  for(int i=0;i<thread_boxes_with_cells.size();i++)
  {
    cout << "i: " << i << " " << thread_boxes_with_cells[i] << endl;
  }
  
#else // this is the radial case
  
    auto bbox = mtboxes.boxes[0];
    //cout <<"found mtboxes" << mtboxes.boxes.size() << " boxes" << endl;
    float multi_box_x_length = (bbox[1][0]-bbox[0][0])*grid.ld.Scale();
    float multi_box_y_length = (bbox[1][1]-bbox[0][1])*grid.ld.Scale();
    float multi_box_z_length = (bbox[1][2]-bbox[0][2])*grid.ld.Scale();
    float diagonal_of_multi = multi_box_x_length*multi_box_x_length+multi_box_y_length*multi_box_y_length+multi_box_z_length*multi_box_z_length;
    diagonal_of_multi = std::sqrt(diagonal_of_multi);
#ifndef NDEBUG
    cout << "multi_box_x_length" << multi_box_x_length<< endl;
    cout << "multi_box_y_length" << multi_box_y_length<< endl;
    cout << "multi_box_z_length" << multi_box_z_length<< endl;
    cout << "diagonal_of_multi" << diagonal_of_multi << endl;
    cout << "max_cell_distance_to_seeding:\t" << max_cell_distance_to_seeding << endl;
#endif

    BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, mtboxes.boxes)
    {
      //get center of multibox
      float multi_box_x_center = (bbox[0][0]*grid.ld.Scale()+0.5*grid.ld.Scale())+origin_of_tum[0];//in vbl coordinates
      float multi_box_y_center = (bbox[0][1]*grid.ld.Scale()+0.5*grid.ld.Scale())+origin_of_tum[1];//in vbl coordinates
      float multi_box_z_center = (bbox[0][2]*grid.ld.Scale()+0.5*grid.ld.Scale())+origin_of_tum[2];//in vbl coordinates
      //get distance from a multibox center to the vbl seeding point
      float radial_distance_of_multibox_to_seeding_point =
                    (multi_box_x_center-seeding_pos_in_vbl_coordinates[0])*(multi_box_x_center-seeding_pos_in_vbl_coordinates[0])+
                    (multi_box_y_center-seeding_pos_in_vbl_coordinates[1])*(multi_box_y_center-seeding_pos_in_vbl_coordinates[1])+
                    (multi_box_z_center-seeding_pos_in_vbl_coordinates[2])*(multi_box_z_center-seeding_pos_in_vbl_coordinates[2]);
        
      radial_distance_of_multibox_to_seeding_point = std::sqrt(radial_distance_of_multibox_to_seeding_point);

#if 0
      cout << "radial_distance_of_multibox_to_seeding_point: \t" << radial_distance_of_multibox_to_seeding_point << endl;
      cout << "multi_box_x_center: " << multi_box_x_center<< endl;
      cout << "multi_box_y_center: " << multi_box_y_center<< endl;
      cout << "multi_box_z_center: " << multi_box_z_center<< endl;
#endif
      //if the discrete box is within the limit of the maximal cell extention, we add the box to the list
      if ( IsInBounds(radial_distance_of_multibox_to_seeding_point, (float)(max_cell_distance_to_seeding-0.5*diagonal_of_multi) ,(float)(max_cell_distance_to_seeding+0.5*diagonal_of_multi)))
        {
          
#ifndef NDEBUG
          
          cout << "bbox: " << bbox << endl;
          cout << "bbox[0]: " << bbox[0] << endl;
          cout << "bbox[1]: " << bbox[1] << endl;
//           cout << "multi_box_x_min" << multi_box_x_min<< endl;
//           cout << "multi_box_x_max" << multi_box_x_max<< endl;
//           cout << "multi_box_y_min" << multi_box_y_min<< endl;
//           cout << "multi_box_y_max" << multi_box_y_max<< endl;
//           cout << "multi_box_z_min" << multi_box_z_min<< endl;
//           cout << "multi_box_z_max" << multi_box_z_max<< endl;
#endif

          thread_boxes_with_cells.push_back(bbox);
        }
    }
  //}

  cout << "found:\t " << thread_boxes_with_cells.size() << " boxes." << endl;
  thread_boxes_with_cells.erase(std::unique( thread_boxes_with_cells.begin(), thread_boxes_with_cells.end() ), thread_boxes_with_cells.end() );
  cout << "found:\t " << thread_boxes_with_cells.size() << " unique boxes." << endl;
  for(int i=0;i<thread_boxes_with_cells.size();i++)
  {
    cout << "i: " << i << " " << thread_boxes_with_cells[i] << endl;
  }
#endif
  
  
  //from continuum-grid.cpp at 263
//   typedef boost::iterator_range<DynArray<DomainDecomposition::ThreadBox>::const_iterator> range_type;
//   boost::array<DynArray<ThreadBox>, 48> by_thread;
//   const DynArray<ThreadBox> &ar = by_thread[omp_get_thread_num()];
  //return range_type(ar.begin(), ar.end());
  
  /** if   loop over cell 
   *            loop over boxes
   * than 
   * when two cell live in the save box--> multiple thread access the same memory
   */
#if 0
  #pragma omp parallel
  {
    //BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, mtboxes.getCurrentThreadRange())
    BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, thread_boxes_with_cells)
    {
      //cout << omp_get_thread_num() << " / " << omp_get_max_threads() << endl;
      for(int i=0; i<tumorcode_pointer_to_currentCellsSystem->Get_ncells();++i)
      {
        Float3 pos(tumorcode_pointer_to_currentCellsSystem->Get_x()[i],tumorcode_pointer_to_currentCellsSystem->Get_y()[i],tumorcode_pointer_to_currentCellsSystem->Get_z()[i]);
        AddSmoothDelta(cell_GFsrc, bbox, grid.ld, grid.dim, pos, (float)(tumorcode_pointer_to_currentCellsSystem->Get_r()[i]/tumorcode_pointer_to_currentCellsSystem->Get_ncells()));
        
        auto this_o2_rate = tumorcode_pointer_to_currentCellsSystem->Get_O2Rate()[i];
        AddSmoothDelta(cell_O2src, bbox, grid.ld, grid.dim, pos, (float) this_o2_rate );
      }
    }
  }
#else 
  for(int i=0; i<tumorcode_pointer_to_currentCellsSystem->Get_ncells();++i)
  {
    //cout << omp_get_thread_num() << " / " << omp_get_max_threads() << endl;
    BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, thread_boxes_with_cells)
    {
      Float3 pos(tumorcode_pointer_to_currentCellsSystem->Get_x()[i],tumorcode_pointer_to_currentCellsSystem->Get_y()[i],tumorcode_pointer_to_currentCellsSystem->Get_z()[i]);
      AddSmoothDelta(cell_GFsrc, bbox, grid.ld, grid.dim, pos, (float)(tumorcode_pointer_to_currentCellsSystem->Get_r()[i]/tumorcode_pointer_to_currentCellsSystem->Get_ncells()));
      
      auto this_o2_rate = tumorcode_pointer_to_currentCellsSystem->Get_O2Rate()[i];
      AddSmoothDelta(cell_O2src, bbox, grid.ld, grid.dim, pos, (float) this_o2_rate );
    }
  }
  
#endif
  gf_model.update(state.gffield, cell_GFsrc);
  o2_uptake_model.update(state.cell_O2_consumption, cell_O2src);
  
#ifdef DEBUG
  cout << "stats: " << state.gffield.valueStatistics() << endl;
#endif
  my::log().pop();
  state.chem_checksum++;
}


/* this has to be tidyed up before proceeding */
void FakeTumMTS::FakeTumorSimMTS::UpdateVesselVolumeFraction()
{
  my::LogScope log_push_(my::log(), "ves:");
  
  // update only when the vessel system state has changed since last update 
  if (!vessel_volume_fraction.empty() && last_vessels_checksum==state.vessels_checksum)
    return;
#ifndef TOTAL_SILENCE
  cout << "volume fraction and o2 sources update!" << endl;
#endif
  
  VesselsInBoxes vessboxes;
  SortVesselsIntoMtBoxGrid(grid.ld, *vl, 2, mtboxes, vessboxes);
  
  // reinit and fill the array
  Array3dOps<float> ops(mtboxes, grid.Box(), grid.dim, 0);
  ops.init(vessel_volume_fraction, true);
//   ops.init(vessel_o2src_clin);
//   ops.init(vessel_o2src_crhs);
//   ops.init(vessel_glucosesrc_clin);
//   ops.init(vessel_glucosesrc_crhs);
  ops.init(cell_GFsrc);
  ops.init(cell_O2src);

#pragma omp parallel
  {
    VesselVolumeGenerator volumegen(*vl, grid.ld, grid.dim, make_ptree("samples_per_cell", 100));
    
    BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, mtboxes.getCurrentThreadRange())
    {
      int dummy_;
      volumegen.Fill(bbox, vessel_volume_fraction, vessboxes[bbox.global_index], dummy_);
    }
  }// end #pragma omp parallel
  last_vessels_checksum = state.vessels_checksum;
}

#if 0
void FakeTumMTS::FakeTumorSimMTS::insertGFCoefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder &mb)
{
  Array3d<float> l_coeff(bbox),
                 rhs(bbox);
                 
  //l_coeff[bbox] += cell_GFsrc_clin[bbox];
  //rhs[bbox] += cell_GFsrc_crhs[bbox];
  
  

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
  mb.AddDiffusion<> (bbox, ConstValueFunctor<float>(1.), -0.005);
  
  FOR_BBOX3(p, bbox)
  {
//     for(int i = 0;i<currentCellsSystem->Get_x().size();++i)
//     {
//       if(sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])<currentCellsSystem->Get_r()[i])
//       {
//         l_coeff(p) = 42;
//       }
//     }
#if 0
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
#endif
    //mb.AddLocally(p, -loc_l_coeff - l_coeff(p), -rhs(p)); /// this was it before
    //mb.AddLocally(p, - 42*sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]) , -rhs(p)); /// this was it before
    mb.AddLocally(p, -l_coeff(p) , -rhs(p));
  }
}
#endif

void FakeTumMTS::FakeTumorSimMTS::WriteCellsSystemHDF_with_nearest_vessel_index( H5::Group &out_cell_group)
{
#ifndef NDEBUG
  cout<< "going to write cells to a hdf file" << endl;
#endif
  int numberOfCells = tumorcode_pointer_to_currentCellsSystem->Get_ncells();
  std::vector<double> x = tumorcode_pointer_to_currentCellsSystem->Get_x();
  std::vector<double> y = tumorcode_pointer_to_currentCellsSystem->Get_y();
  std::vector<double> z = tumorcode_pointer_to_currentCellsSystem->Get_z();
  DynArray<Float3> a(numberOfCells);
  DynArray<int> index_of_nearest_vessel(numberOfCells);
  DynArray<float> min_distances(numberOfCells);
  for( int i = 0; i<numberOfCells ;++i)
  {
    a[i] = Float3(x[i],y[i],z[i]);
    nearest theNearest = vectorOfnearestVessels[i];
    index_of_nearest_vessel[i] = theNearest.indexOfVessel;
    min_distances[i] = theNearest.distance;
  }
  writeDataSetToGroup(out_cell_group, string("cell_center_pos"), a);
  writeDataSetToGroup(out_cell_group, string("distance_to_nearest_vessel"), min_distances);
  writeDataSetToGroup(out_cell_group, string("index_of_nearest_vessel"), index_of_nearest_vessel);
  
  DynArray<float> buffer(numberOfCells);
  DynArray<int> buffer_int(numberOfCells);
  
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = tumorcode_pointer_to_currentCellsSystem->Get_O2Rate()[i];
  }
  writeDataSetToGroup(out_cell_group, string("cell_o2_consumption_rate"), buffer);
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = tumorcode_pointer_to_currentCellsSystem->Get_r()[i];
  }
  writeDataSetToGroup(out_cell_group, string("cell_radii"), buffer);
  // glucose extracellular Get_G_extra()
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = tumorcode_pointer_to_currentCellsSystem->Get_G_extra()[i];
  }
  writeDataSetToGroup(out_cell_group, string("glucose_ex"), buffer);
  
  // ph Get_pH
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = tumorcode_pointer_to_currentCellsSystem->Get_pH()[i];
  }
  writeDataSetToGroup(out_cell_group, string("pH_ex"), buffer);
  
  // oxygen Get_O2
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = tumorcode_pointer_to_currentCellsSystem->Get_O2()[i];
  }
  writeDataSetToGroup(out_cell_group, string("o2"), buffer);
  
  // lactate Get_AcL_extra
  // not needed, could be calculated from pH
//   for( int i = 0; i<numberOfCells; ++i)
//   {
//     buffer[i] = tumorcode_pointer_to_currentCellsSystem->Get_AcL_extra()[i];
//   }
//   writeDataSetToGroup(out_cell_group, string("AcL_ex"), buffer);
  
  // cell phase
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer_int[i] = tumorcode_pointer_to_currentCellsSystem->Get_phase(i);
  }
  writeDataSetToGroup(out_cell_group, string("cell_phase"), buffer_int);
  
  // cell age Get_age
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = tumorcode_pointer_to_currentCellsSystem->Get_age()[i];
  }
  writeDataSetToGroup(out_cell_group, string("cell_age"), buffer);
  
  // cell phage_age Get_phase_age
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = tumorcode_pointer_to_currentCellsSystem->Get_phase_age()[i];
  }
  writeDataSetToGroup(out_cell_group, string("cell_phase_age"), buffer);
  
  // number of neighbours
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer_int[i] = tumorcode_pointer_to_currentCellsSystem->Get_neigh(i);
  }
  writeDataSetToGroup(out_cell_group, string("cell_no_neigh"), buffer_int);
  cout<< "finished writting cells to hdf" << endl;
}


double calculate_distance_from_vessel_to_point_in_space( const Float3 &a, const Float3 &b, ANNpoint &spacePoint)
{
  double tp = 0;
  double nrm = 0;
  double bvd = 0;
  Float3 buffer;
  std::array<double,3> x0 = {0.0,0.0,0.0};
  
  for(int i=0; i<3; i++) // calcolo di t0
  {
      tp += (spacePoint[i]-a[i])*(b[i]-a[i]);
      nrm += (b[i]-a[i])*(b[i]-a[i]);
  }
  tp /= nrm;
  
  // cout << "nrm: " << nrm << endl;
  // cout << "tp: " << tp << endl;
  
  
  if(tp >= 0 && tp <= 1) // Calculating the distance from the cylindrical part
      for(int k=0; k<3; k++)
      {
          x0[k] = tp*(b[k]-a[k]) + a[k];
          bvd += (spacePoint[k] - x0[k])*(spacePoint[k] - x0[k]);
      }
  else if(tp<0)   // Calculating the distance from the spherical shell at the beginning
      for(int k=0; k<3; k++)
      {
          x0[k] = a[k];
          bvd += ( spacePoint[k]-a[k])*( spacePoint[k]-a[k]);
      }
  else if(tp>1)   // Calculating the distance from the spherical shell to the end
      for(int k=0; k<3; k++)
      {
          x0[k] = b[k];
          bvd += ( spacePoint[k]-b[k] )*( spacePoint[k]-b[k] );
      }
  
  bvd = sqrt(bvd);
  
  return(bvd);
}
/** @brief convert units between vbl and tumorcode
 */
float to_vbl_o2_units(float pressure_in_mmhg)
{
  //henrys constant, see milotti paper
  float k_h = 1.3e-3;
  // need pressure in atm
  return pressure_in_mmhg * 0.0013157895567935 * k_h;
}

void FakeTumMTS::FakeTumorSimMTS::findNearestVessel( boost::optional<DetailedPO2::VesselPO2Storage> previous_po2vessels)
{
  /* we will fill the ann structure contnouslsy, so we need that map to restore the 
   * index within the the vessel list
   */
  std::map <uint, uint> ann_to_vl;
  //obviously there is no need for this (so far)
  //std::map <uint, uint> vl_to_ann;
  
  int ecnt      =vl->GetECount();// actual number of data points
  int nPts      = 0; //actual number of data points required for the kdtree
  myAssert(nPts<ANN_maxPts);
  
  //to exploit the full power of ANN, we need to fix the memory and there need to know how many vessels we will consider
  // this could NOT be done in parallel since nPts is a atom variable
  for(int i=0; i<ecnt; ++i)
  {
    Vessel* v = vl->GetEdge(i);
    if( not v->IsCirculated() )
    {
      continue;
    }
    ann_to_vl[nPts]=i;
    //vl_to_ann[i]=nPts;
    nPts++;
  }
  
  // get positions of the cells
  int numberOfCells = tumorcode_pointer_to_currentCellsSystem->Get_ncells();
  cout << "number of cells : " << numberOfCells << endl;
  std::vector<double> x = tumorcode_pointer_to_currentCellsSystem->Get_x();
  std::vector<double> y = tumorcode_pointer_to_currentCellsSystem->Get_y();
  std::vector<double> z = tumorcode_pointer_to_currentCellsSystem->Get_z();
  
  // resultes are stored on the heap by this vector
  vectorOfnearestVessels.resize(numberOfCells);
  
  /* this could be done simultaneously for every cell since 
   * the size of vectorOfnearestVessels is already determinded
   */

  //no, this does not work
  //http://forum.openmp.org/forum/viewtopic.php?f=3&t=757
// #pragma omp parallel
// {
  ANNpointArray    dataPts;         // data points
  ANNpoint         queryPt;         // query point
  queryPt = annAllocPt(ANN_dim);                  // allocate query point
  dataPts = annAllocPts(nPts, ANN_dim);     // allocate data points
  
    // fill kd tree structure
  // we run with the ann structure since the non circulated vessels are gone there
// #pragma for
  for(int i=0; i<nPts; ++i)
  {
    Vessel* v = vl->GetEdge(ann_to_vl[i]);            //get vessel per thread
    Float3 p_a,p_b,a_b,vessel_center;                 //allocate in loop to make sure we are threadsafe
    p_a = vl->Ld().LatticeToWorld(v->NodeA()->lpos);  //read position a
    p_b = vl->Ld().LatticeToWorld(v->NodeB()->lpos);  //read position b
    a_b = p_a-p_b;                                    // difference vector
    vessel_center = p_b+0.5*a_b;                      // vessel center
    for(int k=0;k<ANN_dim;++k)
    {
      dataPts[i][k] = vessel_center[k];
    }
  }// end #pragma
  
  //note the example from the website also works with pointers
  ANNkd_tree *kd_tree_of_vl;        // ann kd tree structurs
  // build search structure with complete vessel tree
  kd_tree_of_vl = new ANNkd_tree(
                                  dataPts,				// the data points
                                  nPts,					// number of points
                                  ANN_dim);				// dimension of space
  //prepare search
  ANNidxArray ANN_nnIdx = new ANNidx[ANN_k];					// allocate near neigh indices
  ANNdistArray ANN_dists = new ANNdist[ANN_k];					// allocate near neighbor dists
// #pragma omp for  
  for( int i = 0; i<numberOfCells ;++i)
  {
    queryPt[0] = x[i];
    queryPt[1] = y[i];
    queryPt[2] = z[i];
    kd_tree_of_vl->annkSearch(			// search
                              queryPt,		// query point
                              ANN_k,		// number of near neighbors
                              ANN_nnIdx,	// nearest neighbors (returned)
                              ANN_dists,	// distance (returned)
                              ANN_eps);		// error bound
//     //check found ANN_dists
//     for( int iii = 0; iii<ANN_k;++iii)
//     {
//       cout << "ANN_dists[" << iii << "]: " << ANN_dists[iii] << endl;
//     }
    
/** at this point the 5 nearset vessel indexes are found for cell number i
 * by ANN, but maybe this is not the best one.
 * Because we filled ANN with the centers of the vessels,
 * now improve the search, by geometric arguments!
 */
    nearest candidate;
    for(int ii=0;ii<ANN_k;ii++)//this could NOT be done in parallel since the candidate is atom!
    {
      uint currentIndex_from_ANN = ANN_nnIdx[ii];
      Vessel* v = vl->GetEdge(ann_to_vl[currentIndex_from_ANN]);
      //cout << "current index: " << currentIndex_from_ANN << endl;
      double distance = calculate_distance_from_vessel_to_point_in_space(
      vl->Ld().LatticeToWorld(v->LPosA()),
      vl->Ld().LatticeToWorld(v->LPosB()),
      queryPt);
      if( distance < candidate.distance )
      {
        candidate.distance = distance;
        candidate.indexOfVessel = ann_to_vl[ANN_nnIdx[ii]];
      }
    }
    vectorOfnearestVessels[i] = candidate;
    //cout << "shortes distance for index "<< candidate.indexOfVessel << ": " << candidate.distance <<endl;
  }
// }//end #pragma omp parallel
  delete kd_tree_of_vl;
  delete [] ANN_nnIdx;
  delete [] ANN_dists;
  annClose();
#ifndef NDEBUG
  cout << "annClose called" << endl;
#endif
  //transfere this nice informations to vbl
  //note this could nicely done in parallel
  tumorcode_pointer_to_currentCellsSystem->clean_BloodVesselVector();
#ifndef NDEBUG
  cout << "cleaned BloodVesselVector" << endl;cout.flush();
#endif

  if(params.useTumorcodeVessels)
  {
    for( int i = 0; i<numberOfCells ;++i)
    {
      Float3 buffer;
      std::array<double,3> bufferToFill;
  #ifndef NDEBUG
      //printf("ecnt: %i, po2Store.size(): %i,  cell_i: %i, ann_to_vl[i]: %i", ecnt, po2Store.size(), i, ann_to_vl[i]);
      if(previous_po2vessels)
      {
        if(ann_to_vl[i]<previous_po2vessels->size())
        {
          printf("ann_to_vl[i]: %i, po2Store.size(): %i \n", ann_to_vl[i], previous_po2vessels->size());
          myAssert(ann_to_vl[i]<previous_po2vessels->size());
        }
      }
  #endif
      const Vessel* v= vl->GetEdge(ann_to_vl[i]);
      vbl::BloodVessel suggestion = vbl::BloodVessel();
      
      /****** topology ****/
      suggestion.SetBloodVesselR(v->r);
      double otherr=suggestion.GetBloodVesselR();
      //we use the Eigen3 library to store array, this is faster
      //pos a
      buffer = vl->Ld().LatticeToWorld(v->LPosA());
      bufferToFill = {buffer[0], buffer[1], buffer[2]};
      suggestion.SetBloodVessela(bufferToFill);
      //pos b
      buffer = vl->Ld().LatticeToWorld(v->LPosB());
      bufferToFill = {buffer[0], buffer[1], buffer[2]};
      suggestion.SetBloodVesselb(bufferToFill);
      
      /****** dynamics ****/
      //cout << "Main: chemical blood vessel variables " << endl; 

    //     suggestion.SetBloodVesselO2start( envO2 );
    //     suggestion.SetBloodVesselO2end( envO2 );
      if(params.useConstO2)
      {
        //do not use the po2Store!!!
        float envO2 = vbl::O2_BV;
        suggestion.SetBloodVesselO2start(envO2);
        suggestion.SetBloodVesselO2end(envO2);
      }
      else
      {
        float o2_a = to_vbl_o2_units((*previous_po2vessels)[v->Index()][0]);
        suggestion.SetBloodVesselO2start( o2_a );
        float o2_b = to_vbl_o2_units((*previous_po2vessels)[v->Index()][1]);
        suggestion.SetBloodVesselO2end( o2_b );
      }

      suggestion.SetBloodVesselCO2start( 0. );
      suggestion.SetBloodVesselCO2end( 0. );
      
      double envG = vbl::G_BV;
      suggestion.SetBloodVesselG( envG );
      
      double envA = vbl::A_BV;
      suggestion.SetBloodVesselA( envA );

      suggestion.SetBloodVesselAcL( 0. );
      
      tumorcode_pointer_to_currentCellsSystem->Add_BloodVesselVector(suggestion);
    }
  }
  else
  {
    //tumorcode vessels not used
  }
#ifndef NDEBUG
  cout << "exit find nearest" << endl;cout.flush();
#endif
}

