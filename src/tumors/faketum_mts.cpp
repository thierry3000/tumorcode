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



#ifdef USE_DETAILED_O2
void update_milotti_vessels(CellsSystem &currentCellSys, VesselList3d &vl, DetailedPO2::VesselPO2Storage &po2Store)
{
  int ecnt = vl.GetECount();
  
  
  /* copy from milotti 
   * 
   * // lettura del file dei vasi
    
    ifstream BloodVesselFile( "BloodVessels.txt" );
    cout << BloodVesselFile.is_open() << endl;
    cout << "\nMain: blood vessels from BloodVessels.txt" << endl;
    int nbv;
    BloodVesselFile >> nbv;
    cout << nbv << " vaso(i) nel file" << endl;
    for(int k=0; k<nbv; k++)
    {
        cout << k;
        BloodVessel NewBV;
        double dbv;
        BloodVesselFile >> dbv;
        NewBV.SetBloodVesselR( dbv );
        cout << "\t R: " << dbv;
        BloodVesselFile >> dbv;
        NewBV.SetBloodVesselvR( dbv );
        cout << "vR: " << dbv << " ( x1 y1 z1 ) = ( ";
        for(int j=0; j<3; j++)
        {
            BloodVesselFile >> dbv;
            NewBV.SetBloodVesselak( dbv, j );
            cout << dbv << " ";
            // cout << "k: " << k << " a[k]: " << dbv << endl;
        }
        cout << "); ( x2 y2 z2 ) = ( ";
        for(int j=0; j<3; j++)
        {
            BloodVesselFile >> dbv;
            NewBV.SetBloodVesselbk( dbv, j );
            cout << dbv << " ";
            // cout << "k: " << k << " b[k]: " << dbv << endl;
        }
        cout << "); ( vx1 vy1 vz1 ) = ( ";
        for(int j=0; j<3; j++)
        {
            BloodVesselFile >> dbv;
            NewBV.SetBloodVesselvak( dbv, j );
            cout << dbv << " ";
        }
        cout << "); ( vx2 vy2 vz2 ) = ( ";
        for(int j=0; j<3; j++)
        {
            BloodVesselFile >> dbv;
            NewBV.SetBloodVesselvbk( dbv, j );
            cout << dbv << " ";
        }
        cout << ")" << endl;
        
        // chemical blood vessel variables are set equal to environmental values
        
        cout << "Main: chemical blood vessel variables " << endl; 
        
        double envO2 = O2_BV;
        NewBV.SetBloodVesselO2start( envO2 );
        NewBV.SetBloodVesselO2end( envO2 );
        
        NewBV.SetBloodVesselCO2start( 0. );
        NewBV.SetBloodVesselCO2end( 0. );
        
        double envG = G_BV;
        NewBV.SetBloodVesselG( envG );
        
        double envA = A_BV;
        NewBV.SetBloodVesselA( envA );

        NewBV.SetBloodVesselAcL( 0. );
        
        cout << endl;
        
        
        CellsSystem.Add_BloodVessel( NewBV );  // qui si copia il vettore dei vasi nel vettore di CellsSystem

    }
    */
  //create new entry
  BloodVessel suggestion;
  Float3 buffer;
  vector<double> bufferToFill;
  
  for(int i = 0; i<ecnt; i++)
  {
    const Vessel* v= vl.GetEdge(i);
    
    if( not v->IsCirculated() )
    {
      cout << "does not make sense to consider uncirculated vessels here! " << endl;
    }
    else
    {
      /****** topology ****/
      suggestion.SetBloodVesselR(v->r);
      //we use the Eigen3 library to store array, this is faster
      //pos a
      buffer = vl.Ld().LatticeToWorld(v->LPosA());
      bufferToFill = {buffer[0], buffer[1], buffer[2]};
      suggestion.SetBloodVessela(bufferToFill);
      //pos b
      buffer = vl.Ld().LatticeToWorld(v->LPosB());
      bufferToFill = {buffer[0], buffer[1], buffer[2]};
      suggestion.SetBloodVesselb(bufferToFill);
      
      /****** dynamics ****/
      //cout << "Main: chemical blood vessel variables " << endl; 
          
      double envO2 = O2_BV;

  //     suggestion.SetBloodVesselO2start( envO2 );
  //     suggestion.SetBloodVesselO2end( envO2 );

      suggestion.SetBloodVesselO2start( po2Store[i][0] );
      suggestion.SetBloodVesselO2end( po2Store[i][1] );

      suggestion.SetBloodVesselCO2start( 0. );
      suggestion.SetBloodVesselCO2end( 0. );
      
      double envG = G_BV;
      suggestion.SetBloodVesselG( envG );
      
      double envA = A_BV;
      suggestion.SetBloodVesselA( envA );

      suggestion.SetBloodVesselAcL( 0. );
      
      currentCellSys.Add_BloodVessel( suggestion );
    }//end if circulated
  }  
}
#else
void update_milotti_vessels(CellsSystem &currentCellSys, VesselList3d &vl)
{
  int ecnt = vl.GetECount();
  
  
  /* copy from milotti 
   * 
   * // lettura del file dei vasi
    
    ifstream BloodVesselFile( "BloodVessels.txt" );
    cout << BloodVesselFile.is_open() << endl;
    cout << "\nMain: blood vessels from BloodVessels.txt" << endl;
    int nbv;
    BloodVesselFile >> nbv;
    cout << nbv << " vaso(i) nel file" << endl;
    for(int k=0; k<nbv; k++)
    {
        cout << k;
        BloodVessel NewBV;
        double dbv;
        BloodVesselFile >> dbv;
        NewBV.SetBloodVesselR( dbv );
        cout << "\t R: " << dbv;
        BloodVesselFile >> dbv;
        NewBV.SetBloodVesselvR( dbv );
        cout << "vR: " << dbv << " ( x1 y1 z1 ) = ( ";
        for(int j=0; j<3; j++)
        {
            BloodVesselFile >> dbv;
            NewBV.SetBloodVesselak( dbv, j );
            cout << dbv << " ";
            // cout << "k: " << k << " a[k]: " << dbv << endl;
        }
        cout << "); ( x2 y2 z2 ) = ( ";
        for(int j=0; j<3; j++)
        {
            BloodVesselFile >> dbv;
            NewBV.SetBloodVesselbk( dbv, j );
            cout << dbv << " ";
            // cout << "k: " << k << " b[k]: " << dbv << endl;
        }
        cout << "); ( vx1 vy1 vz1 ) = ( ";
        for(int j=0; j<3; j++)
        {
            BloodVesselFile >> dbv;
            NewBV.SetBloodVesselvak( dbv, j );
            cout << dbv << " ";
        }
        cout << "); ( vx2 vy2 vz2 ) = ( ";
        for(int j=0; j<3; j++)
        {
            BloodVesselFile >> dbv;
            NewBV.SetBloodVesselvbk( dbv, j );
            cout << dbv << " ";
        }
        cout << ")" << endl;
        
        // chemical blood vessel variables are set equal to environmental values
        
        cout << "Main: chemical blood vessel variables " << endl; 
        
        double envO2 = O2_BV;
        NewBV.SetBloodVesselO2start( envO2 );
        NewBV.SetBloodVesselO2end( envO2 );
        
        NewBV.SetBloodVesselCO2start( 0. );
        NewBV.SetBloodVesselCO2end( 0. );
        
        double envG = G_BV;
        NewBV.SetBloodVesselG( envG );
        
        double envA = A_BV;
        NewBV.SetBloodVesselA( envA );

        NewBV.SetBloodVesselAcL( 0. );
        
        cout << endl;
        
        
        CellsSystem.Add_BloodVessel( NewBV );  // qui si copia il vettore dei vasi nel vettore di CellsSystem

    }
    */
  //create new entry
  BloodVessel suggestion;
  Float3 buffer;
  vector<double> bufferToFill;
  
  for(int i = 0; i<ecnt; i++)
  {
    /****** topology ****/
    const Vessel* v= vl.GetEdge(i);
    suggestion.SetBloodVesselR(v->r);
    //we use the Eigen3 library to store array, this is faster
    //pos a
    buffer = vl.Ld().LatticeToWorld(v->LPosA());
    bufferToFill = {buffer[0], buffer[1], buffer[2]};
    suggestion.SetBloodVessela(bufferToFill);
    //pos b
    buffer = vl.Ld().LatticeToWorld(v->LPosB());
    bufferToFill = {buffer[0], buffer[1], buffer[2]};
    suggestion.SetBloodVesselb(bufferToFill);
    
    /****** dynamics ****/
    //cout << "Main: chemical blood vessel variables " << endl; 
        
    double envO2 = O2_BV;

    suggestion.SetBloodVesselO2start( envO2 );
    suggestion.SetBloodVesselO2end( envO2 );

//     suggestion.SetBloodVesselO2start( po2Store[i][0] );
//     suggestion.SetBloodVesselO2end( po2Store[i][1] );

    suggestion.SetBloodVesselCO2start( 0. );
    suggestion.SetBloodVesselCO2end( 0. );
    
    double envG = G_BV;
    suggestion.SetBloodVesselG( envG );
    
    double envA = A_BV;
    suggestion.SetBloodVesselA( envA );

    suggestion.SetBloodVesselAcL( 0. );
    
    currentCellSys.Add_BloodVessel( suggestion );
  }  
}
#endif

void initMilotti(CellsSystem &currentCellsSystem)
//void initMilotti(CellsSystem &currentCellsSystem, VesselList3d &vl, DetailedPO2::VesselPO2Storage &po2Store)
{
    /**   INIT Milotti   */
#ifndef undo
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
//   currentCellsSystem.Set_Commands( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/commands.txt" );
//   currentCellsSystem.Set_CellTypeFile( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/CellType.txt" );
//   currentCellsSystem.Set_CellTypeFileAlt( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/CellType.txt" );
//   currentCellsSystem.Set_EnvironmentFile( "/home/usersHR/thierry/git_codes/Sim3D-v3/parameters/Environment.txt" );
  currentCellsSystem.Set_Commands( "commands.txt" );
  currentCellsSystem.Set_CellTypeFile( "CellType.txt" );
  currentCellsSystem.Set_CellTypeFileAlt( "CellType.txt" );
  currentCellsSystem.Set_EnvironmentFile( "Environment.txt" );
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
}

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
  std::vector<double> allCellVolumina = currentCellsSystem.Get_volume();
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
#define someInteraction
#ifdef someInteraction // do the super cool new interaction stuff
  float estimated_tum_radius = estimateTumorRadiusFromCells();
  float currentShell = pos.norm();
  return std::max(estimated_tum_radius + params.rGf - currentShell, 0.);
#else // do it like fake tum has done it before
  float r = pos.norm();
  return std::max(tumor_radius + params.rGf - r, 0.);
#endif
}

//float FakeTumMTS::FakeTumorSimMTS::getPress(const Float3 &pos) const
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
  {
    /* setup basic things */
    my::SetNumThreads(params.num_threads);
    /* open hdf5 file containing the vessels */
    h5cpp::File file(params.fn_vessel, "r");
    /* instructions on how to read the vessels */
    ptree pt;
    pt.put("scale subdivide", 10.); //subdivide lattice by factor 10
    /* read vessels 
     * n.b. this file is keep open during the complete simulation
     * an all our vessels are pointer to that file. this keep ram at minimum
     */
    vl = ReadVesselList3d(file.root().open_group("vessels"), pt);
    
    //create kd-tree from vessellist
    fillKdTreeFromVl();
    
    /*
     * this test worked, know we try to find the closest vessel for every cell in a similar fashion
    //test kd tree
    queryPt = annAllocPt(ANN_dim);
    queryPt[0] = 100.0;
    queryPt[1] = 150.0;
    queryPt[2] = 200.0;
    // execute search
    kd_tree_of_vl->annkSearch(						// search
                              queryPt,		// query point
                              ANN_k,			// number of near neighbors
                              ANN_nnIdx,	// nearest neighbors (returned)
                              ANN_dists,	// distance (returned)
                              ANN_eps);		// error bound
                              */
    //print search results
#ifdef DEBUG
#if 0
    cout << "\tNN:\tIndex\tDistance\n";
		for (int i = 0; i < ANN_k; i++) 
    {			// print summary
			ANN_dists[i] = sqrt(ANN_dists[i]);			// unsquare distance
			cout << "\t" << i << "\t" << ANN_nnIdx[i] << "\t" << ANN_dists[i] << "\n";
			cout << "Point: ";
      cout << "(" << dataPts[ANN_nnIdx[i]][0];
      for (int i = 1; i < 3; i++) 
      {
        cout << ", " << dataPts[ANN_nnIdx[i]][i];
      }
      cout << ")\n";
    }
#endif
#endif
    // adjust vessel list ld
    const Float3 c = 0.5 * (vl->Ld().GetWorldBox().max + vl->Ld().GetWorldBox().min);
    vl->SetDomainOrigin(vl->Ld().LatticeToWorld(Int3(0))-c);

    cout << "--------------------"<< endl;
    cout << "Vessel Lattice is: " << endl;
    vl->Ld().print(cout); cout  << endl;
    cout << "--------------------"<< endl;

    params.vesselfile_message = file.root().open_group("parameters").attrs().get<string>("MESSAGE");
    params.vesselfile_ensemble_index = file.root().open_group("parameters").attrs().get<int>("ENSEMBLE_INDEX");
    
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
    model.Init(vl.get(), this->model.params, callbacks);
  }

  /* set initial conditions of the FakeTumMTS simulation */
  tumor_radius = params.tumor_radius;
  time = 0.;
  num_iteration = 0.;
  output_num = 0;
  double next_output_time = 0.;
  double next_adaption_time = 0.;
  
  //for the adaption it could be usefull to have the
  //vessel network after the adaption in the beginning   ---> done automatically since adaption is in tum-only-vessels
//   bool writeVesselsafter_initial_adaption = false;
//   if(model.params.badaption_on_off)
//   {
//     writeVesselsafter_initial_adaption = true;
//   }
  
  // this is needed to pass tumor information to the DetailedPO2 simulation
  boost::optional<h5cpp::Group> lastTumorGroupWrittenByFakeTum;
  std::string lastTumorGroupWrittenByFakeTumName;

  
#ifdef USE_DETAILED_O2
  /* set up detailed oxygen calculation */
  //this should be read from file later on
  //DetailedPO2::Parameters oxy_params;
  double grid_lattice_const = 40;
  double safety_layer_size = 120;
  boost::optional<Int3> grid_lattice_size;
#endif
  /* continum lattice stuff
   * 
   * needed for solving diffusion equations
   */
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
    
  //oxyops.init(mtboxes, grid.Box(), grid.dim, 2);
  //oxyops.init(state.o2field);
  
  /* glucose is not yet used, but could be handeled like this*/
  glucoseOps.init(mtboxes, grid.Box(), grid.dim, 2);
  glucoseOps.init(state.glucoseField);
  all_pt_params.add_child("glucose" , glucoseParams.as_ptree());
  /* vessel volume fractions are need to set the source strengths */
  UpdateVesselVolumeFraction();
  /* init the cell system */
  initMilotti(currentCellsSystem);
  
  last_chem_update = -1;
  last_vessels_checksum = -1;
  
  /* start main loop */
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
//#pragma omp barrier
    if (time >= next_output_time - params.dt * 0.1)
    {
#ifdef USE_DETAILED_O2
      /* right now, the o2 simulation is called at every output time */
      {
        o2_sim.init(o2_params, bfparams,*vl,grid_lattice_const, safety_layer_size, grid_lattice_size, lastTumorGroupWrittenByFakeTum);
        cout << "\nInit O2 completed" << endl;
        o2_sim.run(*vl);
        cout << "\nDetailed O2 completed" << endl;
      }
#else
      /// no detailed o2 model
      
#endif
      /*    o2 data is created
       *  and feed to the cells
       */
#ifdef USE_DETAILED_O2
      // feed milotti structure with detailed o2 simulation result
      update_milotti_vessels(currentCellsSystem, *vl, o2_sim.po2vessels);
#else
      // simple version which set oxygen level of the vessels to a constant value
      update_milotti_vessels(currentCellsSystem, *vl);
#endif
      lastTumorGroupWrittenByFakeTumName = writeOutput();
      h5::File f(params.fn_out + ".h5", "a");
      lastTumorGroupWrittenByFakeTum = f.root().open_group(lastTumorGroupWrittenByFakeTumName+"/tumor");
      
      //provide addition information about the tissue phases to the diffusion solver
      SetupTissuePhases(phases, grid, mtboxes, lastTumorGroupWrittenByFakeTum);//filling
      f.close();
      next_output_time += params.out_intervall;
    }

    if (time > params.tend) break;
    
    /* stop if tumor reaches certain fraction of volume */
    double size_limit = 0.5*maxCoeff(Size(vl->Ld().GetWorldBox())) * params.stopping_radius_fraction; 
    //cout << format("size_limit = %f vs tumor_radius = %f\n") % size_limit % tumor_radius;
    if (tumor_radius >  size_limit) break;

    /* do a vessel model remodeling step */
    cout << boost::format("start vessel remodel step! \n");
    doStep(params.dt);
    fillKdTreeFromVl();
    cout << boost::format("finished vessel remodel step! \n");
    /* increment tumor time */
    time += params.dt;
    /* propergate cells in time until current fake tumor time */
    cout << boost::format("advance milotti until: %f\n") % time;
    //currentCellsSystem.Set_tmax(time);
    doMilottiStep();

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
  CalcFlow(*vl, bfparams);
  /* do not use adaptation stuff for mts*/
#ifdef USE_ADAPTION
  model.DoStep(dt, &adap_params, &bfparams);
#else
  model.DoStep(dt, &bfparams);
#endif
  
  /* this calculates the simple diffusion of substances */
  //calcChemFields();
  tumor_radius += dt * params.tumor_speed;
}

#ifndef undo
void FakeTumMTS::FakeTumorSimMTS::doMilottiStep()
{
  cout << format("start mts at tumor time: %f\n" ) % time;
  //from milotti
  /* the tumor time is given in hours 
   * evolve cells until that
   */
  uint returnValue = currentCellsSystem.runMainLoop( time * 3600 );
  
  // for safety reasons we use both output structures
  // this one the the output of milotti, see WriteCellsSystemHDF for the hdf output
  if ( currentCellsSystem.Get_ready2start() )
  {
    currentCellsSystem.Print2logfile("Cells at the end of the run");
    currentCellsSystem.WriteCellsSystem( );					// dump of the final configuration
  }
  
  cout << format("finished FakeTumMTS::FakeTumorSimMTS::doMilottiStep(),  mts at tumor time: %f\n" ) % time;
  //end milotti
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
    // write the simulation parameters to file
    a.set("MESSAGE",params.message);
    a.set("VESSELTREEFILE",params.fn_vessel);
    a.set("OUTPUT_NAME", params.fn_out);
    a.set("VESSELFILE_MESSAGE", params.vesselfile_message);
    a.set("VESSELFILE_ENSEMBLE_INDEX", params.vesselfile_ensemble_index);
    g = root.create_group("parameters");
    WriteHdfPtree(g.create_group("vessels"), model.params.as_ptree());
    WriteHdfPtree(g.create_group("o2_params"), o2_params.as_ptree());
    WriteHdfPtree(g.create_group("calcflow"), bfparams.as_ptree());
    WriteHdfPtree(g.create_group("faketumorMTS"), params.as_ptree());
  }
  
  /* create time slices */
  std::string tumOutName = str(format("out%04i") % output_num);
  h5::Group gout = root.create_group(tumOutName);
  a = gout.attrs();
  a.set("time", time);
  a.set("OUTPUT_NUM",output_num);
  h5::Group cells_out = gout.create_group("cells");
  /* writes the cell system stuff */
  //WriteCellsSystemHDF(currentCellsSystem, cells_out);
  WriteCellsSystemHDF_with_nearest_vessel_index(currentCellsSystem,kd_tree_of_vl, cells_out);
  /* writes the vessel list */
  WriteVesselList3d(*vl, gout.create_group("vessels"));
  {
    /* write continuous fields */
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
    WriteScalarField(gtum, "fieldGlucose", state.glucoseField, grid.ld, field_ld_group);
//     WriteScalarField(gtum, "fieldOxy", state.o2field, ld, field_ld_group);
    // vessel continuum
    UpdateVesselVolumeFraction();
    WriteScalarField(gtum, "vessel_volume_fraction", vessel_volume_fraction, grid.ld, field_ld_group);
    WriteScalarField(gtum, "oxy_source_lin", vessel_o2src_clin, grid.ld, field_ld_group);
  }
  {
#ifdef USE_DETAILED_O2
    //write the detailed o2 stuff necessary for the cells
    int ecnt = o2_sim.po2vessels.size();
    DynArray<float> po2(2*ecnt);
    for (int i=0; i<ecnt; ++i)
    {
      // note: po2 values could be zero, if vessel is uncirculated, and even when circulated, angiogenesis + calcflow, ---> check that
      po2[2*i+0] = o2_sim.po2vessels[i][0];
      po2[2*i+1] = o2_sim.po2vessels[i][1];
    }
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(gout.create_group("detailedPo2"), "po2_vessels", h5cpp::Dataspace::simple_dims(ecnt,2), &po2[0]);
#else
    int ecnt = currentCellsSystem.Get_nbv();
    vector<BloodVessel> bfvessels = currentCellsSystem.Get_BloodVesselVector();
    DynArray<float> po2(2*ecnt);
    for (int i=0; i<ecnt; ++i)
    {
      // note: po2 values could be zero, if vessel is uncirculated, and even when circulated, angiogenesis + calcflow, ---> check that
      po2[2*i+0] = bfvessels[i].GetBloodVesselO2start();
      po2[2*i+1] = bfvessels[i].GetBloodVesselO2end();
    }
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(gout.create_group("detailedPo2"), "po2_vessels", h5cpp::Dataspace::simple_dims(ecnt,2), &po2[0]);
    //WriteScalarField(gtum, "fieldOxy", state.o2field, grid.ld, field_ld_group);
#endif
  }
  ++output_num;
  f.flush();
  f.close();
  cout << format("files %s flushed and closed")  % params.fn_out << endl;
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
//     StationaryDiffusionSolve(grid, mtboxes, state.glucoseField, boost::bind(&FakeTumorSimMTS::insertGlucoseCoefficients, this, _1, _2, boost::cref(state), _3) , pt_params);
//     #pragma omp parallel
//     {
//       BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
//       {
//         CopyBorder(state.o2field, bbox, grid.Box(), grid.dim, 2);
//       }
//     }
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

/* this has to be tidyed up before proceeding */
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
//   #pragma omp parallel
//   {
//     VesselVolumeGenerator volumegen(*vl, grid.ld, grid.dim, make_ptree("samples_per_cell", 100));
//     
//     BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, mtboxes.getCurrentThreadRange())
//     {
//       int dummy_;
//       volumegen.Fill(bbox, vessel_volume_fraction, vessboxes[bbox.global_index], dummy_);
//       O2Model::AddSourceDistribution( bbox, grid.ld, 
//                                       grid.dim, 
//                                       vessel_o2src_clin, 
//                                       vessel_o2src_crhs, 
//                                       vl->Ld(), 
//                                       vessboxes[bbox.global_index], 
//                                       all_pt_params.get_child("simple_o2"));
//     }
//   }
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

void FakeTumMTS::FakeTumorSimMTS::WriteCellsSystemHDF_with_nearest_vessel_index(CellsSystem &currentCellsSystem, ANNkd_tree *kd_tree_of_vl, h5cpp::Group &out_cell_group)
{
  cout<< "going to write cells to a hdf file" << endl;
  int numberOfCells = currentCellsSystem.Get_ncells();
  std::vector<double> x = currentCellsSystem.Get_x();
  std::vector<double> y = currentCellsSystem.Get_y();
  std::vector<double> z = currentCellsSystem.Get_z();
  DynArray<float> a(3*numberOfCells);
  DynArray<int> index_of_nearest_vessel(numberOfCells);
  for( int i = 0; i<numberOfCells ;++i)
  {
    a[3*i+0] = x[i];
    a[3*i+1] = y[i];
    a[3*i+2] = z[i];
    queryPt[0] = x[i];
    queryPt[1] = y[i];
    queryPt[2] = z[i];
    kd_tree_of_vl->annkSearch(						// search
                              queryPt,		// query point
                              ANN_k,			// number of near neighbors
                              ANN_nnIdx,	// nearest neighbors (returned)
                              ANN_dists,	// distance (returned)
                              ANN_eps);		// error bound
    index_of_nearest_vessel[i] = ANN_nnIdx[0];
  }
  if(!out_cell_group.exists("index_of_nearest_vessel"))
  {
//     h5cpp::Dataset ds = h5cpp::create_dataset<float>(vesselgroup, "cell_center_pos", h5cpp::Dataspace::simple_dims(a.size()/3,3), &a[0]);
    h5cpp::Dataset ds = h5cpp::create_dataset<int>(out_cell_group, "index_of_nearest_vessel", h5cpp::Dataspace::simple_dims(numberOfCells,1), &index_of_nearest_vessel[0]);
  }
  if(!out_cell_group.exists("cell_center_pos"))
  {
//     h5cpp::Dataset ds = h5cpp::create_dataset<float>(vesselgroup, "cell_center_pos", h5cpp::Dataspace::simple_dims(a.size()/3,3), &a[0]);
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "cell_center_pos", h5cpp::Dataspace::simple_dims(numberOfCells,3), &a[0]);
  }
  DynArray<float> buffer(numberOfCells);
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_r()[i];
  }
  if(!out_cell_group.exists("cell_radii"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "cell_radii", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  // glucose extracellular Get_G_extra()
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_G_extra()[i];
  }
  if(!out_cell_group.exists("glucose_ex"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "glucose_ex", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  // ph Get_pH
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_pH()[i];
  }
  if(!out_cell_group.exists("pH_ex"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "pH_ex", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  // oxygen Get_O2
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_O2()[i];
  }
  if(!out_cell_group.exists("o2"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "o2", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  // lactate Get_AcL_extra
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_AcL_extra()[i];
  }
  if(!out_cell_group.exists("AcL_ex"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "AcL_ex", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  cout<< "finished writting cells to hdf" << endl;
}

void FakeTumMTS::FakeTumorSimMTS::WriteCellsSystemHDF(CellsSystem &currentCellsSystem, h5cpp::Group &out_cell_group)
{
  cout<< "going to write cells to a hdf file" << endl;
  int numberOfCells = currentCellsSystem.Get_ncells();
  std::vector<double> x = currentCellsSystem.Get_x();
  std::vector<double> y = currentCellsSystem.Get_y();
  std::vector<double> z = currentCellsSystem.Get_z();
  DynArray<float> a(3*numberOfCells);
  for( int i = 0; i<numberOfCells ;++i)
  {
//     a[3*i+0]= vl.GetNode(i)->worldpos[0]; /// from world part
//     a[i] = x[i];
//     a[i+1] = y[i];
//     a[i+2] = z[i];
    a[3*i+0] = x[i];
    a[3*i+1] = y[i];
    a[3*i+2] = z[i];
  }
  if(!out_cell_group.exists("cell_center_pos"))
  {
//     h5cpp::Dataset ds = h5cpp::create_dataset<float>(vesselgroup, "cell_center_pos", h5cpp::Dataspace::simple_dims(a.size()/3,3), &a[0]);
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "cell_center_pos", h5cpp::Dataspace::simple_dims(numberOfCells,3), &a[0]);
  }
  DynArray<float> buffer(numberOfCells);
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_r()[i];
  }
  if(!out_cell_group.exists("cell_radii"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "cell_radii", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  // glucose extracellular Get_G_extra()
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_G_extra()[i];
  }
  if(!out_cell_group.exists("glucose_ex"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "glucose_ex", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  // ph Get_pH
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_pH()[i];
  }
  if(!out_cell_group.exists("pH_ex"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "pH_ex", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  // oxygen Get_O2
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_O2()[i];
  }
  if(!out_cell_group.exists("o2"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "o2", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  // lactate Get_AcL_extra
  for( int i = 0; i<numberOfCells; ++i)
  {
    buffer[i] = currentCellsSystem.Get_AcL_extra()[i];
  }
  if(!out_cell_group.exists("AcL_ex"))
  {
    h5cpp::Dataset ds = h5cpp::create_dataset<float>(out_cell_group, "AcL_ex", h5cpp::Dataspace::simple_dims(numberOfCells,1), &buffer[0]);
  }
  cout<< "finished writting cells to hdf" << endl;
}
/** @brief
 * This will fill the kd-tree structure with the midpoint of the vessels 
 * stored in the vessel list vl
 */
void FakeTumMTS::FakeTumorSimMTS::fillKdTreeFromVl()
{
  int              nPts       = vl->GetECount();					// actual number of data points
  myAssert(nPts<ANN_maxPts);
  
	queryPt = annAllocPt(ANN_dim);                  // allocate query point
	dataPts = annAllocPts(ANN_maxPts, ANN_dim);     // allocate data points
	
	// fill kd tree structur
  #pragma omp_parallel_for
  for(int i=0; i<nPts; ++i)
  {
    Vessel* v = vl->GetEdge(i);                       //get vessel per thread
    Float3 p_a,p_b,a_b,vessel_center;                 //allocate in loop to make sure we are threadsafe
    p_a = vl->Ld().LatticeToWorld(v->NodeA()->lpos);  //read position a
    p_b = vl->Ld().LatticeToWorld(v->NodeB()->lpos);  //read position b
    a_b = p_a-p_b;                                    // difference vector
    vessel_center = p_b+0.5*a_b;                      // vessel center
    ANNpoint p;
    for(int k=0;k<ANN_dim;++k)
    {
      dataPts[i][k] = vessel_center[k];
    }
  }
  kd_tree_of_vl = new ANNkd_tree(					// build search structure
                                  dataPts,					// the data points
                                  nPts,						// number of points
                                  ANN_dim);						// dimension of space
  //prepare search
  ANN_nnIdx = new ANNidx[ANN_k];						// allocate near neigh indices
  ANN_dists = new ANNdist[ANN_k];						// allocate near neighbor dists
}
