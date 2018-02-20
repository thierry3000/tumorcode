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
#include "iff_drug_app.h"
#include "vesselmodel1.h"
#include "calcflow.h"
#include "bulktissuemodel1_new.h"

IffDrugApp3d::IffDrugApp3d()
{
// #ifdef DEBUG
//   out_intervall = 0.1/60.;
//   tend = 240;
// #else
//   out_intervall = 10./60;
//   tend = 5;
// #endif
  message = "(None)";
  output_number = 0;
}


void ReadInto(const BBox3 &bb, H5::DataSet ds, Array3d<float> arr)
{
  Array3d<float> tmp;
  ReadArray3D(ds, tmp);
  arr[bb].fill(tmp[bb]);
}


// void IffDrugApp3d::InitFieldsAndLd(h5::Group ldgroup)
// {
//   ReadHdfLd(ldgroup, field_ld); // read lattice data
//
//   phi_cells = MakeArray3dWithBorder<float>(field_ld.Box(), dim, 1, SET_BOX_INNER);
//   phi_water = MakeArray3dWithBorder<float>(field_ld.Box(), dim, 1, SET_BOX_INNER);
//   theta_tumor = MakeArray3dWithBorder<float>(field_ld.Box(), dim, 1, SET_BOX_INNER);
// }

void IffDrugApp3d::InitAsSimpleModels(Array3d<float> theta_tumor, Array3d<float> theta_necro, Array3d<float> phi_water, Array3d<float> phi_cells)
{
  SimpleTissueModel *tmp1 = new SimpleTissueModel();
  tmp1->params = &iff_params;
  tmp1->phi_water = phi_water;
  tmp1->theta_necro = theta_necro;
  tmp1->theta_tumor = theta_tumor;
  tissue_f.reset((tmp1));

  SimpleLymphModel *tmp2 = new SimpleLymphModel();
  tmp2->params = &iff_params;
  tmp2->theta_tumor = theta_tumor;
  lymph_sources_f.reset(tmp2);

  phi_cells_f.reset(new VirtualArray3dGridFunctions<float,1>(phi_cells));
  phi_water_f.reset(new VirtualArray3dGridFunctions<float,1>(phi_water));
  phi_cells_wo_necro_f.reset(new VirtualGridFuncCellsWoNecro(phi_cells, theta_necro));

  tissueCompositionCallback.phi_cells = phi_cells;
  tissueCompositionCallback.phi_water = phi_water;
}


bool IffDrugApp3d::InitNewState()
{
  cout << "--- starting new iff run ---\n";
  cout << "--- input tum file: " << fn_tumor <<endl;
  cout << "--- output iff file: " << fn_out <<endl;
  {
  iff_params.assign(all_pt_params.get_child("iff"));
  boost::optional<ptree&> pt_drug = all_pt_params.get_child_optional("drugflow");
  if (pt_drug)
    ift_params.assign(*pt_drug);
  }

  {
  H5::H5File file(fn_tumor, H5F_ACC_RDONLY);
  H5::Group root = file.openGroup("/");
  H5::Group tum_grp;
  bool tumorIsPresent = true;
  try
  {
    tum_grp = root.openGroup(h5_path_tumor);
  }
  catch( H5::Exception )
  {
    tumorIsPresent=false;
  }
  
//   if (root.exists(h5_path_tumor))
//   {
//     tum_grp = root.open_group(h5_path_tumor);
//     h5::Attributes tum_attrs = tum_grp.attrs();
//     if (tumor_type_string.empty())
//       tumor_type_string = tum_attrs.get<string>("TYPE");
//   }
//   else
//   {
//     tumorIsPresent=false;
//   }
  
  H5::Group h5_lattice = root.openGroup(h5_path_lattice);
  std::unique_ptr<polymorphic_latticedata::LatticeData> ldp = polymorphic_latticedata::LatticeData::ReadHdf(h5_lattice);
  //vl.reset(new VesselList3d(ldp));
  vl->Init(*ldp);
  H5::Group vesselgroup = root.openGroup(h5_path_vessel);
  //ReadHdfGraph(vesselgroup, *vesselList);
  ReadHdfGraph(vesselgroup, vl.get());
  { // read node pressure values
  DynArray<float> press;
  //h5::read_dataset(vesselgroup.open_dataset("nodes/pressure"), press);
  //readDataSetFromGroup<DynArray<float>>(vesselgroup.openGroup("nodes"),string("pressure"), &press);
  H5::Group h5_nodes = vesselgroup.openGroup("nodes");
  readDataSetFromGroup(h5_nodes,string("pressure"), press);
  int ncnt = vl->GetNCount();
  for(int i=0; i<ncnt; ++i)
  {
    vl->GetNode(i)->press = press[i];
  }
  // read edge flow values
  DynArray<double> a;
  //h5::read_dataset(vesselgroup.open_dataset("edges/flow"), a);
  //readDataSetFromGroup<DynArray<double>>(vesselgroup.openGroup("edges"),string("flow"), &a);
  H5::Group h5_edges = vesselgroup.openGroup("edges");
  readDataSetFromGroup(h5_edges,string("flow"), a);
  int ecnt = vl->GetECount();
  for(int i=0; i<ecnt; ++i)
  {
    vl->GetEdge(i)->q = a[i];
  }
  //h5::read_dataset(vesselgroup.open_dataset("edges/maturation"), a);
  //readDataSetFromGroup<DynArray<double>>(vesselgroup.openGroup("edges"),string("maturation"), &a);
  readDataSetFromGroup(h5_edges,string("maturation"), a );
  for(int i=0; i<ecnt; ++i) vl->GetEdge(i)->maturation = a[i];
  }

#if 0
  if (tumor_type_string == "3 Phase crystal")
  {
    // read the tumor
    h5::Dataset ds_phi_tumor = tum_grp.open_dataset("tc_density");
    InitFieldsAndLd(root.open_group(ds_phi_tumor.attrs().get<string>("LATTICE_PATH")));

    ReadInto(field_ld.Box(), ds_phi_tumor, theta_tumor);
    ReadInto(field_ld.Box(), tum_grp.open_dataset("dead_density"), phi_water);

    FOR_BBOX3(p, field_ld.Box())
    {
      float tum = theta_tumor(p);
      float dead = phi_water(p);
      // theta_tumor = tum + necro
      theta_tumor(p) = my::cut<float>(tum+dead, 0., 1.);
      // cells = 0.5 - necro
      phi_cells(p) = my::cut<float>(0.5 - dead, 0., 1.);
      // water = 0.3 + necro
      phi_water(p) = my::cut<float>(0.3 + dead, 0., 0.7);
    }
  }
#endif

  if (tumor_type_string == "faketumor" and tumorIsPresent)
  {
    H5::DataSet ds_phi_tumor = tum_grp.openDataSet("tc_density");
    //double radius = tum_grp.attrs().get<double>("TUMOR_RADIUS");
    double radius;
    readAttrFromH5(tum_grp, string("TUMOR_RADIUS"), radius);
    { 
      LatticeDataQuad3d ld;
      SetupFieldLattice(vl->Ld().GetWorldBox(), 3, 30., -30, ld);
      grid = ContinuumGrid(ld, dim); 
    }
    phi_water = MakeArray3dWithBorder<float>(grid.ld.Box(), dim, 3);
    phi_cells = MakeArray3dWithBorder<float>(grid.ld.Box(), dim, 1);
    theta_tumor = MakeArray3dWithBorder<float>(grid.ld.Box(), dim, 1);
    theta_necro = MakeArray3dWithBorder<float>(grid.ld.Box(), dim, 1);

    FOR_BBOX3(p, grid.ld.Box())
    {
      theta_tumor(p) = my::smooth_heaviside_sin(radius - grid.ld.LatticeToWorld(p).norm(), 30.);
      phi_water(p) = 0.2;
      phi_cells(p) = 0.5;
    }

    CopyBorder(phi_cells, dim, 1);
    CopyBorder(phi_water, dim, 3);
    CopyBorder(theta_tumor, dim, 1);
    CopyBorder(theta_necro, dim, 1);

    InitAsSimpleModels(theta_tumor, theta_necro, phi_water, phi_cells);
  }
  else if ((tumor_type_string == "BulkTissue" || tumor_type_string == "BulkTissueFormat1") and tumorIsPresent)
  {
//     h5::Dataset ds_theta_tumor = tum_grp.open_dataset("ptc");
//     h5::Group ld_group = root.open_group(ds_theta_tumor.attrs().get<string>("LATTICE_PATH"));
    H5::DataSet ds_theta_tumor = tum_grp.openDataSet("ptc");
    string lattice_path;
    readAttrFromH5(ds_theta_tumor, "LATTICE_PATH", lattice_path);
    H5::Group ld_group = root.openGroup(lattice_path);
    { 
      LatticeDataQuad3d ld;
      ReadHdfLd(ld_group, ld);
      grid = ContinuumGrid(ld, dim); 
    }
    Array3d<float> phi_obstacle(MakeArray3dWithBorder<float>(grid.ld.Box(), dim, 1));
    phi_water = MakeArray3dWithBorder<float>(grid.ld.Box(), dim, 3);
    phi_cells = MakeArray3dWithBorder<float>(grid.ld.Box(), dim, 1);
    theta_tumor = MakeArray3dWithBorder<float>(grid.ld.Box(), dim, 1);
    theta_necro = MakeArray3dWithBorder<float>(grid.ld.Box(), dim, 1);

    ReadInto(grid.ld.Box(), ds_theta_tumor, theta_tumor);
    ReadInto(grid.ld.Box(), tum_grp.openDataSet("conc"), phi_cells);
    ReadInto(grid.ld.Box(), tum_grp.openDataSet("obstacle"), phi_obstacle);
//     if (tum_grp.exists("necro"))
//       ReadInto(grid.ld.Box(), tum_grp.open_dataset("necro"), theta_necro);
    try
    {
      H5::DataSet necro = tum_grp.openDataSet("necro");
      ReadInto(grid.ld.Box(), necro, theta_necro);
    }
    catch( H5::Exception error)
    {
      //std::cout << error.printError();
    }

    FOR_BBOX3(p, grid.ld.Box())
    {
      float t = phi_obstacle(p) + phi_cells(p);
      phi_water(p) = my::clamp(1. - t, 0., 0.99);
    }

    CopyBorder(phi_cells, dim, 1);
    CopyBorder(phi_water, dim, 3);
    CopyBorder(theta_tumor, dim, 1);
    CopyBorder(theta_necro, dim, 1);

    InitAsSimpleModels(theta_tumor, theta_necro, phi_water, phi_cells);
  }
  else
    throw std::runtime_error("unkown tumor type");

  //T.F. why is this 64 here?
  mtboxes = MakeMtBoxGrid(grid.ld.Box(), Int3(64, 32, 32));
  SortVesselsIntoMtBoxGrid(grid.ld, *vl, 1, mtboxes, vesselsInBoxes);
  VesselFluidExavasationModel *tmp = new VesselFluidExavasationModel();
  tmp->Init(*vl, grid, vesselsInBoxes, iff_params);
  vessel_sources_f.reset(tmp);
  }

#ifdef ANALYZE_FLOW_CHANGE
// compute pressure and flow without IF coupling
  { 
    CalcFlow(vl.get(), 0, 0);
    int ecnt = vl->GetECount();
    int ncnt = vl->GetNCount();
    org_press.resize(ncnt);
    org_flow.resize(ecnt);
    for (int i=0; i<ecnt; ++i)
      org_flow[i] = vl->GetEdge(i)->q;
    for (int i=0; i<ncnt; ++i) {
      org_press[i] = vl->GetNode(i)->press;
      //cout << org_press[i] << endl;
    }
  }
#endif
  cout << "read complete" << endl;

  cout << "vessel ld: " << vl->Ld() << endl;
  cout << "field ld: " << grid.ld << endl;

  iffstate.reset( new IfFlowState() );
  return true;
}


void IffDrugApp3d::MeasureIfFlowState(H5::Group g)
{
#if 0
  g.attrs().set("MESSAGE",message);
  g.attrs().set("INPUTFILE",fn_tumor);
  g.create_group("parameters");
  iff_params.WriteH5(g.create_group("parameters/iff"));
  h5cpp::Datatype disktype = h5cpp::get_disktype<float>();

  const LatticeDataQuad3d &ld = grid.ld;

  h5::Group field_ld_grp = g.create_group("field_ld");
  WriteHdfLd(field_ld_grp, ld);

  h5::Dataset ds;
  Array3d<Float3> flowvel;

  g = g.create_group("iff");

  {
    Array3d<double> cond = MakeArray3dWithBorder<double>(ld.Box(), dim, 1);
    Array3d<double> phi_water = MakeArray3dWithBorder<double>(ld.Box(), dim, 1);
    VirtualGridFunctions<double, 2>::List l1 = {{ cond, phi_water }};
    tissue_f->GetValues(-1, ld.Box(), l1);

    CopyBorder(cond, dim, 1);
    CopyBorder(phi_water, dim, 1);

    WriteScalarField(g, "iff_sources", iffstate->source_field, ld, field_ld_grp, disktype);
    WriteScalarField(g, "iff_sources_vess_out", iffstate->vessel_outflow_field, ld, field_ld_grp, disktype); // this is uptake
    WriteScalarField(g, "iff_sources_vess_in", iffstate->vessel_inflow_field, ld, field_ld_grp, disktype);
    WriteScalarField(g, "iff_sources_lymph_out", iffstate->lymph_outflow_field, ld, field_ld_grp, disktype);
    WriteScalarField(g, "phi_water", phi_water, ld, field_ld_grp, disktype);
    WriteScalarField(g, "phi_cells", phi_cells, ld, field_ld_grp, disktype);
    WriteScalarField(g, "cond", cond, ld, field_ld_grp, disktype);

    WriteScalarField( g,"iff_pressure", iffstate->pfield, ld, field_ld_grp, disktype);
    WriteScalarField(g,"theta_tumor", theta_tumor, ld, field_ld_grp, disktype);
    WriteScalarField(g,"theta_necro", theta_necro, ld, field_ld_grp, disktype);

#if 0
    {
      // check sources
      Array3d<double> div = MakeArray3dWithBorder<double>(ld.Box(), grid.dim, 1);
      // the unit of this field is (micron^3/s)[of fluid volume]/(micron^3)[tissue volume]
      for (int axis=0; axis<grid.dim; ++axis)
      {
        FOR_BBOX3(p, ld.Box())
        {
          Int3 q = add_to_axis(p, axis, -1);
          float loc_water = 0.5 * (phi_water(q) + phi_water(p));
          float flux = loc_water * iffstate->velocities[axis](p) / ld.Scale();
          div(p) -= flux;
          div(q) += flux;
        }
      }
      WriteScalarField(g, "iff_div", div[ld.Box()], ld, field_ld_grp);
      div[ld.Box()] -= iffstate->source_field[ld.Box()];
      WriteScalarField(g, "iff_div_minus_source", div[ld.Box()], ld, field_ld_grp);
      div.clear();
    }
#endif

    flowvel.initFromBox(ld.Box());
    FOR_BBOX3(p, ld.Box())
    {
      Float3 u;
      for (int axis=0; axis<dim; ++axis)
      {
        Int3 pp = p; ++pp[axis];
        u[axis] = 0.5 * (iffstate->velocities[axis](p) + iffstate->velocities[axis](pp));
      }
      flowvel(p) = u;
    }
    ds = WriteVectorField(g, "iff_velocity", flowvel, ld, field_ld_grp);

    //if (iff_params.debugOut)
//     {
//       CalcFlowField(iffstate->pfield, cond, ld, flowvel);
//       ds = WriteVectorField(g, "iff_velocity_check", flowvel, ld, field_ld_grp);
//     }
  }

  h5::Group vessel_group = g.create_group("vessels");
  WriteHdfGraph(vessel_group, *vl);
  vl->Ld().WriteHdf(g.create_group("lattice"));


  { 
    int ecnt = vl->GetECount();
    DynArray<float> tmp1(ecnt);
    for (int i=0; i<ecnt; ++i)
    {
      double coeffa, coeffb, vala, valb;
      vessel_sources_f->GetVesselInfo(vl->GetEdge(i), coeffa, coeffb, vala, valb);
      tmp1[i] = 0.5*(coeffa+coeffb);
    }
    h5::create_dataset(vessel_group, "edges/wall_conductivity", tmp1); 
  }

#ifdef ANALYZE_FLOW_CHANGE
  {
    int ecnt = vl->GetECount();
    int ncnt = vl->GetNCount();
    for (int i=0; i<ecnt; ++i)
      org_flow[i] = vl->GetEdge(i)->q - org_flow[i];
    for (int i=0; i<ncnt; ++i)
      org_press[i] = vl->GetNode(i)->press - org_press[i];
    h5::create_dataset_range(vessel_group,"edges/flow_if_delta", org_flow);
    h5::create_dataset_range(vessel_group,"nodes/pressure_if_delta", org_press);
    org_flow.clear();
    org_press.clear();
  }
#endif
#endif
}



void IffDrugApp3d::DoIffCalc()
{
  my::log().push("iff:");
  UncoupledIfModel model;
  model.Init(grid, mtboxes, *vessel_sources_f, *lymph_sources_f, *tissue_f, iff_params);
  model.calculate();
  model.getOutState(*iffstate, IfFlowState::ALL);
  CopyBorder(iffstate->pfield, dim, 1);
  cout << format("Max Velocity: %f, CFL dt: %f") % iffstate->max_velocity % (grid.ld.Scale()/(2*dim*iffstate->max_velocity)) << endl;
  cout << format("Max Memory Usage: %s") % FormatMemSize(GetMemoryUsage().rss_peak) << endl;
  my::log().pop();
}


void IffDrugApp3d::WriteDrugOutput(double t, const IfDrug::Calculator::State &conc_field, const IfDrug::Calculator& model, const ptree& params)
{
  H5::H5File f(params.get<string>("fn_out"), H5F_ACC_TRUNC);
  try
  {
    H5::Group ift = f.openGroup("parameters/ift");
  }
  catch( H5::Exception error)
  {
    H5::Group out =f.createGroup("parameters/ift");
    ift_params.WriteH5(out);
  }
//   if (output_number == 0 && !f.root().exists("parameters/ift"))
//     ift_params.WriteH5(f.root().create_group("parameters/ift"));

  cout << format("drug hdf output t=%f -> %s") % t % f.getFileName() << endl;
  H5::Group g = f.createGroup((format("out%04i") % output_number).str());

//   g.attrs().set("time", t);
// 
//   g.attrs().set("real_time", (my::Time() - real_start_time).to_s());
//   MemUsage memusage = GetMemoryUsage();
//   g.attrs().set<uint64>("mem_vsize", memusage.vmem_peak);
//   g.attrs().set<uint64>("mem_rss", memusage.rss_peak);

  model.writeH5(f, g, conc_field, t, f.openGroup("field_ld"));

  ++output_number;
}

/** @brief this is called by the user
 */
int  IffDrugApp3d::Main(const ptree &read_params, const string &outfilename, py::object &drug_measurement_function)
{
  { //****** read in parameters *******
    fn_out = outfilename;
    all_pt_params.put_child("iff", IffParams().as_ptree());
    all_pt_params.put_child("ift", IfDrug::Params().as_ptree());
    all_pt_params.put_child("vessels", VesselModel1::Params().as_ptree());
    all_pt_params.put_child("calcflow", BloodFlowParameters().as_ptree());
    all_pt_params.put_child("tumor", NewBulkTissueModel::Params().as_ptree());
    #define DOPT(name) all_pt_params.put(#name, name)
    DOPT(out_intervall);
    DOPT(fn_out);
    DOPT(fn_tumor);
    //DOPT(tend);
    DOPT(message);
    DOPT(h5_path_vessel);
    DOPT(h5_path_lattice);
    DOPT(h5_path_tumor);
    DOPT(parameterset_name);
    #undef DOPT
    //all_pt_params.put("num_threads", 1);

    ptree default_params = all_pt_params;
    { 
      boost::property_tree::ptree p, &iter_array = boost::property_tree::require_child(default_params, "out_times");
      for (float t = 0.; t <= 5.; t += 1.)  // put some default output times in there
      { // i can't have them in the all_pt_params, or they will be mixed up with the read values
	p.put_value(t);
	iter_array.push_back(std::make_pair("", p));
      }
    }
    ptree vararg_params = boost::property_tree::make_ptree("out_times","");

    // copy & pasted from HandleSimulationProgramArguments
    ptree unknown_params = boost::property_tree::subtract(read_params, default_params);
    unknown_params = boost::property_tree::remove(unknown_params, vararg_params);
    // fail if a parameter is not in the default list
    if (unknown_params.begin() != unknown_params.end())
    {
      std::cerr << "--- WARNING: Unknown input parameters detected ---" << endl;
      boost::property_tree::write_info(std::cerr, unknown_params);
      std::cerr << "-------------------------------------------------" << endl;
      throw std::invalid_argument("parameter fail!");
    }
    // update defaults
    boost::property_tree::update(all_pt_params, read_params);

    // read parameters from ptree into variables
    #define DOPT(name) boost::property_tree::get(name, #name, all_pt_params)
    DOPT(out_intervall);
    DOPT(fn_out);
    DOPT(fn_tumor);
    //DOPT(tend);
    DOPT(message);
    DOPT(h5_path_vessel);
    DOPT(h5_path_lattice);
    DOPT(h5_path_tumor);
    DOPT(parameterset_name);
    #undef DOPT
    iff_params.assign(all_pt_params.get_child("iff"));
    ift_params.assign(all_pt_params.get_child("ift"));
    // obtain list of output times
    out_times.clear();
    {  const boost::property_tree::ptree &iter_array = boost::property_tree::require_child(all_pt_params,"out_times");
      BOOST_FOREACH(const ptree::value_type &p, iter_array)
      {
        out_times.push_back(p.second.get_value<float>());
      }
    }
    std::sort(out_times.begin(), out_times.end());
  }

  if(fn_tumor.empty() ) {
    cerr << "Error: no file specified" << endl;
    return -1;
  }
  if(fn_out.empty() || fn_out == fn_tumor) {
    cerr << "Error: no output file specified" << endl;
    return -1;
  }

  //my::SetNumThreads(all_pt_params.get<int>("num_threads", 1)); // set this in python
  //get time
  real_start_time = my::Time();
  /** @brief this initializes a new simulation state
   * checking in parameters, vessel and tumor files
   * 
   * if everthing goes right, InitNewState() returns true and 
   * the simulation can starting
   * 
   * if errors occure, InitNewState() returns false and 
   * the whole simulation exits with state -1 !!!
   */
  if(!InitNewState())
  {
    return -1;
  }
  cout << "init successfully passed" << endl;
  real_start_time = my::Time();
  /** @brief calculating the pressure and flow field
   */
  DoIffCalc();
  {
    //write some output
    H5::H5File file(fn_out, H5F_ACC_RDWR);
    MeasureIfFlowState(file.openGroup("/")); 
  } // don't keep file
  
#ifdef USE_IFDRUGSIM
  /** from here on we calculate drug related stuff */
  if (out_times.size() > 0)
  {
    IfDrug::VesselDrugExavasationModel drugexavasationmodel;
    drugexavasationmodel.Init(*vl, grid, vesselsInBoxes, ift_params);

    IfDrug::Calculator drugcalc;
    drugcalc.Init(grid, mtboxes, *iffstate, tissueCompositionCallback, drugexavasationmodel, ift_params);

    IfDrug::Calculator::State state; 
    state.init(2, grid.ld, grid.dim, 3);

    real_start_time = my::Time();

    //store some more parameters, which is never bad
    ptree pt;
    pt.put("fn_out", fn_out);
    pt.put("out_intervall", out_intervall);
    //pt.put("tend", tend);
    pt.put("save_hdf", true);
    pt.put("hdf_clear_file", false);

    //boost::function4<void, double, const IfDrug::Calculator::State&, const IfDrug::Calculator&, const ptree&> observer = boost::bind(&IffDrugApp3d::WriteDrugOutput, this, _1, _2, _3, _4);
    //::run_model(state, drugcalc, observer, pt);

    // manual loop for variable output time step
    Steppers::StepControl ctrl;
    /* first dimension is the field index 
     * (i.e. tissue compartment or w/e), 
     * higher dimensions are space dimensions*/
    Vec<Py_ssize_t, 4> dim; dim[0] = 2; 
    for (int i=0; i<3; ++i)
    {
      dim[i+1] = grid.ld.Size()[i];
    }
#if BOOST_VERSION>106300
    np::ndarray a_drug = np::zeros(py::make_tuple(dim[0],dim[1],dim[2],dim[3]), np::dtype::get_builtin<float>());
#else
    np::arrayt<float> a_drug = np::zeros(4, dim.data(), np::getItemtype<float>());
#endif
    /**
      * main drug loop
      */
    double next_time = 0, next_time_py_measure = 0;
    while (true)
    {
      //basically write the concentration field
      if (ctrl.t >= next_time - 0.1 * ctrl.dt)
      {
        WriteDrugOutput(ctrl.t, state, drugcalc, pt);
        if (output_number >= out_times.size())
          ctrl.stopflag = true;
        else
          next_time = out_times[output_number];
      }
      if (my::checkAbort())
	ctrl.stopflag = true;
      
      /** @brief later added for more frequent measurement via python
       * see submitIff.py --> class Measure
       */
      if (!drug_measurement_function.is_none() && (next_time_py_measure - 0.1 * ctrl.dt <= ctrl.t))
      {
#ifdef DEBUG
	printf("next_time_py_measure: %f, ctrl.dt: %f, ctrl.t %f, next_time_py_measure - 0.1 * ctrl.dt: %f\n", next_time_py_measure,ctrl.dt,ctrl.t,next_time_py_measure - 0.1 * ctrl.dt);
#endif        
	Int3 org = grid.ld.Box().min;
        FOR_BBOX3(p, grid.ld.Box())
        {
          float v[2] = { state.field[0](p),
                         state.field[1](p), };
          float loc_phi_cells = phi_cells(p);
          for (int i=0; i<2; ++i)
          {
            v[i] /= loc_phi_cells + 1.e-13;
#if BOOST_VERSION>106300
            a_drug[i][p[0]-org[0]][ p[1]-org[1]][p[2]-org[2]] = v[i];//shift to center?
#else
            a_drug(i, p[0]-org[0], p[1]-org[1], p[2]-org[2]) = v[i];//shift to center?
#endif
          }
        }
        py::tuple py_stats = py::make_tuple(ctrl.t, ctrl.dt);
	/* 
	 * tight shit, is this a python function acting inside C++? Yes it is!
	 */
#if BOOST_VERSION>106300
        drug_measurement_function(py_stats, a_drug);
#else
        drug_measurement_function(py_stats, a_drug.getObject());
#endif
        next_time_py_measure += out_intervall;
      }
      
      if (ctrl.stopflag) break;
      ctrl.dt = next_time - ctrl.t;
      /** @brief here the c++ model is iterated!!!!
       */
      ctrl =  drugcalc.doStep(state, ctrl);
    }
  }
#endif
  // if we are here, everything worked smoothly and we gladdly return state 0
  return 0;
}
