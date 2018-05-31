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
#include "bulktissue-with-vessels.h"

BulkTissue::Params::Params()
{
  #define DOPT(name, value) name = value
  DOPT(override_scale, -1);
  DOPT(rGf, 200);
  DOPT(out_intervall, 10);
  DOPT(tend, 100);
  DOPT(message, string());
  DOPT(lattice_scale, 30);
  DOPT(gf_production_threshold, 0.1);
//   DOPT(hematocrit_init, 0.45);
//   DOPT(reference_intercapillary_distance, 80);
//   DOPT(o2_level_normal, 0.6);
//   for (int i=0; i<3; ++i)
//   {
//     o2_range[i] = 100.;
//     o2_cons_coeff[i] = 1./my::sqr(o2_range[i]);
//   }
//   capillary_wall_permeability = O2Model::CalcHomogeneousCoeffOxy(o2_cons_coeff[0], o2_level_normal, 4., reference_intercapillary_distance);
  DOPT(create_single_output_file, true);
  DOPT(vessel_volume_exclusion, false);
  vesselfile_ensemble_index = 0;
  DOPT(paramset_name, "aname");
  #undef DOPT
}
void BulkTissue::Params::update_ptree(ptree& dst, const ptree& src)
{
  boost::property_tree::update(dst, src);
  for (int i=0; i<3; ++i)
  {
    double range = 0., cons_coeff = 0.;
    O2Model::assignRangeParam(range, cons_coeff, "_"+tissue_name[i], src.get_child("simple_o2"));
    dst.put("o2_range_"+tissue_name[i], range);
    dst.put("o2_cons_coeff_"+tissue_name[i], cons_coeff);
  }
  double capillary_wall_permeability = 0.;
  double o2_level_normal = 0.;
  O2Model::assignWallPermeabilityCoefficient(capillary_wall_permeability,
					      o2_level_normal,
					      dst.get<double>("o2_cons_coeff_"+tissue_name[TISSUE]),
					      dst.get<double>("reference_intercapillary_distance"),
					      4.,
					      src.get_child("simple_o2"));
  dst.put("capillary_wall_permeability", capillary_wall_permeability);
  dst.put("o2_level_normal", o2_level_normal);
}

void BulkTissue::Params::assign(const ptree& pt)
{
    #define DOPT(name) boost::property_tree::get(name, #name, pt)
    DOPT(override_scale);
    DOPT(rGf);
    DOPT(out_intervall);
    DOPT(tend);
    DOPT(message);
    DOPT(lattice_scale);
    DOPT(lattice_size);
    DOPT(gf_production_threshold);
//     DOPT(hematocrit_init);
//     DOPT(reference_intercapillary_distance);
//     for (int i=0; i<3; ++i) {
//       boost::property_tree::get(o2_cons_coeff[i], "o2_cons_coeff_"+tissue_name[i], pt);
//       boost::property_tree::get(o2_range[i], "o2_range_"+tissue_name[i], pt);
//     }
//     DOPT(capillary_wall_permeability);
//     DOPT(o2_level_normal);
    DOPT(create_single_output_file);
    DOPT(vessel_volume_exclusion);
    DOPT(fn_out);
    DOPT(fn_vessel);
    DOPT(paramset_name);
    #undef DOPT
    const auto bfparamsPtree = pt.get_child_optional("calcflow");
    if (bfparamsPtree) bfparams.assign(*bfparamsPtree);
}

ptree BulkTissue::Params::as_ptree() const
{
    boost::property_tree::ptree pt;
    #define DOPT(name) pt.put(#name, name)
    DOPT(override_scale);
    DOPT(rGf);
    DOPT(out_intervall);
    DOPT(tend);
    DOPT(message);
    DOPT(lattice_scale);
    DOPT(lattice_size);
    DOPT(gf_production_threshold);
//     DOPT(hematocrit_init);
//     DOPT(reference_intercapillary_distance);
//     DOPT(o2_level_normal);
//     DOPT(capillary_wall_permeability);
//     #define DOPT2(name, i) pt.put(#name"_"+tissue_name[i], name[i])
//     for (int i=0; i<3; ++i)
//     {
//       DOPT2(o2_range, i);
//       DOPT2(o2_cons_coeff, i);
//     }
    DOPT(create_single_output_file);
    DOPT(vessel_volume_exclusion);
    DOPT(fn_out);
    DOPT(fn_vessel);
    DOPT(paramset_name);
    #undef DOPT
    #undef DOPT2
    pt.put_child("calcflow", bfparams.as_ptree());
    return pt;
}

//int NewTumorSim::run(int argc, char **argv)
int BulkTissue::NewTumorSim::run(const ptree &pparams)
{
  // direct cout through log
  cout.rdbuf(my::log().rdbuf());
  // handle parameters
  all_pt_params = BulkTissue::Params().as_ptree();
  //params = BulkTissue::Params();
  //all_pt_params = pparams;
  //all_pt_params.put("num_threads", 1);
  all_pt_params.put("fake_height_2d", 200.);
  all_pt_params.put_child("vessels", VesselModel1::Params().as_ptree());
  all_pt_params.put_child("tumor", NewBulkTissueModel::Params().as_ptree());
  all_pt_params.put_child("simple_o2", O2Model::SimpleO2Params().as_ptree() );
#ifdef USE_ADAPTION
  all_pt_params.put_child("adaption", Adaption::Parameters().as_ptree());
#endif
  { 
    //boost::optional<ptree> read_params = HandleSimulationProgramArguments(all_pt_params, argc, argv);
    //boost::optional<ptree> read_params = all_pt_params;
    //if (!pparams) return 0;
    Params::update_ptree(all_pt_params, pparams);
    this->params.assign(all_pt_params); // this sets some options based on which parameters are supplied
  }
  //HACK2018
  //my::SetNumThreads(all_pt_params.get<int>("num_threads"));

  my::log().push(" init: ");
  {
    H5::H5File readInfile = H5::H5File(params.fn_vessel, H5F_ACC_RDONLY);
    ptree pt;
    //factor by which lattice is subdivided for tumor growth
    pt.put("scale subdivide", 10.);
    pt.put("scale override", params.override_scale);
    //pt.put("filter", true); // does not help, is also filtered in oxygen model
    H5::Group h5_vessels = readInfile.openGroup("/vessels");
    std::unique_ptr<VesselList3d> vl = ReadVesselList3d(h5_vessels,pt);
    // adjust vessel list ld
    const Float3 c = 0.5 * (vl->Ld().GetWorldBox().max + vl->Ld().GetWorldBox().min);
    vl->SetDomainOrigin(vl->Ld().LatticeToWorld(Int3(0))-c);

    cout << "--------------------"<< endl;
    cout << "Vessel Lattice is: " << endl;
    vl->Ld().print(cout); cout  << endl;
    H5::Group h5params = readInfile.openGroup("/parameters");
    string message;
    readAttrFromH5(h5params, string("MESSAGE"),message);
    params.vesselfile_message = message;
    int index;
    readAttrFromH5(h5params, string("ENSEMBLE_INDEX"),index);
    params.vesselfile_ensemble_index = index;
    //params.vesselfile_message = file.openGroup("/parameters").attrs().get<string>("MESSAGE");
    //params.vesselfile_ensemble_index = file.root().open_group("parameters").attrs().get<int>("ENSEMBLE_INDEX");
    
    state.vessels.reset(vl.release());
    last_vessels_checksum = -1;
  }

  {
    Int3 s = params.lattice_size;
    int dim = s[2]<=1 ? (s[1]<=1 ? 1 : 2) : 3;
    LatticeDataQuad3d field_ld;
    Bool3 centering = Bool3::mapIndex([=](int i) { return i<dim; });
    field_ld.Init(params.lattice_size, params.lattice_scale);
    field_ld.SetCellCentering(centering);
    field_ld.SetOriginPosition(-field_ld.GetWorldBox().max.cwiseProduct(centering.cast<float>()) * 0.5); // set origin = lower left side
    //grid.init(field_ld, dim);
    grid = ContinuumGrid(field_ld, dim);
    mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
    
    cout << "--------------------"<< endl;
    cout << format("Tumor Lattice is %i dimensional ") % dim << endl;
    field_ld.print(cout); cout  << endl;
    cout << "--------------------"<< endl;
    cout << "With DomainDecomposition:"<< endl;
    mtboxes.print(cout);
    cout << "--------------------"<< endl;
  }

  {
    //fieldinterp_extrapolate.init(CONT_EXTRAPOLATE);
    //fieldinterp_const.init(CONT_CONST, make_ptree("value", 0.)("gradient", Float3(0.)));
    
    VesselModel1::Callbacks callbacks;
    callbacks.getGf = boost::bind(&NewTumorSim::getGf, boost::ref(*this), _1);
    callbacks.getPress = boost::bind(&NewTumorSim::getPress, boost::ref(*this), _1);
    callbacks.getTumorDens = boost::bind(&NewTumorSim::getTumorDens, boost::ref(*this), _1);
    callbacks.getGfGrad = boost::bind(&NewTumorSim::getGfGrad, boost::ref(*this), _1);

    VesselModel1::Params params_vess; 
    params_vess.assign(all_pt_params.get_child("vessels"));
    
    /* need to compute flow because shearforce must be
     * known and be consistent with current parameters.
     * Shear force is used e.g. in model.Init to initialize
     * f_initial. */
    CalcFlow(*state.vessels, params.bfparams); 
    vessel_model.Init(state.vessels.get(), params_vess, callbacks);
  }

  gf_model.init(grid, mtboxes, all_pt_params);
  gf_model.initField(state.gffield);
  
  oxyops.init(mtboxes, grid.Box(), grid.dim, 2);
  oxyops.init(state.o2field);
  last_chem_update = -1;

  UpdateVesselVolumeFraction();
  std::cout<< "vessel volume fraction updated" << std::endl;

  auto obtain_vessel_volume = [=](const BBox3 &bbox) -> Array3df
  {
    return this->vessel_volume_fraction;
  };

  auto obtain_oxygen_field = [=](const BBox3 &bbox) -> Array3df
  {
    return this->state.o2field;
  };
  
  tumor_model.init(grid, mtboxes, state.tumor,
                      params.vessel_volume_exclusion ? obtain_vessel_volume : NewBulkTissueModel::DataRequestCallback(),
                      obtain_oxygen_field,
                      all_pt_params.get_child("tumor", ptree()));
  tumor_radius_estimate = EstimateTumorRadius(grid, mtboxes, state.tumor.ls);
  
  std::cout<< "tumor_model initialized" << std::endl;
  
  calcChemFields();
  
  std::cout<< "calcChemFields executed" << std::endl;

//   {
//     //store this now for the first evaluation of the vessel model
//     const TumorModel::OutState &ost = tumor_model.getOutState(TumorModel::OutState::THETA | TumorModel::OutState::PRESSURE);
//     tumor_phase = ost.theta; // shallow copy
//     cell_pressure = ost.pressure;
//   }

  my::log().pop();

  output_num = 0;
  iteration_num = 0;
  real_start_time = my::Time();
  failFlag = false;
  checkStop = false;

  auto doStep = [&](NewSteppers::StepControl &ctrl) -> bool
  {
    my::log().push(str(format("%5i: ") % this->iteration_num));
    bool ret = this->doStep(ctrl);
    my::log().pop();
    ++this->iteration_num;
    return ret;
  };

  auto doObserve = [&](double t) {
    writeOutput(t);
  };
  
  //this starts the simulation
  std::cout<< "bulktissue simulation started" << std::endl;
  NewSteppers::run(doStep, doObserve, all_pt_params);

  if (failFlag)
  {
    std::cout<< "failFlag case" << std::endl;
    //intel bug?
    //writeOutput(std::numeric_limits<double>::quiet_NaN());
    writeOutput(std::numeric_limits<double>::max());
  }

  return 1;
}


bool BulkTissue::NewTumorSim::doStep(NewSteppers::StepControl &ctrl)
{
  ctrl.dt = std::min(1., ctrl.dt);
  
  // advance the subsystem which is most back in time
  if (vessel_step_ctrl.t < tumor_step_ctrl.t)
  {
    advanceVesselState();
#ifndef NDEBUG
    std::cout << "adavanceVesselState done" << std::endl;std::cout.flush();
#endif
  }
  else
  {
    advanceTumorState(ctrl.dt);
#ifndef NDEBUG
    std::cout << "adavanceTumorState done" << std::endl;std::cout.flush();
#endif
  }

  double time = my::min(vessel_step_ctrl.t, tumor_step_ctrl.t);
  ctrl.dt = time - ctrl.t;
  ctrl.t  = time;
  
  // chemical fields
  if (last_chem_update + std::min(0.99, 2.*tumor_model.max_dt_vel) < time)
  {
    calcChemFields();
    last_chem_update = time;
  }
#ifndef NDEBUG
  std::cout << "chemicals done" << std::endl;std::cout.flush();
#endif
  bool stopFlag = false;
  if (checkStop)
  {
    tumor_radius_estimate = EstimateTumorRadius(grid, mtboxes, state.tumor.ls);
    double size_limit = 0.5*maxCoeff(Size(grid.ld.GetWorldBox())) * 0.90;
    if (tumor_radius_estimate >  size_limit)
      stopFlag = true;
    checkStop = false;
  }

  cout << format("master t = %f, dt = %f") % ctrl.t % ctrl.dt << endl;
#ifndef NDEBUG
  std::cout << "stopFlag" << stopFlag << std::endl;std::cout.flush();
#endif
  return !stopFlag;
}


void BulkTissue::NewTumorSim::advanceTumorState(double dt)
{
  my::log().push("tum:");
  tumor_step_ctrl.dt = dt;
  // euler_dt remains from last step!
  
  UpdateVesselVolumeFraction();
  
  tumor_model.doPreStep(state.tumor, tumor_step_ctrl);

  auto splitStepCallback = [=](int num, NewBulkTissueModel::State &state, NewSteppers::StepControl &ctrl) -> bool
  {
    return this->tumor_model.doSplitStep(num, state, ctrl);
  };
  auto cloneStateCallback = [=](NewBulkTissueModel::State &dst, const NewBulkTissueModel::State &src)
  {
    this->tumor_model.cloneState(dst, src);
  };
  
  NewSteppers::doStrangeSplittingSteps<NewBulkTissueModel::State>(state.tumor, tumor_step_ctrl, 2, splitStepCallback, cloneStateCallback);

  state.tumor_checksum++;

  //if (this->iteration_num % 10 == 0)
  {
    cout << format("t = %f, euler_dt = %f")  % tumor_step_ctrl.t % tumor_step_ctrl.euler_dt;
    cout << format(", by vel: %f, by diff = %f, by src = %f") % tumor_model.max_dt_vel % tumor_model.max_dt_diff % tumor_model.max_dt_src;
    if (tumor_model.params.surface_tension>0.)
      cout << format(", by stf = %f")  % tumor_model.max_dt_stf;
    cout << format(", tum rad. = %f") % tumor_radius_estimate;
    cout << endl;
  }
  my::log().pop();
}


void BulkTissue::NewTumorSim::advanceVesselState()
{
  my::log().push("ves:");
  
  vessel_step_ctrl.dt = 1.;
  vessel_step_ctrl.euler_dt = my::Smallest<double>(1.);
  
  const NewBulkTissueModel::OutState &ost = tumor_model.getOutState(state.tumor,
    NewBulkTissueModel::OutState::THETA | NewBulkTissueModel::OutState::PRESSURE);
  tumor_phase = ost.theta; // shallow copy
  cell_pressure = ost.pressure;  
  
  // do step
  vessel_model.DoStep(1.,nullptr);
  CalcFlow(*state.vessels, params.bfparams);

  // update data for interaction
  tumor_phase.clear();
  cell_pressure.clear();

  vessel_step_ctrl.update();
  state.vessels_checksum++;

  cout << format("t = %f, dt = %f")  % vessel_step_ctrl.t % vessel_step_ctrl.dt << endl;
  my::log().pop();
}


void BulkTissue::NewTumorSim::UpdateVesselVolumeFraction()
{
  my::LogScope log_push_(my::log(), "ves:");
  
  // update only when the vessel system state has changed since last update 
  if (!vessel_volume_fraction.empty() && last_vessels_checksum==state.vessels_checksum)
    return;

  cout << "volume fraction and o2 sources update!" << endl;
  
  VesselsInBoxes vessboxes;
  SortVesselsIntoMtBoxGrid(grid.ld, *state.vessels, 2, mtboxes, vessboxes);
  
  // reinit and fill the array
  Array3dOps<float> ops(mtboxes, grid.Box(), grid.dim, 0);
  ops.init(vessel_volume_fraction, true);
  ops.init(vessel_o2src_clin);
  ops.init(vessel_o2src_crhs);

  #pragma omp parallel
  {
    VesselVolumeGenerator volumegen(*state.vessels, grid.ld, grid.dim, make_ptree("samples_per_cell", 20));
    
    BOOST_FOREACH(const DomainDecomposition::ThreadBox &bbox, mtboxes.getCurrentThreadRange())
    {
      int dummy_;
      volumegen.Fill(bbox, vessel_volume_fraction, vessboxes[bbox.global_index], dummy_);
      O2Model::AddSourceDistribution( bbox, 
                                      grid.ld, 
                                      grid.dim, 
                                      vessel_o2src_clin, 
                                      vessel_o2src_crhs, 
                                      state.vessels->Ld(), 
                                      vessboxes[bbox.global_index], 
                                      all_pt_params);
    }
  }

  last_vessels_checksum = state.vessels_checksum;
}




void BulkTissue::NewTumorSim::insertO2Coefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder &mb)
{
  Array3d<float> l_coeff(bbox),
                 rhs(bbox);
  //O2Model::AddSourceDistribution(bbox, field_ld, l_coeff, rhs, vl->Ld(), vessboxes[box_index], all_pt_params);
  l_coeff[bbox] += vessel_o2src_clin[bbox];
  rhs[bbox] += vessel_o2src_crhs[bbox];

  Array3df phases[3];
  tie(phases[TISSUE], phases[TCS], phases[DEAD]) = tumor_model.getTissuePhases(bbox, state.tumor);

  float hemostatic_cell_frac_norm[3];
  hemostatic_cell_frac_norm[TISSUE] = 1./tumor_model.params.ncells_norm;
  hemostatic_cell_frac_norm[TCS] = 1./tumor_model.params.ncells_tumor;
  hemostatic_cell_frac_norm[DEAD] = 1./tumor_model.params.ncells_tumor;
  
  mb.AddDiffusion(bbox, ConstValueFunctor<float>(1.), -1.);
  
  FOR_BBOX3(p, bbox)
  {
    const float loc_phases[3] = { phases[0](p), phases[1](p), phases[2](p) };
    // we must ensure that the resulting diffusion distance agrees with the parameters which define it
    // by l = sqrt(c/D). So the effect that the cell volume fraction is about 0.5 must be canceled out.
    float loc_l_coeff = 0.;
    for (int i=0; i<3; ++i)
    {
      loc_l_coeff -= hemostatic_cell_frac_norm[i] * loc_phases[i] * params.o2Params.o2_cons_coeff[i];
      
    }
    mb.AddLocally(p, -loc_l_coeff - l_coeff(p), -rhs(p));
  }

}

 

void BulkTissue::NewTumorSim::calcChemFields()
{
  {
    my::log().push("o2:");

    ptree pt_params = make_ptree("preconditioner","multigrid")("output", 1)("max_iter", 100);
    StationaryDiffusionSolve(grid, mtboxes, state.o2field, boost::bind(&NewTumorSim::insertO2Coefficients, this, _1, _2, boost::cref(state), _3) , pt_params);
    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
      {
        CopyBorder(state.o2field, bbox, grid.Box(), grid.dim, 2);
      }
    }
#ifdef DEBUG
    cout << "stats: " << state.o2field.valueStatistics() << endl;
#endif
    my::log().pop();
  }

  {
    my::log().push("gf:");
    cout << "update" << endl;

    Array3df src;
    Array3dOps<float>(mtboxes, grid.Box(), grid.dim, 0).init(src, false);

    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bb, mtboxes.getCurrentThreadRange())
      {
        Array3df phases[3];
        tie(phases[TISSUE], phases[TCS], phases[DEAD]) = tumor_model.getTissuePhases(bb, state.tumor);
        
        FOR_BBOX3(p, bb)
        {
          float o2 = state.o2field(p);
          src(p) = (o2 < params.gf_production_threshold) ? phases[TCS](p) : 0;
        }
      }
    }
    gf_model.update(state.gffield, src);
#ifdef DEBUG
    cout << "stats: " << state.gffield.valueStatistics() << endl;
#endif
    my::log().pop();
  }
  state.chem_checksum++;
}


float BulkTissue::NewTumorSim::getTumorDens(const Float3 &pos) const
{
  return FieldInterpolate::ValueAveraged(tumor_phase, grid.ld, FieldInterpolate::Extrapolate(), pos);
  //return fieldinterp_extrapolate.value(tumor_phase, grid.ld, pos);
}

float BulkTissue::NewTumorSim::getPress(const Float3 &pos) const
{
  return FieldInterpolate::ValueAveraged(cell_pressure, grid.ld, FieldInterpolate::Extrapolate(), pos);
  //return fieldinterp_extrapolate.value(cell_pressure, grid.ld, pos);
}

float BulkTissue::NewTumorSim::getGf(const Float3 &pos) const
{
  return FieldInterpolate::ValueAveraged(state.gffield, grid.ld, FieldInterpolate::Const(0.f), pos);
  //return fieldinterp_extrapolate.value(state.gffield, grid.ld, pos);
}

Float3 BulkTissue::NewTumorSim::getGfGrad(const Float3 &pos) const
{
  return FieldInterpolate::Gradient(state.gffield, grid.ld, FieldInterpolate::Extrapolate(), pos);
  //return fieldinterp_extrapolate.gradient(state.gffield, grid.ld, _pos);
}


void BulkTissue::NewTumorSim::writeOutput(double time)
{
  my::log().push("out:");
  H5::H5File f;
  H5::Group g, gout;
  //h5::Attributes a;
  if (!params.create_single_output_file)
  {
    const string fn = str(format("%s-%03i.h5") % RemoveExtension(params.fn_out) % output_num);
    cout << format("output %i -> %s") % output_num % fn << endl;
    f = H5::H5File(fn, H5F_ACC_RDWR);
    gout = f.openGroup("/");
  }
  else
  {
    const string fn = RemoveExtension(params.fn_out)+".h5";
    cout << format("output %i +-> %s") % output_num % fn << endl;
    f = H5::H5File(fn, output_num==0 ? H5F_ACC_TRUNC : H5F_ACC_RDWR);
    //f = h5::File(fn, output_num==0 ? "w" : "a");
    //gout = f.root().require_group(str(format("out%04i") % output_num));
    //gout = f.openGroup(str(format("out%04i") % output_num));
    gout = f.createGroup(str(format("out%04i") % output_num));
  }
  if (output_num == 0 || !params.create_single_output_file)
  {
    H5::Group root = f.openGroup("/");
    writeAttrToH5(root, string("MESSAGE"), params.message);
    writeAttrToH5(root, string("VESSELTREEFILE"), params.fn_vessel);
    writeAttrToH5(root, string("OUTPUT_NUM"), output_num);
    writeAttrToH5(root, string("OUTPUT_NAME"), params.fn_out);
    writeAttrToH5(root, string("VESSELFILE_MESSAGE"), params.vesselfile_message);
    writeAttrToH5(root, string("VESSELFILE_ENSEMBLE_INDEX"), params.vesselfile_ensemble_index);
    H5::Group h5_params = f.createGroup("/parameters");
    H5::Group h5_params_vessel = h5_params.createGroup("vessels");
    WriteHdfPtree(h5_params_vessel, vessel_model.params.as_ptree());
    H5::Group h5_params_tumor = h5_params.createGroup("tumor");
    WriteHdfPtree(h5_params_tumor, tumor_model.pt_params);
    WriteHdfPtree(h5_params, params.as_ptree());
//     a= f.root().attrs();
//     a.set("MESSAGE",params.message);
//     a.set("VESSELTREEFILE", params.fn_vessel);
//     a.set("OUTPUT_NUM",output_num);
//     a.set("OUTPUT_NAME", params.fn_out);
//     a.set("VESSELFILE_MESSAGE", params.vesselfile_message);
//     a.set("VESSELFILE_ENSEMBLE_INDEX", params.vesselfile_ensemble_index);
//     // parameters
//     g = f.root().create_group("parameters");
//     WriteHdfPtree(g.create_group("vessels"), vessel_model.params.as_ptree());
//     WriteHdfPtree(g.create_group("tumor"), tumor_model.pt_params);
//     WriteHdfPtree(g, params.as_ptree());
  }
  // snapshot attributes
  MemUsage memusage = GetMemoryUsage();
  writeAttrToH5(gout, string("time"), time);
  writeAttrToH5(gout, string("real_time"),(my::Time() - real_start_time).to_s());
  writeAttrToH5(gout, string("mem_vsize"), (int)memusage.vmem_peak);
  writeAttrToH5(gout, string("mem_rss"), (int)memusage.rss_peak);
  
//   a = gout.attrs();
//   a.set("time", time);
//   a.set("real_time", (my::Time() - real_start_time).to_s());
//   a.set<uint64>("mem_vsize", memusage.vmem_peak);
//   a.set<uint64>("mem_rss", memusage.rss_peak);
  bool has_grp = false;
  // vessels
  H5::Group h5_vessels = gout.createGroup("vessels");
  std::cout << " start Write Vessel list " << std::endl;
  WriteVesselList3d(*state.vessels, h5_vessels);
  std::cout << " done Write Vessel list " << std::endl;
  // tumor
  
  H5::Group ld_group_tum;
  try{
    ld_group_tum = f.openGroup("field_ld");
    has_grp = true;
  }
  catch(H5::Exception e)
  {
    //in the first iteration, there is no group yet!
    ld_group_tum = f.createGroup("field_ld");
    grid.ld.WriteHdfLd(ld_group_tum);
  }
  
  
  //H5::Group ld_group_tum = f.openGroup("field_ld");
//   if (!has_grp)
//     //WriteHdfLd(ld_group_tum, grid.ld);
//     grid.ld.WriteHdfLd(ld_group_tum);
  g = gout.createGroup("tumor");
  tumor_model.writeH5(g, state.tumor, time, ld_group_tum);
  
  // chem fields
  WriteScalarField(gout, "fieldGf", state.gffield, grid.ld, ld_group_tum);
  WriteScalarField(gout, "fieldOxy", state.o2field, grid.ld, ld_group_tum);
  
  // vessel continuum
  UpdateVesselVolumeFraction();
  std::cout << " UpdateVesselVolumeFraction " << std::endl;
  WriteScalarField(gout, "vessel_volume_fraction", vessel_volume_fraction, grid.ld, ld_group_tum);
  WriteScalarField(gout, "oxy_source_lin", vessel_o2src_clin, grid.ld, ld_group_tum);
  ++output_num;
  checkStop = true;
  my::log().pop();
}
