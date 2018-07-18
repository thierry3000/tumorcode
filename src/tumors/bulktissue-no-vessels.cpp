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
#include "bulktissue-no-vessels.h"

/*------------------------------
 *    oxygen
 * -----------------------------*/

class OxySourceModel :  public boost::noncopyable
{
public:
  typedef Array3df State;
  typedef NewSteppers::ImprovedExplicitEuler<OxySourceModel*, Array3dOps<float> > Stepper;

private:
  const ContinuumGrid *grid;
  const DomainDecomposition *mtboxes;
  double decay_rate;
  Stepper stepper;
  typedef boost::function<boost::tuple<Array3df,Array3df,Array3df> (const BBox3 &)> ObtainTumorFractionCallback;
  ObtainTumorFractionCallback obtain_tumor_fraction;

public: 
  void init(const ContinuumGrid &grid_, const DomainDecomposition &mtboxes_,
            State &state, ObtainTumorFractionCallback obtain_tumor_fraction_, const ptree &pt_params)
  {
    grid = &grid_;
    mtboxes = &mtboxes_;
    obtain_tumor_fraction = obtain_tumor_fraction_;
    stepper.set_model(this);
    stepper.ops().init(*mtboxes, grid->Box(), grid->dim, 0);
    
    decay_rate = 1./pt_params.get<double>("o2_source_decay_time");

    stepper.ops().init(state, false);
    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 bbox, mtboxes->getCurrentThreadRange())
      {
        state[bbox].fill(1.);
      }
    }
  }

  void cloneState(State &u, const State &v)
  {
    stepper.ops().initFrom(u, v, ConsMode::AS_COPY);
  }
    
  void calcSlope(Array3df &o2sources, State &slope, NewSteppers::StepControl &ctrl)
  {
    if (slope.empty()) stepper.ops().init(slope, false);

    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 bbox, mtboxes->getCurrentThreadRange())
      {
        Array3df phases[3];
        tie(phases[TISSUE], phases[TCS], phases[DEAD]) = obtain_tumor_fraction(bbox);
        
        FOR_BBOX3(p, bbox)
        {
          const float loc_phases[3] = { phases[0](p), phases[1](p), phases[2](p) };
          float f = (loc_phases[TCS]+loc_phases[DEAD])/(loc_phases[0]+loc_phases[1]+loc_phases[2]);
          slope(p) = -decay_rate * f * o2sources(p);
        }
      }
    }
    ctrl.euler_dt.min(0.9 * 0.5 / decay_rate);
  }

  bool doSplitStep(int num, State &state, NewSteppers::StepControl &ctrl)
  {
    return stepper.doStep(state, ctrl);
  }
};

/**
 * Contains all the things present in bulktissue model
 * like ...
 * paper?
 */
struct Model : public boost::noncopyable
{
  void init(const Int3 &size, double scale, int ndims_, BulkTissueWithoutVessels::State &state, const ptree &pt_params_)
  {
    pt_params = O2Model::SimpleO2Params().as_ptree();
    boost::property_tree::update(pt_params, pt_params_);
    params.assign(pt_params);
    
    my::log().push(" init: ");
    {
      LatticeDataQuad3d ld;
      ld.Init(size, scale);
      ld.SetCellCentering(Bool3(true, ndims_>1, ndims_>2));
      ld.SetOriginPosition(-ld.GetWorldBox().max * 0.5); // set origin = lower left side
      //grid.init(ld, ndims_);
      grid = ContinuumGrid(ld, ndims_);
    }
    //this is for the multithreading
    mtboxes.init(MakeMtBoxGrid(grid.Box(), pt_params.get<Int3>("mtboxsize")));

    mtboxes.print(cout);
    cout << "------------------------------------------" << endl;

    master_state = &state;  

    last_o2_reinit = -std::numeric_limits<double>::max();
    checkStop = false;
    stopFlag  = false;
    failFlag  = false;
    iteration_num = 0;
    tumor_radius_estimate = 0;
    
    if (params.test_obstacle > 0)
    {
      obstacle_volume.initFromBox(grid.Box());
      FloatBBox3 wb = grid.ld.GetWorldBox();
      float width = grid.ld.Scale()*10; 
      float distance = (wb.max-wb.min)[0]*0.25;
      FOR_BBOX3(p, grid.Box())
      {
        Float3 wp = grid.ld.LatticeToWorld(p);
        float f = 1. - std::abs(std::abs(wp[0])-distance)/width;
        // f > 0 within the obstacle
        obstacle_volume(p) = 0.8 * my::smooth_heaviside(f, 0.3f);
      }
    }

    auto obtain_vessel_volume = [=](const BBox3 &bbox) -> Array3df
    {
      return this->obstacle_volume;
    };

    auto obtain_oxygen_field = [=](const BBox3 &bbox) -> Array3df
    {
      return this->master_state->o2field;
    };
    
    bulktissue_model.init(grid, mtboxes, state.tumor_state,
                        params.test_obstacle ? obtain_vessel_volume : NewBulkTissueModel::DataRequestCallback(),
                        obtain_oxygen_field,
                        pt_params.get_child("tumor", ptree()));

    if (params.use_o2_source_decay)
    {
      auto obtain_tumor_fraction = [=](const BBox3 &bbox)
      {
        return this->bulktissue_model.getTissuePhases(bbox, this->master_state->tumor_state);
      };
      oxysource_model.init(grid, mtboxes, state.o2sources, obtain_tumor_fraction, pt_params);
    }

    Array3dOps<float>(mtboxes, grid.Box(), grid.dim, 1).init(state.o2field, false);
    calcO2(state);
    
    my::log().pop();
  }//end init

  void writeH5(H5::H5File f, H5::Group g, BulkTissueWithoutVessels::State &state, double t, int out_num)
  {
    myAssert(&state == master_state);
    try
    {
      H5::Group gparams = f.openGroup("parameters");
    }
    catch(H5::Exception &e)
    {
      H5::Group gparams = f.createGroup("parameters");
      WriteHdfPtree(gparams, pt_params);
    }
//     if (!f.root().exists("parameters"))
//     {
//       h5cpp::Group gparams = f.root().create_group("parameters");
//       WriteHdfPtree(gparams, pt_params);
//     }
    H5::Group ld_group = RequireLatticeDataGroup(f,"field_ld", grid.ld);
    bulktissue_model.writeH5(g, state.tumor_state, t, ld_group);
    writeAttrToH5(g, string("tumor_radius"), estimateRadius(state));
    //g.attrs().set("tumor_radius", estimateRadius(state));
    WriteScalarField(g, "oxy", state.o2field[grid.Box()], grid.ld, ld_group);
    if (params.use_o2_source_decay)
      WriteScalarField(g, "oxy_sources", state.o2sources[grid.Box()], grid.ld, ld_group);
    checkStop = true;
  }

  
  double estimateRadius(const BulkTissueWithoutVessels::State &state)
  {
    return EstimateTumorRadius(grid, mtboxes, state.tumor_state.ls);
  }


  bool doSplitStep(int num, BulkTissueWithoutVessels::State &state, NewSteppers::StepControl &ctrl)
  {
    if (num == 2)
    {
      myAssert(params.use_o2_source_decay);
      return oxysource_model.doSplitStep(0, state.o2sources, ctrl);
    }
    else
    {
      return bulktissue_model.doSplitStep(num, state.tumor_state, ctrl);
    }
  }

  std::pair<bool, string> isStateFucked(BulkTissueWithoutVessels::State &state)
  {
    bool failFlag = false;
    string result;
    tbb::spin_mutex mutex;
    #pragma omp parallel
    {
      string reason;
      Array3df phases[3], other_volume;
      BOOST_FOREACH(const BBox3 bbox, mtboxes.getCurrentThreadRange())
      {
        tie(phases[TISSUE], phases[TCS], phases[DEAD]) = bulktissue_model.getTissuePhases(bbox, state.tumor_state);
        other_volume = bulktissue_model.getObstructionPhase(bbox, state.tumor_state);

        FOR_BBOX3(p, bbox)
        {
          float total = phases[0](p) + phases[1](p) + phases[2](p) + other_volume(p);
          if (total<-0.1 || total>1.1)
          {
            failFlag = true;
            if (reason.empty())
            {
              reason = str(format("at point %s: TISSUE: %f TCS: %f DEAD: %f OBST %f, sum: %f") %
                p % phases[TISSUE](p) % phases[TCS](p) % phases[DEAD](p) % other_volume(p) % total);
            }
          }
        }
      }
      
      if (!reason.empty())
      {
        tbb::spin_mutex::scoped_lock lock(mutex);
        if (!result.empty()) result += '\n';
        result += reason;
      }
    }
    return std::make_pair(failFlag, result);
  }

  bool doStep(BulkTissueWithoutVessels::State &state, NewSteppers::StepControl &ctrl)
  {
    my::LogScope logscope(my::log(), str(format("%5i: ") % (iteration_num++)));

    ctrl.is_rewind_step = false;
    
    // calc o2 if needed
    if ((params.o2_rel_tumor_source_density>=0. && params.o2_rel_tumor_source_density<1.) 
      || (params.o2_range[0]!=params.o2_range[1] || params.o2_range[1]!=params.o2_range[2]))
    {
      if (last_o2_reinit + bulktissue_model.max_dt_vel < ctrl.t)
      {
        calcO2(state);
        last_o2_reinit = ctrl.t;
      }
    }
   
    auto splitStepCallback = [=](int num, BulkTissueWithoutVessels::State &state, NewSteppers::StepControl &ctrl)
    {
      return this->doSplitStep(num, state, ctrl);
    };
    auto cloneStateCallback = [=](BulkTissueWithoutVessels::State &dst, const BulkTissueWithoutVessels::State &src)
    {
      this->bulktissue_model.cloneState(dst.tumor_state, src.tumor_state);
      if (this->params.use_o2_source_decay)
        this->oxysource_model.cloneState(dst.o2sources, src.o2sources);
    };

    bulktissue_model.doPreStep(state.tumor_state, ctrl);
    
    NewSteppers::doStrangeSplittingSteps<BulkTissueWithoutVessels::State>(state, ctrl, params.use_o2_source_decay ? 3 : 2, splitStepCallback, cloneStateCallback);

    if (checkStop)
    {
      tumor_radius_estimate = estimateRadius(state);
      double size_limit = 0.5*maxCoeff(Size(grid.ld.GetWorldBox())) * 0.80;
      if (tumor_radius_estimate >  size_limit)
        stopFlag = true;
      checkStop = false;
    }

    if (this->iteration_num % 10 == 0)
    {
      string failReason;
      boost::tie(failFlag, failReason) = isStateFucked(state);
      if (failFlag)
      {
        stopFlag = true;
        cout << "model failure!" << endl;
        cout << failReason << endl;
      }
    }

    if (this->iteration_num % 10 == 0)
    {
      cout << format("t = %f, euler_dt = %f")  % ctrl.t % ctrl.euler_dt;
      cout << format(", by vel: %f, by diff = %f, by src = %f") % bulktissue_model.max_dt_vel % bulktissue_model.max_dt_diff % bulktissue_model.max_dt_src;
      if (bulktissue_model.params.surface_tension>0.)
        cout << format(", by stf = %f")  % bulktissue_model.max_dt_stf;
      cout << format(", tum rad. = %f") % tumor_radius_estimate;
      cout << endl;
    }
    return !stopFlag;
  }

  /*
   * this is called concurrently from multiple threads
   * it should fill the coefficients for the region defined by bbox
   */
  void insertO2Coefficients(int box_index, const BBox3& bbox, const BulkTissueWithoutVessels::State &state, FiniteVolumeMatrixBuilder &mb)
  {
    //float kcons = params.o2_cons_coeff_normal/bulktissue_model.params.ncells_norm;
    float kvess_norm = params.o2_level_normal * (params.o2_cons_coeff[TISSUE]*0.5) / (1 - params.o2_level_normal);
    float kvess_tum = params.o2_rel_tumor_source_density * kvess_norm;

    Array3df phases[3];
    tie(phases[TISSUE], phases[TCS], phases[DEAD]) = bulktissue_model.getTissuePhases(bbox, state.tumor_state);
    
    mb.AddDiffusion(bbox, ConstValueFunctor<float>(1.), -1.);
    
    FOR_BBOX3(p, bbox)
    {
      //float ncells = ost.phi_cells(p) + 0.2; // 0.2 is a hack to make non-cell filled regions degrade o2
      float q;
      const float loc_phases[3] = { phases[0](p), phases[1](p), phases[2](p) };
      if (params.use_o2_source_decay)
      {
        float f = state.o2sources(p);
        q = (1.-f) * kvess_tum + f * kvess_norm;
      }
      else
      {
        float rel_tum_fraction = (loc_phases[TCS]+loc_phases[DEAD])/(loc_phases[0]+loc_phases[1]+loc_phases[2]);
        q = rel_tum_fraction * kvess_tum + (1.-rel_tum_fraction) * kvess_norm;
      }

      float loc_l_coeff = 0.;
      for (int i=0; i<3; ++i)
      {
        loc_l_coeff -= loc_phases[i] * params.o2_cons_coeff[i];
      }

      float lin_coeff = loc_l_coeff - q;
      float rhs = - q * 1.;

      mb.AddLocally(p, -lin_coeff, -rhs);
    }
  }


  void calcO2(BulkTissueWithoutVessels::State &state)
  {
    /*
     * solve
     * D c - lambda + q (c0 - c) = 0
     */
    my::log().push("o2:");
    ptree pt_params = boost::property_tree::make_ptree<>("output", 1)
                                                        ("max_iter", 200)
                                                        ("max_resid", 1.e-6);
    StationaryDiffusionSolve(grid, mtboxes, state.o2field, boost::bind(&Model::insertO2Coefficients, this, _1, _2, boost::cref(state), _3) , pt_params);

    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
      {
        CopyBorder(state.o2field, bbox, grid.Box(), grid.dim, 1);
      }
    }

    my::log().pop();
  }

  // grid info
  ContinuumGrid grid;
  DomainDecomposition mtboxes;
  // the simulation state
  Array3d<float> obstacle_volume;  // distribution of excluded volume due to vesses, not available for cell population.
  BulkTissueWithoutVessels::State *master_state;             // Cell density (actually cell volume fraction) and levelset function.
  // bookkeeping data
  double last_o2_reinit;           // o2 recomputed only every bulktissue_model.max_dt_vel time units
  bool checkStop, stopFlag, failFlag;  // signal stop and failure state across subsystems
  int iteration_num;
  double tumor_radius_estimate; // last known estimate, mean distance from origin taken over tumor interface.
  // parameters
  O2Model::SimpleO2Params params;
  ptree pt_params;
  // models
  NewBulkTissueModel::Model bulktissue_model; // continuum mechanics model, responsible for propagating State::tumor_state in time + utility functions like determining source strength of cell density
  OxySourceModel oxysource_model;         // calculation of oxygen distribution
};


void BulkTissueWithoutVessels::run(const ptree &params)
{
  Int3 size = params.get<Int3>("lattice_size");
  float scale = params.get<float>("lattice_scale");

  Model model;
  State state;
  NewSteppers::ObserverPde observer(params);
  LatticeDataQuad3d ld;

  model.init(size, scale, size[1]<=1 ? 1 : size[2] <= 1 ? 2 : 3, state, params);
  my::Time t_;

  auto doStep = [&](NewSteppers::StepControl& ctrl) -> bool {
    bool ok = model.doStep(state, ctrl);
    if (model.failFlag) return false;
    ctrl.is_rewind_step = false;
    return ok;
  };
  
  auto doObserve = [&](double t) {
    observer(t, state, model, params);
  };

  NewSteppers::run(doStep, doObserve, params);
  
  if (model.failFlag)
  {
    observer.writeH5("debugout", 0., state, model);
  }
  
  cout << "runtime: " << (my::Time() - t_) << endl;
}

//}

BulkTissueWithoutVessels::SimulationParameters::SimulationParameters()
{
  lattice_size=Int3(200,1,1);
  lattice_scale=30.;
  num_threads = 2;
  tend = 200.;
  out_intervall = 10.;
  save_hdf = true;
  safe_images = false;
  fn_out = "bulktissuetumor";
  fn_vessel = "setme";
  mtboxsize = Int3(64,32,32);
  paramset_name = "aname";
}

void BulkTissueWithoutVessels::SimulationParameters::assign(const ptree& pt)
{
  #define DOPT(name) boost::property_tree::get(name, #name, pt)
  DOPT(num_threads);
  DOPT(lattice_scale);
  DOPT(lattice_size);
  DOPT(tend);
  DOPT(out_intervall);
  DOPT(save_hdf);
  DOPT(safe_images);
  DOPT(fn_out);
  DOPT(fn_vessel);
  DOPT(paramset_name);
  DOPT(mtboxsize);
  #undef DOPT
}

ptree BulkTissueWithoutVessels::SimulationParameters::as_ptree() const
{
  boost::property_tree::ptree pt;
  #define DOPT(name) pt.put(#name, name)
  DOPT(num_threads);
  DOPT(lattice_scale);
  DOPT(lattice_size);
  DOPT(tend);
  DOPT(out_intervall);
  DOPT(save_hdf);
  DOPT(safe_images);
  DOPT(fn_out);
  DOPT(fn_vessel);
  DOPT(paramset_name);
  DOPT(mtboxsize);
  #undef DOPT
  return pt;
}

void BulkTissueWithoutVessels::SimulationParameters::update_ptree(ptree& dst, const ptree& src)
{
  boost::property_tree::update(dst, src);
}

// int main(int argc, char **argv)
// {
//   feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
//   my::MultiprocessingInitializer mpinit(argc, argv);
//   // default parameters
//   SimulationParameters sparams;
//   NewBulkTissueModel::Params pparams;
//   BloodFlowParameters bfparams;
//   VesselModel1::Params vessel_params;
//   O2Model::PrezO2Params prezO2params;
//   Adaption::Parameters adaption_params;
//   
//   
//   ptree all_pt_params;
//   all_pt_params = sparams.as_ptree();
//   all_pt_params.put_child("tumor", pparams.as_ptree());
//   all_pt_params.put_child("calcflow", bfparams.as_ptree());
//   all_pt_params.put_child("vessels", vessel_params.as_ptree());
//   all_pt_params.put_child("prez_o2", prezO2params.as_ptree());
//   all_pt_params.put_child("adaption", adaption_params.as_ptree());
// cout.rdbuf(my::log().rdbuf());
//   {
//   boost::optional<ptree> read_params = HandleSimulationProgramArguments(all_pt_params, argc, argv);
//   if (!read_params) 
//     return 0;
//   SimulationParameters::update_ptree(all_pt_params,*read_params);
//   
//   all_pt_params.put<Int3>("lattice_size", Int3(200,1,1));
//   sparams.assign(all_pt_params);
//   //boost::property_tree::update(params, BulkTissueWithoutVessels::Params().as_ptree());
//   
//   all_pt_params.put_child("tumor", NewBulkTissueModel::Params().as_ptree());
//   { 
// #ifdef DEBUG
//     cout << "read params in main are: ";
//     boost::property_tree::write_info(cout, all_pt_params);
//     cout << endl;
// #endif
//   }
// 
//   }//end cout.buffer
//   
//   my::SetNumThreads(sparams.num_threads);
//   BulkTissueWithoutVessels::run(all_pt_params);
// }

