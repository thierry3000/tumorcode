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
#include <fstream>

#include "remodeler.h"

#include "../calcflow.h"
#include "../hdfio.h"
#include <deque>


//void WriteHdfHistogram( h5::Group f, const string &id, const BasicHistogram1D<float> &h );

/* forward declaration, defined in vessgen_support.cpp */
void DoOutput(H5::H5File &file,
              const VesselList3d &vl,
              const Grower &grower,
              const boost::property_tree::atree &additional_data,
              const ptree &input_pt
             );
H5::Group DebugOutVessels(const Grower &grower, const string &name);


struct QualityAnalysis
{
  std::deque<float> avg_window;
  float maxQuality;
  float windowAverage;
  float lastWindowAverage;
  float changeRate;
  int maxWindowSize;
  float changeRateThreshold; // stop if quality change rate / windowAverageQuality is less than this, (and other conditions are satisfied)

  void init(int maxWindowSize_, float changeRateThreshold_)
  {
    maxWindowSize = maxWindowSize_;
    maxQuality = 0;
    windowAverage = 0;
    lastWindowAverage = 0;
    changeRate = std::numeric_limits<float>::max();
    avg_window.clear();
    changeRateThreshold = changeRateThreshold_;
  }
  
  void insert(float val)
  {
    avg_window.push_back(val);
    if (avg_window.size() > maxWindowSize)
      avg_window.pop_front();
    maxQuality = std::max(maxQuality, val);
    lastWindowAverage = windowAverage;
    windowAverage = 0.;
    for (int i=0; i<avg_window.size(); ++i)
      windowAverage += avg_window[i];
    windowAverage /= avg_window.size();
    changeRate = windowAverage - lastWindowAverage;
  }

  bool testStop() const
  {
     if (avg_window.size() < maxWindowSize) return false;
     return (changeRate / windowAverage) < changeRateThreshold &&
            (maxQuality - windowAverage) / maxQuality < 0.5; // window average is close to max measured quality?
  }

  void print(std::ostream &os) const
  {
    os << format("||Q: avg: %f, rate: %f, best: %f, stop=%i") % windowAverage % (changeRate/windowAverage) % ((maxQuality-windowAverage)/maxQuality) % testStop();
  }
};




struct VessGenApp
{
  VessGenApp();
  void run(const ptree& input_pt_);


  int quality_buff_size;
  int measure_interval;
  string outfilename;
  my::Time real_start_time;
  

  ptree input_pt;

  boost::property_tree::atree iter_data, last_iter_data;

  QualityAnalysis qa;
  float last_saved_quality;
  int finished_countdown;

  bool Callback(const Grower& grower);
  void RecordData(const Grower &grower);
};



VessGenApp::VessGenApp()
{
  // protect from ever beeing deleted
#ifdef DEBUG  
  quality_buff_size = 10;
  measure_interval = 10;
#else
  quality_buff_size = 30;
  measure_interval = 50;
#endif
  outfilename  = string();
  last_saved_quality = 0;
  finished_countdown = 0;
}



void VessGenApp::RecordData(const Grower& grower)
{
  int sites;
  float bloodVolume, frac_sites;

  const VesselList3d &vl = grower.get_vl();
  const VesselList3d::LatticeData &ld = vl.Ld();
  int dummy1;
  bloodVolume = MeasureBloodVolume(vl, 10.0f, dummy1, sites);
  frac_sites = (float(sites)/Volume(ld.Box()));
  iter_data.put("num_iter", grower.iteration_number);
  boost::property_tree::atree p, &iter_array = boost::property_tree::require_child(iter_data, "iters");
  p.put("iter", grower.iteration_number);
  p.put("rBV", bloodVolume);
  p.put("sites", sites);
  p.put("frac_sites", frac_sites);
  last_iter_data = p;
  iter_array.push_back(std::make_pair("", p));
  qa.insert(frac_sites);
}


bool VessGenApp::Callback(const Grower& grower) // returns if the iteration should continue
{
  if (grower.iteration_number_on_level == 0) 
    qa.init(quality_buff_size, input_pt.get<float>("changeRateThreshold", 0.1e-3));
  
  RecordData(grower);
        
  cout << format("hit: %i/%i, it: %i, rBV: %f, sites: %f ") % grower.iteration_number_on_level % grower.hierarchy_level % grower.iteration_number % last_iter_data.get<float>("rBV") % last_iter_data.get<float>("frac_sites");
  qa.print(cout); 
  cout << endl;

  if (qa.testStop() || (grower.iteration_number_on_level >= grower.max_num_iter - (quality_buff_size)))
    ++finished_countdown;
  else
    finished_countdown = 0;
  

  if (finished_countdown > 0 && last_iter_data.get<float>("rBV")>last_saved_quality && grower.hierarchy_level>=grower.max_hierarchy_level)
  {
    iter_data.put<my::Time>("real_start_time", real_start_time);
    try{
	// new way  H5F_ACC_RDWR
       //H5::H5File *readInFile = new H5::H5File(fn, H5F_ACC_RDONLY );
      //H5::H5File file(outfilename+".h5", H5F_ACC_RDWR);H5F_ACC_EXCL
      
      H5::H5File file = H5::H5File(outfilename + ".h5", H5F_ACC_TRUNC);
      DoOutput(file, grower.get_vl(), grower, iter_data, input_pt);
      file.close();
    }
    catch(H5::Exception &e)
    {
      e.printErrorStack();
    }
    //CalcFlow(grower.get_vl(),);
    
  }
#ifndef NDEBUG
  std::cout << "before return in callback" << std::endl;
#endif
  
  return finished_countdown < quality_buff_size;
}



void VessGenApp::run(const ptree &input_pt_)
{
    input_pt = input_pt_;
    outfilename = input_pt.get<string>("out_fn");
    quality_buff_size = input_pt.get<int>("max_num_iter") / 10;
    /* TF:
     * I think vessel gen should be single threaded and
     * I will get rid of this hand make thread management to use the system queueing system
     */
    //my::SetNumThreads(input_pt.get<int>("num_threads", 1));
    //my::SetNumThreads(1);
    
    Grower grower;
    grower.Run(input_pt, boost::bind(&VessGenApp::Callback, this, _1));
    std::cout << "grower.Run finished!" << std::endl;
}


namespace VesselGenerator
{

#define DOPT(name, default) pt.put(#name, default)
#define DOPT2(name, default) calcflow.put(#name, default)
ptree getDefaultParameters()
{
  ptree pt;
  ptree calcflow;
  DOPT2(viscosityPlasma, 1.2e-06);
  DOPT2(rheology, "RheologySecomb2005");
  DOPT2(inletHematocrit, 0.37);
  DOPT2(includePhaseSeparationEffect, true);
  pt.add_child("calcflow", calcflow);
  pt.put("out_fn", string());
  DOPT(max_num_iter, 300);
  DOPT(num_threads, 1);
  pt.put("lattice_size", Int3(100, 100, 100));
  DOPT(lattice_spacing, 150.);
  DOPT(seed, 0);
  DOPT(lattice_type, "FCC");
  DOPT(o2.diffusion_range, 90.);
  DOPT(max_sprout_radius_artery, 5.);
  DOPT(max_sprout_radius_vein, 15.); // was 6 for both
  DOPT(radius_vein, 5.);
  DOPT(radius_capi, 4.);
  DOPT(radius_artery, 4.);
  DOPT(murray_alpha, 3.); // 3  might be better, also 2.55, was 2.4
  DOPT(full_debug_output, false);
  DOPT(num_hierarchical_iterations, 0);
  DOPT(optimization_mode, "v1");
  DOPT(samples_per_cell, 3);
  DOPT(full_debug_output, false);
  DOPT(message, "aVesselName");
  DOPT(capillariesUntilLevel, 0);
  DOPT(changeRateThreshold, 0.001);
  DOPT(murray_alpha_vein, 3.0);
  DOPT(murray_alpha_artery, 3.0);
  DOPT(ensemble_index, 5);
  DOPT(generate_more_capillaries, 0);
  return pt;
}



ptree getVariableArgumentCountParameters()
{
  ptree pt;
  pt.put_child("roots", ptree());
  return pt;
}
#undef DOPT

  
void run(const ptree &parameters)
{
  VessGenApp app;
  app.run(parameters);
  std::cout << "VessGenApp.run finished " << std::endl;
}

}


