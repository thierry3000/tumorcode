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
#ifndef IFF_DRUG_APP_H
#define IFF_DRUG_APP_H

#include "python_krebsutils/python_helpers.h"

#include "calc-ifflow.h"
#include "calc-ifdrug.h"

//#ifdef USE_IFFLOW

#include "common/common.h"
#include "common/hdfio.h"

#include "mwlib/log.h"
#include "mwlib/drawing.h"

#include <fenv.h>
#include <fstream>
#include <boost/property_tree/info_parser.hpp>

#include "common/vessels3d.h"
#include "common/continuum-flow.h"
#include "common/time_stepper_utils.h"

namespace h5 = h5cpp;


class VirtualGridFuncCellsWoNecro : public VirtualGridFunctions<float, 1>
{
  ConstArray3d<float> cells, necro;
public:
  VirtualGridFuncCellsWoNecro(const ConstArray3d<float> &cells, const ConstArray3d<float> &necro) : cells(cells), necro(necro) {}
  virtual void GetValues(int boxindex, const BBox3 &bbox, List &values) const
  {
    FOR_BBOX3(p, bbox)
    {
      values[0](p) = cells(p)-necro(p);
    }
  }
};

struct TissueCompositionCallback
{
  Array3d<float> phi_cells, phi_water;
  void operator()(int boxindex, const BBox3 &bbox, Array3d<float> &phi_cells_, Array3d<float> &phi_water_, Array3d<float> &phi_ecm_)
  {
    FOR_BBOX3(p, bbox)
    {
      phi_cells_(p) = phi_cells(p);
      phi_water_(p) = phi_water(p);
      phi_ecm_(p) = 1. - phi_cells(p) - phi_water(p);
    }
  }
};




class IffDrugApp3d
{
  enum { dim = 3 };
  float out_intervall;
  DynArray<float> out_times;
  string fn_out,
         fn_tumor,
         message,
         h5_path_vessel,
         h5_path_lattice,
         h5_path_tumor,
         tumor_type_string,
	 parameterset_name;
  int output_number;



  //void InitFieldsAndLd(h5::Group ldgroup);
  //LatticeDataQuad3d field_ld;

  ptree all_pt_params;
  my::Time real_start_time;

  std::auto_ptr<VesselFluidExavasationModel> vessel_sources_f;
  std::auto_ptr<VirtualGridFunctions<double,2> > lymph_sources_f;
  std::auto_ptr<VirtualGridFunctions<double,2> > tissue_f;
  std::auto_ptr<VirtualGridFunctions<float,1> > phi_cells_f;
  std::auto_ptr<VirtualGridFunctions<float,1> > phi_cells_wo_necro_f;
  std::auto_ptr<VirtualGridFunctions<float,1> > phi_water_f;
  TissueCompositionCallback tissueCompositionCallback;

  IffParams iff_params;
  IfDrug::Params ift_params;
  ContinuumGrid grid;
  DynArray<BBox3> mtboxes;
  std::auto_ptr<VesselList3d> vl;
  VesselsInBoxes vesselsInBoxes;
  Array3d<float> theta_tumor, theta_necro, phi_cells, phi_water;

  std::auto_ptr<IfFlowState> iffstate;

  DynArray<double> org_press, org_flow;

  bool InitNewState();
  void InitAsSimpleModels(Array3d< float > theta_tumor, Array3d< float > theta_necro, Array3d< float > phi_water, Array3d< float > phi_cells);

  void DoIffCalc();
  void MeasureIfFlowState(h5::Group g);

#ifdef USE_IFDRUGSIM
  void WriteDrugOutput(double t, const IfDrug::Calculator::State &conc_field, const IfDrug::Calculator& model, const ptree& params);
#endif

public:
  IffDrugApp3d();
  int  Main(const ptree &read_params, const string &outfilename, py::object &drug_measurement_function);
};

#endif //#ifndef IFF_DRUG_APP_H
