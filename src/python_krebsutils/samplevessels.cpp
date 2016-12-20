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

/**@brief This file contains mostly functions to obtains samples of various quantities 
   distributed over blood vessel networks, e.g. blood pressure.
**/

#include "python-helpers.h"
#include "numpy.hpp"
#include "shared-objects.h"
#include "lattice-data-polymorphic.h"
#include "continuum-utils.h"
#include "continuum-grid.h"
#include <vesselgen/remodeler.h>
#include <vesselgen/vessgen_support.h>
#include <boost/range/numeric.hpp>


enum Mode {
  DATA_PER_NODE = 1, //@brief Associated with nodes. Should be linearly interpolated.
  DATA_CONST = 2,  //@brief No interpolation. Should not be used together with DATA_PER_NODE.
  DATA_LINEAR = 4, //@brief Interpolate. If not DATA_PER_NODE is not also set, then we expect a 2xn sized data array with values for both ends of each vessel.
};

/**@brief Generate sample data associated with vessels.

pos - world position of nodes
edges - n x 2 node indices
data - data associated with edges or nodes, depending on mode
sample_len - average distance between samples. The samples are taken in regular intervals if this value is smaller than the length of vessels. Otherwise the sampling becomes random.
mode - see Mode enum
*/
template<class T>
np::arraytbase sample_edges(np::arrayt<float> pos, np::arrayt<int> edges, np::arrayt<T> data, float sample_len, int mode)

{
  int cnt = edges.shape()[0];

  int ncomps = data.rank() > 1 ? pos.shape()[1] : 1;
  int itemtype = data.itemtype();
  T c[16][2];

  DynArray<T> tmp(1024, ConsTags::RESERVE);

  int num_total_samples = 0;
  CylinderNetworkSampler sampler;
  sampler.Init(sample_len, ptree());

  for(int i=0; i<cnt; ++i)
  {
    Float3 p0, p1;
    for (int j=0; j<3; ++j)
    {
      p0[j] = pos(edges(i,0), j);
      p1[j] = pos(edges(i,1), j);
    }

    sampler.Set(p0, p1, 0.); // 3rd arg is the radius
    int num_samples = sampler.GenerateLineSamples();

    for (int k=0; k<ncomps; ++k)
    {
      if (mode&DATA_PER_NODE)
      {
        myAssert(mode&DATA_LINEAR);
        c[k][0] = data(edges(i,0), k);
        c[k][1] = data(edges(i,1), k);
      }
      else if (mode&DATA_LINEAR)
      {
        c[k][0] = data(i, k, 0);
        c[k][1] = data(i, k, 1);
      }
      else
      {
        c[k][0] = data(i, k);
      }
    }
    for(int j=0; j<num_samples; ++j)
    {
      CylinderNetworkSampler::Sample ls = sampler.GetSample(j);
      if(mode&DATA_LINEAR)
      {
        float f = ls.fpos[2];
        for(int k=0; k<ncomps; ++k)
        {
          tmp.push_back((1.f-f)*c[k][0] + f*c[k][1]);
        }
      }
      else if(mode&DATA_CONST)
      {
        for(int k=0; k<ncomps; ++k) {
          tmp.push_back(c[k][0]);
        }
      }
    }
    num_total_samples += num_samples;
  }

  np::arrayt<T> acc_res(np::empty(2, Int2(num_total_samples, ncomps).cast<np::ssize_t>().eval().data(), itemtype));
  for(int i=0, k=0; i<num_total_samples; ++i)
  {
    for(int j=0; j<ncomps; ++j,++k)
    {
      acc_res(i,j) = tmp[k];
    }
  }
  return acc_res;
}


py::object sample_edges_weights(nm::array pypos, nm::array pyedges, float sample_len)
{
  np::arrayt<float> pos(pypos);
  np::arrayt<int> edges(pyedges);
  int cnt = edges.shape()[0];

  DynArray<float> tmp(1024, ConsTags::RESERVE);

  CylinderNetworkSampler sampler;
  sampler.Init(sample_len, ptree());

  for(int i=0; i<cnt; ++i)
  {
    Float3 p0, p1;
    for (int j=0; j<3; ++j)
    {
      p0[j] = pos(edges(i,0), j);
      p1[j] = pos(edges(i,1), j);
    }

    sampler.Set(p0, p1, 0.); // 3rd arg is the radius
    int num_samples = sampler.GenerateLineSamples();

    for (int j=0; j<num_samples; ++j)
    {
      tmp.push_back(sampler.weight);
    }
  }
  //np::arrayt<float> acc_res(np::empty(1, &tmp.size(), np::getItemtype<float>()));
  int strides[] = { 1 };
  int dims[] = { (int)tmp.size() };
  return np::copy<float,1>(dims, get_ptr(tmp), strides);
}



template<class T>
np::arraytbase sample_field(const nm::array py_pos, const np::arrayt<T> field, const py::object &py_ld, bool linear_interpolation, bool use_extrapolation_value, const T extrapolation_value)
{
  /*
   * sample values from a grid, where the sample locations are given by the py_pos n x 3 array
    */
  np::arrayt<float> pos(py_pos);

//   FieldInterpolator<T> interpolate;
//   if (use_extrapolation_value)
//     interpolate.init(CONT_CONST, make_ptree("value", extrapolation_value));
//   else
//     interpolate.init(CONT_EXTRAPOLATE);

  np::ssize_t num_samples = pos.shape()[0];
  
  np::arrayt<T> res = np::zeros(1, &num_samples, np::getItemtype<T>());

  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ld);

  Array3d<T> arr3d = Array3dFromPy<T>(const_cast<np::arrayt<T>&>(field));
  arr3d.move(ld.Box().min);

  auto interpolate = [&](const Float3 &p) -> T
  {
    if (use_extrapolation_value)
    {
      if (linear_interpolation)
        return FieldInterpolate::ValueAveraged(arr3d, ld, FieldInterpolate::Const<T>(extrapolation_value), p);
      else
        return FieldInterpolate::Value(arr3d, ld, FieldInterpolate::Const<T>(extrapolation_value), ld.WorldToLattice(p));
    }
    else
    {
      if (linear_interpolation)
        return FieldInterpolate::ValueAveraged(arr3d, ld, FieldInterpolate::Extrapolate(), p);
      else
        return FieldInterpolate::Value(arr3d, ld, FieldInterpolate::Extrapolate(), ld.WorldToLattice(p));
    }
  };
  
  for (int i=0; i<num_samples; ++i)
  {
    Float3 p(pos(i,0),pos(i,1),pos(i,2));
    T r = interpolate(p);
    res(i) = r;
  }

  return res;
}


np::arraytbase make_position_field(const py::object &py_ldobj)
{
  /*
   * this fills a 4d array with world coordinates of its grid points (or cell centers, depending on the configuration of the lattice data)
   */
  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ldobj);
  const auto box = ld.Box();
  Int3 size = Size(box);
  np::ssize_t ndims[4] = { size[0], size[1], size[2], 3 };
  np::arrayt<float> res = np::zeros(4, ndims, np::getItemtype<float>());
  FOR_BBOX3(p, box)
  {
    Float3 w = ld.LatticeToWorld(p);
    for (int i=0; i<3; ++i)
      res(p[0]-box.min[0], p[1]-box.min[1], p[2]-box.min[2], i) = w[i];
  }
  return res;
}


/**@brief Approximately fill a grid with how much of each voxel is occupied by vessels.
*/
np::arraytbase compute_vessel_volume_fraction_field(np::arrayt<float> pos, np::arrayt<int> edges, np::arrayt<float> radius, const py::object &py_ldfield, int samples_per_cell)
{
  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ldfield);

  int cnt = edges.shape()[0];

  Array3d<float> tmp(ld.Box());

  int num_total_samples = 0;
  CylinderNetworkSampler sampler;
  sampler.Init(ld.Scale(), make_ptree("samples_per_cell", samples_per_cell));

  for(int i=0; i<cnt; ++i)
  {
    Float3 p0, p1;
    for (int j=0; j<3; ++j)
    {
      p0[j] = pos(edges(i,0), j);
      p1[j] = pos(edges(i,1), j);
    }
    sampler.Set(p0, p1, radius(i)); // 3rd arg is the radius
    int num_samples = sampler.GenerateVolumeSamples();
    for (int k=0; k<num_samples; ++k)
    {
      AddSmoothDelta(tmp, ld.Box(), ld, 3, sampler.GetSample(k).wpos, sampler.weight_per_volume);
    }
    //cout << sampler.GetSample(0).wpos << endl;
  }
  
  np::arrayt<float> res = np::zeros(3, Cast<np::ssize_t>(::Size(ld.Box())).data(), np::getItemtype<float>());
  FOR_BBOX3(p, ld.Box())
  {
    res(p[0], p[1], p[2]) = std::min<float>(1., tmp(p));
  }

  return res;
}


/**@brief Used to determine the Fractal Dimension of the network. (Or more correctly the box counting dimension.)
*/
double compute_vessel_boxcounts(nm::array pypos, nm::array pyedges, nm::array pyradius, const py::object &py_ldfield, double volume_scaling, double volume_threshold)
{
  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ldfield);

  np::arrayt<float> pos(pypos);
  np::arrayt<int> edges(pyedges);
  np::arrayt<float> radius(pyradius);
  int cnt = edges.shape()[0];

  // first sort vessels into boxes
  DynArray<BBox3> boxes = MakeMtBoxGridLarge(ld.Box(), 160);
  DynArray<DynArray<int> > vesselgrid(boxes.size());
  
  DynArray<FloatBBox3> wboxes(boxes.size());
  for (int i=0; i<boxes.size(); ++i)
  {
    BBox3 bb = boxes[i];
    wboxes[i] = FloatBBox3(ld.LatticeToWorld(bb.min) - Float3(0.5*ld.Scale()),
                           ld.LatticeToWorld(bb.max) + Float3(0.5*ld.Scale()));
  }

  for (int vi=0; vi<cnt; ++vi)
  {
    Float3 wp[2];
    for (int dim=0; dim<3; ++dim)
      for (int i=0; i<2; ++i)
        wp[i][dim] = pos(edges(vi,i), dim);
    
    FloatBBox3 vwbb = FloatBBox3().Add(wp[0]).Add(wp[1]);
    vwbb.Extend(radius(vi));
    for (int i=0; i<boxes.size(); ++i)
    {
      if (!(vwbb.Overlaps(wboxes[i]))) continue;
      vesselgrid[i].push_back(vi);
    }
  }

  int64 boxcount = 0;
  #pragma omp parallel
  {
    CylinderNetworkSampler sampler;
    sampler.Init(ld.Scale(), make_ptree("samples_per_cell", 10));
    Array3d<float> buffer;
    my::Time t_;

    #pragma omp for nowait schedule(dynamic, 1) reduction(+:boxcount)
    for (int ibox=0; ibox<boxes.size(); ++ibox)
    {
      t_ = my::Time();
      int box_num_samples = 0;
      buffer.initFromBox(boxes[ibox]);
      for(int vi=0; vi<vesselgrid[ibox].size(); ++vi)
      {
        int i = vesselgrid[ibox][vi];
        Float3 p0, p1;
        for (int j=0; j<3; ++j)
        {
          p0[j] = pos(edges(i,0), j);
          p1[j] = pos(edges(i,1), j);
        }
        sampler.Set(p0, p1, radius(i));
        int num_samples = sampler.GenerateVolumeSamples();
        for (int k=0; k<num_samples; ++k)
        {
          AddSmoothDelta(buffer, boxes[ibox], ld, 3, sampler.GetSample(k).wpos, sampler.weight_per_volume);
        }
        box_num_samples += num_samples;
      }

//       #pragma omp critical
//       {
//         cout << format("box %i/%i, samples = %i, time = %f ms") % ibox % boxes.size() % box_num_samples % (my::Time()-t_).to_ms() << endl;
//         Image imgdbg; imgdbg.Init(::Size(boxes[ibox])[0], ::Size(boxes[ibox])[1]);
//         FOR_BBOX3(p, boxes[ibox])
//           if (p[2] == (boxes[ibox].max[2]+boxes[ibox].min[2])/2)
//             if (buffer(p)*volume_scaling > volume_threshold)
//               imgdbg.SetPixel(p[0]-boxes[ibox].min[0], p[1]-boxes[ibox].min[1], 255,255,255);
//         imgdbg.Write(str(format("boxcount_h%f_part%04i.png") % ld.Scale() % ibox));
//       }
      
      FOR_BBOX3(p, boxes[ibox])
      {
        if (buffer(p)*volume_scaling > volume_threshold) {
          ++boxcount;
        }
      }
    }
  }
  return boxcount;
}

py::object calculate_within_fake_tumor_lattice_based(const py::str &property_name, py::object &ld_grp_obj, const py::object &vess_grp_obj, const py::object &tumor_range, const bool av)
{
  std::string property = py::extract<std::string>(property_name);
  std::vector<std::string> mylist{"lengths","radii", "pressures", "flows"};
  if (std::find(std::begin(mylist), std::end(mylist), property) != std::end(mylist))
  {
    printf("Property %s is good!" , property.c_str());
  }
  else
  {
    printf("Property %s is not supported!" , property.c_str());
  }
  FpExceptionStateGuard exception_state_guard(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  h5cpp::Group g_vess = PythonToCppGroup(vess_grp_obj);
  float theTumorRange = py::extract<float>(tumor_range);

  std::auto_ptr<VesselList3d> vl = ReadVesselList3d(g_vess, make_ptree("filter", true));
  
  //discriminate unwanted vessels
  //std::vector<bool> in_tumor(vl->GetECount(),ConsTags::RESERVE);
  //std::vector<bool> outside_tumor(vl->GetECount(),ConsTags::RESERVE);
  DynArray<bool> in_tumor(vl->GetECount(),ConsTags::RESERVE);
  DynArray<bool> outside_tumor(vl->GetECount(),ConsTags::RESERVE);
  for(int i=0;i<vl->GetECount();++i)
  {
    Vessel* v = vl->GetEdge(i);
    float normA = vl->Ld().LatticeToWorld(v->GetNode(0)->lpos).transpose().norm();
    float normB = vl->Ld().LatticeToWorld(v->GetNode(1)->lpos).transpose().norm();
    if(normA<theTumorRange and normB<theTumorRange)
    {
      in_tumor[i] = true;
      outside_tumor[i] = false;
    }
    else
    {
      in_tumor[i] = false;
      outside_tumor[i] = true;
    }
  }
  std::printf("%i vessels in the tumor.\n", boost::accumulate(in_tumor, 0));
  std::printf("%i vessels outside the tumor.\n", boost::accumulate(outside_tumor, 0));
  
  py::list ret;//this will be the return value
  
  if(property == "lengths")
  {
    std::vector<BranchDat> lengths_in = MeasureBranchLengths(*vl, in_tumor);
    std::vector<BranchDat> lengths_out = MeasureBranchLengths(*vl, outside_tumor);
  
    np::ssize_t dims_in[2] = {3,(int)lengths_in.size()};
    np::ssize_t dims_out[2] = {3,(int)lengths_out.size()};
    // create numpy array
    np::arrayt<float> buffer_in = np::zeros(2, dims_in, np::getItemtype<float>());
    for (int i=0; i<lengths_in.size(); ++i)
    {
      buffer_in(0, i) = lengths_in[i].r;
      buffer_in(1, i) = lengths_in[i].l;
      buffer_in(2, i) = lengths_in[i].n;
    }
    np::arrayt<float>buffer_out = np::zeros(2, dims_out, np::getItemtype<float>());
    for (int i=0; i<lengths_out.size(); ++i)
    {
      buffer_out(0, i) = lengths_out[i].r;
      buffer_out(1, i) = lengths_out[i].l;
      buffer_out(2, i) = lengths_out[i].n;
    }
    ret.append(buffer_in);
    ret.append(buffer_out);
  }
  if(property == "radii")
  {
    std::vector<float> radii_in;
    std::vector<float> radii_out;
    for( int i=0;i<vl->GetECount();++i)
    {
      Vessel* v=vl->GetEdge(i);
      if(in_tumor[i])
      {
	radii_in.push_back(v->r);
      }
      if(outside_tumor[i])
      {
	radii_out.push_back(v->r);
      }
    }
    np::ssize_t dims_in[2] = {1,(int)radii_in.size()};
    np::ssize_t dims_out[2] = {1,(int)radii_out.size()};
    // create numpy array
    np::arrayt<float> buffer_in = np::zeros(2, dims_in, np::getItemtype<float>());
    for (int i=0; i<radii_in.size(); ++i)
    {
      buffer_in(0, i) = radii_in[i];
    }
    np::arrayt<float>buffer_out = np::zeros(2, dims_out, np::getItemtype<float>());
    for (int i=0; i<radii_out.size(); ++i)
    {
      buffer_out(0, i) = radii_out[i];
    }
    ret.append(buffer_in);
    ret.append(buffer_out);
  }
  if(property == "flows")
  {
    std::vector<float> radii_in;
    std::vector<float> radii_out;
    for( int i=0;i<vl->GetECount();++i)
    {
      Vessel* v=vl->GetEdge(i);
      if(in_tumor[i])
      {
	radii_in.push_back(v->q);
      }
      if(outside_tumor[i])
      {
	radii_out.push_back(v->q);
      }
    }
    np::ssize_t dims_in[2] = {1,(int)radii_in.size()};
    np::ssize_t dims_out[2] = {1,(int)radii_out.size()};
    // create numpy array
    np::arrayt<float>buffer_in = np::zeros(2, dims_in, np::getItemtype<float>());
    for (int i=0; i<radii_in.size(); ++i)
    {
      buffer_in(0, i) = radii_in[i];
    }
    np::arrayt<float>buffer_out = np::zeros(2, dims_out, np::getItemtype<float>());
    for (int i=0; i<radii_out.size(); ++i)
    {
      buffer_out(0, i) = radii_out[i];
    }
    ret.append(buffer_in);
    ret.append(buffer_out);
  }
  if(property == "pressures")
  {
    std::vector<float> pressure_in_a;
    std::vector<float> pressure_in_v;
    std::vector<float> pressure_out_a;
    std::vector<float> pressure_out_v;
    for( int i=0;i<vl->GetECount();++i)
    {
      Vessel* v=vl->GetEdge(i);
      if(in_tumor[i])
      {
	if(v->IsArtery())
	  pressure_in_a.push_back(0.5*(v->NodeA()->press+v->NodeB()->press));
	
	if(v->IsVein())
	  pressure_in_v.push_back(0.5*(v->NodeA()->press+v->NodeB()->press));
      }
      if(outside_tumor[i])
      {
	if(v->IsArtery())
	  pressure_out_a.push_back(0.5*(v->NodeA()->press+v->NodeB()->press));
	  pressure_out_v.push_back(0.5*(v->NodeA()->press+v->NodeB()->press));
      }
    }
    if(not av)
    {
      pressure_in_a.insert(pressure_in_a.end(),pressure_in_v.begin(),pressure_in_v.end());
      pressure_out_a.insert(pressure_out_a.end(),pressure_out_v.begin(),pressure_out_v.end());
      np::ssize_t dims_in[2] = {1,(int)pressure_in_a.size()};
      np::ssize_t dims_out[2] = {1,(int)pressure_out_a.size()};
      // create numpy array
      np::arrayt<float>buffer_in = np::zeros(2, dims_in, np::getItemtype<float>());
      for (int i=0; i<pressure_in_a.size(); ++i)
      {
	buffer_in(0, i) = pressure_in_a[i];
      }
      np::arrayt<float>buffer_out = np::zeros(2, dims_out, np::getItemtype<float>());
      for (int i=0; i<pressure_out_a.size(); ++i)
      {
	buffer_out(0, i) = pressure_out_a[i];
      }
      ret.append(buffer_in);
      ret.append(buffer_out);
    }
    if(av)
    {
      np::ssize_t dims_in_a[2] = {1,(int)pressure_in_a.size()};
      np::ssize_t dims_in_v[2] = {1,(int)pressure_in_v.size()};
      np::ssize_t dims_out_a[2] = {1,(int)pressure_out_a.size()};
      np::ssize_t dims_out_v[2] = {1,(int)pressure_out_v.size()};
      // create numpy array
      // this is quick and dirty!!!
      np::arrayt<float>buffer_in_a = np::zeros(2, dims_in_a, np::getItemtype<float>());
      for (int i=0; i<pressure_in_a.size(); ++i)
      {
	buffer_in_a(0, i) = pressure_in_a[i];
      }
      np::arrayt<float>buffer_in_v = np::zeros(2, dims_in_v, np::getItemtype<float>());
      for (int i=0; i<pressure_in_v.size(); ++i)
      {
	buffer_in_v(0, i) = pressure_in_v[i];
      }
      np::arrayt<float>buffer_out_a = np::zeros(2, dims_out_a, np::getItemtype<float>());
      for (int i=0; i<pressure_out_a.size(); ++i)
      {
	buffer_out_a(0, i) = pressure_out_a[i];
      }
      np::arrayt<float>buffer_out_v = np::zeros(2, dims_out_v, np::getItemtype<float>());
      for (int i=0; i<pressure_out_v.size(); ++i)
      {
	buffer_out_v(0, i) = pressure_out_v[i];
      }
      
      ret.append(buffer_in_a);
      ret.append(buffer_in_v);
      ret.append(buffer_out_a);
      ret.append(buffer_out_v);
    }
  }
  return ret;
}
/* 
 * returns an 3xn array with BranchDat to python
 * there this could be use to easily create histograms
 * TODO make this function also work on world coordinates?
 */
/*
py::object calculate_lengths_lattice_based(const py::object &ld_grp_obj, const py::object &vess_grp_obj, const py::object &tumor_range)
{ 
  FpExceptionStateGuard exception_state_guard(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  my::h5::Group g_ld = PythonToCppGroup(ld_grp_obj);
  my::h5::Group g_vess = PythonToCppGroup(vess_grp_obj);
  float theTumorRange = py::extract<float>(tumor_range);

  std::auto_ptr<VesselList3d> vl = ReadVesselList3d(g_vess, g_ld, make_ptree("filter", true));
  
  //discriminate unwanted vessels
  std::vector<bool> in_tumor(vl->GetECount(),RESERVE);
  std::vector<bool> outside_tumor(vl->GetECount(),RESERVE);
  for(int i=0;i<vl->GetECount();++i)
  {
    Vessel* v = vl->GetEdge(i);
    float normA = vl->Ld().LatticeToWorld(v->GetNode(0)->lpos).transpose().norm();
    float normB = vl->Ld().LatticeToWorld(v->GetNode(1)->lpos).transpose().norm();
    if(normA<theTumorRange and normB<theTumorRange)
    {
      in_tumor[i] = true;
      outside_tumor[i] = false;
    }
    else
    {
      in_tumor[i] = false;
      outside_tumor[i] = true;
    }
  }
  std::printf("%i vessels in the tumor.\n", boost::accumulate(in_tumor, 0));
  std::printf("%i vessels outside the tumor.\n", boost::accumulate(outside_tumor, 0));
  
  std::vector<BranchDat> lengths_in = MeasureBranchLengths(*vl, in_tumor);
  std::vector<BranchDat> lengths_out = MeasureBranchLengths(*vl, outside_tumor);
  
  int dims_in[2] = {3,lengths_in.size()};
  int dims_out[2] = {3,lengths_out.size()};
  // create numpy array
  np::array py_branch_dat_in = np::zeros(2, dims_in, np::get_itemtype<float>());
  np::arrayt<float>buffer_in(py_branch_dat_in)  ;
  for (int i=0; i<lengths_in.size(); ++i)
  {
    buffer_in(0, i) = lengths_in[i].r;
    buffer_in(1, i) = lengths_in[i].l;
    buffer_in(2, i) = lengths_in[i].n;
  }
  np::array py_branch_dat_out = np::zeros(2, dims_out, np::get_itemtype<float>());
  np::arrayt<float>buffer_out(py_branch_dat_out)  ;
  for (int i=0; i<lengths_out.size(); ++i)
  {
    buffer_out(0, i) = lengths_out[i].r;
    buffer_out(1, i) = lengths_out[i].l;
    buffer_out(2, i) = lengths_out[i].n;
  }

  py::list ret;
  ret.append(py_branch_dat_in);
  ret.append(py_branch_dat_out);
  return ret;
  
  
//   py::list abc(lengths);
//   return abc;
  
}

py::object get_radii_within_fake_tumor(const py::object &ld_grp_obj, const py::object &vess_grp_obj, const py::object &tumor_range)
{ 
  FpExceptionStateGuard exception_state_guard(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  my::h5::Group g_ld = PythonToCppGroup(ld_grp_obj);
  my::h5::Group g_vess = PythonToCppGroup(vess_grp_obj);
  float theTumorRange = py::extract<float>(tumor_range);

  std::auto_ptr<VesselList3d> vl = ReadVesselList3d(g_vess, g_ld, make_ptree("filter", true));
  
  //discriminate unwanted vessels
  std::vector<bool> in_tumor(vl->GetECount(),RESERVE);
  std::vector<bool> outside_tumor(vl->GetECount(),RESERVE);
  for(int i=0;i<vl->GetECount();++i)
  {
    Vessel* v = vl->GetEdge(i);
    float normA = vl->Ld().LatticeToWorld(v->GetNode(0)->lpos).transpose().norm();
    float normB = vl->Ld().LatticeToWorld(v->GetNode(1)->lpos).transpose().norm();
    if(normA<theTumorRange and normB<theTumorRange)
    {
      in_tumor[i] = true;
      outside_tumor[i] = false;
    }
    else
    {
      in_tumor[i] = false;
      outside_tumor[i] = true;
    }
  }
  std::printf("%i vessels in the tumor.\n", boost::accumulate(in_tumor, 0));
  std::printf("%i vessels outside the tumor.\n", boost::accumulate(outside_tumor, 0));
  
  std::vector<float> radii_in;
  std::vector<float> radii_out;
  
  for( int i=0;i<vl->GetECount();++i)
  {
    Vessel* v=vl->GetEdge(i);
    if(in_tumor[i])
    {
      radii_in.push_back(v->r);
    }
    if(outside_tumor[i])
    {
      radii_out.push_back(v->r);
    }
  }
  
  int dims_in[2] = {1,radii_in.size()};
  int dims_out[2] = {1,radii_out.size()};
  // create numpy array
  np::array py_radii_in = np::zeros(2, dims_in, np::get_itemtype<float>());
  np::arrayt<float>buffer_in(py_radii_in)  ;
  for (int i=0; i<radii_in.size(); ++i)
  {
    buffer_in(0, i) = radii_in[i];
  }
  np::array py_radii_out = np::zeros(2, dims_out, np::get_itemtype<float>());
  np::arrayt<float>buffer_out(py_radii_out)  ;
  for (int i=0; i<radii_out.size(); ++i)
  {
    buffer_out(0, i) = radii_out[i];
  }

  py::list ret;
  ret.append(py_radii_in);
  ret.append(py_radii_out);
  return ret;
  
  
//   py::list abc(lengths);
//   return abc;
  
}
*/
void export_samplevessels()
{
//   py::def("export_network_for_povray", export_network_for_povray);
  py::enum_<Mode>("VesselSamplingFlags")
    .value("DATA_PER_NODE", DATA_PER_NODE)
    .value("DATA_CONST", DATA_CONST)
    .value("DATA_LINEAR", DATA_LINEAR);
  #define DEFINE_sample_edges_t(T)\
      py::def("sample_edges_"#T, sample_edges<T>);
  DEFINE_sample_edges_t(float)
  DEFINE_sample_edges_t(double)
  DEFINE_sample_edges_t(int)
  DEFINE_sample_edges_t(uint)
  DEFINE_sample_edges_t(char)
  DEFINE_sample_edges_t(uchar)

#define DEFINE_sample_field_t(T)\
  py::def("sample_field_"#T, sample_field<T>);
  DEFINE_sample_field_t(float)
  DEFINE_sample_field_t(double)
  py::def("sample_edges_weights", sample_edges_weights);
  py::def("make_position_field", make_position_field);
  py::def("make_vessel_volume_fraction_field", compute_vessel_volume_fraction_field);
  py::def("calc_vessel_boxcounts", compute_vessel_boxcounts);
//  py::def("calculate_lengths_lattice_based", calculate_lengths_lattice_based);
//  py::def("get_radii_within_fake_tumor", get_radii_within_fake_tumor);
  py::def("calculate_within_fake_tumor_lattice_based", calculate_within_fake_tumor_lattice_based);
}


