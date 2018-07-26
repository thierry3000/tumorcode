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

#include "python_helpers.h"

//#include "numpy.hpp"

#include <fstream>

enum {
  CLIP_NONE = 0,
  CLIP_PARTIAL = 1,
  CLIP_FULL = 2
};

class ClipBase
{
public:
  virtual ~ClipBase() {}
  virtual int clipSphere(const Float3 &pos, float rad) { return CLIP_NONE; }
  virtual int clipCylinder(const Float3 &posa, const Float3 &posb, float rad) { return CLIP_NONE; }
  virtual string asPovRayObject(const string &stylestr) { return string(""); }
};

typedef boost::shared_ptr<ClipBase> CP;


class ClipPlane : public ClipBase
{
  Float3 n;
  float d;
public:
  ClipPlane(const Float3 &n_, float d_) : n(n_), d(d_) {}
  ClipPlane(const Float3 &n_, const Float3 &o_) : n(n_), d(n_.dot(o_)) {}
  virtual int clipSphere(const Float3 &pos, float r)
  {
    float q = pos.dot(n)-d;
    r *= 1.1;
    if (q>r) return CLIP_FULL;
    if (q<-r) return CLIP_NONE;
    return CLIP_PARTIAL;
  }
  virtual int clipCylinder(const Float3 &pos1, const Float3 &pos2, float r)
  {
    float q1 = pos1.dot(n)-d;
    float q2 = pos2.dot(n)-d;
    r *= 1.1;
    if (q1>r && q2>r) return CLIP_FULL;
    if (q1<-r && q2<-r) return CLIP_NONE;
    return CLIP_PARTIAL;
  }
  virtual string asPovRayObject(const string &stylestr)
  {
    const char* plane_template = "plane\n{ <%f,%f,%f>,%f\n%s }\n";
    return str(format(plane_template) % n[0] % n[1] % n[2] % d % stylestr);
  }
};


static const int clip_logical_and[3][3] = {
  { 0, 0, 0 },
  { 0, 1, 1 },
  { 0, 1, 2 }
};

static const int clip_logical_or[3][3] = {
  { 0, 1, 2 },
  { 1, 1, 2 },
  { 2, 2, 2 }
};


class ClipAnd : public ClipBase
{
  boost::shared_ptr<ClipBase> c0, c1;
public:
  ClipAnd(CP c0_, CP c1_) : c0(c0_), c1(c1_) {}
  virtual int clipSphere(const Float3 &pos, float rad)
  {
    return clip_logical_and[c0->clipSphere(pos,rad)]
                           [c1->clipSphere(pos,rad)];
  }
  virtual int clipCylinder(const Float3 &posa, const Float3 &posb, float rad)
  {
    return clip_logical_and[c0->clipCylinder(posa,posb,rad)]
                           [c1->clipCylinder(posa,posb,rad)];
  }
  virtual string asPovRayObject(const string &stylestr)
  {
    return str(format("union {\n %s %s}\n") % c0->asPovRayObject(stylestr)
                                                   % c1->asPovRayObject(stylestr));
  }
};


class ClipOr : public ClipBase
{
  boost::shared_ptr<ClipBase> c0, c1;
public:
  ClipOr(CP c0_, CP c1_) : c0(c0_), c1(c1_) {}
  virtual int clipSphere(const Float3 &pos, float rad)
  {
    return clip_logical_or[c0->clipSphere(pos,rad)]
                          [c1->clipSphere(pos,rad)];
  }
  virtual int clipCylinder(const Float3 &posa, const Float3 &posb, float rad)
  {
    return clip_logical_or[c0->clipCylinder(posa,posb,rad)]
                          [c1->clipCylinder(posa,posb,rad)];
  }
  virtual string asPovRayObject(const string &stylestr)
  {
    return str(format("intersection {\n %s %s}\n") % c0->asPovRayObject(stylestr)
                                            % c1->asPovRayObject(stylestr));
  }
};


enum ClipId {
  CLIP_ID_NONE = 0,
  CLIP_ID_PIE,
  CLIP_ID_SLICE
};


CP clipper_factory(const py::tuple &py_clip_data)
{
 CP cp;
  {
    int clip_id = py::extract<int>(py_clip_data[0]);
    switch (clip_id)
    {
      case CLIP_ID_PIE:
      {
        Float3 n1 = py::extract<Float3>(py_clip_data[1]);
        Float3 n2 = py::extract<Float3>(py_clip_data[2]);
        Float3 o  = py::extract<Float3>(py_clip_data[3]);
        cp.reset(new ClipAnd(
          CP(new ClipPlane(n1, o)),
          CP(new ClipPlane(n2, o))
                     ));
      }
      break;
      case CLIP_ID_SLICE:
      {
        Float3 n1 = py::extract<Float3>(py_clip_data[1]);
        Float3 n2 = py::extract<Float3>(py_clip_data[2]);
        Float3 o1  = py::extract<Float3>(py_clip_data[3]);
        Float3 o2  = py::extract<Float3>(py_clip_data[4]);
        cp.reset(new ClipOr(
          CP(new ClipPlane(n1, o1)),
          CP(new ClipPlane(n2, o2))
                     ));
      }
      break;
    }
  }
  return cp;
}

#if BOOST_VERSION>106300
void export_network_for_povray(const np::ndarray edges,
                               const np::ndarray pos,
                               const np::ndarray rad,
                               const py::object &py_styler,
                               const py::object &py_clip_styler,
                               const py::tuple &py_clip_data,
                               const std::string filename)
{
#ifndef NDEBUG
  cout << "edges.get_nd(): " << edges.get_nd() << endl;
  for(int i=0; i< edges.get_nd(); i++)
  {
    cout << "dim: " << i << "edges.get_shape()[i]: " << edges.get_shape()[i] << endl;
  }
  cout << "rad.get_nd(): " << rad.get_nd() << endl;
  for(int i=0; i< rad.get_nd(); i++)
  {
    cout << "dim: " << i << "rad.get_shape()[i]: " << rad.get_shape()[i] << endl;
  }
#endif
  int num_edges = edges.get_shape()[0];
  int num_nodes = pos.get_shape()[0];

  std::vector<float> noderad(num_nodes);
  for (int i=0; i<num_edges; ++i)
  {
//     noderad[edges[i,0]] = std::max(noderad[(int)edges[i,0]], rad[i]);
//     noderad[edges[i,1]] = std::max(noderad[edges[i,1]], rad[i]);
//     int index_noderad_a = py::extract<int const>(edges[i,0]);
//     noderad[index_noderad_a] = std::max(noderad[edges[index_noderad_a,0]], rad[i]);
//     noderad[edges[i]] = std::max(noderad[edges[i]], rad[i]);
    int a = py::extract<int>(edges[i][0]);
    int b = py::extract<int>(edges[i][1]);
    double radius = py::extract<float>(rad[i][0]);
    noderad[a] = std::fmax(noderad[a], radius);
    noderad[b] = std::fmax(noderad[b], radius);
  }

  CP cp = clipper_factory(py_clip_data);

  std::ofstream os(filename.c_str());

  for (int i=0; i<num_edges; ++i)
  {
    int a = py::extract<int>(edges[i][0]);
    int b = py::extract<int>(edges[i][1]);
    float pos_a_0 = py::extract<float>(pos[a][0]);
    float pos_a_1 = py::extract<float>(pos[a][1]);
    float pos_a_2 = py::extract<float>(pos[a][2]);
    Float3 pa(pos_a_0,pos_a_1,pos_a_2);
    float pos_b_0 = py::extract<float>(pos[b][0]);
    float pos_b_1 = py::extract<float>(pos[b][1]);
    float pos_b_2 = py::extract<float>(pos[b][2]);
    Float3 pb(pos_b_0,pos_b_1,pos_b_2);

    double radius = py::extract<float>(rad[i][0]);
    int intersect = cp ? cp->clipCylinder(pa, pb, radius) : CLIP_NONE;
    if (intersect == CLIP_FULL) continue;

    string stylestr = py::extract<string>(py_styler.attr("edge_style")(i, a, b));
    string objstr = str(format("cylinder {\n<%f,%f,%f>, <%f,%f,%f>, %f\n%s}\n") %
                                pa[0] % pa[1] % pa[2] % pb[0] % pb[1] % pb[2] % radius % stylestr);
    if (intersect == CLIP_PARTIAL)
    {
      string clip_stylestr = py::extract<string>(py_clip_styler.attr("edge_style")(i, a, b));
      string clip_str = cp->asPovRayObject(clip_stylestr);
      objstr = str(format("intersection {\n %s %s}\n") % clip_str % objstr);
    }

    os << objstr;
  }

  for (int i=0; i<num_nodes; ++i)
  {
    //Float3 p(pos(i,0),pos(i,1),pos(i,2));
    float pos_1 = py::extract<float>(pos[i][0]);
    float pos_2 = py::extract<float>(pos[i][1]);
    float pos_3 = py::extract<float>(pos[i][2]);
    Float3 p(pos_1, pos_2, pos_3);

    int intersect = cp ? cp->clipSphere(p, noderad[i]) : CLIP_NONE;
    if (intersect == CLIP_FULL) continue;

    string stylestr = py::extract<string>(py_styler.attr("node_style")(i));
    string objstr = str(format("sphere {\n<%f,%f,%f>, %f\n%s}\n") %
                                p[0] % p[1] % p[2] % noderad[i] % stylestr);
    if (intersect == CLIP_PARTIAL)
    {
      string clip_stylestr = py::extract<string>(py_clip_styler.attr("node_style")(i));
      string clip_str = cp->asPovRayObject(clip_stylestr);
      objstr = str(format("intersection {\n %s %s}\n") % clip_str % objstr);
    }

    os << objstr;
  }
}
#else
void export_network_for_povray(const np::arrayt<int> edges,
                               const np::arrayt<float> pos,
                               const np::arrayt<float> rad,
                               const py::object &py_styler,
                               const py::object &py_clip_styler,
                               const py::tuple &py_clip_data,
                               const std::string filename)
{
  int num_edges = edges.shape()[0];
  int num_nodes = pos.shape()[0];

  std::vector<float> noderad(num_nodes);
  for (int i=0; i<num_edges; ++i)
  {
    noderad[edges(i,0)] = std::max(noderad[edges(i,0)], rad(i));
    noderad[edges(i,1)] = std::max(noderad[edges(i,1)], rad(i));
  }

  CP cp = clipper_factory(py_clip_data);

  std::ofstream os(filename.c_str());

  for (int i=0; i<num_edges; ++i)
  {
    int a = edges(i,0);
    int b = edges(i,1);
    Float3 pa(pos(a,0),pos(a,1),pos(a,2));
    Float3 pb(pos(b,0),pos(b,1),pos(b,2));

    int intersect = cp ? cp->clipCylinder(pa, pb, rad(i)) : CLIP_NONE;
    if (intersect == CLIP_FULL) continue;

    string stylestr = py::extract<string>(py_styler.attr("edge_style")(i, a, b));
    string objstr = str(format("cylinder {\n<%f,%f,%f>, <%f,%f,%f>, %f\n%s}\n") %
                                pa[0] % pa[1] % pa[2] % pb[0] % pb[1] % pb[2] % rad(i) % stylestr);
    if (intersect == CLIP_PARTIAL)
    {
      string clip_stylestr = py::extract<string>(py_clip_styler.attr("edge_style")(i, a, b));
      string clip_str = cp->asPovRayObject(clip_stylestr);
      objstr = str(format("intersection {\n %s %s}\n") % clip_str % objstr);
    }

    os << objstr;
  }

  for (int i=0; i<num_nodes; ++i)
  {
    Float3 p(pos(i,0),pos(i,1),pos(i,2));

    int intersect = cp ? cp->clipSphere(p, noderad[i]) : CLIP_NONE;
    if (intersect == CLIP_FULL) continue;

    string stylestr = py::extract<string>(py_styler.attr("node_style")(i));
    string objstr = str(format("sphere {\n<%f,%f,%f>, %f\n%s}\n") %
                                p[0] % p[1] % p[2] % noderad[i] % stylestr);
    if (intersect == CLIP_PARTIAL)
    {
      string clip_stylestr = py::extract<string>(py_clip_styler.attr("node_style")(i));
      string clip_str = cp->asPovRayObject(clip_stylestr);
      objstr = str(format("intersection {\n %s %s}\n") % clip_str % objstr);
    }

    os << objstr;
  }
}
#endif

py::object pv_clip_object_str(const py::tuple &py_clip_data, const py::object &py_clip_styler)
{
  CP cp = clipper_factory(py_clip_data);
  if (!cp) return py::object();
  string stylestr = py::extract<string>(py_clip_styler);
  return py::str(cp->asPovRayObject(stylestr));
}


void export_povray_export()
{
  py::def("export_network_for_povray", export_network_for_povray);
  py::enum_<ClipId>("ClipShape")
    .value("slice", CLIP_ID_SLICE)
    .value("pie", CLIP_ID_PIE)
    .value("none", CLIP_ID_NONE);
  py::def("povray_clip_object_str", pv_clip_object_str);
}

