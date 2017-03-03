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
#include <boost/foreach.hpp>

#include <boost/property_tree/info_parser.hpp>
#include "mwlib/ptree_ext.h"

#include "common.h"
#include "continuum-utils.h"
#include "continuum-grid.h"

#include "mwlib/math_ext.h"

#include "calcflow.h"
//#include "mwlib/lattice-data.h"
#include "lattice-data-polymorphic.h"
#include "vessels3d.h"
#include <stdio.h>
#include "shared-objects.h"
//#include "mwlib/lattice-data.h"
#include "shared-objects.h"
#include "hdfio.h"

using namespace std;
//typedef VesselList3d::LatticeData LatticeData;
typedef polymorphic_latticedata::LatticeData LatticeData;
#if 0
struct Tester
{
  //create a lattice structure
  //std::auto_ptr<LatticeData> myLatticeData = LatticeData::Make("quad",BBox3(0,0,0,64,64,64), 10.);
  std::auto_ptr<LatticeData> myLatticeData = LatticeData::Make("quad",BBox3(0,0,0,64,64,64), 10.);
  std::auto_ptr<VesselList3d> vl;
  std::auto_ptr<LatticeData> ldp(LatticeData::Make("fcc", BBox3(0,0,0,10,10,10), 30.));
//   Vessel *v, *v2;
//   VesselNode *vc2;
//   VesselList3d vl;
//   vl.Init(*ldp);
  
  LatticeData &ld = *ldp;
  //LatticeDataQuad3d ld;
  void init()
  {
    ld.Init(Vec<int,3>(10, 10, 1), 10.);
    //vl->Init();
    printf("here1\n");
    vl->Init(*myLatticeData);
    printf("here2\n");
    for (int i=0; i<4; ++i)
    {
      Vessel* v  = vl->InsertVessel(
        Vec<int,3>(i,0,0),
        Vec<int,3>(i+1,0,0)
        );
    }
    #if 1
    // a second pipe between same endings
    vl->InsertVessel(Vec<int,3>(0,0,0), Vec<int,3>(0,1,0));
    for (int i=0; i<4; ++i)
    {
      Vessel* v  = vl->InsertVessel(
        Vec<int,3>(i,1,0),
        Vec<int,3>(i+1,1,0)
        );
    }
    vl->InsertVessel(Vec<int,3>(4,0,0), Vec<int,3>(4,1,0));
    #endif
    for (int i=0; i<vl->GetECount(); ++i)
    {
      vl->GetEdge(i)->r = 1.;
      vl->GetEdge(i)->flags |= CIRCULATED;
    }
    VesselNode* na = vl->FindNode(Vec<int,3>(0,0,0));
    VesselNode* nb = vl->FindNode(Vec<int,3>(4,0,0));
    na->press = 0;
    nb->press = 1;
    na->flags |= BOUNDARY;
    nb->flags |= BOUNDARY;
  }
  void solve_flow()
  {
    printf("testing flowsolver\n");
    BloodFlowParameters bfparams;
    printf("bfparams.includePhaseSeparationEffect: %s\n",bfparams.includePhaseSeparationEffect);
    printf("bfparams.inletHematocrit: %f\n",bfparams.inletHematocrit);
    printf("bfparams.rheology: %s\n",bfparams.rheology);
    printf("bfparams.viscosityPlasma: %f\n",bfparams.viscosityPlasma);
    CalcFlow(*vl, bfparams);
  }
  void output()
  {
    printf("viscos (1.,0.45) = %f\n", CalcRelViscosity(1., 0.45,Rheology::RheologyForHuman));
    for (int i=0; i<vl->GetNCount(); ++i)
    {
      const VesselNode* n = vl->GetNode(i);
      cout << boost::format("node: %i press = %f, pos = %s\n") % i % n->press % n->lpos;
    }
    for (int i=0; i<vl->GetECount(); ++i)
    {
      const Vessel* v = vl->GetEdge(i);
      cout << boost::format("edge %i-%i q = %f, f=%f\n") % v->NodeA()->Index() % v->NodeB()->Index() % v->q % v->f;
    }
  }
};
#endif

#if 0 //this used the dealii library
#include <grid/tria.h>
#include <grid/grid_tools.h>

#include <dofs/dof_handler.h>

#include <fe/fe_values.h>
#include <fe/mapping_q1.h>
#include <fe/fe_q.h>

#include <numerics/vectors.h>
#include <boost/foreach.hpp>

#include "common.h"
#include "helpers-vec.h"
#include "math_ext.h"
#include "timer.h"

// for testing
#include <numerics/matrices.h>
#include <grid/grid_generator.h>
#include <numerics/data_out.h>
#include <dofs/dof_tools.h>
#include <fstream>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/sparsity_pattern.h>
#include <base/function.h>

#include "dealii_common.h"
#include "dealii_cell_lookup.h"
#include "dealii_logstream.h"




void test_integrate_lines()
{
  using namespace dealii;

  for (int iteration = 1; iteration < 2; ++iteration)
  {
    my::Time t_;
    Triangulation<2>     triangulation;
    FE_Q<2> fe(1);
    DoFHandler<2> dof_handler(triangulation);

    //GridGenerator::hyper_cube (triangulation, -100, 100);
    //triangulation.refine_global (iteration);
    double spacing = 0.1;
    int base_level, max_level;
    initialize_rect_mesh<2>(dealii::Point<2>(100, 100), spacing, spacing, triangulation, base_level, max_level);
    print_mesh_parameters(cout, triangulation);
    
    dof_handler.distribute_dofs(fe);

    t_ = my::Time();
    CellLookup<2> lookup;
    BBoxd<double,2> bbox = calc_bounding_box<>(triangulation);
    lookup.initialize(to_dealii_point(bbox.min), to_dealii_point(bbox.max), 2 * GridTools::minimal_cell_diameter<>(triangulation), triangulation);
    printf("lookup init time: %f\n", (my::Time() - t_).to_ms());
    #if 0
    {
      std::ofstream f(boost::str(boost::format("lookup_test-%02i.gp") % iteration).c_str());
      lookup.get_tree().write_gnuplot(f);
    }
    #endif
    Double2 offset = Vec<double,2>(0., iteration * 1.);

    std::vector<double> edge_values;
    std::vector<dealii::Point<2> > points;
    std::vector<my::eqpair<int> > edges;
    for (int i=0; i<6; ++i)
    {
      edges.push_back(my::make_eqpair(i*2, i*2+1));
      Double2 p = 50. * Vec<double,2>(std::cos((my::mconst::pi()*i)/6), std::sin((my::mconst::pi()*i)/6));
      points.push_back(to_dealii_point<2>(p + offset));
      points.push_back(to_dealii_point<2>(-p + offset));
      edge_values.push_back(i == 0 ? 5 : 1);
    }

    double length = 0.;
    FOR_EACH2(my::eqpair<int> edge, edges, double value, edge_values, edges.size())
    {
      length += value * (points[edge[0]] - points[edge[1]]).norm();
    }

    Vector<double> values(dof_handler.n_dofs());

    t_ = my::Time();
    add_integral_over_lines(edges, edge_values, points, dof_handler, values, lookup);
    printf("integrate time per cell: %f, ncells: %i\n", (my::Time()-t_).to_ms()/triangulation.n_cells(), triangulation.n_cells());

    t_ = my::Time();
    
    #if 0
    SparsityPattern sparsity_pattern;
    {
      CompressedSparsityPattern csp;
      csp.reinit(dof_handler.n_dofs(), dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern (dof_handler, csp);
      sparsity_pattern.copy_from(csp);
    }
    SparseMatrix<double> mass;
    mass.reinit (sparsity_pattern);

    MatrixCreator::create_mass_matrix(MappingQ1<2>(),
                                      dof_handler,
                                      QGauss<2>(dof_handler.get_fe().degree+1),
                                      mass);

    SolverControl solver_control (500, 1e-12*values.l2_norm());

    PreconditionJacobi<SparseMatrix<double> > preconditioner;
    preconditioner.initialize(mass, 0.6);
    
    Vector<double> solution(dof_handler.n_dofs());

    SolverCG<Vector<double> >   cg (solver_control);
    cg.solve (mass, solution, values, preconditioner);
    #else
    Vector<double> solution(dof_handler.n_dofs());
    build_lumped_mass_vector<>(dof_handler, QGauss<2>(2), solution);
    for (int i=0; i<solution.size(); ++i)
      solution(i) = values(i)/solution(i);
    #endif
    printf("solve time per cell: %f\n", (my::Time()-t_).to_ms()/triangulation.n_cells());

    double integrated = VectorTools::compute_mean_value<>(dof_handler,
                                                          QGauss<2>(2),
                                                          solution,
                                                          0);
    integrated *= cwise_product(bbox.max-bbox.min);

    printf("edgelen: %f, integrated sol: %f\n", length, integrated);

    DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector(values, "lines");
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches (dof_handler.get_fe().degree-1);

    std::ofstream output(boost::str(boost::format("line_integration_test-%02i.vtu") % iteration).c_str());
    data_out.write_vtu(output);
  }
}


#endif



void TestLatticeData()
{
  {
    LatticeDataQuad3d ld, ldcc;
    ld.Init(Int3(10, 100, 100), 0.1);
    ld.SetOriginPosition(Float3(-0.5, 0, 0));
    //ld.SetZeroIndexPosition(Int3(50, 0, 0));
    ld.MoveSiteRange(Int3(50, 0, 0));
    //ld.SetBox(Move(ld.Box(), Int3(50, 0, 0)));
    ldcc = ld;
    ldcc.SetCellCentering(Vec<bool,3>(true));

    cout << "a lattice data" << endl;
    ldcc.print(cout);
    cout << endl;
    cout << "bb: " << ldcc.Box() << endl;

    float xmin = -0.5;
    float dx = 0.01;
    cout << "x    vc  fv    cc  fc" << endl;
    for (int i=0; i<101; ++i)
    {
      float x = xmin + i * dx;
      Int3 ip, ipcc; Float3 fp, fpcc;
      boost::tie(ip, fp) = ld.WorldToLatticeCell(Float3(x, 0, 0));
      boost::tie(ipcc, fpcc) = ldcc.WorldToLatticeCell(Float3(x, 0, 0));
      printf("%f    %i  (%f)        %i   (%f)\n", x, ip[0], fp[0], ipcc[0], fpcc[0]);
    }

    for (int i=ld.Box().min[0]; i<=ld.Box().max[0]; ++i)
    {
      printf("%i    %f     %f     %li\n", i, ld.LatticeToWorld(Int3(i,0,0))[0], ldcc.LatticeToWorld(Int3(i,0,0))[0], ld.LatticeToSite(Int3(i,0,0)));
    }
  }
  {
    typedef LatticeDataFCC::SiteType SiteType;
    BBox3 bbc;
    LatticeDataFCC ld, ldc;
    ld.Init(BBox3(-5,-5,-5,5,5,5), 1.);
    int subdiv = 2;
    //ld.Init(BBox3(0,0,0,5,5,5), 1.);
    FOR_BBOX3(p, ld.Box())
    {
      Int3 q = ld.GetLatticeIndexOnRefinedGrid(p, subdiv);
      bbc.Add(q);
    }
    FOR_BBOX3(p, ld.Box())
    {
      SiteType s = ld.LatticeToSite(p);
      Int3 q = ld.SiteToLattice(s);
      cout << format("%s -> %i -> %s") % p % s % q << endl;
      assert(q == p);
    }
    FOR_BBOX3(p, ld.Box())
    {
      SiteType s = ld.LatticeToSite(p);
      cout << format("%s (%i):") % p % s << endl;
      for (int i=0; i<ld.NbCount(); ++i)
      {
        SiteType snb = ld.NbSite(s, i);
        Int3 snbp = ld.SiteToLattice(snb);
        Int3 nbp = ld.NbLattice(p, i);
        if (!ld.IsInsideLattice(nbp)) continue;
        cout << format("  %s (%i) == %s") % snbp % snb % nbp << endl;
        assert(snbp == nbp);
      }
    }
    ldc.Init(bbc, ld.Scale()/(subdiv+1));
    FOR_BBOX3(p, ld.Box())
    {
      Int3 q = ld.GetLatticeIndexOnRefinedGrid(p, subdiv);
      Float3 wp = ld.LatticeToWorld(p),
             wq = ldc.LatticeToWorld(q);
      bool ok = (wp-wq).norm() < 1.e-3;
      cout << format("%i %s = %s -> %s") % ok % p % wp % wq << endl;
      assert(ok);
    }
  }
}


void TestInterpolation()
{
  LatticeDataQuad3d ld;
  ld.Init(Int3(1, 10, 1), my::mconst::fpi2() * 0.1);
  ld.SetCellCentering(Vec<bool,3>(true));
  FloatBBox3 wbox = ld.GetWorldBox();
  cout << "world box: " << wbox << endl;
  Array3d<float> field; field.initFromBox(ld.Box());
  FOR_BBOX3(p, ld.Box())
  {
    Float3 wp = ld.LatticeToWorld(p);
    field(p) = std::sin(wp[1]);
  }

  ofstream os("interpol-test.dat");

  //FieldInterpolator<float> interpolate(CONT_EXTRAPOLATE);

  double w = 0.25*(wbox.max[1]-wbox.min[1]);
  double x0 = wbox.min[1]-w;
  double x1 = wbox.max[1]+w;
  int n_samples = ld.Size()[1] * 10;
  for (int i=0; i<n_samples; ++i)
  {
    double xf = double(i) / n_samples;
    double x = xf*x1 + (1.-xf)*x0;
    Float3 wp(0.5, x, 0.5);
    Float3 g(0, std::cos(wp[1]), 0);
    double f = std::sin(wp[1]);
    //double int_f = interpolate.value(field, ld, wp);
    double int_f = FieldInterpolate::ValueAveraged(field, ld, FieldInterpolate::Extrapolate(), wp);
    //Float3 int_g = interpolate.gradient(field, ld, wp);
    Float3 int_g = FieldInterpolate::Gradient(field, ld, FieldInterpolate::Extrapolate(), wp);
    os << boost::format("%f %f %f %f %f\n") % x % f % int_f % g[1] % int_g[1];
  }
}



void TestVesselList()
{
  boost::shared_ptr<LatticeData> ldp(LatticeData::Make("fcc", BBox3(0,0,0,10,10,10), 30.));
  Vessel *v, *v2;
  VesselNode *vc2;
  VesselList3d vl(ldp);
  vl.Init(*ldp);
  v = vl.InsertVessel(Int3(0,5,5), Int3(9,5,5));
  vl.IntegrityCheck();
  vl.SplitVessel(v, 5, v2, vc2);
  vl.InsertVessel(Int3(5,5,5), ldp->NbLattice(Int3(5,5,5), 11));
  vl.IntegrityCheck();
}


void TestBoxGrid()
{
  BBox3 boxes[] = {
    BBox3(0,0,0,0,0,0),
    BBox3(-50, 0, 0, 50, 0, 0),
    BBox3(0, -50, 0, 0, 50, 0),
    BBox3(0, 0, -50, 0, 0, 50),
    BBox3(),
    BBox3(-50,-50,0, 50, 50 ,0),
    BBox3(-50,-50,-50, 50, 50 ,50)
  };
  BOOST_FOREACH(const BBox3 &bb, boxes)
  {
    cout << format("box grid for %s is") % bb << endl;
    DynArray<BBox3> bg = MakeMtBoxGrid(bb);
    BOOST_FOREACH(const BBox3 &bl, bg)
    {
      cout << "  " << bl << endl;
    }
  }
}


void TestDeltaAndHeaviside()
{
  double width = 2.;
  int N = 100;
  std::ofstream os("deltaheavitest.txt");

  for (int i=0; i<N; ++i)
  {
    double x = double(i)/(N-1)*10.-5.;
    os << format("%f %f %f\n") % x % my::smooth_delta_cos(x, width) % my::smooth_heaviside_sin(x, width);
  }
}


void TestVesselFieldCoupling()
{
  //LatticeDataQuad3d vessel_ld;
  //create a lattice structure
  boost::shared_ptr<LatticeData> vessel_ld(LatticeData::Make("quad",BBox3(0,0,0,64,64,64), 10.));
  //create a vessel list
  VesselList3d vl(vessel_ld);
  //link the lattice and the vessel
  vl.Init(*vessel_ld);
  //add a vessel
  Vessel* v = vl.InsertVessel(Int3(4,32,32), Int3(60, 32, 32));
  //create a field and initialize it
  LatticeDataQuad3d ld;
  ld.Init(BBox3(0,0,0, 20, 20, 20), 30.);
  ld.SetCellCentering(Bool3(true));

  h5cpp::File f("vesseltofield.h5", "w");

  cout << format("cellvolume=%f") % my::cubed(ld.Scale()) << endl;

  for (int k=0; k<10; ++k)
  {
    v->r = 3. + k*10.;
    v->flags |= CIRCULATED;
    v->NodeA()->flags |= BOUNDARY;
    v->NodeB()->flags |= BOUNDARY;
    v->NodeA()->press = 1;
    v->NodeA()->press = 2;
//     for (int i=0; i<vl->GetECount(); ++i)
//     {
//       vl->GetEdge(i)->r = 1.;
//       vl->GetEdge(i)->flags |= CIRCULATED;
//     }
//     VesselNode* na = vl->FindNode(Vec<int,3>(0,0,0));
//     VesselNode* nb = vl->FindNode(Vec<int,3>(4,0,0));
//     na->press = 0;
//     nb->press = 1;
//     na->flags |= BOUNDARY;
//     nb->flags |= BOUNDARY;
    
    Array3d<float> volume_field(ld.Box());
    Array3d<float> surface_field(ld.Box());
    Array3d<float> length_field(ld.Box());
    //set up a sampler
    CylinderNetworkSampler sampler; sampler.Init(ld.Scale(), make_ptree("seed", 12345)("samples_per_cell", 1));
    sampler.Set(vessel_ld->LatticeToWorld(v->LPosA()), vessel_ld->LatticeToWorld(v->LPosB()), v->r);
    //sample lines
    int line_cnt = sampler.GenerateLineSamples();
    for (int i=0; i<line_cnt; ++i)
      AddSmoothDelta(length_field,ld.Box(), ld, 3, sampler.GetSample(i).wpos, sampler.weight_per_volume);
    sampler.Restart();
    //sample surface
    int surf_cnt = sampler.GenerateSurfaceSamples();
    for (int i=0; i<surf_cnt; ++i)
      AddSmoothDelta(surface_field,ld.Box(), ld, 3, sampler.GetSample(i).wpos, sampler.weight_per_volume);
    sampler.Restart();
    //sample volume
    int vol_cnt = sampler.GenerateVolumeSamples();
    for (int i=0; i<vol_cnt; ++i)
      AddSmoothDelta(volume_field,ld.Box(), ld, 3, sampler.GetSample(i).wpos, sampler.weight_per_volume);

    // compare real stuff with sampled
    float len = v->WorldLength(*vessel_ld);
    float surface = my::mconst::pi2()*v->r * len;
    float volume = my::mconst::pi()*my::sqr(v->r)*len;
    float estimated_volume = volume_field.valueStatistics().Sum() * my::cubed(ld.Scale());
    float estimated_surface = surface_field.valueStatistics().Sum() * my::cubed(ld.Scale());
    float estimated_length = length_field.valueStatistics().Sum() * my::cubed(ld.Scale());
    
    cout << format("-- r=%f, l=%f ---") % v->r % len << endl;
    cout << format("%i=samples, volume=%f, estimated volume=%f, dev=%f%%") % vol_cnt % volume % estimated_volume % (abs(estimated_volume-volume) / volume)  << endl;
    cout << format("%i=samples, surface=%f, estimated surface=%f, dev=%f%%") % surf_cnt % surface % estimated_surface % (abs(estimated_surface-surface) / surface) << endl;
    cout << format("%i=samples, length=%f, estimated length=%f, dev=%f%%") % line_cnt % len % estimated_length % (abs(estimated_length-len) / len) << endl;

    h5cpp::Group g = f.root().create_group(str(format("out%04i") % k));
    h5cpp::Group ld_group = g.create_group("field_ld");
    WriteHdfLd(ld_group, ld);
    WriteScalarField(g, "volume_field", volume_field[ld.Box()], ld, ld_group);
    WriteScalarField(g, "surface_field", surface_field[ld.Box()], ld, ld_group);
    WriteScalarField(g, "length_field", length_field[ld.Box()], ld, ld_group);
    //to do some rendering with povray
    h5cpp::Group vesselgrp = g.create_group("vessels");
    vesselgrp.attrs().set<std::string>("CLASS","GRAPH");
    WriteVesselList3d(vl, vesselgrp, make_ptree("w_all",false)("w_pressure",false)("w_flow",false));
  }
}


void WriteSampling(const VesselList3d &vl, const LatticeDataQuad3d &ld, int dim, h5cpp::Group ld_group,  h5cpp::Group g, const string &name, const ptree &params)
{
  Array3d<float> field(ExtendForDim(ld.Box(), dim, 2));

  ptree res = AddVesselVolumeFraction(vl, ld, dim, field, params);
  int cnt = res.get<int>("num_samples");

  cout << format("vessels=%i, samples=%i") % vl.GetECount() % cnt << endl;

  WriteScalarField(g, name, field[ld.Box()], ld, ld_group);
}


// void TestVesselFieldCoupling2(int argc, char** argv)
// {
//   string fn = argv[1];
//   uint seed = 12345;
//   if (argc > 2)
//   {
//     std::istringstream ss;
//     //(argv[2]);
//     //ss >> seed;
//     //if (argc > 3)
//     {
//       ss.str(argv[2]);
//       int np;
//       ss >> np;
//       my::SetNumThreads(np);
//     }
//   }
//   
//   h5cpp::File f(fn, "r");
//   h5cpp::Group vess_grp = f.root().open_group("vessels"),
//                 ld_grp = f.root().open_group("vessels/lattice");
//   std::auto_ptr<LatticeData> ldp(LatticeData::ReadHdf(ld_grp));
//   LatticeData &ld = *ldp;
//   ld.SetOriginPosition(-ld.GetWorldBox().max*0.5);
//   std::auto_ptr<VesselList3d> vl;
//   vl->Init(ld);
//   //VesselList3d vl;
//   //vl.Init(ld);
//   ReadHdfGraph(vess_grp, vl.get());
// 
//   cout << "read" << endl;
//   ld.print(cout);
//   cout << endl;
// 
//   Int3 vlldsize = Size(vl.Ld().Box());
//   int dim = vlldsize[2]<=1 ? (vlldsize[1]<=1 ? 1 : 2) : 3;
//   LatticeDataQuad3d field_ld;
//   SetupFieldLattice(vl.Ld().GetWorldBox(), dim, 30, 0, field_ld);
//   cout << "field ld: "; field_ld.print(cout); cout << endl;
//   
//   f = h5cpp::File(RemovePath(RemoveExtension(fn))+"-rasterized.h5", "w");
//   h5cpp::Group ld_group = f.root().create_group("field_ld30");
//   WriteHdfLd(ld_group, field_ld);
//   
//   ptree params = make_ptree("seed", seed)("cut_at", 1.)("samples_per_cell", 1);
//   WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("basecase"), params);
// 
//   WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("noclip"), make_ptree(params)("cut_at", -1));
//   WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("samples2"), make_ptree(params)("samples_per_cell", 4));
// 
//   WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("seed1234"), make_ptree(params)("seed", 1234));
//   WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("seed2345"), make_ptree(params)("seed", 2345));
//   WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("seed3456"), make_ptree(params)("seed", 3456));
//   WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("seed4567"), make_ptree(params)("seed", 4567));
// }
/* you need to provide an iff simulation dataset here */
void TestVesselFieldCoupling3(int argc, char** argv)
{
  string fn = argv[1];
  uint seed = 12345;
  if (argc > 2)
  {
    std::istringstream ss;
    //(argv[2]);
    //ss >> seed;
    //if (argc > 3)
    {
      ss.str(argv[2]);
      int np;
      ss >> np;
      my::SetNumThreads(np);
    }
  }
  
  h5cpp::File f(fn, "r");
  h5cpp::Group vess_grp = f.root().open_group("iff/vessels"),
                ld_grp = f.root().open_group("iff/vessels/lattice");
  boost::shared_ptr<LatticeData> ldp(LatticeData::ReadHdf(ld_grp));
  LatticeData &ld = *ldp;
  ld.SetOriginPosition(-ld.GetWorldBox().max*0.5);
  //VesselList3d vl;
  //vl.Init(ld);
  //ReadHdfGraph(vess_grp, vl);
  std::auto_ptr<const VesselList3d> vl;
  vl = ReadVesselList3d(vess_grp, make_ptree("filter", false));

  cout << "read" << endl;
  ld.print(cout);
  cout << endl;

  Int3 vlldsize = Size(vl->Ld().Box());
  int dim = vlldsize[2]<=1 ? (vlldsize[1]<=1 ? 1 : 2) : 3;
  LatticeDataQuad3d field_ld;
  SetupFieldLattice(vl->Ld().GetWorldBox(), dim, 30, 0, field_ld);
  cout << "field ld: "; field_ld.print(cout); cout << endl;
  
  f = h5cpp::File(RemovePath(RemoveExtension(fn))+"-rasterized.h5", "w");
  h5cpp::Group ld_group = f.root().create_group("field_ld30");
  WriteHdfLd(ld_group, field_ld);
  
  ptree params = make_ptree("seed", seed)("cut_at", 1.)("samples_per_cell", 1);
  WriteSampling(*vl, field_ld, dim, ld_group, f.root(), ("basecase"), params);

  //WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("noclip"), make_ptree(params)("cut_at", -1));
  //WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("samples2"), make_ptree(params)("samples_per_cell", 4));

  //WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("seed1234"), make_ptree(params)("seed", 1234));
  //WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("seed2345"), make_ptree(params)("seed", 2345));
  //WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("seed3456"), make_ptree(params)("seed", 3456));
  //WriteSampling(vl, field_ld, dim, ld_group, f.root(), ("seed4567"), make_ptree(params)("seed", 4567));
}


void TestDeltaFunctionAndInterpolation()
{
  LatticeDataQuad3d ld(Int3(5, 5, 1), 3.);
  ld.SetCellCentering(Bool3(true, true, false));
  Array3df field(ld.Box());
  Array3df field2(ld.Box());
  field2(Int3(2, 2, 0)) = field2(Int3(0, 2, 0)) = field2(Int3(4, 2, 0)) = 1./my::sqr(ld.Scale()*2.);
  Float3 wp0(7.5, 7.5, 0.);
  AddSmoothDelta<float>(field, ld.Box(), ld, 2, wp0, 1./my::sqr(ld.Scale()));
  DynArray<Vec<float, 7> > result;
  Float3 p0(0, 7.5, 0.);
  const int N = 200;
  Float3 dp(15./(N-1), 0, 0);
  for (int i = 0; i<N; ++i)
  {
    Float3 wp = p0 + i*dp;
    Vec<float, 7> tmp;
    tmp[0] = wp[0];
    tmp[1] = FieldInterpolate::Value(field, ld, FieldInterpolate::Const<float>(0.),wp);
    tmp[2] = my::smooth_delta_cos<float>((wp-wp0)[0], ld.Scale());
    tmp[3] = FieldInterpolate::Value(field2, ld, FieldInterpolate::Const<float>(0.),wp);
    tmp[4] = field(ld.WorldToLattice(wp));
    tmp[5] = FieldInterpolate::ValueAveraged(field, ld, FieldInterpolate::Const<float>(0.), wp);
    tmp[6] = FieldInterpolate::Gradient(field, ld, FieldInterpolate::Const<float>(0.), wp)[0];
    result.push_back(tmp);
  }
  cout << ld << endl;
  h5cpp::File f("interpolate.h5", "w");
  //h5cpp::create_dataset_range(f.root(), "data", result);
  h5cpp::create_dataset(f.root(), "data", h5cpp::Dataspace::simple_dims(N, 7), result[0].data());
}


void testptree()
{
  ptree a, b;
  { std::ifstream f("pttesta.info");
  boost::property_tree::read_info(f, a); }
  boost::property_tree::write_info(cout, a);

  cout << " - " << endl;


  { std::ifstream f("pttestb.info");
  boost::property_tree::read_info(f, b); }
  boost::property_tree::write_info(cout, b);

  ptree c = boost::property_tree::remove(a, b);

  cout << " = " << endl;
  
  boost::property_tree::write_info(cout, c);
}



void TestConvolution()
{
  BBox3 bb1(5,5,0, 10, 20, 0);
  LatticeDataQuad3d ld(bb1, 1./5);

  Array3d<float> fun = MakeArray3dWithBorder<float>(bb1, 2, 1, MAKE_ARRAY3D_BOX_OUTER);
  FOR_BBOX3(p, fun.getBox())
  {
    const Float3 wp = ld.LatticeToWorld(p);
    fun(p) = sin(wp.x() * my::mconst::pi()) + wp.y()*wp.y();
  }

  Array3d<float> diff[2] = {
    Stencils::Diff1(0, 2, ld.Scale(), Stencils::W2, Stencils::FWD),
    Stencils::Diff1(1, 2, ld.Scale(), Stencils::W1, Stencils::CENTERED)
  };
  diff[0].move(Int3(-1,-1,0));

  Array3d<float> res1(CellRangeToVertexRange(bb1, 2));
  Array3d<float> res2(bb1);
  Convolutor<float, float> diff1_conv(diff[0], fun);
  Convolutor<float, float> diff2_conv(diff[1], fun, true);
  FOR_BBOX3(p, res1.getBox())
  {
    res1(p) = diff1_conv.point_convolve(p);
  }
  FOR_BBOX3(p, res2.getBox())
  {
    res2(p) = diff2_conv.point_convolve(p);
  }

  Image img, bigimg;
  std::vector<Image> images;
  images.push_back(DrawArray(fun[fun.getBox()], DrawArrayOpts().title("fun").outputRange()));
  images.push_back(DrawArray(res1[res1.getBox()], DrawArrayOpts().title("dx").outputRange()));
  images.push_back(DrawArray(res2[res2.getBox()], DrawArrayOpts().title("dy").outputRange()));
  DrawImageGrid(bigimg, images);
  bigimg.Write("convolution_test.png");
}


void TestArray3dRemoveDims()
{
  BBox3 bb(0,0,0,10,0,5);
  Array3d<int> a(bb);
  a.RemoveSingularDims();
  cout << a.getBox() << endl;
}


void TestVecDebug()
{
  Vec<float, 3> v;
  cout << v << endl;
}


void TestH5Overwrite()
{
  {
    h5cpp::File f("test.h5","a");
    h5cpp::Group g = f.root().create_group("testgroup");
    f.root().create_group("testgroup2");
    g.attrs().set("fail", 1);
  }
  {
    h5cpp::File f("test.h5","w");
    h5cpp::Group g = f.root().create_group("testgroup");
    g.attrs().set("fail", 2);
    f.root().create_group("testgroup3");
  }
  // works all right
}
//#endif

void TestMath()
{
  cout << format("iceil(%f) = %i\n") % 2.3 % my::iceil(2.3);
  cout << format("iceil(%f) = %i\n") % -2.3 % my::iceil(-2.3);
  cout << format("iceil(%f) = %i\n") % 0. % my::iceil(0.);
  cout << format("iceil(%f) = %i\n") % 0.3 % my::iceil(0.3);
  cout << format("iceil(%f) = %i\n") % -0.3 % my::iceil(-0.3);
  cout << format("iceil(%f) = %i\n") % -5. % my::iceil(-5.);
  cout << format("iceil(%f) = %i\n") % 5. % my::iceil(5.);
  cout << format("ifloor(%f) = %i\n") % 2.3 % my::ifloor(2.3);
  cout << format("ifloor(%f) = %i\n") % -2.3 % my::ifloor(-2.3);
  cout << format("ifloor(%f) = %i\n") % 0. % my::ifloor(0.);
  cout << format("ifloor(%f) = %i\n") % 0.3 % my::ifloor(0.3);
  cout << format("ifloor(%f) = %i\n") % -0.3 % my::ifloor(-0.3);
  cout << format("ifloor(%f) = %i\n") % -5.0 % my::ifloor(-5.);
  cout << format("ifloor(%f) = %i\n") % 5. % my::ifloor(5.);
}
void TestBB()
{
  BBoxd<float, 2> bb;
  bb.Add(Vec<float,2>(-1, -1));
  bb.Add(Vec<float,2>(-0.5, 1));
  cout << bb << endl;
}
#if 0
void TestTester()
{
  Tester test;
  test.init();
  test.solve_flow();
  test.output();
}
#endif
int main(int argc, char **argv)
{
  TestMath();
  TestBB();
#if 0
  TestTester();
#endif
  my::MultiprocessingInitializer mpinit(argc, argv, 4);
  //TestArray3dRemoveDims();
  //TestConvolution();
  //testptree();
  //TestDeltaAndHeaviside();
  //TestVesselFieldCoupling();
  //TestVesselFieldCoupling2(argc, argv);
  //TestVesselFieldCoupling3(argc,argv);
  //TestVesselList();
  //TestLatticeData();
  //TestVecDebug();
  //TestInterpolation();
  //TestH5Overwrite();
  //TestDeltaFunctionAndInterpolation();
  
  #if 0
  test_integrate_lines();
  #endif
  //TestLatticeData();
  //TestInterpolation();
  //TestBoxGrid();
  #if 0
  my::Log log(cout);
  Float3 v(-1.);
  //log.push();
  log << v << " " << Abs(v) << v.cwise().abs().eval() << endl;
  //log.pop();
  log << "testing" << endl;
  log.push();
  log << "test" << " new " << " line " << endl;
  log.pop();
  #endif
#if 1
  
#endif
  return 0;
}
