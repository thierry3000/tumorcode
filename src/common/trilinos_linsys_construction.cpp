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
#include "trilinos_linsys_construction.h"
#include "shared-objects.h"
#include "mwlib/timer.h"

//it is good to know the epetra settings before we set up mpi and multithreading
//#include <Epetra_ConfigDefs.h>
//everything done in trilinos_linsys_construction.h like it is ususal

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Version.hpp>
#include <BelosTpetraAdapter.hpp>
//#include <Ifpack2_Factory.hpp>
#include <AztecOO.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_VectorOut.h>
#include <Ifpack.h>
#include <Ifpack_IC.h>
#include <Ifpack_ILU.h>

/**
 * @brief Set up system with 7 point stencil
 * 
 * Adds void functions to the FiniteVolumeMatrixBuilder structure.
 */
void FiniteVolumeMatrixBuilder::FiniteVolumeMatrixBuilder::Init7Point(const LatticeDataQuad3d& ld_, int dim_)
{
  ld = ld_;
  dim = dim_;
  const int d = 3; // dimension
  const BBox3 bb(ld.Box());
  int num_dof = (bb.max - bb.min + BBox3::Vec(1)).prod();

  DynArray<int> num_entries(num_dof);
  FOR_BBOX3(p,bb)
  {
    int n = 1; // one for the diagonal element
    // cell has more neighbors if not at the boundary
    for(int i=0; i<3; ++i)
    {
      if(p[i]>bb.min[i]) n += 1;
      if(p[i]<bb.max[i]) n += 1;
    }
    int site = ld.LatticeToSite(p);
    num_entries[site] = n;
  }

#ifdef EPETRA_MPI
    #warning "Compiling with MPI Enabled"
    Epetra_MpiComm epetra_comm(MPI_COMM_SELF);
#else
    #warning "Compiling without MPI"
    Epetra_SerialComm epetra_comm;
#endif
  
  Epetra_Map epetra_map(num_dof, 0, epetra_comm);
  Epetra_CrsGraph graph(Copy, epetra_map, get_ptr(num_entries), true); // static profile, memory alloction is fixed, map object is copied

  DynArray<BBox3> mtboxes = MakeMtBoxGrid(bb);
  #pragma omp parallel
  {
    DynArray<int> columns(2 * d + 1, ConsTags::RESERVE); // temp buffer for column indices; each thread has its own.
    #pragma omp for nowait
    for (int i=0; i<mtboxes.size(); ++i)
    {
      FOR_BBOX3(p, mtboxes[i])
      {
        int site = ld.LatticeToSite(p);
        columns.push_back(site);
        for(int d=0; d<LatticeDataQuad3d::DIR_CNT; ++d)
        {
          const Int3 nbp = ld.NbLattice(p,d);
          if(!ld.IsInsideLattice(nbp)) continue;
          columns.push_back(ld.LatticeToSite(nbp));
        }
        graph.InsertGlobalIndices(site, columns.size(), get_ptr(columns)); // this function should better be thread safe, at least if the memory layout has been determined in advance!
        columns.remove_all();
      }
    }
  }
  
  graph.FillComplete();
  graph.OptimizeStorage();

  m.reset(new Epetra_CrsMatrix(Copy, graph)); // static profile, memory alloction is fixed, map object is copied
  rhs.reset(new Epetra_Vector(epetra_map));
}


/**
 * @brief Set up system with 9 point stencil
 * 
 * Adds void functions to the FiniteVolumeMatrixBuilder structure.
 */
void FiniteVolumeMatrixBuilder::FiniteVolumeMatrixBuilder::Init27Point(const LatticeDataQuad3d& ld_, int dim_)
{
  ld = ld_;
  dim = dim_;
  const BBox3 bb(ld.Box());
  int num_dof = Size(bb).prod();

  DynArray<int> num_entries(num_dof);
  FOR_BBOX3(p,bb)
  {
    int n = 1; // one for the diagonal element
    // cell has more neighbors if not at the boundary
    Int3 stencil_size(1);
    for(int i=0; i<3; ++i)
    {
      if(p[i]>bb.min[i]) stencil_size[i] += 1;
      if(p[i]<bb.max[i]) stencil_size[i] += 1;
    }
    int site = ld.LatticeToSite(p);
    num_entries[site] = stencil_size.prod();
  }
#ifdef EPETRA_MPI
    #warning "Compiling with MPI Enabled"
    Epetra_MpiComm epetra_comm(MPI_COMM_SELF);
#else
    #warning "Compiling without MPI"
    Epetra_SerialComm epetra_comm;
#endif
  Epetra_Map epetra_map(num_dof, 0, epetra_comm);
  Epetra_CrsGraph graph(Copy, epetra_map, get_ptr(num_entries), true); // static profile, memory alloction is fixed, map object is copied

  DynArray<BBox3> mtboxes = MakeMtBoxGrid(bb);
  #pragma omp parallel
  {
    DynArray<int> columns(my::ipow(3, dim), ConsTags::RESERVE); // temp buffer for column indices; each thread has its own.
    #pragma omp for nowait
    for (int i=0; i<mtboxes.size(); ++i)
    {
      FOR_BBOX3(p, mtboxes[i])
      {
        int site = ld.LatticeToSite(p);
        BBox3 b;
        int axis=0;
        for (; axis<dim; ++axis)
        {
          b.min[axis] = p[axis]-1;
          b.max[axis] = p[axis]+1;
        }
        for (; axis<3; ++axis)
          b.min[axis] = b.max[axis] = 0;
        b.Intersection(bb);
        FOR_BBOX3(q, b)
        {
          columns.push_back(ld.LatticeToSite(q));
        }
        graph.InsertGlobalIndices(site, columns.size(), get_ptr(columns)); // this function should better be thread safe, at least if the memory layout has been determined in advance!
        columns.remove_all();
      }
    }
  }

  graph.FillComplete();
  graph.OptimizeStorage();

  m.reset(new Epetra_CrsMatrix(Copy, graph)); // static profile, memory alloction is fixed, map object is copied
  rhs.reset(new Epetra_Vector(epetra_map));
}



void FiniteVolumeMatrixBuilder::SetDirichletBoundaryConditions(const BBox3 &bb, int bc, double boundary_value)
{
  //cout << "SetDirichletBoundaryConditions:" << endl;
  const bool on_x = bc & FiniteVolumeMatrixBuilder::DIRICHLET_X;
  const bool on_yz = bc & FiniteVolumeMatrixBuilder::DIRICHLET_YZ;
  Int3 ext((dim > 0 && on_x) ? -1 : 0,
           (dim > 1 && on_yz) ? -1 : 0,
           (dim > 2 && on_yz) ? -1 : 0);
  const BBox3 ldbox = ld.Box();
  const BBox3 bbox_inner_full = Extend(ldbox, ext);
  const BBox3 bbox_inner = Intersection(bb, bbox_inner_full);
  //BBox3 dirichletBBoxes[8]; dirichletBBoxes[0] = bb;
  //const int num_boxes = SplitBBox3(bbox_inner, dirichletBBoxes);
  DynArray<BBox3> dirichletBoxes;
  SplitBBoxByBox(bb, bbox_inner, dirichletBoxes, false, BBOXA);
  const int num_boxes = dirichletBoxes.size();

  //cout << "ldbox = " << bb << " bbox_inner = " << bbox_inner << endl;
  //for (int box_idx = 0; box_idx < num_boxes; ++box_idx)
  //{
  //  cout << "  box " << dirichletBoxes[box_idx] << endl;
  //}
  
  double  *values;
  int     *indices;
  int     num_entries;

  for (int box_idx = 0; box_idx < num_boxes; ++box_idx)
  {
    //cout << "  box " << dirichletBoxes[box_idx] << endl;
    FOR_BBOX3(p, dirichletBoxes[box_idx])
    {
      //cout << p << " nbfixes ";
      int site = ld.LatticeToSite(p);
      m->ExtractMyRowView(site, num_entries, values, indices);
      for (int i=0; i<num_entries; ++i)
      {
        if (indices[i] == site)
        {
          values[i] = 1.;
          (*rhs)[site] = boundary_value;
        }
        else
        {
          values[i] = 0.;
          // fix the neighbor
          int neighbor_site = indices[i];
          //cout << ld.SiteToLattice(neighbor_site) << ": ";
          double *nb_values = NULL;
          int    *nb_indices = NULL;
          int    nb_num_entries = 0;
          m->ExtractMyRowView(neighbor_site, nb_num_entries, nb_values, nb_indices);
          double diag_addition = 0;
          for (int j=0; j<nb_num_entries; ++j)
          {
            if (nb_indices[j] == site)
            {
              (*rhs)[neighbor_site] -= nb_values[j] * boundary_value;
              nb_values[j] = 0.;
              //cout << "nbidx=" << j << ",";
            }
          }
          //cout << " ";
        }
        //cout << endl;
      }
    }
  }
}


#if 0
void FiniteVolumeMatrixBuilder::InitFromStencil(const LatticeDataQuad3d& ld_, int dim, const ConstArray3d< bool >& pattern)
{
  myAssert(pattern.isContiguous() && pattern.getBox().min == Int3(0));

  ld = ld_;
  bb = ld.Box();
  int num_dof = product_of_components(Size(bb));

  const BBox3 bb_pattern = pattern.getBox();
  myAssert(bb_pattern.max == -bb_pattern.min);

  bb_inner = Extend(bb, -bb_pattern.max);
  boxes[0] = bb;
  int num_boxes = SplitBBox3(bb_inner, boxes);

  stencil_points.remove_all();
  FOR_BBOX3(q, bb_pattern)
  {
    if (pattern(q))
      stencil_points.push_back(q);
  }

  DynArray<int> num_entries(num_dof);
  for (int boxid = 0; boxid < num_boxes; ++boxid)
  {
    FOR_BBOX3(p, boxes[boxid])
    {
      int n = 0;
      if (boxid == 0)
      {
        // all stencil points contribute within the inner volume
        n = stencil_points.size();
      }
      else
      {
        // check each stencil point if it is inside the lattice
        BOOST_FOREACH(Int3 q, stencil_points)
        {
          for (int axis=0; axis<dim; ++axis)
          {
            if (q[axis]+p[axis]<bb.min[axis] || q[axis]+p[axis]>bb.max[axis]) continue;
            ++n;
          }
        }
      }
      num_entries[ld.LatticeToSite(p)] = n;
    }
  }

  Epetra_SerialComm epetra_comm;
  Epetra_Map epetra_map(num_dof, 0, epetra_comm);
  Epetra_CrsGraph graph(Copy, epetra_map, get_ptr(num_entries), true); // static profile, memory alloction is fixed, map object is copied

  DynArray<int> columns(stencil_points.size());

  for (int boxid = 0; boxid < num_boxes; ++boxid)
  {
    FOR_BBOX3(p, boxes[boxid])
    {
      int k=0;
      BOOST_FOREACH(Int3 q, stencil_points)
      {
        if (boxid == 0 || ld.IsInsideLattice(p+q))
        {
          columns[k++] = ld.LatticeToSite(p+q);
        }
      }
      graph.InsertGlobalIndices(ld.LatticeToSite(p), k, get_ptr(columns));
    }
  }

  graph.FillComplete();
  graph.OptimizeStorage();

  m.reset(new Epetra_CrsMatrix(Copy, graph)); // static profile, memory alloction is fixed, map object is copied
  rhs.reset(new Epetra_Vector(epetra_map));
}
#endif

enum EllipticSolveOutputLevel
{
  OUT_SILENT = 0,
  OUT_NORMAL = 1,
  OUT_FULL = 2,
  OUT_NORMAL_IN_DEBUG = 3
};

#if 1
EllipticEquationSolver::EllipticEquationSolver(RCP<Epetra_CrsMatrix> &_matrix, RCP<Epetra_Vector> &_rhs, const ptree& _params):
sys_matrix(_matrix), rhs(_rhs), params(_params)
{
  int success = init(sys_matrix,rhs,params);
}
int EllipticEquationSolver::init(RCP<Epetra_CrsMatrix> &_matrix, RCP<Epetra_Vector> &_rhs, const ptree& _params)
{
#ifdef EPETRA_MPI
  int MyPID = 0;
  int isMPIinitialized;
  int error = MPI_Initialized(&isMPIinitialized);
  if( isMPIinitialized == 1)
  {
    Epetra_MpiComm Comm (MPI_COMM_WORLD);
    MyPID = Comm.MyPID ();
  }
  else
  {
    //MPI_Init();
    printf("did you use the MPI wrapper?\n");
  }
#else
  Epetra_SerialComm Comm;
#endif
  
  ifpackList = Teuchos::rcp (new Teuchos::ParameterList ("my Teuchos List"));
  
  bool verbose = false;
  bool success = true;
  
  try {
    bool proc_verbose = true;
    int frequency = 10;        // frequency of status test output.
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxiters = 1000;         // maximum number of iterations allowed per linear 
    MT tol = 1.0e-5;           // relative residual tolerance

    sys_matrix->OptimizeStorage ();
    proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */

    // allocates an IFPACK factory. No data is associated
    // to this object (only method Create()).
    Ifpack Factory;

    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    
    //std::string PrecType = "point relaxation stand-alone"; // incomplete Cholesky
    std::string PrecType = params.get<string>("preconditioner","jacobi"); // default jacobi
    keep_preconditioner = params.get<bool>("keep_preconditioner", false);
    string sol_out_lvl = params.get<string>("verbosity","silent");
    
    int OverlapLevel = 0; // must be >= 0. If Comm.NumProc() == 1,
    // it is ignored.
    // create ifpack preconditioner
    if(sys_matrix.get() && PrecType=="multigrid" && (!keep_preconditioner || !ml_prec.get()))
    {
#ifdef DEBUG
      cout<<"hey"<<endl;
#endif
      Teuchos::ParameterList mllist;
      ML_Epetra::SetDefaults("SA",mllist);
      
      mllist.set("max levels", params.get<int>("max_levels", 10));
      mllist.set("cycle applications", 2);
      
      if (!params.get<bool>("use_smoothed_aggregation", false))
      mllist.set("aggregation: damping factor", 0.);
      
      mllist.set("smoother: type","Chebyshev"); // <---- this preconditioner ist much more effective than jacobi for complicated tumor vessel networks!!!
      
      mllist.set("smoother: pre or post", "both");
      mllist.set("smoother: sweeps", 2);
      
      mllist.set("coarse: max size",10000);
      mllist.set("coarse: type", "Amesos-UMFPACK");
      
      int ml_output = OUT_SILENT;
      if (sol_out_lvl == "max" || sol_out_lvl == "full")
        ml_output = 10;
      mllist.set("ML output", ml_output);
      
      try {
        ml_prec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*sys_matrix,mllist,true));
    } catch (std::exception& e) {
      std::ostringstream os;
      os << "ML preconditioner construction threw an exception: "
         << e.what ();
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str ());
    }
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp (new Belos::EpetraPrecOp (ml_prec));
    this->belos_Prec = belosPrec;
    }
//     else if(sys_matrix && prec_name=="ic" && (!keep_preconditioner || !ifpack_prec.get()))
//     {
//       ifpack_Prec.reset(new Ifpack_IC(const_cast<Epetra_RowMatrix*>(sys_matrix)));
//       ifpack_Prec->Initialize();
//       ifpack_Prec->Compute();
//     }
//     else if(sys_matrix && prec_name=="ilu" && (!keep_preconditioner || !ifpack_prec.get()))
//     {
//       ifpack_Prec.reset(new Ifpack_ILU(const_cast<Epetra_RowMatrix*>(sys_matrix)));
//       ifpack_Prec->Initialize();
//       ifpack_Prec->Compute();
//     }
    try {
      ifpack_Prec = Teuchos::rcp (Factory.Create (PrecType, sys_matrix.get(), OverlapLevel));
    } catch (std::exception& e) {
      std::ostringstream os;
      os << "Ifpack preconditioner construction threw an exception: "
         << e.what ();
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (ifpack_Prec.is_null (), std::runtime_error, "Failed to create Ifpack "
       "preconditioner!");

    try {
      ifpackList->set ("relaxation: type", "Jacobi");
    } catch( std::exception& e)
    {
      std::ostringstream os;
      os << "Teuchos parameter list construction threw an exception: "
         << e.what ();
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (ifpackList.is_null (), std::runtime_error, "Failed to create Ifpack "
       "list!");
    
    //
    // Create parameter list for the Belos solver
    //
    const int NumGlobalElements = _rhs->GlobalLength ();
    if (maxiters == -1) {
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
    }
    
    
    //Will be the arguments for the belos Solver
    Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::rcp (new Teuchos::ParameterList ("Belos"));
    belosList->set ("Maximum Iterations", maxiters);
    belosList->set ("Convergence Tolerance", tol);
    if (numrhs > 1) {
      // Show only the maximum residual norm
      belosList->set ("Show Maximum Residual Norm Only", true);
    }
    if (verbose) {
      belosList->set ("Verbosity", Belos::Errors + Belos::Warnings +
                     Belos::TimingDetails + Belos::StatusTestDetails);
      if (frequency > 0) {
        belosList->set ("Output Frequency", frequency);
      }
    }
    else {
      belosList->set ("Verbosity", Belos::Errors + Belos::Warnings +
                      Belos::FinalSummary);
    }
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
}

int EllipticEquationSolver::solve ( RCP<Epetra_Vector> _lhs )
{
  bool verbose = false;
  bool success = true;
  try {
    //
    // Construct a preconditioned linear problem
    // i think this should move to solve
    //
    bool leftprec = false;     // left preconditioning or right.
    bool proc_verbose = true;
    RCP<Belos::LinearProblem<double,MV,OP> > problem
      = rcp (new Belos::LinearProblem<double,MV,OP> (sys_matrix, _lhs, rhs));
    if (leftprec) {
      problem->setLeftPrec (belos_Prec);
    }
    else {
      problem->setRightPrec (belos_Prec);
    }
    bool set = problem->setProblem ();
    if (! set) {
      if (proc_verbose) {
        cout << endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << endl;
      }
      return -1;
    }

    // Create a Belos solver.
    RCP<Belos::SolverManager<double,MV,OP> > solver
      = rcp (new Belos::BiCGStabSolMgr<double,MV,OP> (problem, belosList));

//     if (proc_verbose) {
//       cout << endl << endl;
//       cout << "Dimension of matrix: " << NumGlobalElements << endl;
//       cout << "Number of right-hand sides: " << numrhs << endl;
//       cout << "Max number of CG iterations: " << maxiters << endl;
//       cout << "Relative residual tolerance: " << tol << endl;
//       cout << endl;
//     }
    // Ask Belos to solve the linear system.
    Belos::ReturnType ret = solver->solve();
    //vector_print(*_lhs);
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)

}
#endif


#if 1
int SolveEllipticEquation(Teuchos::RCP<Epetra_CrsMatrix> &matrix, Teuchos::RCP<Epetra_Vector> &rhs, Teuchos::RCP<Epetra_Vector> &lhs, const boost::property_tree::ptree &params)
{
  EllipticEquationSolver solver(matrix, rhs, params);
  solver.solve(lhs);
}
#endif

#if 0
//note: this is good for experimenting with trilinos
//do NOT delete
// void SolveEllipticEquation(const Epetra_CrsMatrix &matrix, const Epetra_Vector &rhs, Epetra_Vector &lhs, const boost::property_tree::ptree &params)
int SolveEllipticEquation(Teuchos::RCP<Epetra_CrsMatrix> _matrix, Teuchos::RCP<Epetra_Vector> _rhs, Teuchos::RCP<Epetra_Vector> _lhs, const boost::property_tree::ptree &params)
{
  cout<<"BELOS?????"<<endl;
  
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;
  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
  
    int MyPID = 0;
#ifdef EPETRA_MPI
  MPI_Init (&argc, &argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
  MyPID = Comm.MyPID ();
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;
  bool success = true;
  
  // This "try" relates to TEUCHOS_STANDARD_CATCH_STATEMENTS near the
  // bottom of main().  That macro has the corresponding "catch".
  
  
//   int output = params.get<int>("output", OUT_NORMAL_IN_DEBUG);
//   if (output == OUT_NORMAL_IN_DEBUG)
// #ifdef DEBUG
//     output = OUT_NORMAL;
// #else
//     output = OUT_SILENT;
// #endif
//   if (params.get<bool>("output_matrix", false))
//   {
//     std::ofstream f("matrix.txt");
//     matrix.Print(f);
//   }
// 
//   my::Time _t;
//   
//   // Set up Tpetra typedefs.
//   typedef double scalar_type;
//   typedef int local_ordinal_type;
//   typedef long global_ordinal_type;
//   typedef KokkosClassic::DefaultNode::DefaultNodeType node_type;
//   
//   typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> matrix_type;
//   typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> op_type;
//   typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> vec_type;
//   
//   // An Ifpack2::Preconditioner is-a Tpetra::Operator.  Ifpack2
//   // creates a Preconditioner object, but users of iterative methods
//   // want a Tpetra::Operator.  That's why create() returns an Operator
//   // instead of a Preconditioner.
//   typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, 
//                                   global_ordinal_type, node_type> prec_type;
//   
//   Teuchos::RCP<prec_type> prec;
//   Ifpack2::Factory factory2;
//   //prec = factory2.create("ILUT", &matrix);
//   Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> _ml_prec;
//   Teuchos::RCP<Ifpack_Preconditioner> _ifpack_preconditioner;
//   
//   Teuchos::oblackholestream blackHole;
//   Teuchos::RCP<const Teuchos::Comm<int>> comm=
//    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
//   
//   Teuchos::ParameterList mllist;
//   mllist.set ("Num Blocks", 40);
//   if (params.get<bool>("use_multigrid", true))
//   {
//     ML_Epetra::SetDefaults("SA",mllist);
// 
//     mllist.set("max levels", params.get<int>("max_levels", 10));
//     mllist.set("cycle applications", 2);
//     //mllist.set("prec type", "MGW");
// 
//     //mllist.set("increasing or decreasing","increasing");
//     //mllist.set("aggregation: type", "Uncoupled");
//     //mllist.set("aggregation: threshold", 0.1);
//     if (!params.get<bool>("use_smoothed_aggregation", false))
//       mllist.set("aggregation: damping factor", 0.);
//     //mllist.set("eigen-analysis: iterations", 3);
// 
//     //mllist.set("smoother: type","Chebyshev");
//     mllist.set("smoother: type","Jacobi");
// 
//     //mllist.set("smoother: pre or post", "both");
//     //mllist.set("smoother: sweeps", 2);
// 
//     //mllist.set("coarse: max size",128*128);
//     //mllist.set("coarse: type", "Chebyshev");
// 
//     mllist.set("coarse: max size",5000);
//     mllist.set("coarse: type", "Amesos-UMFPACK");
//     
//     // ML_Set_SpectralNormScheme_PowerMethod
//     // ML_Set_SpectralNormScheme_Calc
//     // ML_Aggregate_Set_DampingFactor
//     // ML_Aggregate_Set_NullSpace).
// 
//     int ml_output = output;
//     if (output == OUT_FULL)
//       ml_output = 10;
//     mllist.set("ML output", ml_output);
//     
//     //_ml_prec.reset(new  ML_Epetra::MultiLevelPreconditioner(matrix, mllist, true));
//   }
//   //_ml_prec->SetParameterList(mllist);
//   my::Time t_pc = my::Time() - _t;
// 
//   //Epetra_LinearProblem problem(const_cast<Epetra_CrsMatrix*>(&matrix), &lhs, const_cast<Epetra_Vector*>(&rhs));
//   //AztecOO solver(problem);
//   //AztecOO solver(const_cast<Epetra_CrsMatrix*>(&matrix), &lhs, const_cast<Epetra_Vector*>(&rhs));
//   typedef Belos::LinearProblem<BelosScalarType, BelosMultiVector, BelosOperator> BelosLinearProblem;
//   Teuchos::RCP<BelosLinearProblem> linear_problem = Teuchos::rcp(
//     //new BelosLinearProblem(Teuchos::rcp(const_cast<Epetra_CrsMatrix*>(&matrix),false),Teuchos::rcp(&lhs,false),Teuchos::rcp(&rhs,false)));
//     new BelosLinearProblem(Teuchos::rcp(&matrix,false),Teuchos::rcp(&lhs,false),Teuchos::rcp(&rhs,false)));
//   const bool success = linear_problem->setProblem();
//   
//   Teuchos::RCP<Teuchos::ParameterList> belosList(new Teuchos::ParameterList());
//   belosList->set("Num Blocks", 40);
// //   Teuchos::ParameterList belosList;
// //   belosList.set("Verbosity", Belos::TimingDetails);
// //   belosList.set("Num Blocks", 40);
//   //int verbosity = Belos::Error + Belos::Warnings;
//   //solver.SetOutputStream(cout);
//   
//   // Look up the Belos name of the method in _methods. This is a
//   // little complicated since std::maps<> don't have const lookup.
//   //std::map<std::string, std::string>::const_iterator it = _methods.find(_method);
// 
//   //if (it == _methods.end())
// 
//   //_ml_prec->set(*linear_problem, &lhs);
//   // set-up linear solver
//   Belos::SolverFactory<BelosScalarType, BelosMultiVector, BelosOperator> factory;
//   //Teuchos::RCP<Teuchos::ParameterList>
//   //Teuchos::RCP<Belos::SolverManager<BelosScalarType,BelosMultiVector,BelosOperator> > solver = factory.create("cg", Teuchos::rcp(&belosList, false));
//   Teuchos::RCP<Belos::SolverManager<BelosScalarType,BelosMultiVector,BelosOperator> > solver = factory.create("cg",belosList);
//   //linear_problem->setRightPrec(_ml_prec);
//   linear_problem->setProblem();
//   solver->setProblem(linear_problem);
//   
//   // Start solve
//   Belos::ReturnType ret = solver->solve();
//   
// 
//   //solver.SetAztecOption(AZ_precond, AZ_none);
// //   if (prec.get())
// //     solver.SetPrecOperator(prec.get());
// //   else
// //     solver.SetAztecOption(AZ_precond, AZ_Jacobi);
// //     //solver.SetAztecOption(AZ_precond, AZ_ls);
// //   const string solver_str = params.get("solver", "cg");
// //   if (solver_str == "bicgstab")
// //     solver.SetAztecOption(AZ_solver, AZ_bicgstab);
// //   else if (solver_str == "cg")
// //     solver.SetAztecOption(AZ_solver, AZ_cg);
// //   else
// //     throw std::runtime_error(str(format("wtf is solver %s") % solver_str));
// //   int az_output = AZ_none;
// //   switch(output)
// //   {
// //     case OUT_NORMAL: az_output = AZ_warnings; break;
// //     case OUT_FULL: az_output = AZ_all; break;
// //   }
// //   solver.SetAztecOption(AZ_output, az_output);
// //   solver.SetAztecOption(AZ_conv, AZ_noscaled); //AZ_Anorm);
// //   solver.Iterate(params.get<int>("max_iter", 50), params.get<double>("max_resid", 1.e-9));
// //   const double *status = solver.GetAztecStatus();
// //   if (status[AZ_why] != AZ_normal && status[AZ_why] != AZ_loss)
// //   {
// //     if (params.get<bool>("output_matrix_on_failure", false))
// //     {
// //       std::ofstream f("failmatrix.txt");
// //       matrix.Print(f);
// //     }
// //     //throw ConvergenceFailureException("linear system solve did not converge");
// //   }
// //   if (output >= OUT_NORMAL)
// //   {
// //     cout << boost::format("ElEq: time: %s (pc %s), iters: %i, residual: %e") % (my::Time() - _t) % t_pc % status[AZ_its] % status[AZ_r] << endl;
// //   }
}
#endif


/*------------------------------------------------------
------------------------------------------------------*/


template<class T>
void StationaryDiffusionSolve(const LatticeDataQuad3d &ld,
                              const DynArray<BBox3> &mtboxes,
                              int dim,
                              Array3d<T> result,
                              DiffSolveBuildFuncType buildfunc,
                              const ptree &pt_params)
{
  cout << format("stationary diffusion solve called!\n");
  my::Time t_;
  FiniteVolumeMatrixBuilder mb;
  mb.Init7Point(ld, dim);

  #pragma omp parallel for schedule(dynamic, 1)
  for (int i=0; i<mtboxes.size(); ++i)
  {
    buildfunc(i, mtboxes[i], mb);
  }

  optional<double> scaling = pt_params.get_optional<double>("system_scaling_factor");
  if (scaling)
  {
    mb.m->Scale(*scaling);
    mb.rhs->Scale(*scaling);
  }

  //Epetra_Vector lhs(mb.rhs->Map());
  Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(mb.rhs->Map()));
  

//   if (pt_params.get<bool>("try_without_multigrid", false))
//   {
//     #pragma omp parallel for schedule(dynamic, 1)
//     for (int i=0; i<mtboxes.size(); ++i)
//     {
//       FOR_BBOX3(p, mtboxes[i])
//       {
//         lhs[ld.LatticeToSite(p)] = result(p);
//       }
//     }
//     ptree pt = pt_params;
//     pt.put<bool>("use_multigrid", false);
//     try {
//       SolveEllipticEquation(*mb.m, *mb.rhs, lhs, pt);
//     }
//     catch(const ConvergenceFailureException &e)
//     {
//       pt.put<bool>("use_multigrid", true);
//       SolveEllipticEquation(*mb.m, *mb.rhs, lhs, pt);
//     }
//   }
//   else
//   SolveEllipticEquation(*mb.m, *mb.rhs, lhs, pt_params);
  SolveEllipticEquation(mb.m, mb.rhs, lhs, pt_params);

  #pragma omp parallel for schedule(dynamic, 1)
  for (int i=0; i<mtboxes.size(); ++i)
  {
    FOR_BBOX3(p, mtboxes[i])
    {
      //result(p) = lhs[ld.LatticeToSite(p)];
    }
  }
  cout << format("stationary diffusion solve total time: %f ms") % (my::Time()-t_).to_ms() << endl;
}


template<class T>
void StationaryDiffusionSolve(const ContinuumGrid &grid,
                              const DomainDecomposition &mtboxes,
                              Array3d<T> result,
                              DiffSolveBuildFuncType buildfunc,
                              const ptree &pt_params)
{
  my::Time t_;
  FiniteVolumeMatrixBuilder mb;
  mb.Init7Point(grid.ld, grid.dim);

  #pragma omp parallel
  {
    BOOST_FOREACH(const DomainDecomposition::ThreadBox bbox, mtboxes.getCurrentThreadRange())
    {
      buildfunc(bbox.global_index, bbox, mb);
    }
  }

  optional<double> scaling = pt_params.get_optional<double>("system_scaling_factor");
  if (scaling)
  {
    mb.m->Scale(*scaling);
    mb.rhs->Scale(*scaling);
  }

//   Epetra_Vector lhs(mb.rhs->Map());
  Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(mb.rhs->Map()));

//   SolveEllipticEquation(*mb.m, *mb.rhs, lhs, pt_params);
  SolveEllipticEquation(mb.m, mb.rhs, lhs, pt_params);
  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        //result(p) = *lhs[grid.ld.LatticeToSite(p)];
      }
    }
  }
  cout << format("stationary diffusion solve total time: %f ms") % (my::Time()-t_).to_ms() << endl;
}


#define INSTANTIATE_DIFFSOLVE(T) \
  template void StationaryDiffusionSolve<T>(const LatticeDataQuad3d &ld, const DynArray<BBox3> &mtboxes, int dim, Array3d<T> result, DiffSolveBuildFuncType buildfunc, const ptree &pt_params); \
  template void StationaryDiffusionSolve<T>(const ContinuumGrid &grid, const DomainDecomposition &mtboxes, Array3d<T> result, DiffSolveBuildFuncType buildfunc, const ptree &pt_params); 
  
INSTANTIATE_DIFFSOLVE(float)
INSTANTIATE_DIFFSOLVE(double)
/*------------------------------------------------------
------------------------------------------------------*/


