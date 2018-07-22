# include <iostream>
//# include <H5Cpp.h> //independent of tumorcode

#define with_tumorcode_support

#ifdef with_tumorcode_support
# include "../../common/hdfio.h"
# include "../../mwlib/dynamicarray.h"
#endif

int main()
{
  H5::H5File f;
  f = H5::H5File("dbg_stack.h5", H5F_ACC_TRUNC);
  
  DynArray<double> a(1100001);
  //DynArray<double> c(1e10);
  DynArray<double> *c = new DynArray<double>(1e10);
  //1100001 was not working on the stack, when not using boost multi array
  DynArray<double> *b = new DynArray<double>(1100001);
  double foo[20];
  H5::Group root;
  root = f.openGroup("/");
  writeDataSetToGroup(root, string("test_set"), a);
  std::cout << "created file" << std::endl;
}
//   //h5::File f("dbgvessels.h5", number==0 ? "w" : "a");
//   H5::Group grp = f.createGroup(str(format("%s") % name));
//   H5::Group vesselgrp = grp.createGroup("vessels");
//   H5::Group h5_nodes = vesselgrp.openGroup("nodes");
//   //since world coordinates this is also needed for proper hdf output
//   //vesselgrp.attrs().set<std::string>("CLASS","GRAPH");
//   writeAttrToH5(vesselgrp, string("CLASS"), string("GRAPH"));
//   WriteVesselList3d(vl, vesselgrp, make_ptree("w_all",false)("w_pressure",true)("w_flow",true));
// #if GFFIELD_ENABLE
//   {
//     DynArray<float> gf;
//     grower.GetGfAtNodes(gf);
//     //h5::create_dataset(vesselgrp.open_group("nodes"), "gf", gf);
//     
//     writeDataSetToGroup(h5_nodes, string("gf"), gf);
//     //h5::create_dataset(vesselgrp.open_group("nodes"), "gf", gf);
//   }
//   H5::Group field_ld_grp = grp.createGroup("field_ld");
//   //WriteHdfLd(field_ld_grp, grower.get_field_ld());
//   //grower.get_field_ld().Lattice2Hdf(field_ld_grp);
//   grower.get_field_ld().WriteHdfLd(field_ld_grp);
//   WriteScalarField(grp, "gf", grower.GetGf(), grower.get_field_ld(), field_ld_grp);
// #endif
//   {
//     DynArray<uchar> tmp2(vl.GetNCount());
//     for (int i=0; i<vl.GetNCount(); ++i)
//     {
//       tmp2[i] = vl.GetNode(i)->flags;
//     }
//     //h5::create_dataset(vesselgrp.open_group("nodes"), "nodeflags", tmp2);
//     writeDataSetToGroup(h5_nodes, string("nodeflags"), tmp2);
// #ifdef WRITE_REMODELING_ACTIONS    
//     for (int i=0; i<vl.GetNCount(); ++i)
//     {
//       tmp2[i] = grower.last_remodeling_action(vl.GetNode(i));
//     }
//     //h5::create_dataset(vesselgrp.open_group("nodes"), "action", tmp2);
//     //writeDataSetToGroup<DynArray<uchar>>(vesselgrp.openGroup("nodes"), string("action"), tmp2);
//     writeDataSetToGroup(h5_nodes, string("action"), tmp2);
// #endif
//   }
//   ++number;
//   return vesselgrp;