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

// #include "common/hdfio.h"
// #include "hdf_wrapper_array3d.h"


// template<class LD>
// void WriteHdfLdGenericPart_(H5::Group g, const LD &ld)
// {
//   writeAttrToH5(g,"SIZEX", ld.Size()[0]);
//   writeAttrToH5(g,"SIZEY", ld.Size()[1]);
//   writeAttrToH5(g,"SIZEZ", ld.Size()[2]);
//   writeAttrToH5(g, "SIZE", ld.Size());
//   //H5::Attribute attrs = g.attrs();
// //   attrs.set("SIZEX",ld.Size()[0]);
// //   attrs.set("SIZEY",ld.Size()[1]);
// //   attrs.set("SIZEZ",ld.Size()[2]);
// //   set_array(attrs, "SIZE", ld.Size());
//   BBox3 bb = ld.Box();
//   // box is stored in a 6 component vector, xxyyzz, must match python code
//   Int6 bv;
//   //Int6 bv;
//   for (int i=0; i<3; ++i)
//   {
//     bv[i*2  ] = bb.min[i];
//     bv[i*2+1] = bb.max[i];
//   }
//   writeAttrToH5(g, string("BOX"), bv);
//   writeAttrToH5(g, string("WORLD_OFFSET"),ld.GetOriginPosition() );
//   //set_array(attrs, "BOX", bv);
//   //attrs.set("SCALE",ld.Scale());
//   //set_array(attrs, "WORLD_OFFSET", ld.GetOriginPosition());
// }
// 
// void WriteHdfLd( H5::Group g, const LatticeDataQuad3d &ld )
// {
//   //writeAttrToH5<Bool3>(g, "CENTERING", ld.GetCellCentering().cast<int>());
//   writeAttrToH5(g, "CENTERING", ld.GetCellCentering());
//   writeAttrToH5(g, "TYPE", string("QUAD3D"));
// //   h5cpp::Attributes attrs = g.attrs();
// //   attrs.set("TYPE","QUAD3D");
// //   set_array<int,3>(attrs, "CENTERING", ld.GetCellCentering().cast<int>());
//   WriteHdfLdGenericPart_(g, ld);
// }
// 
// void WriteHdfLd( H5::Group g, const LatticeDataFCC &ld )
// {
// //   h5cpp::Attributes attrs = g.attrs();
// //   attrs.set("TYPE","FCC");
//   writeAttrToH5(g, "TYPE", string("FCC"));
//   WriteHdfLdGenericPart_(g, ld);
// }
// 
// 
// 
// 
// template<class LD>
// void ReadHdfLdGenericPart_(H5::Group g, LD &ld)
// {
//   //h5cpp::Attributes attrs = g.attrs();
//   BBox3 bb; float scale;
//   try
//   {
//     // box is stored in a 6 component vector, xxyyzz
//     Int6 bv;
//     readAttrFromH5(g, string("BOX"), bv);
//     //Vec<int, 6> bv = get_array<int, 6>(attrs, "BOX");
//     for (int i=0; i<3; ++i)
//     {
//       // this must match the python code!!
//       bb.min[i] = bv[i*2  ];
//       bb.max[i] = bv[i*2+1];
//     }
//   }
//   catch (const H5::AttributeIException &e)
//   {
//     // legacy code :[
//     readAttrFromH5(g,string("SIZEX"),bb.max[0]);
//     readAttrFromH5(g,string("SIZEY"),bb.max[1]);
//     readAttrFromH5(g,string("SIZEZ"),bb.max[2]);
// //     attrs.get("SIZEX", bb.max[0]);
// //     attrs.get("SIZEY", bb.max[1]);
// //     attrs.get("SIZEZ", bb.max[2]);
//     for (int i=0; i<3; ++i)
//     {
//       bb.min[i] = 0;
//       bb.max[i] -= 1;
//     }
//   }
//   readAttrFromH5(g, string("SCALE"), scale);
//   //attrs.get("SCALE", scale);
//   ld.Init(bb, (double)scale);
//   try 
//   {
//     Float3 world_offset;
//     readAttrFromH5(g, "WORLD_OFFSET", world_offset);
//     //ld.SetOriginPosition(get_array<float,3>(attrs, "WORLD_OFFSET"));
//     ld.SetOriginPosition(world_offset);
//   } 
//   catch (const H5::AttributeIException &e) {}
// }
// 
// 
// // void ReadHdfLd( h5cpp::Group g, LatticeDataQuad3d &ld )
// // {
// //   h5cpp::Attributes attrs = g.attrs();
// //   if(attrs.get<string>("TYPE")!="QUAD3D") throw std::runtime_error("LatticeDataQuad3d from hdf5 mismatch");
// //   ReadHdfLdGenericPart_(g, ld);
// //   try {
// //     ld.SetCellCentering(get_array<int,3>(attrs, "CENTERING").cast<bool>());
// //   } catch (const h5cpp::NameLookupError &e) {}
// // }
// // 
// // 
// // void ReadHdfLd( h5cpp::Group g, LatticeDataFCC &ld )
// // {
// //   h5cpp::Attributes attrs = g.attrs();
// //   if(attrs.get<string>("TYPE")!="FCC") throw std::runtime_error("LatticeDataFCC from hdf5 mismatch");
// //   ReadHdfLdGenericPart_(g, ld);
// // }
// 
// void ReadHdfLd( H5::Group g, LatticeDataQuad3d &ld )
// {
//   //h5cpp::Attributes attrs = g.attrs();
//   
//   // WARNING HERE IS A BUG!
//   string type;
//   //readAttrFromH5(g, string("TYPE"), type);
//   if(type!="QUAD3D") throw std::runtime_error("LatticeDataQuad3d from hdf5 mismatch");
//   ReadHdfLdGenericPart_(g, ld);
//   //do this properly!!!
// //   try {
// //     ld.SetCellCentering(get_array<int,3>(attrs, "CENTERING").cast<bool>());
// //   } catch (const h5cpp::NameLookupError &e) {}
// }
// 
// 
// void ReadHdfLd( H5::Group g, LatticeDataFCC &ld )
// {
//   //h5cpp::Attributes attrs = g.attrs();
//   // __cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEvT_RKS7_RT0_
//   // __cxx1112basic_stringIcSt11char_traitsIcESaIcEEERT0_
//   // IN2H55GroupENSt7
//   // IN2H55GroupEiEvT_RKNSt7
//   // WARNING HERE IS A BUG!
//   // maybe from inherence of differnt string types e.g. H5string std::string etc...
//   string type;
//   readAttrFromH5(g, string("TYPE"), type);
//   if(type!="FCC") throw std::runtime_error("LatticeDataFCC from hdf5 mismatch");
//   ReadHdfLdGenericPart_(g, ld);
// }


