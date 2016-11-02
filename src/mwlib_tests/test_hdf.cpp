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

#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <list>

#include "hdf_wrapper.h"

using namespace std;
namespace h5 = h5cpp;

void WriteFile()
{
  cout << "=== Writing ===" << endl;
  cout << "-- type caching and ref. counting sanity.--" << endl;
  {
    h5::Datatype dt_int1 = h5::get_disktype<int>();
    h5::Datatype dt_int2 = h5::get_disktype<int>(); // refer to the same type object
    cout << "dt_int1 id = " << dt_int1.get_id() << " vs. dt_int2 id = " << dt_int2.get_id() << endl;
    assert(dt_int1.get_id() == dt_int2.get_id());
  }
  {
    h5::Datatype dt = h5::get_disktype<int>();
    cout << "dt ref count = " << dt.get_ref() << endl;
    assert (dt.get_ref() == 2);  // one for cached instance, and another reference accounting for the dt variable here.
  }
  
  cout << "-- group management and attribute R/W tests --" << endl;
  
  /* like in h5py:
    w = create or truncate existing file
    a = append to file or create new file
    r = read only; file must exist
    w- = new file; file must not already exist
    r+ = read/write; file must exist
  */
  h5::File file("test.h5", "w");
  h5::Group root = file.root();  // because File does not have group functionality here. This is reserved for future work.
  
  // Attributes class inspired by h5py
  h5::Attributes a = root.attrs();
  a.create<char>("a_char",'q'); // create scalar attribute
  a.set<double>("set double test", 2.); // replace attribute with new data. Can be of different type, because old attribute will just be deleted if neccessary.
  a.create<float>("to_be_deleted", 9000.0);
  a.create<int>("attrib_overwrite_test1", 1);
  a.create<int>("attrib_overwrite_test2", 2);
  a.create<string>("a_string_attr","teststring");
  a.create("c_str_attr", "it's a c string");
  try
  {
      a.create<float>("to_be_deleted", 9002.0); 
  }
  catch (h5::Exception &e)
  {
    cout << "Attempt to create two attributes under the same name correctly failed with error message:" << endl;
    cout << e.what() << endl;
  }
  a.remove("to_be_deleted");
  a.create("to_be_deleted", "really gone? (it is okay if you see this)");
  a.set<string>("attrib_overwrite_test1", "erased and recreated!");
  a.set<int>("attrib_overwrite_test2", 3); // should just write into the existing space
  
  // try to make a group and set some attributes
  h5::Group g = root.create_group("testing_the_group");
  assert(g.is_valid());
  g.attrs().create<double>("a1", 5.);
  g.attrs().create<int>("second_attribute", 1);

  bool group_already_there;
  h5::Group g2 = root.require_group("testing_the_group", &group_already_there);
  assert(group_already_there && g2.is_valid());
  
  {
    h5::Group g3 = root.require_group("second_group");
    assert(g3.is_valid());
  }

  h5::Group g3 = root.open_group("second_group").create_group("third_group");

  // dataspaces
  h5::Dataspace sp_scalar = h5::Dataspace::scalar();
  assert(sp_scalar.is_valid());
  h5::Dataspace sp1       = h5::Dataspace::simple_dims(10,20,30); // three dimensional dataset with dimensions 10 x 20 x 30
  assert(sp1.is_valid());
  hsize_t dims[3] = { 10, 20, 30};
  h5::Dataspace sp2       = h5::Dataspace::simple(3, dims);
  assert(sp2.is_valid());
  #if 0
  vector<int> dims_vec; 
  dims_vec.push_back(10);
  dims_vec.push_back(20);
  dims_vec.push_back(30);
  // following two create function are found in hdf_wrapper_stl.h
  h5::Dataspace sp5       = h5::Dataspace::(dims_vec); // calls dims_vec.begin(), dims_vec.end() 
  assert(sp5.is_valid());
  h5::Dataspace sp6       = h5::create_dataspace_from_iter((hsize_t*)dims, dims+3); // good old pointer arithmetic
  assert(sp6.is_valid());
  // TODO: check if these dataspaces are okay
#endif
  // array attributes
  int static_ints[6] = { 1, 2, 3, 4, 5, 6 };
  g3.attrs().create("ints", h5::Dataspace::simple_dims(2,3), static_ints);

  // variable length string attributes
  const char* strings1[] = {
    "string1",
    "string2 long",
    "string3 very long"
  };
  g3.attrs().create("strings1", h5::Dataspace::simple_dims(3), strings1); 
  
  // using the wrapper for stl constainers.
  vector<string> strings2;
  strings2.push_back("test");
  strings2.push_back("string");
  strings2.push_back("array attrib");
  h5::set_array(g3.attrs(), "strings2", strings2); // should use the overloaded function for std::vector
  
  cout << "-- datasets --" << endl;
  
  vector<float> data(10);
  for(int i=0; i<data.size(); ++i)
    data[i] = (float)i*i;
  
  h5::Dataset ds = h5::create_dataset(g, "testds", h5::Dataspace::simple_dims(data.size()), &data[0]);
  ds.attrs().create("someAttr", 7);
  
  // make nice data in memory
  vector<double> bigdata(dims[0]*dims[1]*dims[2]);
  int i=0;
  for (int z=0; z<dims[2]; ++z)  for (int y=0; y<dims[1]; ++y) for (int x=0; x<dims[0]; ++x)
  {
    bigdata[i++] = cos(x * M_PI/10.) * cos(y * M_PI/10.) * cos(z * M_PI/10.);
  }
  
  // most general data writing
  h5::Dataset ds2 = h5::Dataset::create(g, "bigdata", 
					h5::get_disktype<double>(), 
                                        sp2,
                                        h5::Dataset::create_creation_properties(sp2, h5::CREATE_DS_COMPRESSED));
  ds2.write(sp2, sp2, &bigdata[0]);

  
  h5::create_dataset(root, "dataset_from_range_vector", data); // automatic determination if contiguous, according to h5cpp::internal::contiguous_storage_traits class

  std::list<int> listdata;
  for (int i = 0; i < 10; ++i)
    listdata.push_back(i*i*i);
  h5::create_dataset(root, "dataset_from_list", std::vector<int>(listdata.begin(), listdata.end()));

  cout << "-- hyperslab selection --" << endl;
  { // selection test
    h5::Dataspace sp = h5::Dataspace::simple_dims(5, 2);
    h5::Dataset ds = h5::Dataset::create<int>(root, "select_test_ds", sp);
    hsize_t offset[2] = { 0, 0 };
    hsize_t stride[2] = { 1, 1 };
    hsize_t block[2] = { 1, 1 };
    hsize_t count[2] = { 4, 1 };
    sp.select_hyperslab(offset, stride, count, block); // writes to [0 .. 4] x [0 .. 0]
    
    int memdata[] = { 1, 2, 3, 4 };

    h5::Dataspace memsp = h5::Dataspace::simple_dims(4);
    ds.write(memsp, sp, &memdata[0]);

    offset[0] = 1;
    offset[1] = 1;
    sp.select_hyperslab(offset, stride, count, block);
    ds.write(memsp, sp, &memdata[0]);  // writes to [1 .. 5] x [1 .. 1]
  }
}

#if 1
void ReadFile()
{
  cout << "=== Reading ===" << endl;
  
  h5::File file("test.h5", "r");
  h5::Group root = file.root();
  h5::Dataset ds;

  cout << "-- group management and attribute R/W tests --" << endl;

  assert(root.is_valid());
  cout << root.get_name() << " has " << root.size() << " children." << endl;
  cout << "children of " << root.get_name() << endl;
  for (int i = 0; i < root.size(); ++i)
  {
    cout << "\t" << root.get_link_name(i) << endl;
  }
  cout << "and again ..." << endl;

  const auto end = root.end(); // takes a little bit of time
  for (auto it = root.begin(); it != end; ++it)
  {
    cout << "\t" << *it << endl;
  }
 
  string s = root.attrs().get<string>("c_str_attr");
  cout << "string attribute: " << s << endl;
  
  h5::Attributes attrs = root.open_group("second_group/third_group").attrs();
  h5::Attribute a = attrs.open("ints");
  assert(a.get_dataspace().get_npoints() == 6);
  int static_ints[6];
  a.read(static_ints);
  cout << "static ints";
  for (int i=0; i<6; ++i)
    cout << " " << static_ints[i];
  cout << endl;
  
  // read n-dimensional attribute into vector
  vector<string> string_vec; h5::get_array(attrs, "strings1", string_vec);

  cout << "string_vec (" << string_vec.size() << ") ";
  for (int i=0; i<string_vec.size(); ++i)
    cout << " " << string_vec[i];
  cout << endl;

  cout << "-- datasets --" << endl;

  {
    ds = file.root().open_group("testing_the_group").open_dataset("testds");

    // check reading the dataspace dimensions first
    h5::Dataspace sp = ds.get_dataspace();
    hsize_t dims[H5S_MAX_RANK];
    int rank = sp.get_dims(dims);
    cout << "dataset rank = " << rank << ", dims = ";
    for (int i = 0; i < rank; ++i)
    {
      cout << dims[i] << " ";
    }
    cout << endl;

    { // check reading to memory block pointed to by a c-style pointer
      vector<float> data(sp.get_npoints());
      ds.read<float>(&data[0]);

      cout << "dataset = ";
      for (auto it = data.begin(); it != data.end(); ++it)
      {
        cout << *it << ", ";
      }
      cout << endl;
    }

    { // check reading to vector
      ds = root.open_dataset("dataset_from_list");
      vector<int> data; h5::read_dataset<int>(ds, data);
    }
  }

  cout << "-- hyperslab selection --" << endl;
  /* for reference: "Dataspace dimensions are numbered from 1 to rank. 
  HDF5 uses C storage conventions, assuming that the last listed dimension 
  is the fastest-changing dimension and the first-listed dimension is the 
  slowest chang­ing."
  */
  int expected_value[5][2] = {
    { 1, 0 }, { 2, 1 }, { 3, 2 }, { 4, 3 }, { 0, 4 }
  };
  ds = file.root().open_dataset("select_test_ds");
  hsize_t dims[H5S_MAX_RANK];
  ds.get_dataspace().get_dims(dims);
  cout << "dims = (" << dims[0] << ", " << dims[1] << ")" << endl;
  vector<int> data; h5::read_dataset<int>(ds, data);
  const int NX = 5, NY = 2;
  for (int y = 0; y < NY; ++y)
  for (int x = 0; x < NX; ++x)
    {
      int value = data[x*NY + y];
      printf("(%i, %i) = %i\n", x, y, value);
      assert(value == expected_value[x][y]);
    }
}
#endif


int main(int argc, char **argv)
{
  h5::disableAutoErrorReporting();  // don't print to stderr
  WriteFile();
  ReadFile();
  cin.get();
  return 0;
}
