
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/format.hpp>

#include "../src/common/vessels3d.h"
#include "../src/common/shared-objects.h"
#include "hdf_wrapper.h"
#include <BloodVessel.h>

namespace po = boost::program_options;
using namespace std;


void wrap_vessels_witho2(VesselList3d &vl, std::vector<BloodVessel> &bloodVesselVector, DynArray<float> &a)
{
  int ecnt = vl.GetECount();
  
  //create new entry
  BloodVessel suggestion;
  Float3 buffer;
  vector<double> bufferToFill;
  
  for(int i = 0; i<ecnt; i++)
  {
    const Vessel* v= vl.GetEdge(i);
    suggestion.SetBloodVesselR(v->r);
    //we use the Eigen3 library to store array, this is faster
    //pos a
    buffer = vl.Ld().LatticeToWorld(v->LPosA());
    bufferToFill = {buffer[0], buffer[1], buffer[2]};
    suggestion.SetBloodVessela(bufferToFill);
    //pos b
    buffer = vl.Ld().LatticeToWorld(v->LPosB());
    bufferToFill = {buffer[0], buffer[1], buffer[2]};
    suggestion.SetBloodVesselb(bufferToFill);
    
    /********* O2 ***********/
    suggestion.SetBloodVesselO2start(a[2*i+0]);
    suggestion.SetBloodVesselO2end(a[2*i+1]);
    
    bloodVesselVector.push_back(suggestion);
  }
  
}
void wrap_vessels(VesselList3d &vl, std::vector<BloodVessel> &bloodVesselVector)
{
  int ecnt = vl.GetECount();
  
  //create new entry
  BloodVessel suggestion;
  Float3 buffer;
  vector<double> bufferToFill;
  
  for(int i = 0; i<ecnt; i++)
  {
    const Vessel* v= vl.GetEdge(i);
    suggestion.SetBloodVesselR(v->r);
    //we use the Eigen3 library to store array, this is faster
    //pos a
    buffer = vl.Ld().LatticeToWorld(v->LPosA());
    bufferToFill = {buffer[0], buffer[1], buffer[2]};
    suggestion.SetBloodVessela(bufferToFill);
    //pos b
    buffer = vl.Ld().LatticeToWorld(v->LPosB());
    bufferToFill = {buffer[0], buffer[1], buffer[2]};
    suggestion.SetBloodVesselb(bufferToFill);
    
    
    bloodVesselVector.push_back(suggestion);
  }
  
}

int main(int argc, char* argv[])
{
  std::cout << "creating blood vessel structure" << std::endl;
  
  std::string fn, fno2;
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("h5", po::value<std::string>(), "set input vessel filename")
    ("o2", po::value<std::string>(), "set input filename for o2 simulation")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
      cout << desc << "\n";
      return 1;
  }

  h5cpp::File file;
  bool withO2 = false;
  
  if (vm.count("h5")) 
  {
    fn = vm["h5"].as<string>();
    cout << "input vessel filename was set to: "  << fn << ".\n";
    file = h5cpp::File(fn, "r");
  } 
  else 
  {
    cout << "vessel file not set.\n";
  }
  
  if (vm.count("o2")) 
  {
    fno2 = vm["o2"].as<string>();
    cout << "o2 filename was set to: "  << fn << ".\n";
    file = h5cpp::File(fno2, "r");
    withO2 = true;
  } 
  else 
  {
    cout << "vessel file not set.\n";
  }

  std::auto_ptr<VesselList3d> vl;
  DynArray<float> a;
  std::vector<BloodVessel> milottiVesselList;
  boost::property_tree::ptree pt;
  if( !withO2 )
  {
    vl = ReadVesselList3d(file.root().open_group("vessels"), pt);
    wrap_vessels(*vl, milottiVesselList);
  }
  else
  {
    vl = ReadVesselList3d(file.root().open_group("recomputed_flow/vessels"), pt);
    const int ecnt = vl->GetECount();
    h5cpp::Group theO2Group = file.root().open_group("po2/vessels");
    h5cpp::Dataset thePO2ofVessels = theO2Group.open_dataset("po2vessels");
    
    h5cpp::Dataspace sp = thePO2ofVessels.get_dataspace();
    hsize_t dims[H5S_MAX_RANK];
    int rank = sp.get_dims(dims);
    h5cpp::read_dataset<float>(thePO2ofVessels,a);
    
    wrap_vessels_witho2(*vl, milottiVesselList, a);
  }
  
  for( auto vessel: milottiVesselList)
  {
    cout<< boost::format("r: %f \n") % vessel.GetBloodVesselR();
    cout<< "position a:\n" << endl;
    cout<< "x\t y\t z\n" << endl;
    for( auto entry: vessel.GetBloodVessela() )
    {
      cout<< boost::format("%f\t ") % entry; 
    }
    cout<< "o2 start \n" << endl;
    cout<< boost::format("%f\t ") % vessel.GetBloodVesselO2start();
    
    cout<< "o2 end \n" << endl;
    cout<< boost::format("%f\t ") % vessel.GetBloodVesselO2end(); 
    

  }
}
