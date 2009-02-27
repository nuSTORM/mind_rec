#include <mind/fitter.h>
#include <mind/MINDplotter.h>
#include <bhep/EventManager2.h>
#include <bhep/gstore.h>
#include <bhep/sreader.h>
#include <recpack/Measurement.h>
#include <recpack/Ring.h>

using namespace std;

/*


 */

int main(int argc, char* argv[]){
    
  //---- set parameters file ----//
    
  string param_file="param/fit_cats.param";
  
  int nevents=1;
  bool fitOk;
  if (argc==2) param_file = argv[1];
  
  else if(argc==3) {param_file = argv[1];nevents = atoi(argv[2]);}
  
  else{cout<<
	 "\n++ Usage: ./fit_tracks param_file number_of_events\n"<<endl;
  exit(1);}
  
  // generate stores for analysis parameters 
  bhep::gstore data_store,ana_store,run_store;
  
  // read analysis parameters 
  bhep::sreader reader(ana_store);
  reader.file(param_file); 
  reader.group("ANA");reader.read();
  bhep::sreader reader2(ana_store);
  reader2.file(param_file); 
  reader2.group("RUN");reader2.read();
  //
  
  // read run properties
  bhep::sreader preader(run_store);
  preader.file(param_file);
  preader.group("RUN");
  preader.read();
  //
  
  // read input/ouput data
  bhep::sreader data_reader(data_store);
  data_reader.file(param_file);
  data_reader.group("DATA");
  data_reader.read();
  //
  
  EventManager2* eman = new EventManager2(data_store,bhep::MUTE);
  
  fitter* fit = new fitter(ana_store,bhep::MUTE);
  
  MINDplotter* plot = new MINDplotter();
  
  eman->initialize();
  
  fit->initialize(eman->get_dst_properties());

  
  plot->initialize("/home/alaing/ntuples/MINDpulls.root",bhep::MUTE);
  
  //add run properties to output dst header
  
  eman->add_run_property("MINDfit","1");
  
  eman->add_run_properties(run_store);
  
  //Vectors and Histograms required for pull calculation.
  EVector v;
  EMatrix M;  

  bool patR = ana_store.fetch_istore("patRec");
  
  for(int i=0; i < nevents; i++) {
    
    bool ok = eman->status();
    
    if (!ok) break; // checks eof
    
    if (i%100==0) cout<< "Number of events read "<<i<<endl;
    
    bhep::event& e = eman->read(); 
    
    // loop over particles
    
    vector<bhep::particle*> parts = e.digi_particles(); 
    cout <<"There are " << parts.size() << " digis in event " << e.event_number() <<endl;
    if (parts.size() != 0) {
      for (size_t part=0; part<parts.size();part++){
	
	bhep::particle& p = *parts[part];
	
	if (p.name()=="void") continue;
	
	fitOk = fit->execute(p,e.event_number());
	
	plot->execute(*fit, e, fitOk, patR);
      }
    }
      //save event containing fit info 
      eman->write(e);
  }
  
  fit->finalize();
  
  eman->finalize();

  plot->finalize();

  //outPut->Write();
  
  return 0;
  
}
