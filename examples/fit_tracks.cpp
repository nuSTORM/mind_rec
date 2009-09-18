#include <mind/fitter.h>
#include <mind/MINDplotter.h>
//#include <bhep/EventManager2.h>
#include <bhep/gstore.h>
#include <bhep/sreader.h>
#include <bhep/reader_root.h>

#include <recpack/Measurement.h>
#include <recpack/Ring.h>
#include <recpack/stc_tools.h>

using namespace std;

/*


 */

int main(int argc, char* argv[]){
    
  //---- set parameters file ----//
    
  string param_file="param/fit_cats.param";
  
  int nevents=1;
  bool fitOk, catchOk;
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
  
  //EventManager2* eman = new EventManager2(data_store,bhep::MUTE);
  bhep::reader_root inDst;
  inDst.open( data_store.fetch_sstore("idst_file") );

  fitter* fit = new fitter(ana_store,bhep::MUTE);
  
  MINDplotter* plot = new MINDplotter();
  
  //catchOk = eman->initialize();
  
  catchOk = fit->initialize();//eman->get_dst_properties());

  catchOk = plot->initialize(run_store.fetch_sstore("out_file"),bhep::MUTE);
  
  //add run properties to output dst header
  
  //will something be missing without these steps??
  //catchOk = eman->add_run_property("MINDfit","1");
  
  //catchOk = eman->add_run_properties(run_store);
  
  bool patR = ana_store.fetch_istore("patRec");
  
  for(int i=0; i < nevents; i++) {
    
    //bool ok = eman->status();
    
    //if (!ok) break; // checks eof
    
    if (i%100==0) cout<< "Number of events read "<<i<<endl;
    
    //bhep::event& e = eman->read(); 
    bhep::event& e = inDst.read_event( i );
    
    // loop over particles
    vector<bhep::particle*> parts = e.digi_particles(); 
    
    cout <<"There are " << parts.size() << " digis in event " << e.event_number() <<endl;
    if (parts.size() != 0) {
      for (size_t part=0; part<parts.size();part++){
	
	if (parts[part]->name()=="void") continue;
	
	fitOk = fit->execute(*parts[part],e.event_number());
	
	catchOk = plot->execute(*fit, e, fitOk, patR);
      }
    }
    
    //save event containing fit info 
    //eman->write(e);
    
    parts.clear();
    e.clear();
  }
  
  catchOk = fit->finalize();
  
  //catchOk = eman->finalize();
  inDst.close();

  catchOk = plot->finalize();
  
  return 0;
  
}
