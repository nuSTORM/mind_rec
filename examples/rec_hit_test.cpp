//Test for voxel creation.

#include <mind/hit_clusterer.h>
#include <mind/fitter.h>

#include <bhep/gstore.h>
#include <bhep/hit.h>
#include <bhep/sreader.h>
#include <bhep/reader_root.h>

#include <recpack/Measurement.h>
#include <recpack/stc_tools.h>

#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>

int main(){

  int nhit1, nhit2, evt;
  double x1[1500], y1[1500], z1[1500];
  double x2[1500], y2[1500], z2[1500];
  TFile *f1 = new TFile("test2.root","recreate");
  TTree *tree = new TTree("h1","stuff");
  tree->Branch("Event", &evt,"evtno/I");
  tree->Branch("bhepHit",&nhit1,"nhit1/I");
  tree->Branch("recHit",&nhit2,"nhit2/I");
  tree->Branch("bhepX",&x1,"x1[nhit1]/D");
  tree->Branch("bhepY",&y1,"y1[nhit1]/D");
  tree->Branch("bhepZ",&z1,"z1[nhit1]/D");
  tree->Branch("recX",&x2,"x2[nhit2]/D");
  tree->Branch("recY",&y2,"y2[nhit2]/D");
  tree->Branch("recZ",&z2,"z2[nhit2]/D");
  int nevents=100;

  string param_file = "examples/param/fit_cats.param";

  bhep::gstore store;
  bhep::sreader reader(store);
  reader.file(param_file);
  reader.group("RUN");
  reader.read();

  hit_clusterer* hct = new hit_clusterer( store );

  bhep::reader_root inDst;

  inDst.open("../digi_out/muCC_test20k_digi.dst.root");

  //measurement_vector meas;
  std::vector<cluster*> meas;
  size_t countification;
  for (int i=0;i<nevents;i++){
    
    bhep::event& e = inDst.read_event( i );
    
    std::vector<bhep::particle*> parts = e.digi_particles();
    
    for (size_t j=0;j < parts.size();j++){
      countification = 0;  
      std::vector<bhep::hit*> hits = parts[j]->hits("tracking");
      
      nhit1 = (int)hits.size();
      std::cout << "Voxels before: " << nhit1 << std::endl;
      for (int k=0;k<nhit1;k++){
	x1[k] = hits[k]->x()[0];
	y1[k] = hits[k]->x()[1];
	z1[k] = hits[k]->x()[2];
      }
      hct->execute(hits, meas);
      std::cout << "Before: " << meas[0]->position()[2] <<", " << meas[meas.size()-1]->position()[2] << ", " << meas[0]->get_mu_prop() << std::endl;
      sort( meas.begin(), meas.end(), reverseSorter() );
      std::cout << "After: " << meas[0]->position()[2] <<", " << meas[meas.size()-1]->position()[2]<< ", " << meas[0]->get_mu_prop() << std::endl;
      sort( meas.begin(), meas.end(), forwardSorter() );
      std::cout << "After2: " << meas[0]->position()[2] <<", " << meas[meas.size()-1]->position()[2]<< ", " << meas[0]->get_mu_prop() << std::endl;
      nhit2 = (int)meas.size();
      for (int l=0;l<nhit2;l++){
	x2[l] = meas[l]->position()[0];
	y2[l] = meas[l]->position()[1];
	z2[l] = meas[l]->position()[2];
      }
      for (int ii=0;ii<nhit2;ii++){
	countification += meas[ii]->get_nVox();
	std::cout << "MuProp meas["<<ii<<"] = " << meas[ii]->get_mu_prop() << std::endl;
      }
      std::cout << "Total Vox in clusts: " << countification << std::endl;
      evt = i;
      tree->Fill();
    }
    
    stc_tools::destroy(meas);
    parts.clear();
    e.clear();

  }

  inDst.close();
  
  f1->Write();
  f1->Close();

  delete hct;

  return 0;
}
