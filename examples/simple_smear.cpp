/******************************************************************
 * Program to create a digital particle for fitting from the hits *
 * in the true particles in the bhep output of GEANT4 MIND.       *
 *                                                                *
 * Temporary measure until full digitization is available.        *
 *                                                                *
 * Execute as:                                                    *
 *       ./simple_smear <param_file_name>                         *
 *                                                                *
 * Author: Andrew Laing.                                          *
 ******************************************************************/

#include <mind/root2dst.h>

#include <bhep/EventManager2.h>
#include <bhep/gstore.h>
#include <bhep/sreader.h>

using namespace std;

int main(int argc, char* argv[]) {

  string param_file;
  long rndmSeed;

  if (argc == 2) param_file = argv[1];
  else {

    std::cout << "Execute as: ./simple_smear <paramater file>" << std::endl;

    return -1;
  }

  //Stores for information from parameter file.
  bhep::gstore run_store, data_store;

  //Read run parameters.
  bhep::sreader reader1(run_store);
  reader1.file(param_file);
  reader1.group("RUN");
  reader1.read();
  //
  //Read files to be processed.
  bhep::sreader reader2(data_store);
  reader2.file(param_file);
  reader2.group("DATA");
  reader2.read();
  //

  int nEvents;
  double smearRes[2];

  if ( !run_store.find_istore("nEvents") ) {
    std::cout << "Parameter file must contain number of events "
	      << "to be processed in group RUN" << std::endl;
    return -1;
  }
  else nEvents = run_store.fetch_istore("nEvents");

  if ( !run_store.find_dstore("Gaus_Sigma") ) {
    std::cout << "Parameter file must contain smear sigma in cm "
	      << "as double Gaus_Sigma in group RUN" << std::endl;
    return -1;
  }
  else smearRes[0] = run_store.fetch_dstore("Gaus_Sigma");

  if ( !run_store.find_dstore("Eng_Res") ) {
    std::cout << "Parameter file must contain Energy resolution"
	      << "as double Eng_Res in group RUN" << std::endl;
    return -1;
  }
  else smearRes[1] = run_store.fetch_dstore("Eng_Res");

  if ( !run_store.find_dstore("Gen_seed") ) {
    std::cout << "Parameter file must contain generator seed "
	      << "as double Gen_seed in group RUN" << std::endl;
    return -1;
  }
  else rndmSeed = (long)run_store.fetch_dstore("Gen_seed");

  EventManager2* eman = new EventManager2(data_store, bhep::NORMAL);

  root2dst* cvt = new root2dst(bhep::NORMAL);

  eman->initialize();

  cvt->initialize( smearRes, rndmSeed);

  for (int iEvent = 0;iEvent < nEvents;iEvent++) {

    bool ok = eman->status();

    if (!ok) break;

    bhep::event& e = eman->read();

    //Get Vector of all true particles from event.
    vector<bhep::particle*> particles = e.true_particles();
    //

    particle* digi_part = cvt->create_digital_representation( particles );

    e.add_digi_particle( digi_part );

    eman->write( e );

    delete digi_part;

  }

  eman->finalize();

  return 0;
}
