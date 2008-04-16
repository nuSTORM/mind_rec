
#include <mind/root2dst.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;


/*****************************************************************
 * Root tree to bHEP dst converter. Takes in a root file or list *
 * of root files reads the information in the trees and creates  *
 * a gz file with the relevant events as objects.                *
 *                                                               *
 * Execution via calls of the form:                              *
 *  ./root_2_dst <Input Filename> <Output Filename> <No. events> *
 *                                                               *
 * Authors: Andrew Laing, Pau Novella.                           *
 *****************************************************************/

int main(int argc, char* argv[]){
    
  if (argc<2){
    
    cout << "Execution requires 1 or 2 arguments." << endl;
    cout << "Call with ./root_2_dst <InFile> " << endl;
    cout << "and optional  <No. events of interest>" << endl;

    return -1;
  }


  Char_t *inFileName;
  Int_t nEvents;

  inFileName = argv[1];

  TString outFileName = inFileName;
  outFileName.Replace(outFileName.Length()-4, 4, "gz");

  //Retreive Tree from file.
  TTree *Data = NULL;

  TFile In(inFileName);
  In.GetObject("h10", Data);

  //Set number of events to be read.
  if (argc==3) nEvents = atoi(argv[2]);
  else nEvents = (Int_t)Data->GetEntries();

  if (nEvents>(Int_t)Data->GetEntries()){
    cout << "Input no. events greater than available in file.\n"
	 << "Replacing with max entries" << endl;
    nEvents = (Int_t)Data->GetEntries();
  }
  
  bhep::prlevel c = bhep::VERBOSE;
    
  root2dst* cvt = new root2dst(c);
    
  cvt->initialize(Data, outFileName);
  
  cout << "Starting Loop over events" << endl;
  for(int i=1; i < nEvents+1; i++) {
    
    if (i%100==0) cout<< "Number of events read "<<i<<endl;
    
    cvt->execute();
    
  }
  
  cvt->finalize();
  
  return 0;
  
}
