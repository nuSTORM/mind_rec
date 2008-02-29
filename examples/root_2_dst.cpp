
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
    
  if (argc<3){
    
    cout << "execution requires 2 or 3 arguments." << endl;
    cout << "Call with ./root_2_dst <InFile> <OutFile> " << endl;
    cout << "and optional  <No. events of interest>" << endl;

    return -1;
  }

  Char_t *inFileName;
  Char_t *outFileName;
  Int_t nEvents;

  inFileName = argv[1];
  outFileName = argv[2];

  //Retreive Tree from file.
  TTree *Data = NULL;

  TFile In(inFileName);
  In.GetObject("h10", Data);

  //Set number of events to be read.
  if (argc==4) nEvents = atoi(argv[3]);
  else nEvents = (Int_t)Data->GetEntries();
  
  bhep::prlevel c = bhep::NORMAL;
    
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
