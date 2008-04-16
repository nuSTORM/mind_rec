/* -*- mode: c++ -*- */
#ifndef _root_dst___
#define _root_dst___

#include <bhep/event.h>
#include <bhep/particle.h>
#include <bhep/messenger.h>
#include <bhep/bprint.h>
#include <bhep/bhep_svc.h>
#include <bhep/ray.h>
#include <bhep/hit.h>

#include <TTree.h>
#include <Riostream.h>
#include <sstream>

using namespace bhep;

//! root2dst Class
/*!
  Converts ROOT tree into bhep DST file
*/

#define MAXEVENTS 50000

class root2dst{
  
public:
  
  root2dst(bhep::prlevel);
  
  ~root2dst(){};
  
  bool initialize(TTree *InPutTree, TString OutFileName);
  bool execute();
  bool finalize();

  //Function to get names of Branches of interest
  void fillBranchVector();

  //Functions to create bHEP object from the ntuple
  void createEvent();
  void make_particles();
  particle* define_lead_particle();
  vector<hit*> append_hits();
  particle* define_hadron();
  particle* hits_to_fit();

  //Handy int/float to string converters.
  TString ToString(Int_t num){
    ostringstream start;
    start<<num;
    TString start1=start.str();
    return start1;
    
  }
  
 TString ToString(Float_t num){
    ostringstream start;
    start<<num;
    TString start1=start.str();
    return start1;
    
  }
  
protected:
  
  //verbosity level
  bhep::prlevel level;
  
  //messenger
  bhep::messenger m;
  
  //event counter
  size_t nevt;
  
private:

  //The Event
  event *nuEvent[MAXEVENTS];

  //file to write dst to.
  writer_gz outgz;

  //Root tree to be read.
  TTree *dataIn;

  //Vector of Branch Names.
  vector<TString> branches;

};





#endif
