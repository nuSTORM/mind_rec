
#include <root2dst.h>

using namespace std;
using namespace bhep;


//*************************************************************
root2dst::root2dst(bhep::prlevel vlevel){
//*************************************************************
  
  level = vlevel;
  
  m = bhep::messenger(level);

  dataIn = NULL;
  
  m.message("+++ root2dst Constructor +++",bhep::VERBOSE); 
  
}


//*************************************************************
bool root2dst::initialize(TTree *InPutTree, TString OutFileName) {
//*************************************************************
    
  m.message("+++ root2dst init  function ++++",bhep::NORMAL);
  
  outgz.open((string)OutFileName);

  dataIn = InPutTree;
  dataIn->SetBranchStatus("*",0);

  fillBranchVector();

  nevt=0;
  
  return true;
}


//*************************************************************
bool root2dst::execute(){
//*************************************************************
    
  /*
    Take a ROOT tree event and convert it into bhep DST event
   */

  nuEvent[nevt] = NULL;
  
  m.message("+++ root2dst execute function ++++",bhep::VERBOSE);

  createEvent();

  make_particles();
  
  outgz.write(*nuEvent[nevt], nevt);

  delete nuEvent[nevt];
  
  nevt++;
 
  return true;

}



//*************************************************************
bool root2dst::finalize() {
//*************************************************************
  
  outgz.close();
  delete dataIn;

  m.message("+++ root2dst finalize function ++++",bhep::NORMAL);
    
  m.message("++ Number of analyzed events: ",nevt,bhep::NORMAL);

  return true;
}


//*************************************************************
void root2dst::fillBranchVector() {
//*************************************************************

/*This function takes in an existing file which is a list of 
  branch names and enters these names in the branches vector.
  If a line starts with # only the # will be entered as a 
  marker for the event creator function */


  TString value;
  Int_t count = 0;

  ifstream branchIn;
  branchIn.open("/home/alaing/mind/examples/Branches.txt");

  while (branchIn.good()){

    branchIn >> value;

    if (value[0] == 35) branches.push_back(value[0]);
    else branches.push_back(value);
    
    count++;

  }

  branchIn.close();

}

//*************************************************************
void root2dst::createEvent() {
//*************************************************************

/* Defines an event and adds the event specific information. */

  Int_t integProps[4][2];
  Float_t floatProps[6][3];

  nuEvent[nevt] = new event((Int_t)nevt);

  for (Int_t iBran = 0;iBran<10;iBran++){

    if (branches[iBran] != "#"){
      dataIn->SetBranchStatus(branches[iBran], 1);
      if (iBran < 4) dataIn->SetBranchAddress(branches[iBran], &integProps[iBran]);
      else dataIn->SetBranchAddress(branches[iBran], &floatProps[iBran-4]);
    }

  }

  dataIn->GetEntry((Int_t)nevt);

  nuEvent[nevt]->set_vertex(floatProps[5][0] * cm, floatProps[5][1] * cm, floatProps[5][2] * cm);

  //Add aditional properties
  for (Int_t iProp = 0;iProp<9;iProp++){

    if (branches[iProp] != "#"){
      if (iProp<2 || iProp==3) nuEvent[nevt]->add_property((string)branches[iProp], (string)ToString(integProps[iProp][0]));
      else if (iProp==2) {
	TString n1, n2;
	n1 = branches[iProp]+"Seed1";
	n2 = branches[iProp]+"Seed2";

	nuEvent[nevt]->add_property((string)n1, (string)ToString(integProps[iProp][0]));
	nuEvent[nevt]->add_property((string)n2, (string)ToString(integProps[iProp][1]));
      }
      else nuEvent[nevt]->add_property((string)branches[iProp], (string)ToString(floatProps[iProp-4][0]));
    }

  }

  dataIn->SetBranchStatus("*", 0);

  cout <<"Event Defined" << endl;
}

//***************************************************************
void root2dst::make_particles() {
//***************************************************************

/* Function to add relevant particles to event */
  particle *par1, *par2, *digPar;

  par1 = define_lead_particle();

  vector<hit*> par1Hits = append_hits();

  if (par1Hits.size() != 0)
    for (Int_t iHit=0;iHit<(Int_t)par1Hits.size();iHit++){
      par1Hits[iHit]->set_mother_particle(*par1);
      par1->add_hit("MIND", par1Hits[iHit]);
    }

  digPar = hits_to_fit();
  
  par2 = define_hadron();
  
  nuEvent[nevt]->add_true_particle(par1);
  nuEvent[nevt]->add_true_particle(par2);
  nuEvent[nevt]->add_digi_particle(digPar);

  cout << "All particles defined and appended to event" << endl;

}

//***************************************************************
particle* root2dst::define_lead_particle() {
//***************************************************************

/* Function to define the properties of the lead particle from
   the event. e.g. the mu- from a numucc interaction */
  
  Int_t ID;
  Float_t moment[5];
  
  ptype pT = TRUTH;
  string pName;
  
  for (Int_t iBran = 10;iBran < 12;iBran++){
    dataIn->SetBranchStatus(branches[iBran], 1);
    if (iBran==10) dataIn->SetBranchAddress(branches[iBran], &ID);
    else dataIn->SetBranchAddress(branches[iBran], &moment);
  }
  
  dataIn->GetEntry((Int_t)nevt);
  
  Point3D pos = nuEvent[nevt]->vertex();
  //multiply by 1000 because ntuple in GeV and DST expects MeV.
  Vector3D mom(moment[2] * GeV, moment[1] * GeV, moment[0] * GeV);
  
  ray *leadstart = new ray(pos, mom);
  
  if (ID == 13) pName = "mu-";
  else pName = "mu+";
  
  particle *lead = new particle(pT, pName, *leadstart);

  dataIn->SetBranchStatus("*", 0);

  delete leadstart;
  
  return lead;
  
}

//**************************************************************
particle* root2dst::define_hadron() {
//**************************************************************

/* Add particle to represent the hadronic shower associated 
   with the event */
  cout << "Defining Hadron" << endl;
  ptype pT = TRUTH;
  string pName = "Hadronic_vector";

  Int_t nHad;
  Float_t moment[5];

  for (Int_t iBran = 38;iBran < 40;iBran++){
    dataIn->SetBranchStatus(branches[iBran], 1);
    if (iBran==38) dataIn->SetBranchAddress(branches[iBran], &nHad);
    else dataIn->SetBranchAddress(branches[iBran], &moment);
  }

  dataIn->GetEntry((Int_t)nevt);
  
  Point3D pos = nuEvent[nevt]->vertex();
  
  Vector3D mom(moment[0] * GeV, moment[1] * GeV, moment[2] * GeV);
  
  ray *hadstart = new ray(pos, mom);
  
  particle *had = new particle(pT, pName, *hadstart);
  cout << "Particle Defined" << endl;
  had->add_property("No_Hadrons_in_Shower", (string)ToString(nHad));

  dataIn->SetBranchStatus("*", 0);

  delete hadstart;

  cout << "Hadronic shower  vector defined" << endl;

  return had;
}


//**************************************************************
vector<hit*> root2dst::append_hits() {
//**************************************************************

/* Add vector of hit locations associated with particle par */

  vector<hit*> track;
  const Int_t maxHits = 400;

  Int_t Nhits;
  Float_t hitVec[3][maxHits];

  for (Int_t iBran = 0;iBran < 4;iBran++){
    
    dataIn->SetBranchStatus(branches[iBran+46], 1);
    if (iBran==0) dataIn->SetBranchAddress(branches[iBran+46], &Nhits);
    else dataIn->SetBranchAddress(branches[iBran+46], &hitVec[iBran-1]);

  }

  dataIn->GetEntry((Int_t)nevt);
  
  if (Nhits==0) return track;

  for (Int_t iHit = 0;iHit < Nhits;iHit++){

    Point3D hitPos(hitVec[0][iHit] * cm, hitVec[1][iHit] * cm, hitVec[2][iHit] * cm);

    hit *xyz = new hit("MIND");
    xyz->set_point(hitPos);

    track.push_back(xyz);

  }

  dataIn->SetBranchStatus("*", 0);

  return track;
}


//***************************************************************
particle* root2dst::hits_to_fit() {
//***************************************************************

/*Makes a digital particle with only the hits, to be used in
  the fitting algorithm */


  ptype pT = DIGI;

  particle *postrack = new particle(pT, "unknown");

  const Int_t maxHits = 400;

  Int_t Nhits;
  Float_t hitVec[3][maxHits];

  for (Int_t iBran = 0;iBran < 4;iBran++){
    
    dataIn->SetBranchStatus(branches[iBran+46], 1);
    if (iBran==0) dataIn->SetBranchAddress(branches[iBran+46], &Nhits);
    else dataIn->SetBranchAddress(branches[iBran+46], &hitVec[iBran-1]);

  }

  dataIn->GetEntry((Int_t)nevt);

  for (Int_t iHit = 0;iHit < Nhits;iHit++){

    Point3D hitPos(hitVec[0][iHit] * cm, hitVec[1][iHit] * cm, hitVec[2][iHit] * cm);

    hit *xyz = new hit(*postrack,"MIND");
    xyz->set_point(hitPos);

    postrack->add_hit("MIND", xyz);

  }

  dataIn->SetBranchStatus("*", 0);

  return postrack;
}
