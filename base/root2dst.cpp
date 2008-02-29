
#include <root2dst.h>

using namespace std;
using namespace bhep;


//*************************************************************
root2dst::root2dst(bhep::prlevel vlevel){
//*************************************************************
  
  level = vlevel;
  
  m = bhep::messenger(level);

  dataIn = NULL;

  nuEvent = NULL;
  
  m.message("+++ root2dst Constructor +++",bhep::VERBOSE); 
  
}


//*************************************************************
bool root2dst::initialize(TTree *InPutTree, Char_t *OutFileName) {
//*************************************************************
    
  m.message("+++ root2dst init  function ++++",bhep::NORMAL);
  
  outgz.open(OutFileName);

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
  
  m.message("+++ root2dst execute function ++++",bhep::VERBOSE);

  createEvent();

  make_particles();
  
  outgz.write(*nuEvent, nevt);
  
  nevt++;
 
  return true;  

}



//*************************************************************
bool root2dst::finalize() {
//*************************************************************
  
  outgz.close();
  delete nuEvent;
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
  branchIn.open("Branches.txt");

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

  nuEvent = new event((Int_t)nevt);

  for (Int_t iBran = 0;iBran<10;iBran++){

    if (branches[iBran] != "#"){
      dataIn->SetBranchStatus(branches[iBran], 1);
      if (iBran < 4) dataIn->SetBranchAddress(branches[iBran], &integProps[iBran]);
      else dataIn->SetBranchAddress(branches[iBran], &floatProps[iBran-4]);
    }

  }

  dataIn->GetEntry((Int_t)nevt);

  nuEvent->set_vertex(floatProps[5][0], floatProps[5][1], floatProps[5][2]);

  //Add aditional properties
  for (Int_t iProp = 0;iProp<9;iProp++){

    if (branches[iProp] != "#"){
      if (iProp<2 || iProp==3) nuEvent->add_property((string)branches[iProp], (string)ToString(integProps[iProp][0]));
      else if (iProp==2) {
	TString n1, n2;
	n1 = branches[iProp]+"Seed1";
	n2 = branches[iProp]+"Seed2";

	nuEvent->add_property((string)n1, (string)ToString(integProps[iProp][0]));
	nuEvent->add_property((string)n2, (string)ToString(integProps[iProp][1]));
      }
      else nuEvent->add_property((string)branches[iProp], (string)ToString(floatProps[iProp-4][0]));
    }

  }

  dataIn->SetBranchStatus("*", 0);

  cout <<"Event Defined" << endl;

}

//***************************************************************
void root2dst::make_particles() {
//***************************************************************

/* Function to add relevant particles to event */
  particle *par1, *par2;

  par1 = define_lead_particle();

  append_hits(par1);

  par2 = define_hadron();

  nuEvent->add_true_particle(par1);
  nuEvent->add_true_particle(par2);

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
  
  Point3D pos = nuEvent->vertex();

  Vector3D mom(moment[3], moment[2], moment[1]);
  
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
  
  Point3D pos = nuEvent->vertex();
  
  Vector3D mom(moment[0], moment[1], moment[2]);
  
  ray *hadstart = new ray(pos, mom);
  
  particle *had = new particle(pT, pName, *hadstart);
  cout << "Particle Defined" << endl;
  had->add_property("No.Hadrons in Shower", (string)ToString(nHad));

  dataIn->SetBranchStatus("*", 0);

  delete hadstart;

  cout << "Hadronic shower  vector defined" << endl;

  return had;
}


//**************************************************************
void root2dst::append_hits(particle *par) {
//**************************************************************

/* Add vector of hit locations associated with particle par */

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

    Point3D hitPos(hitVec[0][iHit], hitVec[1][iHit], hitVec[2][iHit]);

    hit *xyz = new hit("MIND");
    xyz->set_point(hitPos);

    par->add_hit("MIND", xyz);

  }

  dataIn->SetBranchStatus("*", 0);
}
