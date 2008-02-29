
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
bool root2dst::initialize(TTree *InPutTree, Char_t *OutFileName) {
//*************************************************************
    
  m.message("+++ root2dst init  function ++++",bhep::NORMAL);
  
  outgz.open(OutFileName);

  dataIn = InPutTree;
  dataIn->SetBranchStatus("*",0);

  fillBranchVector();

  nevt=0;
  nuEvent = NULL;
  
  return true;
}


//*************************************************************
bool root2dst::execute(){
//*************************************************************
    
  /*
    Take a ROOT tree event and convert it into bhep DST event
   */
  
  m.message("+++ root2dst execute function ++++",bhep::VERBOSE);

  nuEvent = createEvent();

  make_particles();
  cout << "Back in Execute Function" << endl;
  outgz.write(*nuEvent, nevt);
  cout <<"Event written to file" <<endl;
  nevt++;
 
  return true;  

}



//*************************************************************
bool root2dst::finalize() {
//*************************************************************
  
  outgz.close();

  delete nuEvent;

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
event* root2dst::createEvent() {
//*************************************************************

/* Defines an event and adds the event specific information. */

  Int_t integProps[4][2];
  Float_t floatProps[6][3];

  event *evt = new event((Int_t)nevt);

  for (Int_t iBran = 0;iBran<10;iBran++){

    if (branches[iBran] != "#"){
      dataIn->SetBranchStatus(branches[iBran], 1);
      if (iBran < 4) dataIn->SetBranchAddress(branches[iBran], &integProps[iBran]);
      else dataIn->SetBranchAddress(branches[iBran], &floatProps[iBran-4]);
    }

  }

  dataIn->GetEntry((Int_t)nevt);

  evt->set_vertex(floatProps[5][0], floatProps[5][1], floatProps[5][2]);

  //Add aditional properties
  for (Int_t iProp = 0;iProp<9;iProp++){

    if (branches[iProp] != "#"){
      if (iProp<2 || iProp==3) evt->add_property((string)branches[iProp], (string)ToString(integProps[iProp][0]));
      else if (iProp==2) {
	TString n1, n2;
	n1 = branches[iProp]+"Seed1";
	n2 = branches[iProp]+"Seed2";

	evt->add_property((string)n1, (string)ToString(integProps[iProp][0]));
	evt->add_property((string)n2, (string)ToString(integProps[iProp][1]));
      }
      else evt->add_property((string)branches[iProp], (string)ToString(floatProps[iProp-4][0]));
    }

  }

  dataIn->SetBranchStatus("*", 0);

  cout <<"Event Defined" << endl;

  return evt;

}

//***************************************************************
void root2dst::make_particles() {
//***************************************************************

/* Function to add relevant particles to event */
  particle *par;

  par = define_lead_particle();

  append_hits(par);

  nuEvent->add_true_particle(par);

  cout << "All particles defined and appended to event" << endl;

}

//***************************************************************
particle* root2dst::define_lead_particle() {
//***************************************************************

/* Function to define the properties of the lead particle from
   the event. e.g. the mu- from a numucc interaction */
  cout << "Defining Lead Particle" << endl;
  Int_t ID;
  Float_t vectors[2][3];
  
  ptype pT = TRUTH;
  string pName;
  
  for (Int_t iBran = 9;iBran < 12;iBran++){
    dataIn->SetBranchStatus(branches[iBran], 1);
    if (iBran==10) dataIn->SetBranchAddress(branches[iBran], &ID);
    else dataIn->SetBranchAddress(branches[iBran], &vectors[iBran-9]);
  }
  
  dataIn->GetEntry((Int_t)nevt);
  
  point pos(vectors[0][0], vectors[0][1], vectors[0][2]);
  
  Vector3D mom(vectors[2][2], vectors[2][1], vectors[2][0]);
  
  ray *leadstart = new ray(pos, mom);
  
  if (ID == 13) pName = "mu-";
  else pName = "mu+";
  
  particle *lead = new particle(pT, pName, *leadstart);

  dataIn->SetBranchStatus("*", 0);

  delete leadstart;

  cout << "Completed Particle def. and deleted ray" << endl;
  
  return lead;
  
}


//**************************************************************
void root2dst::append_hits(particle *par) {
//**************************************************************

/* Add vector of hit locations associated with particle par */

  Int_t Nhit;
  const Int_t maxHits = 400;

  dataIn->SetBranchStatus(branches[46], 1);
  dataIn->SetBranchAddress(branches[46], &Nhit);

  dataIn->GetEntry((Int_t)nevt);
  dataIn->SetBranchStatus(branches[46], 0);

  Float_t hitVec[3][maxHits];

  for (Int_t iBran = 47;iBran < 50;iBran++){
    
    dataIn->SetBranchStatus(branches[iBran], 1);
    dataIn->SetBranchAddress(branches[iBran], &hitVec[iBran]);

  }

  dataIn->GetEntry((Int_t)nevt);

  for (Int_t iHit = 0;iHit < Nhit;iHit++){

    Point3D hitPos(hitVec[0][iHit], hitVec[1][iHit], hitVec[2][iHit]);

    hit *xyz = new hit("MIND");
    xyz->set_point(hitPos);

    par->add_hit("MIND", xyz);

  }

}
