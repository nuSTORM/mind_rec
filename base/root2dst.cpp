
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
bool root2dst::initialize(TTree *InPutTree, TString OutFileName,
			  double res) {
//*************************************************************
    
  m.message("+++ root2dst init  function ++++",bhep::NORMAL);
  
  outgz.open((string)OutFileName);

  dataIn = InPutTree;
  dataIn->SetBranchStatus("*",0);

  long seed = 178263094;
  ranGen = RanluxEngine(seed, 4);

  sigMa = res;

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
  
  outgz.write(*nuEvent[nevt], nevt);
  
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
void root2dst::createEvent() {
//*************************************************************

/* Defines an event and adds the event specific information. */

  Int_t integs[2][2];
  Float_t floas[6][3];

  event *evt = new event((Int_t)nevt);
  nuEvent.push_back(evt);

  dataIn->SetBranchStatus("Random",1); dataIn->SetBranchAddress("Random",&integs[0]);
  dataIn->SetBranchStatus("Nutype",1); dataIn->SetBranchAddress("Nutype",&integs[1]);
  dataIn->SetBranchStatus("Enu",1); dataIn->SetBranchAddress("Enu",&floas[0]);
  dataIn->SetBranchStatus("Xf",1); dataIn->SetBranchAddress("Xf",&floas[1]);
  dataIn->SetBranchStatus("Yf",1); dataIn->SetBranchAddress("Yf",&floas[2]);
  dataIn->SetBranchStatus("Q2f",1); dataIn->SetBranchAddress("Q2f",&floas[3]);
  dataIn->SetBranchStatus("Wf",1); dataIn->SetBranchAddress("Wf",&floas[4]);
  dataIn->SetBranchStatus("Vertex",1); dataIn->SetBranchAddress("Vertex",&floas[5]);
  
  dataIn->GetEntry((Int_t)nevt);

  nuEvent[nevt]->set_vertex(floas[5][0] * cm, floas[5][1] * cm, floas[5][2] * cm);

  //Add aditional properties
  nuEvent[nevt]->add_property("RandomSeed1", (string)ToString(integs[0][0]));
  nuEvent[nevt]->add_property("RandomSeed2", (string)ToString(integs[0][1]));
  nuEvent[nevt]->add_property("Nutype", (string)ToString(integs[1][0]));
  nuEvent[nevt]->add_property("Enu", (string)ToString(floas[0][0]));
  nuEvent[nevt]->add_property("Xf", (string)ToString(floas[1][0]));
  nuEvent[nevt]->add_property("Yf", (string)ToString(floas[2][0]));
  nuEvent[nevt]->add_property("Q2f", (string)ToString(floas[3][0]));
  nuEvent[nevt]->add_property("Wf", (string)ToString(floas[4][0]));
  
  dataIn->SetBranchStatus("*", 0);

  cout <<"Event Defined" << endl;
}

//***************************************************************
void root2dst::make_particles() {
//***************************************************************

/* Function to add relevant particles to event */
  particle *muon, *hadron, *digPar;

  muon = define_lead_particle();
  
  hadron = define_hadron();

  vector<hit*> muHit, hadHit;
  bool ok = hits_fromFile(muHit, hadHit);

  if (ok) {

    if (muHit.size() != 0)
      for (Int_t iHit=0;iHit<(Int_t)muHit.size();iHit++){
	muHit[iHit]->set_mother_particle(*muon);
	muon->add_hit("MIND", muHit[iHit]);
      }
    if (hadHit.size() != 0)
      for (Int_t jHit=0;jHit<(Int_t)hadHit.size();jHit++){
	hadHit[jHit]->set_mother_particle(*hadron);
	hadron->add_hit("MIND", hadHit[jHit]);
      }
    
    digPar = create_digital_representation(*muon, *hadron, muHit, hadHit);
    
    nuEvent[nevt]->add_true_particle(muon);
    nuEvent[nevt]->add_true_particle(hadron);
    nuEvent[nevt]->add_digi_particle(digPar);

  }
  else cout << "No. hits in event" << endl;

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

  dataIn->SetBranchStatus("Idlead", 1); dataIn->SetBranchAddress("Idlead", &ID);
  dataIn->SetBranchStatus("Plead", 1); dataIn->SetBranchAddress("Plead", &moment);

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

  dataIn->SetBranchStatus("Nhad", 1); dataIn->SetBranchAddress("Nhad", &nHad);
  dataIn->SetBranchStatus("Phadg", 1); dataIn->SetBranchAddress("Phadg", &moment);

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


//**************************************************************************
bool root2dst::hits_fromFile(vector<hit*>& muHit, vector<hit*>& hadHit) {
//**************************************************************************

/* Add vector of hit locations associated with particle par */

  const Int_t maxHits = 400;

  Int_t NmuHits, NhadHits;
  Float_t muHitVec[3][maxHits], hadHitVec[3][maxHits];

  dataIn->SetBranchStatus("Nhits", 1); dataIn->SetBranchAddress("Nhits", &NmuHits);
  dataIn->SetBranchStatus("Xhit", 1); dataIn->SetBranchAddress("Xhit", &muHitVec[0]);
  dataIn->SetBranchStatus("Yhit", 1); dataIn->SetBranchAddress("Yhit", &muHitVec[1]);
  dataIn->SetBranchStatus("Zhit", 1); dataIn->SetBranchAddress("Zhit", &muHitVec[2]);
  dataIn->SetBranchStatus("Nhhits", 1); dataIn->SetBranchAddress("Nhhits", &NhadHits);
  dataIn->SetBranchStatus("Xhhit", 1); dataIn->SetBranchAddress("Xhhit", &hadHitVec[0]);
  dataIn->SetBranchStatus("Yhhit", 1); dataIn->SetBranchAddress("Yhhit", &hadHitVec[1]);
  dataIn->SetBranchStatus("Zhhit", 1); dataIn->SetBranchAddress("Zhhit", &hadHitVec[2]);

  dataIn->GetEntry((Int_t)nevt);
  
  if (NmuHits==0 && NhadHits==0) return false;

  if (NmuHits!=0)
    for (Int_t iHit = 0;iHit < NmuHits;iHit++){
      
      Point3D hitPos(muHitVec[0][iHit] * cm, muHitVec[1][iHit] * cm, muHitVec[2][iHit] * cm);
      
      hit *xyz = new hit("MIND");
      xyz->set_point(hitPos);
      
      muHit.push_back(xyz);
      
    }

  if (NhadHits!=0)
    for (Int_t jHit = 0;jHit < NhadHits;jHit++){
      
      Point3D hitPos(hadHitVec[0][jHit] * cm, hadHitVec[1][jHit] * cm, hadHitVec[2][jHit] * cm);
      
      hit *xyz = new hit("MIND");
      xyz->set_point(hitPos);
      
      hadHit.push_back(xyz);
      
    }

  dataIn->SetBranchStatus("*", 0);

  return true;
}


//***************************************************************
particle* root2dst::create_digital_representation(particle& mu, particle& had,
						  const vector<hit*>& muHit,
						  const vector<hit*>& hadHit) {
//***************************************************************

/*Makes a digital particle with only the hits, to be used in
  the fitting algorithm */


  ptype pT = DIGI;

  particle *hitMap = new particle(pT, "unknown");

  double X, Y, Z;

  for (Int_t iHit = 0;iHit < (Int_t)muHit.size();iHit++) {
    
    X = muHit[iHit]->x().x()/cm + RandGauss::shoot(&ranGen, 0, sigMa);
    Y = muHit[iHit]->x().y()/cm + RandGauss::shoot(&ranGen, 0, sigMa);
    Z = muHit[iHit]->x().z()/cm;
    
    Point3D hitPos(X * cm,Y * cm,Z * cm);
    hit* digHit = new hit("MIND");
    digHit->set_point(hitPos);

    digHit->set_mother_particle(mu);
    hitMap->add_hit("MIND", digHit);

  }
  for (Int_t jHit = 0;jHit < (Int_t)hadHit.size();jHit++) {
    
    X = hadHit[jHit]->x().x()/cm + RandGauss::shoot(&ranGen, 0, sigMa);
    Y = hadHit[jHit]->x().y()/cm + RandGauss::shoot(&ranGen, 0, sigMa);
    Z = hadHit[jHit]->x().z()/cm;

    Point3D hitPos(X * cm,Y * cm,Z * cm);
    hit* digHit = new hit("MIND");
    digHit->set_point(hitPos);

    digHit->set_mother_particle(had);
    hitMap->add_hit("MIND", digHit);

  }

  return hitMap;
}
