
#include <MINDplotter.h>

//*************************************************************************************
MINDplotter::MINDplotter() {
//*************************************************************************************

}

//*************************************************************************************
MINDplotter::~MINDplotter() {
//*************************************************************************************

}

//*************************************************************************************
bool MINDplotter::initialize(TString outFileName, bhep::prlevel vlevel) {
//*************************************************************************************

  bool ok = true;

  level = vlevel;

  m = bhep::messenger(level);

  counterlo = 0;
  counterhi = 0;
  
  FitX = NULL; truX = NULL; xPull = NULL;
  FitY = NULL; truY = NULL; yPull = NULL;
  FitP = NULL; truMom = NULL; pPull = NULL;
  FitTx = NULL; truTx = NULL; tXpull = NULL;
  FitTy = NULL; truTy = NULL; tYpull = NULL;
  pSpec = NULL;
  misIDp = NULL;
  hitSpec = NULL;
  misIDhit = NULL;
  locChi = NULL; locVp = NULL; trajVp = NULL;
  Eff = NULL; Purit = NULL;

  part = NULL;

  m.message("++Creating output root file++",bhep::VERBOSE);

  outFile = new TFile(outFileName, "recreate");

  m.message("++plotter initialized",bhep::VERBOSE);

  return ok;
}

//*************************************************************************************
bool MINDplotter::execute(const EVector& V, 
			  const EMatrix& M, const bhep::event& evt) {
//*************************************************************************************

//Creates output file to take plots.
  truP = false;

  truP = extract_true_particle(V,M,evt);

  return truP;
}

//*************************************************************************************
bool MINDplotter::finalize() {
//*************************************************************************************

  bool ok = true;

  m.message("++Finalizing Output++",bhep::VERBOSE);

  outFile->Write();
  outFile->Close();

  return ok;
}

//*************************************************************************************
bool MINDplotter::extrap_to_vertex(const Trajectory& traj, 
				   const bhep::Point3D& vertexLoc,
				   fitter* fitObj, State& ste) {
//*************************************************************************************

  m.message("++Extrapolation function, Finding best fit to vertex++",bhep::NORMAL);

  ste = traj.node(traj.first_fitted_node()).state();

  EVector pos(3,0); pos[2] = vertexLoc.z();
  EVector axis(3,0); axis[2] = 1;

  double R1 = 1000000000;
  double R2 = 0;
  double l;

  Surface* surf = new Ring(pos, axis, R1, R2);

  fitObj->man().geometry_svc().setup().add_surface("Detector","vertex",surf);
  bool ok = fitObj->man().navigation_svc().propagate(*surf,ste,l);
  fitObj->man().geometry_svc().setup().remove_surface("vertex");

  return ok;
}

//**************************************************************************************
void MINDplotter::position_pulls(const EVector& V, const EMatrix& M, 
				 const bhep::event& evt) {
//**************************************************************************************

//Function to calculate the pull for a measurement
  m.message("++Calculating position pulls++",bhep::VERBOSE);

  if (xPull == NULL){
    xPull = new TH1F("Xpull", "Pulls for vertex X component",120,-15,15);
    FitX = new TH1F("FitX", "Fitted X position",1400,-700,700);
    truX = new TH1F("truX", "True X position",1400,-700,700);
  }
  if (yPull == NULL){
    yPull = new TH1F("Ypull", "Pulls for vertex Y component",120,-15,15);
    FitY = new TH1F("FitY", "Fitted Y position",1400,-700,700);
    truY = new TH1F("truY", "True Y position",1400,-700,700);
  }
  double pull_x, pull_y;

  pull_x = (V[0]-evt.vertex().x())/sqrt(M[0][0]);

  pull_y = (V[1]-evt.vertex().y())/sqrt(M[1][1]);
  
  xPull->Fill(pull_x);
  FitX->Fill(V[0]);
  truX->Fill(evt.vertex().x());
  yPull->Fill(pull_y);
  FitY->Fill(V[1]);
  truY->Fill(evt.vertex().y());

}

//**************************************************************************************
void MINDplotter::momentum_pulls(const EVector& V, const EMatrix& M) {
//**************************************************************************************

///Function to calculate momentum pulls.
  m.message("++Calculating momentum pulls++",bhep::VERBOSE);

  if (pPull==NULL) {
    pPull = new TH1F("pPull", "Pulls for vertex q/p",120,-15,15);
    FitP = new TH1F("FitP", "Fitted q/P at vertex",1000,-1,1);
    truMom = new TH1F("truP", "True q/P at vertex",1000,-1,1);
  }
  double pull_p;
  
  
  pull_p = (V[5]-(1/part->p()))/sqrt(M[5][5]);
  cout << "q/P: "<< V[5] <<"-"<<1/part->p()<<endl;
  pPull->Fill(pull_p);
  FitP->Fill(V[5]);
  truMom->Fill(tru_q_/part->p());

}

//**************************************************************************************
void MINDplotter::direction_pulls(const EVector& V, const EMatrix& M) {
//**************************************************************************************

  if(tXpull==NULL) {
    tXpull = new TH1F("tXpull", "Pulls for vertex #theta_{x} component",120,-15,15);
    FitTx = new TH1F("FitTx", "Fitted slope in x at vertex",1000,-2,2);
    //truTx = new TH1F("truTx", "True slope in x at vertex",1000,-2,2);
  }
  if(tYpull==NULL) {
    tYpull = new TH1F("tYpull", "Pulls for vertex #theta_{y} component",120,-15,15);
    FitTy = new TH1F("FitTy", "Fitted slope in y at vertex",1000,-2,2);
    //truTy = new TH1F("truTy", "True slope in y at vertex",1000,-2,2);
  }

  double pull_x, pull_y;

  double tetaX = part->px()/part->pz();
  double tetaY = part->py()/part->pz();

  pull_x = (V[3]-tetaX)/sqrt(M[3][3]);

  pull_y = (V[4]-tetaY)/sqrt(M[4][4]);

  tXpull->Fill(pull_x);
  FitTx->Fill(V[3]-tetaX);
  //truTx->Fill(tetaX);
  tYpull->Fill(pull_y);
  FitTy->Fill(V[4]-tetaY);
  //truTy->Fill(tetaY);

}

//**************************************************************************************
void MINDplotter::momentum_efficiency(const EVector& V, const EMatrix& M) {
//**************************************************************************************

  m.message("++Momentum efficiency function++",bhep::VERBOSE);

  if (pSpec==NULL) pSpec = new TH1F("pSpec", "Fitted track true momentum spectrum",50,0,50);
  if (misIDp==NULL) misIDp = new TH1F("pMis", "Charge mis-id true momentum spectrum",50,0,50);

  pSpec->Fill(part->p() / GeV);
  
  double q;
  if (V[5] != 0) q = V[5]/fabs(V[5]);
  if (q != tru_q_)
    misIDp->Fill(part->p() / GeV);
  
}

//*************************************************************************************
void MINDplotter::hit_efficiency(const EVector& V, const EMatrix& M) {
//*************************************************************************************

  m.message("++Hit Efficiency function++",bhep::VERBOSE);

  if (hitSpec==NULL) hitSpec = new TH1F("hits", "Fitted track hit spectrum",400,0,400);
  if (misIDhit==NULL) misIDhit = new TH1F("hitMis", "Charge mis-id hit spectrum",400,0,400);

  hitSpec->Fill((int)part->hits("MIND").size());
  
  double q = 1;
  if (V[5] != 0) q = V[5]/fabs(V[5]);
  if (q != tru_q_)
    misIDhit->Fill((int)part->hits("MIND").size());

}

//*************************************************************************************
bool MINDplotter::extract_true_particle(const EVector& V, const EMatrix& M, 
						   const bhep::event& evt) {
//*************************************************************************************

/* sets true particle momentum for calculation and returns a reference
   to the particle */

  const vector<bhep::particle*> Pospart = evt.true_particles();
 
  p_ = pow(fabs(V[5]), -1);
  d_p_ = sqrt(M[5][5])/fabs(V[5]);
  
  int count = 0;
  for (int iParts=0;iParts<(int)Pospart.size();iParts++){
    if (Pospart[iParts]->name().compare("mu-")==0){
      tru_q_ = -1;
      part = Pospart[iParts];
      count++;
    } 
    else if (Pospart[iParts]->name().compare("mu+")==0){
      tru_q_ = 1;
      part = Pospart[iParts];
      count++;
    } 
  }
  if (count == 0) {
    cout << "No particles of muon or antimuon type in file" << endl;
    return false;
  }
  
  return true;
}

//*************************************************************************************
void MINDplotter::max_local_chi2(const Trajectory& traj, double maxChi, const EVector& V) {
//*************************************************************************************

  m.message("++Finding trajectory local chi max++",bhep::VERBOSE);

  if (locChi==NULL) {
    locChi = new TH1F("locChi","Max local #chi^{2} per trajectory",(int)maxChi*2,0,maxChi);
    locVp = new TProfile("locVp","Max local #chi^{2} vs. true particle momentum",50,0,50);
    trajVp = new TProfile("trajVp","Trajectory #chi^{2} vs. true particle momentum",50,0,50);
  }

  size_t nNodes = traj.size();
  //const int trajSize = nNodes;
  //double chi[trajSize];
  double trajMax = 0;

  for (size_t iNode = 0;iNode < nNodes;iNode++){

    trajMax = TMath::Max(trajMax, traj.node(iNode).quality() );
    //chi[iNode] = traj.node(iNode).quality();

  }

  //trajMax = TMath::MaxElement((Long64_t)nNodes, chi);
  locChi->Fill(trajMax);
  locVp->Fill(part->p()/GeV, trajMax);
  trajVp->Fill(part->p()/GeV, traj.quality());

  double q;
  if (V[5] != 0) q = V[5]/fabs(V[5]);

  if (trajMax < 0.2*maxChi && q != tru_q_ && counterlo<10)
    lowChi_MisID_track(traj,trajMax);

  if (trajMax > 0.8*maxChi && q == tru_q_ && counterhi<10)
    highChi_ID_track(traj,trajMax);

}

//****************************************************************************************
void MINDplotter::lowChi_MisID_track(const Trajectory& traj, double trajMax) {
//****************************************************************************************

  m.message("++ Output of low Chi2 mis-ID track ++",bhep::VERBOSE);

  const int nMeas = (int)traj.nmeas();
  double X[nMeas], Z[nMeas];

  TString plotName, plotTitle;
  plotName = "loMisID"+to_string(counterlo);
  plotTitle = "Trajectory of charge Mis-ID particle with low max local #chi^{2} = "+to_string(trajMax);

  for (int iMeas = 0;iMeas < nMeas;iMeas++){

    X[iMeas] = traj.measurement(iMeas).vector()[0]/cm;
    Z[iMeas] = traj.measurement(iMeas).position()[2]/cm;

  }

  TGraph* plot = new TGraph(nMeas,Z,X);
  plot->SetName(plotName);
  plot->SetTitle(plotTitle);
  plot->GetXaxis()->SetTitle("Z position (cm)");
  plot->GetYaxis()->SetTitle("X position (cm)");
  plot->Write();

  counterlo++;

}

//****************************************************************************************
void MINDplotter::highChi_ID_track(const Trajectory& traj, double trajMax) {
//****************************************************************************************

  m.message("++ Output of high Chi2 track ++",bhep::VERBOSE);

  const int nMeas = (int)traj.nmeas();
  double X[nMeas], Z[nMeas];
  
  TString plotName, plotTitle;
  plotName = "hiChiID"+to_string(counterhi);
  plotTitle = "Trajectory of correctly ID'd particle with high max local #chi^{2} = "+to_string(trajMax);
  
  for (int iMeas = 0;iMeas < nMeas;iMeas++){

    X[iMeas] = traj.measurement(iMeas).vector()[0]/cm;
    Z[iMeas] = traj.measurement(iMeas).position()[2]/cm;
    
  }
  
  
  TGraph* plot2 = new TGraph(nMeas,Z,X);
  plot2->SetName(plotName);
  plot2->SetTitle(plotTitle);
  plot2->GetXaxis()->SetTitle("Z position (cm)");
  plot2->GetYaxis()->SetTitle("X position (cm)");
  plot2->Write();

  counterhi++;
  
}

//****************************************************************************************
void MINDplotter::patternStats(const EVector& vec) {
//****************************************************************************************
  if (Eff == NULL){
    Eff = new TProfile("eff","Pattern recognition efficiency as func muon true p",50,0,50);
    Purit = new TProfile("purit","Pattern recognition purity as func muon true p",50,0,50);
  }

  double efficiency;
  double purity;
  efficiency = vec[1]/vec[2];
  purity = vec[1]/vec[0];

  Eff->Fill(part->p()/GeV, efficiency);
  Purit->Fill(part->p()/GeV, purity);

}
