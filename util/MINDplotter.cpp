
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

  m.message("++Creating output root file++",bhep::VERBOSE);

  outFile = new TFile(outFileName, "recreate");

  statTree = new TTree("tree", "Tree with pattern rec and fit data for MIND");

  define_tree_branches();

  m.message("++plotter initialized",bhep::VERBOSE);

  return ok;
}

//*************************************************************************************
bool MINDplotter::execute(fitter& Fit, const bhep::event& evt,
			  bool success, bool patRec) {
//*************************************************************************************
  
  bool ok;

  _evNo = evt.event_number();
  _Fit = success;
  _fail = Fit.get_fail_type();
  
  for (int i = 0;i<4;i++)
    _hitType[i] = 0;

  track_fit( Fit );

  ok = extract_true_particle(evt, Fit, patRec);
  
  if (success) {
    
    State ste;
    ok = extrap_to_vertex(Fit.get_traj(), evt.vertex(), Fit, ste);
    
    if (ok) {
      max_local_chi2( Fit.get_traj() );
      position_pulls();
      direction_pulls();
      momentum_pulls();
    }

    //hadron_direction(Fit);

  }

  //If fit not successful set rec values to zero.
  if (!success){
    _X[1][0] = 0; _X[1][1] = 0;
    _X[2][0] = 0; _X[2][1] = 0;

    _Th[1][0] = 0; _Th[1][1] = 0;
    _Th[2][0] = 0; _Th[2][1] = 0;

    _qP[1] = 0; _qP[2] = 0;

    _Chi[0] = 0; _Chi[1] = 0;

    _Q[1] = 0; _Q[2] = 0;

    _haddot = 99;
    _hadE[1] = -99;
  }
  
  if (patRec)
    patternStats( Fit );
  
  //Fill tree event with the values.
  statTree->Fill();

  return ok;
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
void MINDplotter::define_tree_branches() {
//*************************************************************************************

  statTree->Branch("Evt", &_evNo, "EventNo/I");
  statTree->Branch("Fitted", &_Fit, "success/B");
  statTree->Branch("Fail", &_fail, "FailType/I");
  statTree->Branch("fitTypes", &_fitTracker, "predP/D:firstfit/D:reseed/D:predP2/D:succ/D");
  statTree->Branch("NeuEng", &_nuEng, "NuEng/D");
  statTree->Branch("Position", &_X, "truPos[2]/D:recPos[2]/D:ErrPos[2]/D");
  statTree->Branch("Direction", &_Th, "truTh[2]/D:recTh[2]/D:ErrTh[2]/D");
  statTree->Branch("Momentum", &_qP, "truqP/D:recqP/D:ErrqP/D");
  statTree->Branch("Charge", &_Q, "truQ/I:recQ/I:ID/B");
  statTree->Branch("FitChiInfo", &_Chi, "trajChi/D:MaxLoc/D");
  statTree->Branch("hadronMom", &_hadP, "hadP[3]/D");
  statTree->Branch("hadEng", &_hadE, "truE/D:recE/D");
  statTree->Branch("hadDir", &_haddot, "dotProd/D");
  statTree->Branch("NoHits", &_nhits, "nhits/I");
  statTree->Branch("HitBreakDown", &_hitType, "nTruMu/I:nInMu/I:nMuInMu/I:nFitN/I");
  statTree->Branch("XPositions", &_XPos, "X[nhits]/D");
  statTree->Branch("YPositions", &_YPos, "Y[nhits]/D");
  statTree->Branch("ZPositions", &_ZPos, "Z[nhits]/D");
  statTree->Branch("MuHits", &_mus, "truMu[nhits]/B");
  statTree->Branch("CandHits", &_cand, "inMu[nhits]/B");
  statTree->Branch("FittedNodes",&_node,"fitNode[nhits]/B");
  statTree->Branch("PatRecChi", &_pChi, "maxChiMu/D:MinChiHad/D:MaxConsecHol/D");

}

//*************************************************************************************
void MINDplotter::track_fit(fitter& fit) {
//*************************************************************************************

  if (_fail != 1 || _fail != 7) {
    vector<double> info = fit.get_fit_tracker();

    size_t nEnt = info.size();

    _fitTracker[0] = info[0];
    if (nEnt > 1) _fitTracker[1] = info[1];
    else _fitTracker[1] = 0;
    if (nEnt > 2) _fitTracker[2] = info[2];
    else _fitTracker[2] = 0;
    if (nEnt > 3) _fitTracker[3] = info[3];
    else _fitTracker[3] = 0;
    if (nEnt > 4) _fitTracker[4] = info[4];
    else _fitTracker[4] = 0;
  }
  else
    for (int i=0;i<5;i++)
      _fitTracker[i] = 0;

}

//*************************************************************************************
bool MINDplotter::extrap_to_vertex(const Trajectory& traj, 
				   const bhep::Point3D& vertexLoc,
				   fitter& fitObj, State& ste) {
//*************************************************************************************

  m.message("++Extrapolation function, Finding best fit to vertex++",bhep::VERBOSE);

  ste = traj.node(traj.first_fitted_node()).state();

  EVector pos(3,0); pos[2] = vertexLoc.z();
  EVector axis(3,0); axis[2] = 1;

  double R1 = 1000000000;
  double R2 = 0;
  double l;

  Surface* surf = new Ring(pos, axis, R1, R2);

  fitObj.man().geometry_svc().setup().add_surface("Detector","vertex",surf);
  bool ok = fitObj.man().navigation_svc().propagate(*surf,ste,l);
  fitObj.man().geometry_svc().setup().remove_surface("vertex");

  //Convert to slopes representation.
  fitObj.man().model_svc().model(RP::particle_helix)
	  .representation().convert(ste, RP::slopes_z);

  //Grab fitted vertex information.
  vert = ste.hv().vector();
  vertMat = ste.hv().matrix();

  return ok;
}

//**************************************************************************************
void MINDplotter::position_pulls() {
//**************************************************************************************

//Function to calculate the pull for a measurement
  m.message("++Calculating position pulls++",bhep::VERBOSE);
  
  //Reconstructed x,y position.
  _X[1][0] = vert[0]; _X[1][1] = vert[2];

  //Corresponding Error.
  if (vertMat[0][0]>0)
    _X[2][0] = sqrt(vertMat[0][0]);
  if (vertMat[1][1]>0)
    _X[2][1] = sqrt(vertMat[1][1]);

}

//**************************************************************************************
void MINDplotter::momentum_pulls() {
//**************************************************************************************

///Function to calculate momentum pulls.
  m.message("++Calculating momentum pulls++",bhep::VERBOSE);

  //Reconstructed q/P.
  if (vert[5] !=0) _Q[1] = (int)( vert[5]/fabs(vert[5]) );
  _qP[1] = vert[5];

  //Corresponding Error.
  if (vertMat[5][5]>0)
    _qP[2] = sqrt(vertMat[5][5]);

  //Correctly ID'd charge?.
  if (_Q[0] == _Q[1]) _Q[2] = true;
  else _Q[2] = false;

}

//**************************************************************************************
void MINDplotter::direction_pulls() {
//**************************************************************************************

  //Reconstructed direction.
  _Th[1][0] = vert[3]; _Th[1][1] = vert[4];

  //Corresponding error.
  if (vertMat[3][3]>0)
    _Th[2][0] = sqrt(vertMat[3][3]);
  if (vertMat[4][4]>0)
    _Th[2][1] = sqrt(vertMat[4][4]);

}

//*************************************************************************************
bool MINDplotter::extract_true_particle(const bhep::event& evt, fitter& Fit,
					bool patRec) {
//*************************************************************************************

/* sets true particle momentum for calculation and returns a reference
   to the particle */

  _nuEng = atof( evt.fetch_property("Enu").c_str() ) * GeV;

  const vector<bhep::particle*> Pospart = evt.true_particles();

  int count = 0;
  for (int iParts=0;iParts < (int)Pospart.size();iParts++){
    if (Pospart[iParts]->name().compare("mu-")==0){
      _Q[0] = -1;
      _truPart = Pospart[iParts];
      count++;
    } 
    else if (Pospart[iParts]->name().compare("mu+")==0){
      _Q[0] = 1;
      _truPart = Pospart[iParts];
      count++;
    } 
    if (Pospart[iParts]->name().compare("Hadronic_vector")==0){
      _hadP[0] = Pospart[iParts]->px();
      _hadP[1] = Pospart[iParts]->py();
      _hadP[2] = Pospart[iParts]->pz();
      _hadE[0] = atof( Pospart[iParts]->fetch_property("HadE").c_str() ) * GeV;
    }
  }
  
  if (count == 0) {
    cout << "No particles of muon or antimuon type in file" << endl;
    return false;
  }
  
  //Set true values of muon mom. etc.
  //True q/P.
  _qP[0] = _Q[0]/_truPart->p();
  //True direction.
  _Th[0][0]= _truPart->px()/_truPart->pz();
  _Th[0][1]= _truPart->py()/_truPart->pz();
  //True x,y position.
  _X[0][0] = evt.vertex().x(); _X[0][1] = evt.vertex().y();

  _nhits = Fit.get_nMeas();

  for (int iHits = 0;iHits < _nhits;iHits++){

    _XPos[iHits] = Fit.get_meas(iHits)->vector()[0];
    _YPos[iHits] = Fit.get_meas(iHits)->vector()[1];
    _ZPos[iHits] = Fit.get_meas(iHits)->surface().position()[2];

    if (!patRec)
      if ( Fit.get_traj().node(iHits).status("fitted") )
	_hitType[3]++;
  }
  
  return true;
}

//*************************************************************************************
void MINDplotter::hadron_direction(fitter& fit) {
//*************************************************************************************
  
  double normal;
  EVector fitunit = fit.get_had_unit();
  normal = sqrt(pow(_hadP[0],2)+pow(_hadP[1],2)+pow(_hadP[2],2));
  //cout << "Plotter: "<<fit.get_had_eng()<<endl;
  if (fitunit[0]==0 && fitunit[1]==0) {_haddot = 99; _hadE[1] = -99;}
  else {_haddot = fitunit[0]*(_hadP[0]/normal)+fitunit[1]*(_hadP[1]/normal)+fitunit[2]*(_hadP[2]/normal);
  _hadE[1] = fit.get_had_eng();}

}

//*************************************************************************************
void MINDplotter::max_local_chi2(const Trajectory& traj) {
//*************************************************************************************

  m.message("++Finding trajectory local chi max++",bhep::VERBOSE);

  size_t nNodes = traj.size();
  double trajMax = 0;
  
  for (size_t iNode = 0;iNode < nNodes;iNode++){
    
    if ( traj.node(iNode).qualitymap().has_key("predicted") )
      trajMax = TMath::Max(trajMax, traj.node(iNode).quality("predicted") );

  }

  _Chi[0] = traj.quality();
  _Chi[1] = trajMax;

}

//****************************************************************************************
void MINDplotter::patternStats(fitter& Fit) {
//****************************************************************************************
  //Event classifier version.
  _nhits = Fit.get_nMeas();
  const dict::Key candHit = "inMu";
  int nNode = 0;
  bool isMu;
  
  for (int iHits = 0;iHits < _nhits;iHits++){
    
    if (Fit.get_meas(iHits)->name("MotherParticle").compare("mu+")==0
	|| Fit.get_meas(iHits)->name("MotherParticle").compare("mu-")==0){
      isMu = true;
      _mus[iHits] = true;
      _hitType[0]++;
    }
    else {_mus[iHits] = false; isMu = false;}
    
    if ( _fail != 7 && Fit.get_meas(iHits)->names().has_key(candHit) ){
      if ( Fit.get_meas(iHits)->name(candHit).compare("True")==0 ){//has_key(candHit) ){
	_cand[iHits] = true;
	_hitType[1]++;
	if ( isMu ) _hitType[2]++;
	
	if ( Fit.get_traj().node(nNode).status("fitted") ){	
	  _node[iHits] = true; _hitType[3]++; }
	else _node[iHits] = false;
	nNode++;
      }
      else { _node[iHits] = false; _cand[iHits] = false; }
    } else if ( _fail != 7) { _node[iHits] = false; _cand[iHits] = false; }

  }
  
  if (_fail != 7){
    _pChi[0] = Fit.get_classifier().get_PatRec_Chis()[0];
    _pChi[1] = Fit.get_classifier().get_PatRec_Chis()[1];
    _pChi[2] = Fit.get_classifier().get_PatRec_Chis()[2];
  }
  
  for (int iclear = _nhits;iclear<300;iclear++){
    _mus[iclear] = false; _cand[iclear] = false;
    _node[iclear] = false;
  }

}

// //****************************************************************************************
// void MINDplotter::patternStats(fitter& Fit) {
// //****************************************************************************************
// //Original version. No event classifier.
//   _nhits = Fit.get_nMeas();
//   const dict::Key candHit = "inMu";
//   int nNode;
//   if (_fail != 7) nNode = (int)Fit.get_traj().size()-1;
//   bool isMu;
  
//   for (int iHits = 0;iHits < _nhits;iHits++){
    
//     if (Fit.get_meas(iHits)->name("MotherParticle").compare("mu+")==0
// 	|| Fit.get_meas(iHits)->name("MotherParticle").compare("mu-")==0){
//       isMu = true;
//       _mus[iHits] = true;
//       _hitType[0]++;
//     }
//     else {_mus[iHits] = false; isMu = false;}
    
//     if ( _fail != 7 && Fit.get_meas(iHits)->names().has_key(candHit) ){
//       if ( Fit.get_meas(iHits)->name(candHit).compare("True")==0 ){//has_key(candHit) ){
// 	_cand[iHits] = true;
// 	_hitType[1]++;
// 	if ( isMu ) _hitType[2]++;
	
// 	if ( Fit.get_traj().node(nNode).status("fitted") ){	
// 	  _node[iHits] = true; _hitType[3]++; }
// 	else _node[iHits] = false;
// 	nNode--;
//       }
//       else { _node[iHits] = false; _cand[iHits] = false; }
//     } else if ( _fail != 7) { _node[iHits] = false; _cand[iHits] = false; }

//   }
  
//   if (_fail != 7){
//     _pChi[0] = Fit.get_PatRec_Chis()[0];
//     _pChi[1] = Fit.get_PatRec_Chis()[1];
//     _pChi[2] = Fit.get_PatRec_Chis()[2];
//   }
  
//   for (int iclear = _nhits;iclear<300;iclear++){
//     _mus[iclear] = false; _cand[iclear] = false;
//     _node[iclear] = false;
//   }

// }
