
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
bool MINDplotter::initialize(string outFileName, bhep::prlevel vlevel) {
//*************************************************************************************

  bool ok = true;

  level = vlevel;

  m = bhep::messenger(level);

  m.message("++Creating output root file++",bhep::VERBOSE);

  outFile = new TFile(outFileName.c_str(), "recreate");

  statTree = new TTree("tree", "Tree with pattern rec and fit data for MIND");

  define_tree_branches();

  m.message("++plotter initialized",bhep::VERBOSE);

  return ok;
}

//*************************************************************************************
bool MINDplotter::execute(fitter& Fit, const bhep::event& evt,
			  bool success, bool patRec) {
//*************************************************************************************
  
  bool ok1, ok2;

  _evNo = evt.event_number();
  _Fit = success;
  _fail = Fit.get_fail_type();
  _reFit = Fit.check_reseed();
  if ( _fail != 7 ){
    _visEng = Fit.get_classifier().get_vis_eng();
    Fit.get_classifier().get_planes( _plns );
  } else { _visEng = 0; _plns[0] = 0; _plns[1] = -1; }

  if (_fail != 7)
    _intType = Fit.get_classifier().get_int_type();
  else _intType = 7;
  
  for (int i = 0;i<4;i++)
    _hitType[i] = 0;
  _engTraj = 0;

  ok1 = extract_true_particle(evt, Fit, patRec);

  if (success) {
    
    State ste;
    ok2 = extrap_to_vertex(Fit.get_traj(), evt.vertex(), Fit, ste);

    if (_reFit) _leng = -Fit.get_traj().length();
    else _leng = Fit.get_traj().length();

    if (_leng !=0) {
      Fit.calculate_len_mom( _leng, _rangP );
    } else { _rangP[0] = 0; _rangP[1] = -99; }
    
    if (ok2) {
      max_local_chi2( Fit.get_traj() );
      position_pulls();
      direction_pulls();
      momentum_pulls();
    }

    hadron_direction(Fit);

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
    _leng = 0;
    _rangP[0] = 0;
    _rangP[1] = -99;

    _haddot = 99;
    _hadE[1] = -99;
  }
  
  if (patRec)
    patternStats( Fit );
  
  //Fill tree event with the values.
  statTree->Fill();
  
  return ok1;
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
  statTree->Branch("backFit",&_reFit,"backFit/B");
  statTree->Branch("Fail", &_fail, "FailType/I");
  statTree->Branch("interaction",&_intType,"Inter/I");
  statTree->Branch("NeuEng", &_nuEng, "NuEng/D");
  statTree->Branch("visibleEng", &_visEng, "visEng/D");
  statTree->Branch("visEngTraj",&_engTraj, "engTraj/D");
  statTree->Branch("trajEngVar",&_engvar,"engVar/D");
  statTree->Branch("Position", &_X, "truPos[2]/D:recPos[2]/D:ErrPos[2]/D");
  statTree->Branch("Direction", &_Th, "truTh[2]/D:recTh[2]/D:ErrTh[2]/D");
  statTree->Branch("Momentum", &_qP, "truqP/D:recqP/D:ErrqP/D");
  statTree->Branch("Charge", &_Q, "truQ/I:recQ/I:ID/B");
  statTree->Branch("length", &_leng,"lenTraj/D");
  statTree->Branch("RangeMomentum", &_rangP,"rangP/D:rangErr/D");
  statTree->Branch("FitChiInfo", &_Chi, "trajChi/D:MaxLoc/D");
  statTree->Branch("hadronMom", &_hadP, "hadP[3]/D");
  statTree->Branch("hadEng", &_hadE, "truE/D:recE/D");
  statTree->Branch("hadDir", &_haddot, "dotProd/D");
  statTree->Branch("NoPlanes", &_plns, "nplanes/I:freeplanes/I");
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
bool MINDplotter::extrap_to_vertex(const Trajectory& traj, 
				   const bhep::Point3D& vertexLoc,
				   fitter& fitObj, State& ste) {
//*************************************************************************************

  m.message("++Extrapolation function, Finding best fit to vertex++",bhep::VERBOSE);
  if ( fitObj.check_reseed() ) ste = traj.node(traj.last_fitted_node()).state();
  else ste = traj.node(traj.first_fitted_node()).state();

  EVector pos(3,0); pos[2] = vertexLoc.z();
  EVector axis(3,0); axis[2] = 1;

  double R1 = 1000000000;
  double R2 = 0;
  double l;

  Ring surf(pos, axis, R1, R2);

  fitObj.man().geometry_svc().setup().add_surface("Detector","vertex",&surf);
  bool ok = fitObj.man().navigation_svc().propagate(surf,ste,l);
  fitObj.man().geometry_svc().setup().remove_surface("vertex");

  //Convert to slopes representation.
  fitObj.man().model_svc().model(RP::particle_helix)
	  .representation().convert(ste, RP::slopes_curv_z);

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
  _X[1][0] = vert[0]; _X[1][1] = vert[1];

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
  
  _nhits = Fit.get_nMeas();

  for (int iHits = 0;iHits < _nhits;iHits++){

    _XPos[iHits] = Fit.get_meas(iHits)->vector()[0];
    _YPos[iHits] = Fit.get_meas(iHits)->vector()[1];
    _ZPos[iHits] = Fit.get_meas(iHits)->position()[2];

    if (!patRec)
      if ( Fit.get_traj().node(iHits).status("fitted") )
	_hitType[3]++;
  }

  if (count == 0) {
    cout << "No particles of muon or antimuon type in file" << endl;
    _Q[0] = 0;
    _qP[0] = 0;
    _Th[0][0] = _Th[0][1]= 0;
    _X[0][0] = _X[0][1] = 0;
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
  
  return true;
}

//*************************************************************************************
void MINDplotter::hadron_direction(fitter& fit) {
//*************************************************************************************
  
  double normal;
  EVector fitunit;
  normal = sqrt(pow(_hadP[0],2)+pow(_hadP[1],2)+pow(_hadP[2],2));
  
  if ( _nhits >= 2 ){
    fitunit = fit.get_had_unit();
    
    _haddot = fitunit[0]*(_hadP[0]/normal)
      +fitunit[1]*(_hadP[1]/normal)
      +fitunit[2]*(_hadP[2]/normal);
    
  } else _haddot = 99;
  
  _hadE[1] = fit.get_had_eng();
  
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
  const dict::Key engDep = "E_dep";
  int nNode = 0;
  if ( Fit.check_reseed() ) nNode = (int)Fit.get_traj().size()-1;
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
	_engTraj += bhep::double_from_string( Fit.get_meas(iHits)->name(engDep) ) * GeV;
	if ( isMu ) _hitType[2]++;
	
	if ( Fit.get_traj().node(nNode).status("fitted") && _fail!=1 && _fail<4 ){	
	  _node[iHits] = true; _hitType[3]++; }
	else _node[iHits] = false;
	if ( Fit.check_reseed() ) nNode--;
	else nNode++;
      }
      else { _node[iHits] = false; _cand[iHits] = false; }
    } else if ( _fail != 7) { _node[iHits] = false; _cand[iHits] = false; }

  }
  
  if (_fail != 7){
    _pChi[0] = Fit.get_classifier().get_PatRec_Chis()[0];
    _pChi[1] = Fit.get_classifier().get_PatRec_Chis()[1];
    _pChi[2] = Fit.get_classifier().get_PatRec_Chis()[2];

    _engvar = 0;
    for(int ii=0;ii<_hitType[1];ii++)
      _engvar += pow( bhep::double_from_string( Fit.get_traj().node(ii).measurement().name(engDep) ) * GeV - _engTraj/_hitType[1], 2);
    _engvar /= (_hitType[1]-1);
  } else _engvar = -1;
  
  for (int iclear = _nhits;iclear<300;iclear++){
    _mus[iclear] = false; _cand[iclear] = false;
    _node[iclear] = false;
  }

}
