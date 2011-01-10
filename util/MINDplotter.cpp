
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
void MINDplotter::initialize(string outFileName, bool patRec, bool clust,
			     bhep::prlevel vlevel) {
//*************************************************************************************

  //bool ok = true;
  _patR = patRec;
  _clu = clust;

  level = vlevel;

  m = bhep::messenger(level);

  m.message("++Creating output root file++",bhep::VERBOSE);

  outFile = new TFile(outFileName.c_str(), "recreate");

  statTree = new TTree("tree", "Tree with pattern rec and fit data for MIND");

  define_tree_branches();

  m.message("++plotter initialized",bhep::VERBOSE);

  //return ok;
}

//*************************************************************************************
void MINDplotter::execute(fitter& Fit, const bhep::event& evt, bool success) {
//*************************************************************************************
  
  bool ok1, ok2;
  
  _evNo = evt.event_number();
  bhep::Point3D evVert = evt.vertex();
  _vert[0] = evVert.x(); _vert[1] = evVert.y(); _vert[2] = evVert.z();
  if (_clu){
    get_tru_int_type( evt );
    if ( evt.find_iproperty("Charm") )
      _charm[0] = evt.fetch_iproperty("Charm");
    else _charm[0] = 0;
    if ( evt.find_iproperty("CharmHad") )
      _charm[1] = evt.fetch_iproperty("CharmHad");
    else _charm[1] = 0;
    if ( evt.find_dproperty("Q2") )
      _Q2 = evt.fetch_dproperty("Q2");
    else _Q2 = 0.;
    if ( evt.find_dproperty("EngTrans") )
      _engTrans = evt.fetch_dproperty("EngTrans");
  }
  
  _Fit = success;
  _fail = Fit.get_fail_type();
  _reFit = Fit.check_reseed();
  if ( _fail != 7 ){
    _visEng = Fit.get_classifier().get_vis_eng();
    _plns[0] = Fit.get_classifier().get_planes();
    _plns[1] = Fit.get_classifier().get_free_planes();
  } else { _visEng = 0; _plns[0] = 0; _plns[1] = -1; }
  
  if (_fail != 7)
    _intType = Fit.get_classifier().get_int_type();
  else _intType = 7;
  
  for (int i = 0;i<5;i++){
    _hitType[i] = 0;
    // if ( i < 4 ) { _chadP[i] = 0.; _nhadP[i] = 0.; }
//     //else _hadE[0] = 0.;
    if ( i < 3 ) _hadP[i] = 0.;
    else _hadE[0] = 0.;
  }
  _engTraj = 0;
  _hadE[1] = 0;
  // _nhad[0] = 0;
//   _nhad[1] = 0;

  if ( _clu )
    ok1 = extract_true_particle2(evt, Fit);
  else
    ok1 = extract_true_particle1(evt, Fit);
  
  if (success) {
    
    State ste;
    ok2 = extrap_to_vertex(Fit.get_traj(), evt.vertex(), Fit, ste);
    
    //if (_reFit) _leng = -Fit.get_traj().length();
    //else _leng = Fit.get_traj().length();
    
    //if (_leng !=0) {
    //Fit.calculate_len_mom( _leng, _rangP );
    //} else { _rangP[0] = 0; _rangP[1] = -99; }
    //_rangP[0] = 0; _rangP[1] = -99;

    if (ok2) {
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
    //_leng = 0;
    //_rangP[0] = 0;
    //_rangP[1] = -99;

    //_haddot = 99;
    //_hadE[1] = -99;
  }
  
  if (_patR){
    if ( _clu )
      patternStats2( Fit );
    else
      patternStats1( Fit );
  }
  
  //Fill tree event with the values.
  int fillcatch = statTree->Fill();
  
  //return ok1;
}

void MINDplotter::get_tru_int_type(const bhep::event& evt){
  
  string intName = evt.fetch_sproperty("IntType");

  if ( intName=="CCQE" || intName=="<QES - Weak[CC]>" )
    _truInt = 1;
  else if ( intName=="NCQE" || intName=="<QES - Weak[NC]>" )
    _truInt = 2;
  else if ( intName=="CCDIS" || intName=="<DIS - Weak[CC]>" )
    _truInt = 3;
  else if ( intName=="NCDIS" || intName=="<DIS - Weak[NC]>" )
    _truInt = 4;
  else if ( intName=="1piRes" )
    _truInt = 5;
  else if ( intName=="miscRes" || intName=="<RES - Weak[CC]>" || intName=="<RES - Weak[NC]>" )
    _truInt = 6;
  else if ( intName=="eEl" || intName=="<NuEEL - Weak[NC]>" || intName=="<NuEEL - Weak[CC]>" )
    _truInt = 7;
  else if ( intName=="muINVe" || intName=="<IMD - Weak[CC]>" )
    _truInt = 7;
  else
    _truInt = 8;
  
}

//*************************************************************************************
void MINDplotter::finalize() {
//*************************************************************************************

  //bool ok = true;

  m.message("++Finalizing Output++",bhep::VERBOSE);
  
  outFile->Write();
  outFile->Close();

  //return ok;
}

//*************************************************************************************
void MINDplotter::define_tree_branches() {
//*************************************************************************************

  statTree->Branch("Evt", &_evNo, "EventNo/I");
  statTree->Branch("Fitted", &_Fit, "success/B");
  if (_clu) statTree->Branch("TrueInteraction",&_truInt,"truInt/I");
  statTree->Branch("evVertex", &_vert, "vert[3]/D");
  statTree->Branch("backFit",&_reFit,"backFit/B");
  statTree->Branch("Fail", &_fail, "FailType/I");
  statTree->Branch("interaction",&_intType,"Inter/I");
  statTree->Branch("NeuEng", &_nuEng, "NuEng/D");
  if (_clu) {
    statTree->Branch("Charm", &_charm, "isCharm/I:pdg/I");
    statTree->Branch("InteractionQ2", &_Q2, "Q2/D");
  }
  statTree->Branch("visibleEng", &_visEng, "visEng/D");
  statTree->Branch("visEngTraj",&_engTraj, "engTraj/D");
  statTree->Branch("trajEngVar",&_engvar,"meanDep/D:engVar/D");
  statTree->Branch("Position", &_X, "truPos[2]/D:recPos[2]/D:ErrPos[2]/D");
  statTree->Branch("Direction", &_Th, "truTh[2]/D:recTh[2]/D:ErrTh[2]/D");
  statTree->Branch("Momentum", &_qP, "truqP/D:recqP/D:ErrqP/D");
  statTree->Branch("Charge", &_Q, "truQ/I:recQ/I:ID/B");
  //statTree->Branch("length", &_leng,"lenTraj/D");
  //statTree->Branch("RangeMomentum", &_rangP,"rangP/D:rangErr/D");
  statTree->Branch("FitChiInfo", &_Chi, "trajChi/D:MaxLoc/D");
  statTree->Branch("NChadnNhad", &_nhad, "nChad/I:nNhad/I");
  statTree->Branch("hadronMom", &_hadP, "hadP[3]/D");
  statTree->Branch("EnergyTransfer", &_engTrans, "nu/D");
  statTree->Branch("ChadMom", &_chadP, "chadP[4]/D");
  statTree->Branch("NhadMom", &_nhadP, "nhadP[4]/D");
  statTree->Branch("hadEng", &_hadE, "truE/D:recE/D");
  //statTree->Branch("hadDir", &_haddot, "dotProd/D");
  statTree->Branch("NoPlanes", &_plns, "nplanes/I:freeplanes/I");
  statTree->Branch("NoHits", &_nhits, "nhits/I");
  statTree->Branch("HitBreakDown", &_hitType, "nTruMu/I:nInMu/I:nMuInMu/I:nFitN/I:nhad/I");
  statTree->Branch("XPositions", &_XPos, "X[nhits]/D");
  statTree->Branch("YPositions", &_YPos, "Y[nhits]/D");
  statTree->Branch("ZPositions", &_ZPos, "Z[nhits]/D");
  statTree->Branch("EngDeposit", &_Edep, "E[nhits]/D");
  statTree->Branch("MuHits", &_mus, "truMu[nhits]/B");
  statTree->Branch("CandHits", &_cand, "inMu[nhits]/B");
  statTree->Branch("FittedNodes",&_node,"fitNode[nhits]/B");
  statTree->Branch("HadronNodes",&_had,"hadN[nhits]/B");
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
bool MINDplotter::extract_true_particle1(const bhep::event& evt, fitter& Fit) {
//*************************************************************************************

/* sets true particle momentum for calculation and returns a reference
   to the particle */

  //_nuEng = atof( evt.fetch_property("Enu").c_str() ) * GeV;
  _nuEng = evt.fetch_dproperty("nuEnergy") * MeV;
  // for (int jj=0;jj<4;jj++)
//     _hadP[jj] = 0.;

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
      // _hadP[3] = Pospart[iParts]->fetch_dproperty("HadE") * GeV;
    }
  }
  //STUFF ABOVE FOR REDESIGN.!!!!!
  _nhits = Fit.get_nMeas();
  
  for (int iHits = 0;iHits < _nhits;iHits++){

    _XPos[iHits] = Fit.get_meas(iHits)->vector()[0];
    _YPos[iHits] = Fit.get_meas(iHits)->vector()[1];
    _ZPos[iHits] = Fit.get_meas(iHits)->position()[2];

    if (!_patR)
      if ( Fit.get_traj().node(iHits).status("fitted") )
	_hitType[3]++;
  }

  if (count == 0) {
    //cout << "No particles of muon or antimuon type in file" << endl;
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
bool MINDplotter::extract_true_particle2(const bhep::event& evt, fitter& Fit) {
//*************************************************************************************

  bool primNu = false;
  _nuEng = evt.fetch_dproperty("nuEnergy") * MeV;
  
  const vector<bhep::particle*> Pospart = evt.true_particles();

  int count = 0;
  for (int iParts=(int)Pospart.size()-1;iParts >= 0;iParts--){
    
    if ( Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
	 Pospart[iParts]->find_sproperty("PrimaryLepton") ){
      
      if ( Pospart[iParts]->name() == "mu+" ){
	_truPart = Pospart[iParts];
	_Q[0] = 1;
	count++;
      } else if ( Pospart[iParts]->name() == "mu-" ){
	_truPart = Pospart[iParts];
	_Q[0] = -1;
	count++;
      }
    } else if ( Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
		!Pospart[iParts]->find_sproperty("PrimaryLepton") )
      add_to_hads( *Pospart[iParts] );
    // if ( Pospart[iParts]->name() == "mu+" &&
// 	 Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
// 	 evt.fetch_iproperty("nuType") == -14 && count == 0 && !primNu ){
//       _Q[0] = 1;
//       _truPart = Pospart[iParts];
//       count++;
//     } else if ( Pospart[iParts]->name() == "mu-" &&
// 		Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
// 		evt.fetch_iproperty("nuType") == 14 && count == 0 && !primNu ){
//       _Q[0] = -1;
//       _truPart = Pospart[iParts];
//       count++;
//     } else if ( Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" ){

//       if ( Pospart[iParts]->name() == "nu_e" || Pospart[iParts]->name() == "anti_nu_e"
// 		|| Pospart[iParts]->name() == "nu_mu" || Pospart[iParts]->name() == "anti_nu_mu" ){
// 	if ( !primNu ) primNu = true;
// 	else add_to_hads( *Pospart[iParts] );
//       } else add_to_hads( *Pospart[iParts] );
      
//     }
  }

  std::vector<double> hadInf = evt.fetch_dvproperty("had4vec");
  for (int itn = 0;itn < 4;itn++)
    _hadP[itn] = hadInf[itn];

  _nhits = Fit.get_nMeas();
  for (int iHits = 0;iHits < _nhits;iHits++){
    
    _XPos[iHits] = Fit.get_meas(iHits)->vector()[0];
    _YPos[iHits] = Fit.get_meas(iHits)->vector()[1];
    _ZPos[iHits] = Fit.get_meas(iHits)->position()[2];
    _Edep[iHits] = Fit.get_meas(iHits)->get_eng();
    
    if (!_patR)
      if ( Fit.get_traj().node(iHits).status("fitted") )
	_hitType[3]++;
  }
  if (count == 0) {
    //cout << "No particles of muon or antimuon type in file" << endl;
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

void MINDplotter::add_to_hads(const bhep::particle& part){

  if ( part.charge() != 0 ){
    _chadP[0] += part.px();
    _chadP[1] += part.py();
    _chadP[2] += part.pz();
    _chadP[3] += part.e();
    
    _nhad[0]++;
  } else {
    _nhadP[0] += part.px();
    _nhadP[1] += part.py();
    _nhadP[2] += part.pz();
    _nhadP[3] += part.e();
    
    _nhad[1]++;
  }
  // _hadP[0] += part.px();
//   _hadP[1] += part.py();
//   _hadP[2] += part.pz();

//   _hadE[0] += part.e();

}

//*************************************************************************************
void MINDplotter::hadron_direction(fitter& fit) {
//*************************************************************************************
  
  // double normal;
//   EVector fitunit;
//   normal = sqrt(pow(_hadP[0],2)+pow(_hadP[1],2)+pow(_hadP[2],2));
  
//   if ( _nhits >= 2 ){
//     fitunit = fit.get_had_unit();
    
//     _haddot = fitunit[0]*(_hadP[0]/normal)
//       +fitunit[1]*(_hadP[1]/normal)
//       +fitunit[2]*(_hadP[2]/normal);
    
//   } else _haddot = 99;
  
  //_hadE[1] = fit.get_had_eng();
  
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
void MINDplotter::patternStats1(fitter& Fit) {
//****************************************************************************************
  //Event classifier version.
  //_nhits = Fit.get_nMeas();
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

    _engvar[1] = 0;
    for(int ii=0;ii<_hitType[1];ii++)
      _engvar[1] += pow( bhep::double_from_string( Fit.get_traj().node(ii).measurement().name(engDep) ) * GeV - _engTraj/_hitType[1], 2);
    _engvar[1] /= (_hitType[1]-1);
    _engvar[0] = _engTraj/_hitType[1];
  } else {_engvar[1] = -1; _engvar[0] = -1;}
 
  // for (int iclear = _nhits;iclear<300;iclear++){
//     _mus[iclear] = false; _cand[iclear] = false;
//     _node[iclear] = false;
//   }

}

void MINDplotter::patternStats2(fitter& Fit) {
  
  const dict::Key candHit = "inMu";
  const dict::Key hadHit = "inhad";
  int nNode = 0;
  if ( Fit.check_reseed() ) nNode = (int)Fit.get_traj().size()-1;
  bool isMu;
  
  std::vector<cluster*> hits = Fit.get_meas_vec();
  std::vector<cluster*> inMuC;

  for (int iHits = 0;iHits < _nhits;iHits++){
    
    if ( hits[iHits]->get_mu_prop() > 0.8 ){//good number??
      isMu = true;
      _mus[iHits] = true;
      _hitType[0]++;
    }
    else {_mus[iHits] = false; isMu = false;}

    if ( _fail != 7 && hits[iHits]->names().has_key(candHit) ){

      if ( hits[iHits]->name(candHit).compare("True")==0 ){
	_cand[iHits] = true;
	_had[iHits] = false;
	_hitType[1]++;
	_engTraj += hits[iHits]->get_eng() * MeV;
	inMuC.push_back( hits[iHits] );
	if ( isMu ) _hitType[2]++;
	
	if ( Fit.get_traj().node(nNode).status("fitted") && _fail!=1 && _fail<4 ){	
	  _node[iHits] = true; _hitType[3]++; }
	else _node[iHits] = false;
	if ( Fit.check_reseed() ) nNode--;
	else nNode++;
      } else {
	_node[iHits] = false; _cand[iHits] = false;
	if ( hits[iHits]->names().has_key(hadHit) ){
	  _hitType[4]++;
	  _had[iHits] = true;
	} else _had[iHits] = false;
      }

    } else if ( _fail != 7) {
      _node[iHits] = false; _cand[iHits] = false;
      if ( hits[iHits]->names().has_key(hadHit) ){
	_hitType[4]++;
	_had[iHits] = true;
      } else _had[iHits] = false;
    }
  }

  if (_fail != 7){
    _pChi[0] = Fit.get_classifier().get_PatRec_Chis()[0];
    _pChi[1] = Fit.get_classifier().get_PatRec_Chis()[1];
    _pChi[2] = Fit.get_classifier().get_PatRec_Chis()[2];

    _engvar[1] = 0;
    for(int ii=0;ii<_hitType[1];ii++)
      _engvar[1] += pow( inMuC[ii]->get_eng()*MeV - _engTraj/_hitType[1], 2);
    _engvar[1] /= (_hitType[1]-1);
    _engvar[0] = _engTraj/_hitType[1];
  } else {_engvar[1] = -1; _engvar[0] = -1;}

  hits.clear();
  inMuC.clear();
  
}
