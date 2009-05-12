/***********************************************************************
 * Implementation of neutrino event classification class.              *
 * The functions here will take in the measurements associated to an   *
 * event and attempt to reject neutral current and electron nu events  *
 * while identify CC interaction types and performing appropriate muon *
 * extraction.                                                         *
 *                                                                     *
 * Execute will return an integer code based on what type the event    *
 * has been classified as: 0 - NC/e, 1 - CC DIS, 2 - CC Quasi etc.     *
 *                                                                     *
 * Author: Andrew Laing, Feb. 2009                                     *
 ***********************************************************************/


#include <event_classif.h>
#include <recpack/stc_tools.h>
#include <TMath.h>

//**********************************************************************
event_classif::event_classif() {
//**********************************************************************

  

}

//***********************************************************************
event_classif::~event_classif() {
//***********************************************************************


}

//***********************************************************************
bool event_classif::initialize(const bhep::gstore& pstore, bhep::prlevel vlevel,
			       Setup& det, double wFe) {
//***********************************************************************

  m = bhep::messenger( vlevel );
  m.message("++++ Classifier  init  function ++++",bhep::VERBOSE);

  _infoStore = pstore;
  readParam();

  FeWeight = wFe;
  
  //define parameters for cellular automaton.
  man().matching_svc().retrieve_trajectory_finder<CellularAutomaton>("cat")
    .set_max_distance( _max_sep );
  man().matching_svc().retrieve_trajectory_finder<CellularAutomaton>("cat")
    .set_plane_tolerance( _tolerance );
  man().matching_svc().retrieve_trajectory_finder<CellularAutomaton>("cat")
    .set_max_trajs( _max_traj );
  //
  if ( _outLike ){

    _outFileEv = new TFile("LikeOut.root", "recreate");
    _likeTree = new TTree("h1", "Possible Likelihood parameters");

    set_branches();

  }
  //Temp functions to understand cell auto.
  _outFileCell = new TFile("CellOut.root", "recreate" );
  _cellTree = new TTree("h2", "Cellular automaton results");
  set_cell_branches();
  //
  set_extract_properties( det );

  return true;
}

//***********************************************************************
bool event_classif::execute(measurement_vector& hits,
			    Trajectory& muontraj, measurement_vector& hads) {
//***********************************************************************

  m.message("++++ Classifier Execute Function ++++", bhep::VERBOSE);

  bool ok;
  _intType = 0;
  //reset.
  reset();
  
  //Occupancy.
  ok = get_plane_occupancy( hits );
  
  if ( ok && _outLike)
    output_liklihood_info( hits );
  
  /* Code to discriminate between different event types */
  //if identified as CC
  if ( ok )
    ok = chargeCurrent_analysis(hits, muontraj, hads);
  
  return ok;
}

//***********************************************************************
bool event_classif::finalize() {
//***********************************************************************

  if ( _outLike ){
    _outFileEv->Write();
    _outFileEv->Close();
  }

  _outFileCell->Write();
  _outFileCell->Close();

  return true;
}

//***********************************************************************
void event_classif::readParam() {
//***********************************************************************
 
  m.message("++++ readParam function of classifier ++++",bhep::VERBOSE);
  
  model="particle/helix"; 
  
  if ( _infoStore.find_sstore("fitter") )
    kfitter = _infoStore.fetch_sstore("kfitter");
  else kfitter="kalman";

  if ( _infoStore.find_istore("likeli") )
    _outLike = _infoStore.fetch_istore("likeli");
  else _outLike = false;
  
  patRec_maxChi = _infoStore.fetch_dstore("pat_rec_max_chi");
  patRec_max_outliers = _infoStore.fetch_istore("pat_rec_max_outliers");
  max_consec_missed_planes = _infoStore.fetch_istore("max_consec_missed_planes");
  min_seed_hits = _infoStore.fetch_istore("min_seed_hits");
  min_check =  _infoStore.fetch_istore("min_check_nodes");
  min_hits = _infoStore.fetch_istore("low_Pass_hits");
  
  _tolerance = _infoStore.fetch_dstore("pos_res") * cm;
  _max_sep = _infoStore.fetch_dstore("max_sep") * cm;
  _max_traj = _infoStore.fetch_istore("max_traj");
  chi2_max = _infoStore.fetch_dstore("accept_chi");
  max_coincedence = _infoStore.fetch_dstore("max_coincidence");
  
  vfit = _infoStore.fetch_istore("vfit");
  vnav = _infoStore.fetch_istore("vnav");
  vmod = _infoStore.fetch_istore("vmod");
  vmat = _infoStore.fetch_istore("vmat");
  
}

//***********************************************************************
void event_classif::set_extract_properties(Setup& det) {
//***********************************************************************
  
  std::string info[4]={"MUTE","NORMAL","VERBOSE","VVERBOSE"};
  Messenger::Level l0 = Messenger::str(info[vfit]);
  Messenger::Level l1 = Messenger::str(info[vnav]);
  Messenger::Level l2 = Messenger::str(info[vmod]);
  Messenger::Level l3 = Messenger::str(info[vmat]);

  man().fitting_svc().fitter(model).set_verbosity(l0);
  man().fitting_svc().select_fitter(kfitter);

  man().navigation_svc().set_verbosity(l1);
  man().model_svc().model(model).equation().set_verbosity(l2);
  man().model_svc().model(model).propagator().set_verbosity(l2);
  man().model_svc().model(model).tool("noiser/ms").set_verbosity(l2);

  man().geometry_svc().set_zero_length(1e-5 * mm);
  man().geometry_svc().set_infinite_length(1e12 * mm);

  // add the setup to the geometry service
  man().geometry_svc().add_setup("main", det);
  
  // select the setup to be used by the geometry service
  man().geometry_svc().select_setup("main");
  
  man().navigation_svc().navigator(model).set_unique_surface(true);
  
  man().fitting_svc().retrieve_fitter<KalmanFitter>(kfitter,model).
    set_max_local_chi2ndf(patRec_maxChi);

  man().fitting_svc().retrieve_fitter<KalmanFitter>(kfitter,model).
    set_number_allowed_outliers(patRec_max_outliers);

  //Cell auto verbosity.
  man().matching_svc().retrieve_trajectory_finder<CellularAutomaton>("cat")
    .set_verbosity(l3);

}

//***********************************************************************
void event_classif::reset(){
//***********************************************************************
//Resets relevant members at the start of each execute.

  _nplanes = 0;
  _meanOcc = 0;
  _hitsPerPlane.clear();
  _energyPerPlane.clear();
  _planeZ.clear();

}

//***********************************************************************
bool event_classif::get_plane_occupancy(measurement_vector& hits){
//***********************************************************************
//Gets plane occupancies and total plane energies.
//Needs hits in increasing z order.
  m.message("++++ Calculating plane energies and occupancies ++++",bhep::VERBOSE);

  bool ok = true;

  size_t nHits = hits.size();
  
  int count = 0;
  double EngPlane = 0, testZ, curZ;
  size_t hits_used = 0, imeas = 0;
  const dict::Key Edep = "E_dep";
  
  do {
    
    EngPlane = bhep::double_from_string( hits[imeas]->name( Edep ) ) * GeV;
    testZ = hits[imeas]->position()[2];
    count++;
    hits_used++;

    for (size_t i = hits_used;i < nHits;i++) {
      curZ = hits[i]->position()[2];

      if (curZ <= testZ + _tolerance) {

	EngPlane += bhep::double_from_string( hits[i]->name( Edep ) ) * GeV;
	testZ = hits[i]->position()[2];
	count++;
	hits_used++;
	
      } else break;
    }

    _hitsPerPlane.push_back( count );
    _energyPerPlane.push_back( EngPlane );
    _planeZ.push_back( testZ );

    imeas += count;
    _meanOcc += (double)count;

    count = 0;

  } while (hits_used != nHits);

  _nplanes = (int)_hitsPerPlane.size();

  if ( _nplanes == 0 ) return false;
  _meanOcc /= (double)_nplanes;
  
  return ok;
}

//***********************************************************************
bool event_classif::chargeCurrent_analysis(measurement_vector& hits,
					   Trajectory& muontraj, measurement_vector& hads){
//***********************************************************************
  m.message("++++ Performing CC reconstruction ++++",bhep::VERBOSE);

  bool ok = true;
  _recChi = EVector(3,0);
  _recChi[1] = 100000; //SetLarge dummy value for minChiHadron bit.

  _hitIt = hits.end() - 1;

  _vertGuess = exclude_backwards_particle();
  
  //Maybe another which looks for kinks very early (backwards proton quasi?)
  for (_planeIt = _hitsPerPlane.end()-1;_planeIt>=_hitsPerPlane.begin()+min_check;_planeIt--, _hitIt--){
    if ( (*_planeIt) == 1 ){
      muontraj.add_measurement( *(*_hitIt) );
      const dict::Key candHit = "inMu";
      const dict::Key hit_in = "True";
      (*_hitIt)->set_name(candHit, hit_in);
    } else break;
  }
  
  _lastIso = (int)muontraj.nmeas();

  if ( _meanOcc == 1 ){
    
    _intType = 2;
    if ( muontraj.size() !=0 )
      ok = muon_extraction( hits, muontraj, hads);
    else ok = false;
    
  } else {
    
    if ( (int)muontraj.nmeas() < min_seed_hits ) {
      
      ok = invoke_cell_auto( hits, muontraj, hads);
      if ( !ok ) _failType = 4;
      _intType = 5;
      
    } else
      ok = muon_extraction( hits, muontraj, hads);

  }
  cout << "here3?"<<endl;
  return ok;
}

//***********************************************************************
int event_classif::exclude_backwards_particle(){
//***********************************************************************
//Try to exclude any backwards had by looking for increases in
//occupancy from the start of the event.
  vector<int>::iterator measIt;
  int excluded_hits = 0;

  for (measIt = _hitsPerPlane.begin();measIt != _hitsPerPlane.end();measIt++){

    if (measIt == _hitsPerPlane.begin()) {
      excluded_hits += (*measIt);
      continue;
    }

    if ( (*measIt) > (*(measIt - 1)) ) break;
    else excluded_hits += (*measIt);

  }

  if (measIt == _hitsPerPlane.end())
    excluded_hits = 0;

  return excluded_hits;
}

//***********************************************************************
bool event_classif::muon_extraction(measurement_vector& hits,
				    Trajectory& muontraj, measurement_vector& hads) {
//***********************************************************************
//Decide on appropriate pattern recognition and perform it.
//Are there ways to avoid doing full extraction algorithm.

  bool ok;
  
  if (_vertGuess != 0)
    for (int i = 0;i < _vertGuess;i++)
      hads.push_back( hits[i] );

  State patternSeed;
  
  ok = get_patternRec_seed( patternSeed, muontraj, hits);
  
  if ( ok )
    ok = perform_muon_extraction( patternSeed, hits, muontraj, hads);
  
  if ( ok )
    _seedState = patternSeed;

  return ok;
}

//***********************************************************************
bool event_classif::get_patternRec_seed(State& seed, Trajectory& muontraj,
					measurement_vector& hits) {
//***********************************************************************

  EVector V(6,0); EVector V2(1,0);
  EMatrix M(6,6,0); EMatrix M2(1,1,0);
  double Xtent;

  V[0] = muontraj.nodes()[0]->measurement().vector()[0];
  V[1] = muontraj.nodes()[0]->measurement().vector()[1];
  V[2] = muontraj.nodes()[0]->measurement().position()[2];
  
  //direction
  fit_parabola( V, muontraj);
  
  //Momentum. Estimate from empirical extent function.
  if ( hits.size() != 0 ){

    Xtent = hits[hits.size()-1]->position()[2]
      - hits[_vertGuess]->position()[2];

  } else {
    
    Xtent = muontraj.nodes()[(int)muontraj.nmeas()-1]->measurement().position()[2]
      - muontraj.nodes()[0]->measurement().position()[2];

  }
  
  double pSeed = 668 + 1.06*Xtent; //estimate in MeV, log for fit.

  set_de_dx( pSeed/GeV );

  V[5] = 1./pSeed;

  //Errors
  M[0][0] = M[1][1] = 15.*cm*cm;
  M[2][2] = EGeo::zero_cov()/2;
  M[3][3] = M[4][4] = 1.5;
  M[5][5] = pow(V[5],2)*4;

  //Sense
  V2[0] = 1;
  
  //Seedstate fit properties
  seed.set_name(RP::particle_helix);
  seed.set_name(RP::representation,RP::slopes_z); 
  seed.set_hv(RP::sense,HyperVector(V2,M2));
  seed.set_hv(HyperVector(V,M));
  
  man().model_svc().model(RP::particle_helix).representation(RP::slopes_z)
    .convert(seed,RP::default_rep);
  cout << "here5?"<<endl;
  bool ok = perform_kalman_fit( seed, muontraj);
  if ( !ok )
    _failType = 5;
  cout << "here9q?"<<endl;
  return ok;
}

//***********************************************************************
void event_classif::fit_parabola(EVector& vec, Trajectory& track) {
//***********************************************************************

  size_t nMeas = track.nmeas();

  if (nMeas > 3) nMeas = 3;

  double x[(const int)nMeas], y[(const int)nMeas], z[(const int)nMeas];

  for (int iMeas = nMeas-1;iMeas >= 0;iMeas--){

    x[iMeas] = track.nodes()[iMeas]->measurement().vector()[0];
    y[iMeas] = track.nodes()[iMeas]->measurement().vector()[1];
    z[iMeas] = track.nodes()[iMeas]->measurement().position()[2];

  }
  
  TGraph *gr1 = new TGraph((const int)nMeas, z, x);
  TGraph *gr2 = new TGraph((const int)nMeas, z, y);

  TF1 *fun = new TF1("parfit","[0]+[1]*x",-3,3);
  fun->SetParameters(0.,0.001);
  
  gr1->Fit("parfit", "QN");
  vec[3] = fun->GetParameter(1);

  fun->SetParameters(0.,0.001);
  gr2->Fit("parfit", "QN");
  vec[4] = fun->GetParameter(1);
  
  delete gr1;
  delete gr2;
  delete fun;

}

//***********************************************************************
void event_classif::set_de_dx(double mom){
//***********************************************************************
  
  double de_dx = -( 12.37 * FeWeight * pow( mom, 0.099) 
		    + (1 - FeWeight) * 2.16 * pow( mom, 0.075) );
  de_dx *= MeV/cm;
  
  man().geometry_svc().setup().set_volume_property_to_sons("mother","de_dx",de_dx);

}

//***********************************************************************
bool event_classif::perform_kalman_fit(State& seed, Trajectory& track) {
//***********************************************************************
  cout << "here6?"<<endl;
  bool ok = man().fitting_svc().fit(seed, track);
  cout << "here7?"<<endl;
  if (ok)
    seed = track.state(track.first_fitted_node());
  cout << "here8?"<<endl;
  return ok;

}

//***********************************************************************
bool event_classif::perform_muon_extraction(const State& seed, measurement_vector& hits,
					    Trajectory& muontraj, measurement_vector& hads) {
//***********************************************************************
//Loop through multiple occupancy planes finding the best match to the muon
//traj and filtering it into the trajectory.
  bool ok;
  long ChiMin;
  int nConsecHole = 0;

  while (_hitIt >= hits.begin() + _vertGuess) {

    double Chi2[(const int)(*_planeIt)];

    for (int iMat = (*_planeIt)-1;iMat >= 0;iMat--, _hitIt--){
      ok = man().matching_svc().match_trajectory_measurement(muontraj, (*(*_hitIt)), Chi2[iMat]);
      if ( !ok )
	Chi2[iMat] = 999999999;
    }

    ChiMin = TMath::LocMin( (const int)(*_planeIt), Chi2);

    for (int iFil = 0;iFil < (*_planeIt);iFil++){

      if ( iFil == (int)ChiMin) {

	ok = man().fitting_svc().filter(*(*(_hitIt+iFil+1)), seed, muontraj);

	if ( ok ) {
	  const dict::Key candHit = "inMu";
	  const dict::Key hit_in = "True";
	  (*(_hitIt+iFil+1))->set_name(candHit, hit_in);

	  if ( (*(_hitIt+iFil+1))->name("MotherParticle").compare("mu+")==0 ||
	       (*(_hitIt+iFil+1))->name("MotherParticle").compare("mu-")==0 )
	    _recChi[0] = TMath::Max(Chi2[iFil], _recChi[0]);
	  else
	    _recChi[1] = TMath::Min(Chi2[iFil], _recChi[1]);

	  if ( nConsecHole > _recChi[2] )
	    _recChi[2] = nConsecHole;

	  nConsecHole = 0;

	} else nConsecHole++;

      } else {

	hads.push_back( (*(_hitIt+iFil+1)) );

      }
    }

    _planeIt--;

    if ( nConsecHole > max_consec_missed_planes ) {
      _failType = 6;
      return false;
    }

  }

  return true;
}

//***********************************************************************
bool event_classif::invoke_cell_auto(measurement_vector& hits,
				     Trajectory& muontraj, measurement_vector& hads){
//***********************************************************************
//uses cellular automaton to try and retrieve more complicated events.
  std::cout << "Performing Cellular Automaton analysis " <<_nplanes<< std::endl;
  bool ok;

  if ( _nplanes < min_hits ) return false;

  std::vector<Trajectory*> trajs;
  
  ok = man().matching_svc().find_trajectories( hits, trajs);
  
  if ( !ok || trajs.size() == 0) return false;
   
  output_results_tree( hits, trajs );
  
  if ( trajs.size() == 1 ) {
  
    muontraj.reset();

    muontraj = *trajs[0];

    sort_hits( hits, muontraj, hads);

    ok = true;

  } else {

    ok = sort_trajs( muontraj, trajs);

    if ( ok )
      sort_hits( hits, muontraj, hads);

  }

  if ( ok )
    delete_bad_trajs( muontraj, trajs );
  else stc_tools::destroy( trajs );

  std::cout << "End of Cellular automaton" << std::endl;
  return ok;
}

//***********************************************************************
void event_classif::sort_hits(measurement_vector& hits,
			      Trajectory& muontraj, measurement_vector& hads){
//***********************************************************************
  
  bool inTraj;
  vector<Node*>::iterator inTrIt;
  inTrIt = muontraj.nodes().begin();

  const dict::Key candHit = "inMu";
  const dict::Key not_in = "False";
  const dict::Key hit_in = "True";

  for (_hitIt = hits.begin();_hitIt!=hits.end();_hitIt++){

    inTraj = false;
    if ( (*_hitIt)->names().has_key(candHit) )
      (*_hitIt)->set_name(candHit, not_in);

    for (inTrIt = muontraj.nodes().begin();inTrIt != muontraj.nodes().end();inTrIt++){

      if ( (*_hitIt)->position() == (*inTrIt)->measurement().position() ) {

	(*_hitIt)->set_name(candHit, hit_in);

	inTraj = true;
      }

    }

    if ( !inTraj )
      hads.push_back( (*_hitIt) );

  }

}

//***********************************************************************
void event_classif::delete_bad_trajs(const Trajectory& muontraj,
				     vector<Trajectory*>& trajs){
//***********************************************************************

  vector<Trajectory*>::iterator dIt = trajs.begin();
  bool found = false;
  const dict::Key traj_index = "traj_no";

  while ( !found && dIt != trajs.end() ){
    if ( muontraj.quality( traj_index ) == (*dIt)->quality( traj_index ) ){
      found = true;

      trajs.erase( dIt );
    }
    dIt++;
  }
  stc_tools::destroy( trajs );  

}

//***********************************************************************
bool event_classif::sort_trajs(Trajectory& muontraj, vector<Trajectory*>& trajs){
//***********************************************************************
//Reject possible trajectories based on number of hits, chi2, hits in common.
  
  bool ok;
  std::vector<Trajectory*> trajs2;
  
  ok = reject_small( trajs, trajs2);

  if ( ok ){
    if ( trajs2.size() == 1 ){
      
      muontraj.reset();
      
      muontraj = *trajs2[0];
      
      return ok;

    } else {

      std::vector<Trajectory*> trajs3;
      std::vector<double> Chis;

      ok = reject_high( trajs2, trajs3);

      trajs2.clear();

      if ( ok ){
	if ( trajs2.size() == 1 ){
  
	  muontraj.reset();
      
	  muontraj = *trajs2[0];
      
	  return ok;

	} else {

	  ok = reject_final( trajs3, muontraj);

	}

      }

    }

  }

  return ok;
}

//***********************************************************************
bool event_classif::reject_small(vector<Trajectory*>& trajs,
				 vector<Trajectory*>& trajs2){
//***********************************************************************
  
  vector<Trajectory*>::iterator it1;

  for (it1 = trajs.begin();it1 != trajs.end();it1++){

    if ( (int)(*it1)->nmeas() >= min_hits )
      trajs2.push_back( (*it1) );

  }

  if ( trajs2.size() == 0 ) return false;

  return true;
}

//***********************************************************************
bool event_classif::reject_high(vector<Trajectory*>& trajs,
				vector<Trajectory*>& trajs2){
//***********************************************************************
  
  bool ok1, ok2 = false;
  double traj_no;
  const dict::Key traj_index = "traj_no";

  vector<Trajectory*>::iterator it1;
  State temp_seed;
  measurement_vector dummy_hits;
  double momErr;
  
  for (it1 = trajs.begin();it1 != trajs.end();it1++){
    
    Trajectory& temp_traj = *(*it1);
    traj_no = temp_traj.quality( traj_index );
    
    temp_traj.sort_nodes( -1 );
    
    ok1 = get_patternRec_seed( temp_seed, temp_traj, dummy_hits);
    temp_traj.set_quality( traj_index, traj_no );
    
    if ( ok1 && temp_traj.quality() < chi2_max ){

      (*it1)->set_quality( temp_traj.quality() );
      
      const dict::Key relErr = "relErr";
      momErr = temp_seed.matrix()[6][6] / temp_seed.vector()[6];
      (*it1)->set_quality( relErr, momErr);

      trajs2.push_back( (*it1) );

      ok2 = true;
    }
    
  }

  sort( trajs2.begin(), trajs2.end(), chiSorter() );

  return ok2;
}

//***********************************************************************
bool event_classif::reject_final(vector<Trajectory*>& trajs,
				 Trajectory& muontraj){
//***********************************************************************
//Accept lowest local chi and then any others with few coincidences
//with this trajectory then choose the longest/lowest error track.
  
  std::vector<Trajectory*> temp_trajs;
  temp_trajs.push_back( trajs[0] );

  vector<Trajectory*>::iterator it1, it2;

  double coincidence;

  for (it1 = trajs.begin()+1;it1 != trajs.end();it1++){
    
    vector<Node*>& node2 = (*it1)->nodes();
    it2 = temp_trajs.begin();

    while ( it2 != temp_trajs.end() && coincidence < max_coincedence ){

      vector<Node*>& node1 = (*it2)->nodes();

      coincidence = compare_nodes( node1, node2);

      it2++;
    }
    
    if ( coincidence < max_coincedence )
      temp_trajs.push_back( (*it1) );

  }
  trajs.clear();

  if ( temp_trajs.size() == 1 ) {
    
    muontraj.reset();

    muontraj = *temp_trajs[0];

  } else {

    select_trajectory( temp_trajs, muontraj);

  }
  temp_trajs.clear();

  return true;
}

//***********************************************************************
double event_classif::compare_nodes(const vector<Node*>& n1,
				    const vector<Node*>& n2){
//***********************************************************************
  
  std::vector<Node*>::const_iterator nIt1, nIt2;
  double counter = 0;

  for (nIt1 = n1.begin();nIt1 != n1.end();nIt1++){
    for (nIt2 = n2.begin();nIt2 != n2.end();nIt2++){

      if ( (*nIt1)->measurement().position() == (*nIt2)->measurement().position() )
	counter++;

    }
  }

  return counter / (double)n2.size();
}

//***********************************************************************
void event_classif::select_trajectory(vector<Trajectory*>& trajs,
				      Trajectory& muontraj){
//***********************************************************************
  
  int length[(const int)trajs.size()];
  double Err[(const int)trajs.size()];
  
  muontraj.reset();

  const dict::Key relErr = "relErr";
  
  for (int it1 = 0;it1 < (int)trajs.size();it1++){

    length[it1] = (int)trajs[it1]->nmeas();
    Err[it1] = trajs[it1]->quality( relErr );
    
  }

  int longest = (int)TMath::LocMax( (int)trajs.size(), length);
  int safest = (int)TMath::LocMin( (int)trajs.size(), Err);

  if ( longest == safest )
    muontraj = *trajs[longest];
  else {

    if ( length[longest] - length[safest] > 5)
      muontraj = *trajs[longest];
    else muontraj = *trajs[safest];

  }

}

//***********************************************************************
void event_classif::set_branches(){
//***********************************************************************
  
  _likeTree->Branch("nhits", &_nhit, "nhits/I");
  _likeTree->Branch("VisEng", &_visEng, "visEng/D");
  _likeTree->Branch("nplanes", &_nplanes, "nplanes/I");
  _likeTree->Branch("freePlanes", &_freeplanes, "freeplanes/I");
  _likeTree->Branch("occupancy", &_hitsPerPlane, "occ[nplanes]/I");
  _likeTree->Branch("planeEnergy", &_energyPerPlane, "plEng[nplanes]/D");
  _likeTree->Branch("planeZ", &_planeZ, "planeZ[nplanes]/D");
  _likeTree->Branch("nhitTraj", &_trajhit, "nhitTraj/I");
  _likeTree->Branch("TrajPur", &_trajpur, "trajpur/D");
  _likeTree->Branch("EngTraj", &_trajEng, "engTraj/D");
  _likeTree->Branch("EngTrajPlane", &_trajEngPlan, "engTrajPlane[nplanes]/D");

}

//***********************************************************************
void event_classif::output_liklihood_info(const measurement_vector& hits){
//***********************************************************************
  
  bool multFound = false;
  const dict::Key Edep = "E_dep";
  const dict::Key candHit = "inMu";

  _nhit = 0;
  _freeplanes = 0;
  _visEng = 0;
  _trajhit = 0;
  _trajpur = 0;
  _trajEng = 0;

  vector<int>::iterator itLike;
  vector<double>::iterator itEngL = _energyPerPlane.end()-1;
  measurement_vector::const_iterator hitIt3;

  int counter = _nplanes - 1;

  for (itLike = _hitsPerPlane.end()-1;itLike >= _hitsPerPlane.begin();itLike--, itEngL--){

    _nhit += (*itLike);
    _visEng += (*itEngL);

    if ( !multFound && (*itLike) == 1)
      _freeplanes++;
    else multFound = true;

    counter--;

  }
  //????
  for (hitIt3 = hits.begin();hitIt3 != hits.end();hitIt3++){
    if ( (*hitIt3)->names().has_key( candHit ) ){
      if ( (*hitIt3)->name( candHit ).compare("True") == 0 ){
	_trajhit++;
	_trajEng += bhep::double_from_string( (*hitIt3)->name( Edep ) ) * GeV;
	if ( (*hitIt3)->name("MotherParticle").compare("mu+") == 0
	     || (*hitIt3)->name("MotherParticle").compare("mu-") == 0 )
	  _trajpur++;
      }
    }
  }
  _trajpur /= _trajhit;
  
  _likeTree->Fill();

}
//***********************************************************************
void event_classif::set_cell_branches(){
//***********************************************************************

  _cellTree->Branch("Nhits",&_nhit, "nhits/I");
  _cellTree->Branch("XPositions", &_XPos, "X[nhits]/D");
  _cellTree->Branch("YPositions", &_YPos, "Y[nhits]/D");
  _cellTree->Branch("ZPositions", &_ZPos, "Z[nhits]/D");
  _cellTree->Branch("Ntraj", &_ntraj,"ntraj/I");
  _cellTree->Branch("TrajInfo", &_trInfo,"ntrajhit[ntraj]/D");
  _cellTree->Branch("TrajInfo2", &_trhigh,"highneigh[ntraj]/D");
  _cellTree->Branch("TrajInfo3", &_trInd,"highindex[ntraj]/D");
  _cellTree->Branch("trajHits",&_trajHit,"hitintraj[ntraj][nhit]/I");

}

//***********************************************************************
void event_classif::output_results_tree(const measurement_vector& hits,
					const std::vector<Trajectory*>& trajs){
//***********************************************************************
//Temporary functions to study cellular automaton.
  
  _nhit = (int)hits.size();
  _ntraj = (int)trajs.size();

  const dict::Key max_neigh_index = "max_neigh_pos";
  const dict::Key max_neigh_val = "max_neigh";
  
  measurement_vector::const_iterator hitIt2;
  std::vector<Trajectory*>::const_iterator trIt;
  std::vector<Node*>::iterator noIt;
  
  int hitcount = 0, trajCount = 0;
  
  for (hitIt2 = hits.begin();hitIt2 != hits.end();hitIt2++){

    _XPos[hitcount] = (*hitIt2)->vector()[0];
    _YPos[hitcount] = (*hitIt2)->vector()[1];
    _ZPos[hitcount] = (*hitIt2)->position()[2];
    
    hitcount++;
  }

  for (trIt = trajs.begin();trIt != trajs.end();trIt++){
    
    _trInfo[trajCount] = (*trIt)->nmeas();
    _trhigh[trajCount] = (*trIt)->quality( max_neigh_val );
    _trInd[trajCount] = (*trIt)->quality( max_neigh_index );

    vector<Node*> nodesV = (*trIt)->nodes();
    hitcount = 0;
    
    for (hitIt2 = hits.begin();hitIt2 != hits.end();hitIt2++){
      for (noIt = nodesV.begin();noIt != nodesV.end();noIt++){
	if ( (*hitIt2)->position() == (*noIt)->measurement().position() ){
	  _trajHit[trajCount][hitcount] = 1;
	} else {
	  _trajHit[trajCount][hitcount] = 0;
	}
      }
      hitcount++;
    }
    trajCount++;
  }

  _cellTree->Fill();

}
