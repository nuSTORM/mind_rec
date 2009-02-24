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
#include <TMath.h>

//**********************************************************************
event_classif::event_classif(const bhep::gstore& pstore,
			     bhep::prlevel vlevel) {
//**********************************************************************

  _infoStore = pstore;
  
  m = bhep::messenger( vlevel );

}

//***********************************************************************
event_classif::~event_classif() {
//***********************************************************************


}

//***********************************************************************
bool event_classif::initialize(Setup& det) {
//***********************************************************************

  m.message("++++ Classifier  init  function ++++",bhep::VERBOSE);

  readParam();

  set_extract_properties( det );

  return true;
}

//***********************************************************************
bool event_classif::execute(measurement_vector& hits,
			    Trajectory& muontraj, measurement_vector& hads) {
//***********************************************************************

  m.message("++++ Classifier Execute Function ++++", bhep::VERBOSE);

  bool ok;

  //reset.
  reset();

  //Occupancy.
  ok = get_plane_occupancy( hits );

  /* Code to discriminate between different event types */
  //if identified as CC (should this be a switch?)
  if ( ok )
    ok = chargeCurrent_analysis(hits, muontraj, hads);

  return ok;
}

//***********************************************************************
bool event_classif::finalize() {
//***********************************************************************

  return true;
}

//***********************************************************************
void event_classif::readParam() {
//***********************************************************************
 
  m.message("++++ readParam function of classifier ++++",bhep::VERBOSE);
        
    model="particle/helix"; 
    
    if (_infoStore.find_sstore("fitter"))
      kfitter = _infoStore.fetch_sstore("kfitter");
    else kfitter="kalman";
    
    patRec_maxChi = _infoStore.fetch_dstore("pat_rec_max_chi");
    patRec_max_outliers = _infoStore.fetch_istore("pat_rec_max_outliers");
    max_consec_missed_planes = _infoStore.fetch_istore("max_consec_missed_planes");

    _tolerance = _infoStore.fetch_dstore("pos_res") * cm;
    
    vfit = _infoStore.fetch_istore("vfit");
    vnav = _infoStore.fetch_istore("vnav");
    vmod = _infoStore.fetch_istore("vmod");

}

//***********************************************************************
void event_classif::set_extract_properties(Setup& det) {
//***********************************************************************

  std::string info[4]={"MUTE","NORMAL","VERBOSE","VVERBOSE"};
  Messenger::Level l0 = Messenger::str(info[vfit]);
  Messenger::Level l1 = Messenger::str(info[vnav]);
  Messenger::Level l2 = Messenger::str(info[vmod]);

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

}

//***********************************************************************
void event_classif::reset(){
//***********************************************************************
//Resets relevant members at the start of each execute.

  _nplanes = 0;
  _meanOcc = 0;
  _hitsPerPlane.clear();
  _energyPerPlane.clear();

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

    EngPlane += bhep::double_from_string( hits[imeas]->name( Edep ) ) * GeV;
    testZ = hits[imeas]->surface().position()[2];
    count++;
    hits_used++;

    for (size_t i = hits_used;i < nHits;i++) {
      curZ = hits[i]->surface().position()[2];

      if (curZ <= testZ + _tolerance) {

	EngPlane += bhep::double_from_string( hits[i]->name( Edep ) ) * GeV;
	testZ = hits[i]->surface().position()[2];
	count++;
	hits_used++;
	
      } else break;
    }

    _hitsPerPlane.push_back( count );
    _energyPerPlane.push_back( EngPlane );

    imeas += count;
    _meanOcc += (double)count;

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

  if ( _meanOcc == 1 ){
    _intType = 2; return ok; }//free muon if positive CC ident, no pat rec required.
  //could have backwards proton, how to exclude?

  _hitIt = hits.end() - 1;

  _vertGuess = exclude_backwards_particle();
  //Maybe another which looks for kinks very early (backwards proton quasi?)
  for (_planeIt = _hitsPerPlane.end()-1;_planeIt>=_hitsPerPlane.begin();_planeIt--, _hitIt--){
    if ( (*_planeIt) == 1 ){
      muontraj.add_measurement( *(*_hitIt) );
      const dict::Key candHit = "inMu";
      const dict::Key hit_in = "True";
      (*_hitIt)->set_name(candHit, hit_in);
    } else break;
  }
 
  ok = muon_extraction( hits, muontraj, hads);

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
  

  return ok;
}

//***********************************************************************
bool event_classif::get_patternRec_seed(State& seed, Trajectory& muontraj,
					measurement_vector& hits) {
//***********************************************************************

  EVector V(6,0); EVector V2(1,0);
  EMatrix M(6,6,0); EMatrix M2(1,1,0);

  V[0] = muontraj.nodes()[0]->measurement().vector()[0];
  V[1] = muontraj.nodes()[0]->measurement().vector()[1];
  V[2] = muontraj.nodes()[0]->measurement().surface().position()[2];

  //direction
  fit_parabola( V, muontraj);

  //Momentum. Estimate from empirical extent function.
  double Xtent = hits[hits.size()-1]->surface().position()[2]
    - hits[_vertGuess]->surface().position()[2];

  double pSeed = 668 + 1.06*Xtent; //estimate in MeV, log for fit.
  double de_dx = -(2.08 + 10.3*pow(pSeed/GeV, 0.116)) * MeV/cm;

  V[5] = 1./pSeed;
  man().geometry_svc().setup().set_volume_property_to_sons("mother","de_dx",de_dx);

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

  bool ok = perform_kalman_fit( seed, muontraj);

  return ok;
}

//***********************************************************************
void event_classif::fit_parabola(EVector& vec, Trajectory& track) {
//***********************************************************************

  size_t nMeas = track.nmeas();

  if (nMeas > 10) nMeas = 10;

  double x[(const int)nMeas], y[(const int)nMeas], z[(const int)nMeas];

  for (int iMeas = nMeas-1;iMeas >= 0;iMeas--){

    x[iMeas] = track.nodes()[iMeas]->measurement().vector()[0];
    y[iMeas] = track.nodes()[iMeas]->measurement().vector()[1];
    z[iMeas] = track.nodes()[iMeas]->measurement().surface().position()[2];

  }

  TGraph *gr1 = new TGraph((const int)nMeas, z, x);
  TGraph *gr2 = new TGraph((const int)nMeas, z, y);

  TF1 *fun = new TF1("parfit","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)",-3,3);
  fun->SetParameters(0.,0.,0.001,0.0001,0.0001);

  gr1->Fit("parfit", "QN");
  vec[3] = fun->GetParameter(1) + 2*fun->GetParameter(2)*vec[0]
    + 3*fun->GetParameter(3)*pow(vec[0],2) + 4*fun->GetParameter(4)*pow(vec[0],3);

  gr2->Fit("parfit", "QN");
  vec[4] = fun->GetParameter(1) + 2*fun->GetParameter(2)*vec[1]
    + 3*fun->GetParameter(3)*pow(vec[1],2) + 4*fun->GetParameter(4)*pow(vec[1],3);

}

//***********************************************************************
bool event_classif::perform_kalman_fit(State& seed, Trajectory& track) {
//***********************************************************************

  bool ok = man().fitting_svc().fit(seed, track);
  
  if (ok)
    seed = track.state(track.first_fitted_node());

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
	} else cout << "Filter failed"<<endl;

      } else {

	hads.push_back( (*(_hitIt+iFil+1)) );

      }
    }

    _planeIt--;
  }

  return true;
}
