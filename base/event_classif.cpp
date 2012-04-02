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
// void event_classif::initialize(const bhep::gstore& pstore, bhep::prlevel vlevel,
// 			       Setup& det, double wFe) {
void event_classif::initialize(const bhep::gstore& pstore, bhep::prlevel vlevel, double wFe) {
//***********************************************************************

  m = bhep::messenger( vlevel );
  m.message("++++ Classifier  init  function ++++",bhep::VERBOSE);

  _infoStore = pstore;
  readParam();

  _voxEdge = _infoStore.fetch_dstore("rec_boxX") * bhep::cm;

  FeWeight = wFe;

  if ( _outLike ){
    string likeFile = pstore.fetch_sstore("like_file");
    _outFileEv = new TFile(likeFile.c_str(), "recreate");
    _likeTree = new TTree("h1", "Possible Likelihood parameters");

    set_branches();

  }

}

//***********************************************************************
bool event_classif::execute(vector<cluster*>& hits,
			    Trajectory& muontraj, vector<cluster*>& hads) {
//***********************************************************************

  m.message("++++ Classifier Execute Function ++++", bhep::VERBOSE);

  bool ok;
  _intType = 0;
  //reset.
  reset();
  
  //Occupancy.
  ok = get_plane_occupancy( hits );
  
  if ( _outLike )
    output_liklihood_info( hits );
  
  /* Code to discriminate between different event types */
  //if identified as CC
  if ( ok )
    ok = chargeCurrent_analysis(hits, muontraj, hads);
  
  if ( ok && _outLike)
    traj_like( hits, muontraj );

  if ( _outLike ) out_like();
  
  return ok;
}

//***********************************************************************
void event_classif::finalize() {
//***********************************************************************

  if ( _outLike ){
    _outFileEv->Write();
    _outFileEv->Close();
  }

}

//***********************************************************************
void event_classif::readParam() {
//***********************************************************************
 
  m.message("++++ readParam function of classifier ++++",bhep::VERBOSE);
  
  if ( _infoStore.find_istore("likeli") )
    _outLike = _infoStore.fetch_istore("likeli");
  else _outLike = false;
  
  max_consec_missed_planes = _infoStore.fetch_istore("max_consec_missed_planes");
  min_plane_prop = _infoStore.fetch_dstore("min_use_prop");
  min_seed_hits = _infoStore.fetch_istore("min_seed_hits");
  min_check =  _infoStore.fetch_istore("min_check_nodes");
  min_hits = _infoStore.fetch_istore("low_Pass_hits");
  
  _tolerance = _infoStore.fetch_dstore("pos_res") * bhep::cm;
  chi2_max = _infoStore.fetch_dstore("accept_chi");
  max_coincedence = _infoStore.fetch_dstore("max_coincidence");

  _pieceLength = _infoStore.fetch_dstore("widthI") * bhep::cm
    + _infoStore.fetch_dstore("widthS") * _infoStore.fetch_istore("nplane") * bhep::cm
    + _infoStore.fetch_dstore("widthA") * (_infoStore.fetch_istore("nplane")+1) * bhep::cm;
  
  _maxBlobSkip = _infoStore.fetch_dstore("maxBlobSkip");
  _minBlobOcc = _infoStore.fetch_dstore("minBlobOcc");

  _detX = _infoStore.fetch_dstore("MIND_x") * bhep::m;
  _detY = _infoStore.fetch_dstore("MIND_y") * bhep::m;

  if(_infoStore.find_dstore("WLSatten"))
    _WLSAtten = _infoStore.fetch_dstore("WLSatten");
  else						
    _WLSAtten = 5000.;


}

//***********************************************************************
void event_classif::reset(){
//***********************************************************************
//Resets relevant members at the start of each execute.

  _nplanes = 0;
  _meanOcc = 0;
  badplanes = 0;
  _longestSingle = 0;
  _endProj = false;
  _hitsPerPlane.clear();
  _energyPerPlane.clear();
  _planeZ.clear();

}

//***********************************************************************
bool event_classif::get_plane_occupancy(vector<cluster*>& hits){
//***********************************************************************
//Gets plane occupancies and total plane energies.
//Needs hits in increasing z order.
  m.message("++++ Calculating plane energies and occupancies ++++",bhep::VERBOSE);

  bool ok = true;

  size_t nHits = hits.size();
  
  int count = 0, single_count = 0;
  double EngPlane = 0, testZ, curZ;
  size_t hits_used = 0, imeas = 0;
  _longestSingle = 0;

  if ( !_outLike ) _visEng = 0;
  
  do {
    
    //EngPlane = bhep::double_from_string( hits[imeas]->name( Edep ) ) * bhep::GeV;
    EngPlane = correctEdep(hits[imeas]->get_eng()*bhep::MeV, 
			   hits[imeas]->position()[0], hits[imeas]->position()[1]);
    testZ = hits[imeas]->position()[2];
    count++;
    hits_used++;

    for (size_t i = hits_used;i < nHits;i++) {
      curZ = hits[i]->position()[2];

      if (curZ <= testZ + _tolerance) {
	double x = hits[i]->position()[0];
	double y = hits[i]->position()[1];
	EngPlane += correctEdep(hits[i]->get_eng() * bhep::MeV, x, y);
	testZ = hits[i]->position()[2];
	count++;
	hits_used++;
	
      } else break;
    }

    _hitsPerPlane.push_back( count );
    _energyPerPlane.push_back( EngPlane );
    _planeZ.push_back( testZ );

    if ( count == 1 )
      single_count++;
    else {
      if ( single_count >= _longestSingle ){
	_longestSingle = single_count;
	_endLongSing = imeas - 1;
	_endLongPlane = (int)_hitsPerPlane.size()-2;
      }
      single_count = 0;
    }
    if ( imeas == nHits-1 && single_count != 0 && count == 1 ){
      if ( single_count >= _longestSingle ){
	_longestSingle = single_count;
	_endLongSing = imeas;
	_endLongPlane = (int)_hitsPerPlane.size()-1;
      }
    }
    if ( !_outLike ) _visEng += EngPlane;

    imeas += count;
    _meanOcc += (double)count;

    count = 0;

  } while (hits_used != nHits);

  _nplanes = (int)_hitsPerPlane.size();

  if ( _nplanes == 0 ) return false;
  _meanOcc /= (double)_nplanes;

  if ( _longestSingle >= min_seed_hits && _endLongPlane > 0.5*(double)_nplanes )
    assess_event( hits );
  
  return ok;
}

//***********************************************************************
void event_classif::assess_event(vector<cluster*>& hits)
//***********************************************************************
{
  //Assess the event 'free' sections to look for a viable
  //muon and tidy up the endpoint if neccessary.
  double dispTrans[2], endMeanOcc = 0, skipSize = 0;
  unsigned int evEnd = hits.size() - 1;
  vector<int>::iterator aIt;
  const dict::Key hit_in = "True";
  const dict::Key skipped = "skipped";
  
  if ( _endLongSing == (int)evEnd ){
    
    dispTrans[0] = sqrt( pow( hits[evEnd]->vector()[0]-hits[evEnd-1]->vector()[0], 2)
		      + pow( hits[evEnd]->vector()[1]-hits[evEnd-1]->vector()[1], 2) );
    dispTrans[1] = sqrt( pow( hits[evEnd-1]->vector()[0]-hits[evEnd-2]->vector()[0], 2)
			 + pow( hits[evEnd-1]->vector()[1]-hits[evEnd-2]->vector()[1], 2) );
    //Just in case there is a bad hit at the endpoint.
    if ( dispTrans[0] > _voxEdge*10 && dispTrans[1] <= _voxEdge*10 ){
      hits[evEnd]->set_name(skipped,hit_in);
      badplanes = 1;
      _hitIt = hits.end() - 2;
      _planeIt = _hitsPerPlane.end() - 2;
      
    } else if ( dispTrans[1] > _voxEdge*10 ){
      hits[evEnd]->set_name(skipped,hit_in);
      hits[evEnd-1]->set_name(skipped,hit_in);
      badplanes = 2;
      _hitIt = hits.end() - 3;
      _planeIt = _hitsPerPlane.end() - 3;
      
    } else {
      _hitIt = hits.end() - 1;
      _planeIt = _hitsPerPlane.end() - 1;
    }
    
  } else {
    //A more complicated endpoint.
    for ( aIt = _hitsPerPlane.end()-1;aIt != _hitsPerPlane.begin()+_endLongPlane;aIt--){
      endMeanOcc += (double)(*aIt);
      skipSize++;
    }
    endMeanOcc /= (double)skipSize;
    if ( endMeanOcc > _minBlobOcc ){
      
      if ( skipSize/(double)_hitsPerPlane.size() < _maxBlobSkip
	   || _longestSingle >= 2*min_seed_hits ){
	
	_hitIt = hits.end() - 1;
	while ( _hitIt != hits.begin()+_endLongSing ){
	  (*_hitIt)->set_name(skipped,hit_in);
	  _hitIt--;
	}
	_planeIt = _hitsPerPlane.begin()+_endLongPlane;
	badplanes = (int)skipSize;
      } else _longestSingle = 0; //Brute force way to reject.
    } else {
      _endProj = true;
      _hitIt = hits.begin() + _endLongSing;
      _planeIt = _hitsPerPlane.begin() + _endLongPlane;
    }
      
  }

}

//***********************************************************************
bool event_classif::chargeCurrent_analysis(vector<cluster*>& hits,
					   Trajectory& muontraj, vector<cluster*>& hads){
//***********************************************************************
  m.message("++++ Performing CC reconstruction ++++",bhep::VERBOSE);

  bool ok = true;
  // bool fixed = false;
  _recChi = EVector(3,0);
  _recChi[1] = 100000; //SetLarge dummy value for minChiHadron bit.

  // _hitIt = hits.end() - 1;
  
  _vertGuess = exclude_backwards_particle();
  
  const dict::Key candHit = "inMu";
  const dict::Key hit_in = "True";
  // const dict::Key skipped = "skipped";
  const dict::Key not_in = "False";

  if ( _longestSingle >= min_seed_hits && _endLongPlane > 0.5*(double)_nplanes ){
    while ( _planeIt >= _hitsPerPlane.begin()+_exclPlanes+min_check && (*_planeIt) == 1 ){
      muontraj.add_measurement( *(*_hitIt) );
      (*_hitIt)->set_name(candHit, hit_in);
      _hitIt--;
      _planeIt--;

    }
  }
  
  //Maybe another which looks for kinks very early (backwards proton quasi?)
  // for (_planeIt = _hitsPerPlane.end()-1;_planeIt>=_hitsPerPlane.begin()+_exclPlanes+min_check;_planeIt--, _hitIt--){
//     if ( (*_planeIt) == 1 ){
//       muontraj.add_measurement( *(*_hitIt) );
//       (*_hitIt)->set_name(candHit, hit_in);
//       //  } else if ( _planeIt == _hitsPerPlane.begin() ) {
//       //       break;
//     // } else if ( _planeIt == _hitsPerPlane.end()-1 && (*_planeIt) < 4 && (*(_planeIt-1)) == 1 ){
// //       //Extra logic to avoid unneccessary use of call. auto. on high curve events.
// //       use_mini_cellAuto( (*_planeIt), muontraj );
      
//     } else if ( (int)muontraj.nmeas() < min_seed_hits
// 		&& _planeIt >= _hitsPerPlane.end()-(int)(0.2*(double)_nplanes) ){
      
//       int skiphit = (*_planeIt);
//       badplanes = 1;
      
//       int goodplane, iC1 = badplanes, iC2;
//       bool undone;
      
//       do {
// 	if ( (*(_planeIt-iC1)) != 1 ) { badplanes++; skiphit += (*(_planeIt-iC1)); }
// 	else {
// 	  iC2 = 1; undone = false; goodplane = 1;
// 	  while ( iC2 < min_seed_hits && !undone ){
// 	    if ( (*(_planeIt-iC1-iC2)) == 1 ) goodplane++;
// 	    else undone = true;
// 	    iC2++;
// 	  }
// 	  if ( goodplane == min_seed_hits ){//-(int)muontraj.nmeas() ){
// 	    fixed = true;
// 	    vector<cluster*>::iterator currentIt = _hitIt;
      
// 	    if ( (int)muontraj.nmeas() < min_seed_hits-2 && muontraj.nmeas() != 0 ){
// 	      for (int irem=1;irem<=(int)muontraj.nmeas();irem++){
// 		(*(_hitIt+irem))->set_name(candHit,not_in);
// 		(*(_hitIt+irem))->set_name(skipped,hit_in);}
// 	      muontraj.reset();
// 	    }
// 	    while ( _hitIt > currentIt-skiphit+1 ){
	      
// 	      (*_hitIt)->set_name(skipped,hit_in);
// 	      _hitIt--;
// 	    }
// 	    // for (int ifix=0;ifix<=iC1-1;ifix++)
// // 	      if ( ifix == iC1-1 ) _hitIt -= (*(_planeIt-ifix))-1;
// // 	      else _hitIt -= (*(_planeIt-ifix));
// 	    _planeIt -= iC1-1;
	    
// 	  } else if ( undone ){ badplanes++; skiphit += (*(_planeIt-iC1)); }
// 	}
// 	iC1++;
//       } while ( iC1 < (int)(0.2*(double)_nplanes) && !fixed );
//       if ( !fixed ) { badplanes = 0; break; }
//     } else break;
//   }
  
  _lastIso = (int)muontraj.nmeas();
  
  if ( !_outLike ) _freeplanes = (int)muontraj.nmeas();
  
  //set to reconstruction mode.
  if ( !MINDfitman::instance().in_rec_mode() )
    MINDfitman::instance().rec_mode();

  if ( _meanOcc == 1 ){
    
    _intType = 2;
    if ( muontraj.size() !=0 )
      ok = muon_extraction( hits, muontraj, hads);
    else ok = false;
    
  } else {
    
    if ( (int)muontraj.nmeas() < min_seed_hits ) {
      
      _intType = 5;
      ok = invoke_cell_auto( hits, muontraj, hads);
      if ( !ok ) _failType = 4;
      // ok = false; _failType = 4;
      
    } else {
      if ( badplanes !=0 ) _intType = 3;
      ok = muon_extraction( hits, muontraj, hads);
      // if ( !ok ) _failType = 4; 

      //Now I am going to create a new trajectory out of the hadron
      //hits and see if there is another potential trajectory
      // Trajectory newtraj;
      // vector<cluster*> hads2;
      // track_from_hads( hads, muontraj, hads2); 
      
    }
  }
  
  return ok;
}

//***********************************************************************
int event_classif::exclude_backwards_particle(){
//***********************************************************************
//Try to exclude some occurencies of  backwards had by looking
//for spaces between the early planes
  //occupancy from the start of the event.
  vector<int>::iterator measIt = _hitsPerPlane.begin();
  vector<double>::iterator zIt;
  int check_planes;
  if ( _nplanes >= 50 ) check_planes = 21;
  else check_planes = (int)( 0.4 * _nplanes ) + 1;
  int excluded_hits = 0;
  _exclPlanes = 0;
  bool found = false;

  //for (measIt = _hitsPerPlane.begin()+1;measIt != _hitsPerPlane.end();measIt++){
  for (zIt = _planeZ.begin()+1;zIt != _planeZ.begin()+check_planes;zIt++,measIt++){

    excluded_hits += (*measIt);
    _exclPlanes++;

    if ( abs( (*zIt) - (*(zIt-1)) ) >  2*_pieceLength ){
      found = true;
      break;
    }
    // if (measIt == _hitsPerPlane.begin()) {
//       excluded_hits += (*measIt);
//       continue;
//     }

    // if ( (*measIt) > (*(measIt - 1)) ) break;
//     else excluded_hits += (*measIt);
    

  }

  // if (measIt == _hitsPerPlane.end())
  if ( !found ){
    excluded_hits = 0;
    _exclPlanes = 0;
  }

  return excluded_hits;
}
//
void event_classif::use_mini_cellAuto(const int occ, Trajectory& muontraj){
  //
  double disp[occ];
  long minPos;
  double nextX = (*(_hitIt-occ))->vector()[0];
  double nextY = (*(_hitIt-occ))->vector()[1];
  //calculate x,y displacement to hit in next plane.
  for (int i=0;i<occ;i++){

    disp[i] = sqrt( pow( (*(_hitIt-i))->vector()[0] - nextX, 2) +
		    pow( (*(_hitIt-i))->vector()[1] - nextY, 2) );

  }
  //find the smallest displacement.
  minPos = TMath::LocMin( occ, disp );

  //Add selected hit to trajectory.
  muontraj.add_measurement( *(*(_hitIt-minPos)) );
  const dict::Key candHit = "inMu";
  const dict::Key hit_in = "True";
  (*(_hitIt-minPos))->set_name(candHit, hit_in);

  //move iterator so next plane hit is the next in line.
  _hitIt -= occ-1;

}

//***********************************************************************
bool event_classif::muon_extraction(vector<cluster*>& hits,
				    Trajectory& muontraj, vector<cluster*>& hads) {
//***********************************************************************
//Decide on appropriate pattern recognition and perform it.
//Are there ways to avoid doing full extraction algorithm.

  bool ok;
  
  if (_vertGuess != 0)
    for (int i = 0;i < _vertGuess;i++){
      const dict::Key candHit = "inhad";
      const dict::Key hit_in = "True";
      hits[i]->set_name(candHit, hit_in);
      hads.push_back( hits[i] );
    }
  
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
					vector<cluster*>& hits) {
//***********************************************************************

  EVector V(6,0); EVector V2(1,0);
  EMatrix M(6,6,0); EMatrix M2(1,1,0);
  double Xtent;

  V[0] = muontraj.nodes()[0]->measurement().vector()[0];
  V[1] = muontraj.nodes()[0]->measurement().vector()[1];
  V[2] = muontraj.nodes()[0]->measurement().position()[2];
  
  //direction
  double dqtot = fit_parabola( V, muontraj);
  
  //Momentum. Estimate from empirical extent function.
  if ( hits.size() != 0 ){

    // Xtent = hits[hits.size()-1]->position()[2]
//       - hits[_vertGuess]->position()[2];
    Xtent = hits[_endLongSing]->position()[2]
      - hits[_vertGuess]->position()[2];

  } else {
    
    Xtent = muontraj.nodes()[(int)muontraj.nmeas()-1]->measurement().position()[2]
      - muontraj.nodes()[0]->measurement().position()[2];

  }
  
  //double pSeed = 668 + 1.06*Xtent; //estimate in bhep::MeV, log for fit.
  double pSeed = (9180-6610*FeWeight) + (-2.76+4.01*FeWeight)*Xtent; //best for 3/2 seg
  // double pSeed = 574 + 0.66*Xtent;//test for 1/2 seg
  set_de_dx( pSeed/bhep::GeV );
  //Use slope in bending plane to try and get correct sign for momentum seed.
  // pSeed *= -fabs(V[3])/V[3];//Pos. gradient in X is ~negative curvature.
  // Assume that the R-Z plane is the bending plane for the toroidal field
  double VRad =  dqtot/fabs(dqtot);
  // (V[3]*V[0] + V[4]*V[1])/sqrt(V[0] * V[0] + V[1] * V[1]);
  if( VRad != 0 ) 
    pSeed *= fabs(VRad)/VRad;
  
  
  V[5] = 1./pSeed;
  
  //Errors
  M[0][0] = M[1][1] = 15.*bhep::cm*bhep::cm;
  M[2][2] = EGeo::zero_cov()/2;
  M[3][3] = M[4][4] = 1.5;
  M[5][5] = pow(V[5],2)*4;

  //Sense
  V2[0] = 1;
  
  //Seedstate fit properties
  seed.set_name(RP::particle_helix);
  seed.set_name(RP::representation,RP::slopes_curv_z); 
  seed.set_hv(RP::sense,HyperVector(V2,M2));
  seed.set_hv(HyperVector(V,M));
  
 //  man().model_svc().model(RP::particle_helix).representation(RP::slopes_curv_z)
//     .convert(seed,RP::default_rep);
  
  bool ok = perform_kalman_fit( seed, muontraj);
  
  if ( !ok )
    _failType = 5;
  
  return ok;
}

//***********************************************************************
double event_classif::fit_parabola(EVector& vec, Trajectory& track) {
//***********************************************************************
// No longer well named

  int fitcatcher;
  size_t nMeas = track.nmeas();

  if (nMeas > 4) nMeas = 4;

  EVector pos(3,0);
  EVector Z(3,0); Z[2] = 1;
  double x[(const int)nMeas], y[(const int)nMeas], 
    z[(const int)nMeas], u[(const int)nMeas], r[(const int)nMeas];
  int minindex = nMeas;
  double minR = 99999.999, pdR = 0.0, sumdq=0;
  double pathlength=0, tminR=9999.9;
  

  pos[0] = x[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[0];
  pos[1] = y[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[1];
  z[nMeas-1] = track.nodes()[nMeas-1]->measurement().position()[2];
  pos[2] = 0.0;

  for (int iMeas = nMeas-2;iMeas >= 0;iMeas--){
    x[iMeas] = track.nodes()[iMeas]->measurement().vector()[0];
    y[iMeas] = track.nodes()[iMeas]->measurement().vector()[1];
    z[iMeas] = track.nodes()[iMeas]->measurement().position()[2];
    // get the b-field from the previous step
    EVector B = geom.getBField(pos);
    u[iMeas] = dot(pos, crossprod(Z,B))/crossprod(Z,B).norm();
    pos[0] = x[iMeas];  pos[1] = y[iMeas];   pos[2] = z[iMeas];
    EVector dR = EVector(3,0);
    dR[0] = x[iMeas+1] - x[iMeas];
    dR[1] = y[iMeas+1] - y[iMeas];
    dR[2] = z[iMeas+1] - z[iMeas];
    double dq = 0.;
    pathlength += dR.norm();
    if(iMeas < nMeas-3){
      EVector dR1 = EVector(3,0);
      dR1[0] = x[iMeas+2] - x[iMeas+1];
      dR1[1] = y[iMeas+2] - y[iMeas+1];
      dR1[2] = z[iMeas+2] - z[iMeas+1];
      EVector ddR = dR1 + dR;
      dq = dot(ddR, crossprod(Z,B))/ (crossprod(Z,B).norm());
    }
    if(pdR != 0)
      if(// minR < R && dR/fabs(dR) == -pdR/fabs(pdR) && 
	 (x[iMeas]/fabs(x[iMeas]) != x[iMeas+1]/fabs(x[iMeas+1]) ||
	  y[iMeas]/fabs(y[iMeas]) != y[iMeas+1]/fabs(y[iMeas+1]))){ 
	// Sign change and minimum momentum
	minR = pos.norm();
	minindex = iMeas;
	sumdq = 0.0;
      }
    pdR = dR.norm();
    sumdq += dq;
  }
  
  
  double wFe = geom.get_Fe_prop();
  double p = RangeMomentum(pathlength);

  TGraph *gr1 = new TGraph((const int)minindex, z, x);
  TGraph *gr2 = new TGraph((const int)minindex, z, y);
  TGraph *gr3 = new TGraph((const int)minindex, z, u);

  TF1 *fun = new TF1("parfit","[0]+[1]*x",-3,3);
  fun->SetParameters(0.,0.001);
  TF1 *fun2 = new TF1("parfit2","[0]+[1]*x+[2]*x*x",-3,3);
  fun2->SetParameters(0.,0.001,0.001);
  
  fitcatcher = gr1->Fit("parfit", "QN");
  vec[3] = fun->GetParameter(1);

  fun->SetParameters(0.,0.001);
  fitcatcher = gr2->Fit("parfit", "QN");
  vec[4] = fun->GetParameter(1);

  // fun->SetParameters(0.,0.001);
  fitcatcher = gr3->Fit("parfit2", "QN");
  double qtilde = -0.3*pow(1+fun2->GetParameter(1),3./2.)/
    (2*fun2->GetParameter(2));
  
  // int meansign = sumdq/fabs(sumdq);
  int meansign = qtilde/fabs(qtilde);

  vec[5] = meansign/p;

  delete gr1;
  delete gr2;
  delete fun;
  // return sumdq;
  return qtilde;

}

//***********************************************************************
void event_classif::set_de_dx(double mom){
//***********************************************************************
  
  double de_dx = -( 12.37 * FeWeight * pow( mom, 0.099) 
		    + (1 - FeWeight) * 2.16 * pow( mom, 0.075) );
  de_dx *= bhep::MeV/bhep::cm;
  
  man().geometry_svc().setup().set_volume_property_to_sons("mother","de_dx",de_dx);

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
bool event_classif::perform_muon_extraction(const State& seed, vector<cluster*>& hits,
					    Trajectory& muontraj, vector<cluster*>& hads) {
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

	  if ( (*(_hitIt+iFil+1))->get_mu_prop() > 0.8 )
	    _recChi[0] = TMath::Max(Chi2[iFil], _recChi[0]);
	  else
	    _recChi[1] = TMath::Min(Chi2[iFil], _recChi[1]);

	  if ( nConsecHole > _recChi[2] )
	    _recChi[2] = nConsecHole;

	  nConsecHole = 0;

	} else nConsecHole++;

      } else {
	const dict::Key hadHit = "inhad";
	const dict::Key had_in = "True";
	(*(_hitIt+iFil+1))->set_name(hadHit, had_in);
	hads.push_back( (*(_hitIt+iFil+1)) );

      }
    }

    _planeIt--;
    
    if ( nConsecHole > max_consec_missed_planes ) {
      if ( _endProj ) check_forwards( seed, hits, muontraj );
      if ( muontraj.nmeas() < min_plane_prop*(double)(_nplanes-badplanes) ){
	_failType = 6;
	return false;
      } else {
	sort_hits( hits, muontraj, hads );
	return true;
      }
    }

  }

  if ( _endProj ) check_forwards( seed, hits, muontraj );
  
  return true;
}

//***********************************************************************
void event_classif::check_forwards(const State& seed, vector<cluster*>& hits,
				   Trajectory& muontraj)
//***********************************************************************
{
  //Those not in 'blob' class get checked if there are any more hits
  //that can be added to the trajectory.
  //Reset the iterator positions.
  _hitIt = hits.begin() + _endLongSing + 1;
  _planeIt = _hitsPerPlane.begin() + _endLongPlane + 1;
  
  bool ok;
  int counter = 0;
  long ChiMin;
  double chi2[2];
  const dict::Key candHit = "inMu";
  const dict::Key hit_in = "True";
  const dict::Key skipped = "skipped";
  
  while ( (*_planeIt) == 2 && _planeIt != _hitsPerPlane.end() ){

    for (int j=0;j<2;j++){
      
      ok = man().matching_svc().match_trajectory_measurement(muontraj, (*(*_hitIt)), chi2[j]);
      _hitIt++;

    }
    
    ChiMin = TMath::LocMin( 2, chi2);

    if ( ChiMin == 0 ){
      ok = man().fitting_svc().filter(*(*(_hitIt-2)), seed, muontraj);
      if ( !ok ) (*(_hitIt-2))->set_name(skipped,hit_in);
      else (*(_hitIt-2))->set_name(candHit, hit_in);
      (*(_hitIt-1))->set_name(skipped,hit_in);
    } else {
      ok = man().fitting_svc().filter(*(*(_hitIt-1)), seed, muontraj);
      if ( !ok ) (*(_hitIt-1))->set_name(skipped,hit_in);
      else (*(_hitIt-1))->set_name(candHit, hit_in);
      (*(_hitIt-2))->set_name(skipped,hit_in);
    }
    _planeIt++;
  }
  
  while ( _hitIt != hits.end() ){
    (*_hitIt)->set_name(skipped,hit_in);
    _hitIt++;
  }
  
  while ( _planeIt != _hitsPerPlane.end() ){
    counter++;
    _planeIt++;
  }
  badplanes = counter;
  
}

//***********************************************************************
bool event_classif::invoke_cell_auto(vector<cluster*>& hits,
				     Trajectory& muontraj, vector<cluster*>& hads){
//***********************************************************************
//uses cellular automaton to try and retrieve more complicated events.

  m.message("+++ Performing Cellular Automaton analysis +++",bhep::VERBOSE);
  
  bool ok;

  if ( _nplanes-_exclPlanes < min_hits ) return false;

  std::vector<Trajectory*> trajs;
  //TEST!!
  measurement_vector hit_meas;
  get_cluster_meas( hits, hit_meas );
  //
  //ok = man().matching_svc().find_trajectories( hits, trajs);
  ok = man().matching_svc().find_trajectories( hit_meas, trajs );
  
  if ( !ok || trajs.size() == 0) return false;
  
  // output_results_tree( hits, trajs );
  
  if ( trajs.size() == 1 ) {
    
    if ( trajs[0]->nmeas() >= min_hits ){
      muontraj.reset();
      
      muontraj = *trajs[0];
      
      sort_hits( hits, muontraj, hads);
      
      ok = true;
    } else ok = false;

  } else {

    ok = sort_trajs( muontraj, trajs);

    if ( ok )
      sort_hits( hits, muontraj, hads);

  }

  if ( ok )
    delete_bad_trajs( muontraj, trajs );
  else stc_tools::destroy( trajs );

  m.message("+++ End of Cellular automaton +++",bhep::VERBOSE);
  
  return ok;
}

void event_classif::get_cluster_meas(const vector<cluster*>& hits,
				     measurement_vector& meas)
{

  std::vector<cluster*>::const_iterator cIt;
  for (cIt = hits.begin()+_vertGuess;cIt != hits.end();cIt++)
    meas.push_back( (*cIt) );

}

//***********************************************************************
void event_classif::sort_hits(vector<cluster*>& hits,
			      Trajectory& muontraj, vector<cluster*>& hads){
//***********************************************************************
  
  bool inTraj;
  vector<Node*>::iterator inTrIt;
  inTrIt = muontraj.nodes().begin();

  const dict::Key candHit = "inMu";
  const dict::Key not_in = "False";
  const dict::Key hit_in = "True";
  const dict::Key skipped = "skipped";
  
  //Protection against double entries in hads in max_consec_planes fail case.
  if ( _intType != 5 ) hads.clear();

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

    if ( !inTraj && !(*_hitIt)->names().has_key(skipped) ){
      const dict::Key hadHit = "inhad";
      const dict::Key had_in = "True";
      (*_hitIt)->set_name(hadHit, had_in);
      hads.push_back( (*_hitIt) );
    }

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
      
      //!!!!
      trajs2.clear();

      if ( ok ){
	if ( trajs3.size() == 1 ){
	  
	  muontraj.reset();
      
	  muontraj = *trajs3[0];
      
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
  vector<cluster*> dummy_hits;
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
      momErr = (double)(temp_seed.matrix()[5][5] / temp_seed.vector()[5]);
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

  double coincidence = 0;

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

    length[it1] = trajs[it1]->length();
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
  
  _likeTree->Branch("truInt",&_truInt,"Int/I");
  _likeTree->Branch("nhits", &_nhit, "nhits/I");
  _likeTree->Branch("VisEng", &_visEng, "visEng/D");
  _likeTree->Branch("nplanes", &_nplanes, "nplanes/I");
  _likeTree->Branch("freePlanes", &_freeplanes, "freeplanes/I");
  _likeTree->Branch("occupancy", &_occ, "occ[nplanes]/I");
  _likeTree->Branch("planeEnergy", &_plEng, "plEng[nplanes]/D");
  _likeTree->Branch("nhitTraj", &_trajhit, "nhitTraj/I");
  _likeTree->Branch("TrajPur", &_trajpur, "trajpur/D");
  _likeTree->Branch("EngTraj", &_trajEng, "engTraj/D");
  _likeTree->Branch("EngTrajPlane", &_trajEngPlan, "engTrajPlane[nplanes]/D");
  _likeTree->Branch("HitsInTrajClust",&_trclusthits,"nhittrajclust[nplanes]/I");

}

//***********************************************************************
void event_classif::output_liklihood_info(const vector<cluster*>& hits){
//***********************************************************************
  
  bool multFound = false;

  _nhit = 0;
  _freeplanes = 0;
  _visEng = 0;
  _trajhit = 0;
  _trajpur = 0;
  _trajEng = 0;
  for (int ipl = 0;ipl<_nplanes;ipl++){
    _trajEngPlan[ipl] = 0;
    _trclusthits[ipl] = 0; }
  
  vector<int>::iterator itLike;
  vector<double>::iterator itEngL = _energyPerPlane.end()-1;
  
  int counter = _nplanes - 1;

  for (itLike = _hitsPerPlane.end()-1;itLike >= _hitsPerPlane.begin();itLike--, itEngL--){

    _nhit += (*itLike);
    _visEng += (*itEngL);

    _occ[counter] = (*itLike);
    _plEng[counter] = (*itEngL);
    
    if ( !multFound && (*itLike) == 1 )
      _freeplanes++;
    else multFound = true;
    
    counter--;

  }
  
}

//***********************************************************************
void event_classif::traj_like(const vector<cluster*>& hits, const Trajectory& muontraj){
//***********************************************************************
  
  const dict::Key Edep = "E_dep";
  const dict::Key candHit = "inMu";
  vector<cluster*>::const_iterator hitIt3;
  //vector<Node*>::const_iterator trIt1 = muontraj.nodes().begin();
  vector<double>::iterator planIt;

  int counter = _nplanes - 1;
  
  for (hitIt3 = hits.begin();hitIt3 != hits.end();hitIt3++){
    if ( (*hitIt3)->names().has_key( candHit ) ){
      if ( (*hitIt3)->name( candHit ).compare("True") == 0 ){
	_trajhit++;
	//_trajEng += bhep::double_from_string( (*hitIt3)->name( Edep ) ) * bhep::GeV;
	_trajEng += (*hitIt3)->get_eng() * bhep::MeV;
	_trajEngPlan[counter] = (*hitIt3)->get_eng() * bhep::MeV;
	// if ( (*hitIt3)->name("MotherParticle").compare("mu+") == 0
// 	     || (*hitIt3)->name("MotherParticle").compare("mu-") == 0 )
	if ( (*hitIt3)->get_mu_prop() > 0.8 )//still to be decided.
	  _trajpur++;
	_trclusthits[counter] = (*hitIt3)->get_nVox();
	counter--;
      }
    }
  }
  _trajpur /= _trajhit;
  
  // for (planIt = _planeZ.end()-1;planIt >= _planeZ.begin();planIt--){
//     if ( trIt1 != muontraj.nodes().end() )
//       if ( fabs( (*planIt) - (*trIt1)->measurement().position()[2]) < _tolerance ){
	
// 	_trajEngPlan[counter] = bhep::double_from_string( (*trIt1)->measurement().name( Edep ) ) * bhep::GeV;
// 	trIt1++;
	
//       }
//     counter--;
//   }
  
}

//***********************************************************************
void event_classif::out_like(){
//***********************************************************************
  
  int fillcatch = _likeTree->Fill();
  
}

//Temp function for likelihoods with tru interaction in tree.
void event_classif::set_int_type(const string name){

  if ( name=="CCQE" )
    _truInt = 1;
  else if ( name=="NCQE" )
    _truInt = 2;
  else if ( name=="CCDIS" )
    _truInt = 3;
  else if ( name=="NCDIS" )
    _truInt = 4;
  else if ( name=="1piRes" )
    _truInt = 5;
  else if ( name=="miscRes" )
    _truInt = 6;
  else if ( name=="eEl" )
    _truInt = 7;
  else if ( name=="muINVe" )
    _truInt = 7;
  else if ( name=="unknown" )
    _truInt = -1;
  else
    _truInt = 8;

}


double event_classif::correctEdep(double edep, double X, double Y)
{

  double sum1 = 0;
  
  double slope = (_detY - _detX*tan(atan(1)/2.))/
    (_detY*tan(atan(1)/2.) - _detX);

  //need to take into account drift distance to closest and furthest edge.
  if((fabs(X) < _detY*tan(atan(1)/2.)/2 &&
      fabs(Y) < _detX*tan(atan(1)/2.)/2 ) ){
    sum1 =  exp( -(_detX/2 - fabs(X))/_WLSAtten );
    sum1 += exp( -(_detX/2 + fabs(X))/_WLSAtten );
    sum1 += exp( -(_detY/2 - fabs(Y))/_WLSAtten );
    sum1 += exp( -(_detY/2 + fabs(Y))/_WLSAtten );
  }
  else if(fabs(X) > _detY*tan(atan(1)/2.)/2 &&
	  fabs(Y) < _detX*tan(atan(1)/2.)/2  ){
    double xedge = _detX/2 + (fabs(Y) - 
				   _detX/2.*tan(atan(1)/2.))*slope;
    sum1 =  exp( -(xedge  - fabs(X))/_WLSAtten );
    sum1 += exp( -(xedge + fabs(X))/_WLSAtten );
    sum1 += exp( -(_detY/2 - fabs(Y))/_WLSAtten );
    sum1 += exp( -(_detY/2 + fabs(Y))/_WLSAtten );
  }
  else if(fabs(X) < _detY*tan(atan(1)/2.)/2 &&
	  fabs(Y) > _detX*tan(atan(1)/2.)/2  ){
    double yedge = _detY/2 + (fabs(X) - 
				   _detY/2.*tan(atan(1)/2.))*slope;
    sum1 =  exp( -(_detX/2 - fabs(X))/_WLSAtten );
    sum1 += exp( -(_detX/2 + fabs(X))/_WLSAtten );
    sum1 += exp( -(yedge - fabs(Y))/_WLSAtten );
    sum1 += exp( -(yedge + fabs(Y))/_WLSAtten );
  }
  else if(fabs(X) > _detY*tan(atan(1)/2.)/2 &&
	  fabs(Y) > _detX*tan(atan(1)/2.)/2  ){
    double xedge = _detX/2 + (fabs(Y) - 
				   _detX/2.*tan(atan(1)/2.))*slope;
    double yedge = _detY/2 + (fabs(X) - 
				   _detY/2.*tan(atan(1)/2.))*slope;
    sum1 =  exp( -(xedge - fabs(X))/_WLSAtten );
    sum1 += exp( -(xedge + fabs(X))/_WLSAtten );
    sum1 += exp( -(yedge - fabs(Y))/_WLSAtten );
    sum1 += exp( -(yedge + fabs(Y))/_WLSAtten );
  }
  
  double corrEng = 4*edep/sum1;
  return corrEng;
}


double event_classif::RangeMomentum(double length){
  // Computing the momentum of a track based on the range
  // Parameters calculated from "Atomic Data and Nuclear Data Tables 78, 183-356 (2001)"
  double p = (FeWeight*(0.011844*bhep::GeV * pow((length/bhep::cm),1.03235))
	      + (1- FeWeight)*(0.0023705067*bhep::GeV* pow(length/bhep::cm,1.00711577)));
  return p;
  // Should be compared to G4 simulation directly.
}
