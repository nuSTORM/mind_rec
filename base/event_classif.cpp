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



using namespace bhep;


class sortTrajByHits{
public:
  bool  operator()(const Trajectory* t1,const Trajectory* t2 ){
    if (t1->nmeas() > t2->nmeas()) return true;
    return false;
  }
  
};
class sortTrajByLength{
public:
  bool  operator()(const Trajectory* t1,const Trajectory* t2 ){
    if (t1->length() > t2->length()) return true;
    return false;
  }

};


class sortHitsByR{
public:
  bool  operator()(cluster* hit1, cluster* hit2 ){
    double r1= sqrt(pow( hit1->position()[0],2) + pow( hit1->position()[1],2));
    double r2= sqrt(pow( hit2->position()[0],2) + pow( hit2->position()[1],2));
    if (r1 > r2) return true;
    return false;
  }
  
};
class sortHitsByZ{
public:
  bool  operator()(cluster* hit1, cluster* hit2 ){
    //double r1= sqrt(pow( hit1->position()[0],2) + pow( hit1->position()[1],2));
    // double r2= sqrt(pow( hit2->position()[0],2) + pow( hit2->position()[1],2));
    if (hit2->position()[2] > hit1->position()[2]) return true;
    return false;
  }
  
};


//**********************************************************************
event_classif::event_classif() {
  //**********************************************************************

  

}

//***********************************************************************
event_classif::~event_classif() {
  //***********************************************************************


}

//Initialization of the class
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

  _FeWeight = wFe;

  if ( _outLike ){
    string likeFile = pstore.fetch_sstore("like_file");
    _outFileEv = new TFile(likeFile.c_str(), "recreate");
    _likeTree = new TTree("h1", "Possible Likelihood parameters");
    
    
    set_branches();
    
  }
  
}

//Execution of the class
//***********************************************************************
bool event_classif::execute(vector<cluster*>& hits,
			    vector<Trajectory*> &vmuontrajs, vector<cluster*>& hads) {
  //***********************************************************************
 
  m.message("++++ Classifier Execute Function ++++", bhep::VERBOSE);
  
  ///reset variable for each event
  reset();
  
  ///destray all vector otrajectories
  stc_tools::destroy( vmuontrajs ); 
 
  ///start looking for trajectories
  bool ok;
  while (hits.size() > 5) {  
    
    ///create trajectory     
    Trajectory* muontraj = new Trajectory();
            
    ///sort hits in increasing Z positions
    sort( hits.begin(), hits.end(), sortHitsByZ() );
    
    ///calculate the number of planes containing single hit and arrange the z-position, energy of the hits in increasing z 
    //Occupancy.
    ok = get_plane_occupancy( hits );
        
    ///called if liklihood output has to generate
    if ( _outLike )
      output_liklihood_info( hits );
        
    //get free section along Z

    ///perform muon extraction
    /* if ( ok )
      ok = muon_extraction_through_PatterRec(hits, *muontraj, hads);
    if(!ok) cout<<"muon_extraction_through_PatterRec not ok"<<endl;



     ///if PR is ok and nmeas < 5 then check for CA
    if ( ok && _freeplanes < min_seed_hits) {
      _intType = 5;
      ok = invoke_cell_auto( hits, *muontraj, hads);
      // if(!ok) cout<<" invoke_cell_auto not ok"<<endl;
      if ( !ok ) _failType = 4;
      } */

    
    /// CA and PR both
    ok = chargeCurrent_analysis(hits, *muontraj, hads);
      
    cout<<"nmeas in traj="<<muontraj->nmeas()<<"  and hadrons left="<<hads.size()<<endl;


    ///for liklhood
    if ( ok && _outLike)
      traj_like( hits, *muontraj );
    
    if ( _outLike ) out_like();


    //    cout<<"inside classifier traj nmeas ="<<muontraj->nmeas()<<"  and traj is="<<*muontraj<<endl;
   

 
    //fill the vector of trajectories
    if(ok) vmuontrajs.push_back(muontraj);
    
    
    
    ///set information of trajectory
    fill_traj_info( *muontraj); 
    
    // cout<<"  nmeas="<<muontraj->nmeas()<<endl;  
    
    
    ///clear vector of clusters
    hits.clear();
   
    
    ///if hadrons > 5 then loop
    if(hads.size() >= _min_hits) hits = hads;
    else break;
    
    ///clear hadron container
    hads.clear();
   
   
  }
  cout<<"eventclass::PR, size = "<< _vPR_seed.size()<<" vtraj size ="<<vmuontrajs.size()<<endl;

  
  return ok;
    
}

//***********************************************************************
void event_classif::finalize() {
  //***********************************************************************

  if ( _outLike ){
    _outFileEv->Write();
    _outFileEv->Close();

  }
  // _CAFile->Write();
  // _CAFile->Close();

}

//***********************************************************************
void event_classif::readParam() {
  //***********************************************************************
 
  m.message("++++ readParam function of classifier ++++",bhep::VERBOSE);
  
  if ( _infoStore.find_istore("likeli") )
    _outLike = _infoStore.fetch_istore("likeli");
  else _outLike = false;
  
  _max_consec_missed_planes = _infoStore.fetch_istore("max_consec_missed_planes");
  _min_plane_prop = _infoStore.fetch_dstore("min_use_prop");
  _min_seed_hits = _infoStore.fetch_istore("min_seed_hits");
  _min_check =  _infoStore.fetch_istore("min_check_nodes");
  _min_hits = _infoStore.fetch_istore("low_Pass_hits");
  _max_hits = _infoStore.fetch_istore("plane_occupancy");  
  _max_nmult = _infoStore.fetch_istore("max_multOcc_plane");
  _tolerance = _infoStore.fetch_dstore("pos_res") * bhep::cm;
  _chi2_max = _infoStore.fetch_dstore("accept_chi");
  _max_coincedence = _infoStore.fetch_dstore("max_coincidence");

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
  _badplanes = 0;
  _longestSingle = 0;
  _endProj = false;
  _hitsPerPlane.clear();
  _energyPerPlane.clear();
  _planeZ.clear();
  _vRadCount.clear();
  _radialLongest = 0;
  _Xtent = 0;
  _intType = 0;
  _endLongSing = 0;
  _endLongPlane = 0;
  _lastIso = 0;
  ///25062012
  _vertGuess = 0;
  _vertGuessZ = 0;
  _failType = 0;

  
  
  _vPR_seed.clear();
  
  ///

}

//***********************************************************************
bool event_classif::get_plane_occupancy(vector<cluster*>& hits){
  //***********************************************************************
  //Gets plane occupancies and total plane energies.
  //Needs hits in increasing z order.
  /*calculate the number of planes(_hitsPerPlane) containing single hit i.e., _nplanes and arrange the 
    z-position(_planeZ), energy(engPlane) of the hits in increasing z*/ 
  
  //   std::cout<<"+++++I am in get_plane_occupancy"<<std::endl;
  m.message("++++ Calculating plane energies and occupancies ++++",bhep::VERBOSE);

  bool ok = true;

  /// size of vector<cluster> hits
  size_t nHits = hits.size();
  cout<<"nHits="<<nHits<<endl;
  int count = 0, single_count = 0;
  double EngPlane = 0, testZ, curZ;
  size_t hits_used = 0, imeas = 0;
  _longestSingle = 0;

  if ( !_outLike ) _visEng = 0;
  
  do {
        
    EngPlane = correctEdep(hits[imeas]->get_eng()*bhep::MeV, 
			   hits[imeas]->position()[0], hits[imeas]->position()[1]);

    testZ = hits[imeas]->position()[2];
    count++;
    hits_used++;
    
    ///calculate the z position which is the current z for hits 1 -> total no of hits in the cluster 
    for (size_t i = hits_used;i <nHits;i++) {
      curZ = hits[i]->position()[2];
                 
      if (curZ <= testZ + _tolerance) {

	//	EngPlane += hits[i]->get_eng() * bhep::MeV
	double x = hits[i]->position()[0];
	double y = hits[i]->position()[1];
	EngPlane += correctEdep(hits[i]->get_eng() * bhep::MeV, x, y);

	testZ = hits[i]->position()[2];
	count++;
	hits_used++;
      } else break;      
    }



    ///fill the vectors of single hitsPerPlane, energy and arrange them in increasing z
    _hitsPerPlane.push_back( count );
    _energyPerPlane.push_back( EngPlane );
    _planeZ.push_back( testZ );
    
    ///if single hit plane then increase the single_count, otherwise assign it 0
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
    

    if ( imeas == nHits -1 && single_count != 0 && count == 1 ){
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


  ///total no of planes
  _nplanes = (int)_hitsPerPlane.size();
 
  if ( _nplanes == 0 ) return false;
  _meanOcc /= (double)_nplanes;


  ///if no of single-hits in the longest part is > 5 and free part is 50% of the hitsperplane then 
  ///assess the eventś free section to look for viable muon

  if ( _longestSingle >= _min_seed_hits && _endLongPlane > 0.5*(double)_nplanes )
    assess_event( hits );
  

  return ok;
}


//***********************************************************************
bool event_classif::get_radial_occupancy(vector<cluster*>& hits)
  //***********************************************************************
{
   m.message("++++ Calculating Radial occupancies ++++",bhep::VERBOSE);
 
  
  bool ok = true;
  
 
  int count = 0, single_count = 0;
  double  testR, curR; /// EngPlane = 0, 
  size_t hits_used = 0, imeas = 0;
  _radialLongest = 0;
  _radialFree = 0;
  
  vector<cluster*> hits_temp = hits;
  size_t nHits = hits_temp.size();
  
  
  ///sort hits by R independent of Z  
  // sort(hits_temp.begin(), hits_temp.end(), sortHitsByR());
  // cout<<"  radial_hits end Z="<<hits_temp.back()->position()[2]<<endl;  


  ///start count radial occupancy
  do {
       
    testR = sqrt(pow( hits_temp[imeas]->position()[0],2) + pow( hits_temp[imeas]->position()[1],2));
    count++;
    hits_used++;
    
    ///calculate the z position which is the current z for hits 1 -> total no of hits in the cluster 
    for (size_t i = hits_used;i <nHits;i++) {
      curR = sqrt(pow( hits_temp[i]->position()[0],2) + pow( hits_temp[i]->position()[1],2));
      
      	
      ///tolerance 5mm in R
      if (testR <= curR + 5*bhep::mm) {
	testR = sqrt(pow( hits_temp[i]->position()[0],2) + pow( hits_temp[i]->position()[1],2));
	count++;
	hits_used++;
      } else break;
	
    }
      
    _vRadCount.push_back(count);
      
 
    ///if single hit plane then increase the single_count, otherwise assign it 0
    if ( count == 1 )
      single_count++;
    else {
      if ( single_count >= _radialLongest ){
	_radialLongest = single_count;
	
      }
      single_count = 0;
    }
    if ( imeas == nHits -1 && single_count != 0 && count == 1 ){
      if ( single_count >= _radialLongest ){
	_radialLongest = single_count;
	
      }
    }
    
  
    ///increase the measurement count and also mean0cc
    imeas += count;
      
    count = 0;
  } while (hits_used != nHits);
    
 
  return ok;
}





//***********************************************************************
void event_classif::assess_event(vector<cluster*>& hits)
  //***********************************************************************
{
   m.message("++++ Assess event for complicated end point ++++",bhep::VERBOSE);
  //cout<<"inside assess_event"<<endl;
  //Assess the event 'free' sections to look for a viable
  //muon and tidy up the endpoint if neccessary.
  

  double dispTrans[2], endMeanOcc = 0, skipSize = 0;
  unsigned int evEnd = hits.size() - 1;

  vector<int>::iterator aIt;
  const dict::Key hit_in = "True";
  const dict::Key skipped = "skipped";


  ///At the last hit of the longest free section  
  if ( _endLongSing == (int)evEnd ){
    
    dispTrans[0] = sqrt( pow( hits[evEnd]->vector()[0]-hits[evEnd-1]->vector()[0], 2)
			 + pow( hits[evEnd]->vector()[1]-hits[evEnd-1]->vector()[1], 2) );
    dispTrans[1] = sqrt( pow( hits[evEnd-1]->vector()[0]-hits[evEnd-2]->vector()[0], 2)
			 + pow( hits[evEnd-1]->vector()[1]-hits[evEnd-2]->vector()[1], 2) );
   

    //Just in case there is a bad hit at the endpoint.
    ///_voxEdge=rec_boxX D 3.5 bhep::cm
    if ( dispTrans[0] > _voxEdge*10 && dispTrans[1] <= _voxEdge*10 ){
      hits[evEnd]->set_name(skipped,hit_in);
      _badplanes = 1;
      _hitIt = hits.end() - 2;
      _planeIt = _hitsPerPlane.end() - 2;
     
    } else if ( dispTrans[1] > _voxEdge*10 ){
      hits[evEnd]->set_name(skipped,hit_in);
      hits[evEnd-1]->set_name(skipped,hit_in);
      _badplanes = 2;
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
	   || _longestSingle >= 2*_min_seed_hits ){
	
	_hitIt = hits.end() - 1;
	while ( _hitIt != hits.begin()+_endLongSing ){
	  (*_hitIt)->set_name(skipped,hit_in);
	  _hitIt--;
	}
	_planeIt = _hitsPerPlane.begin()+_endLongPlane;
	_badplanes = (int)skipSize;
      } else _longestSingle = 0; //Brute force way to reject.
    } else {
      
      _endProj = true;
      _hitIt = hits.begin() + _endLongSing;
      _planeIt = _hitsPerPlane.begin() + _endLongPlane;
    }
    
  }
  
}



//***********************************************************************
bool event_classif::muon_extraction_through_PatterRec(vector<cluster*>& hits,
						      Trajectory& muontraj, vector<cluster*>& hads){
  //***********************************************************************
  m.message("++++ Performing muon reconstruction ++++",bhep::VERBOSE);
  //cout<<"inside muon_extraction_through_PatterRec"<<endl;
  
  muontraj.reset();
  bool ok = true;
  
  _recChi = EVector(3,0);
  _recChi[1] = 100000; //SetLarge dummy value for minChiHadron bit.
  
  
  
  /*//looking into beginning of planes to guess vertex hitno, to exclude backward particle hits and if not found then  excluded_hits = 0; _exclPlanes = 0; i.e, vertGuess =0*/
  _vertGuess = exclude_backwards_particle();
  
  
  ///Zpos of vertex
  _vertGuessZ = hits[_vertGuess]->position()[2];
  
  //cout<<"_vertGuess ="<<_vertGuess<<"  "<<hits[_vertGuess]->position()[2]<<endl;
  
  const dict::Key candHit = "inMu";
  const dict::Key hit_in = "True";
  const dict::Key not_in = "False";
  
  
  
  //add measurements to muontraj 
  if ( _longestSingle >= _min_seed_hits && _endLongPlane > 0.5*(double)_nplanes ){
    while ( _planeIt >= _hitsPerPlane.begin()+_exclPlanes+_min_check && (*_planeIt) == 1 ){
      muontraj.add_measurement( *(*_hitIt) );
      (*_hitIt)->set_name(candHit, hit_in);
      _hitIt--;
      _planeIt--;
      
    }
  }
  
  //cout<<" after Z search ::traj.nmeas()="<<muontraj.nmeas()<<endl; 
  
  /*//radial search
    int kount = 0;
    if( (int)muontraj.nmeas() < _min_seed_hits ) {
  
    vector<cluster*> temp_hits = hits;
    vector<cluster*>::iterator temp_hitIt;
  
  
    ///sort hits by R independent of Z  
    sort(temp_hits.begin(), temp_hits.end(), sortHitsByR());
    
    ///calculate radial occupancy 
    ok = get_radial_occupancy( temp_hits );
    //cout<<"after get_radial_occ:: radialOcc ="<<_radialLongest <<endl;
    
    
    //look for radial free section
    m.message("vertex Z="<<_planeZ[0]<<"   hits begin="<<hits[_vertGuess]->position()[2],bhep::DETAILED);
    cout<<"vertex Z="<<_planeZ[0]<<"   hits begin="<<hits[_vertGuess]->position()[2]<<endl;
    if (temp_hits.back()->position()[2] > hits[_vertGuess]->position()[2] ){
    _RadIt = _vRadCount.end()-1;
    temp_hitIt = temp_hits.end() -1;
      
    while ( temp_hitIt >= temp_hits.begin() && (*_RadIt) == 1 ){
    muontraj.add_measurement( *(*temp_hitIt) );
    (*temp_hitIt)->set_name(candHit, hit_in);
    temp_hitIt--;
    _RadIt--;
	
    ///

    if (muontraj.nmeas()>=1) m.message(<<" 1:temp_hits[muontraj.nmeas()]->position()[2] "<<temp_hits[muontraj.nmeas()]->position()[2]<<"  temp_hits[muontraj.nmeas()-1]->position()[2]="<<temp_hits[muontraj.nmeas()-1]->position()[2], bhep::DETAILED);

    if (muontraj.nmeas()>=1) cout<<" 1:temp_hits[muontraj.nmeas()]->position()[2] "<<temp_hits[muontraj.nmeas()]->position()[2]<<"  temp_hits[muontraj.nmeas()-1]->position()[2]="<<temp_hits[muontraj.nmeas()-1]->position()[2]<<endl;
    if (muontraj.nmeas()>=1 && temp_hits[muontraj.nmeas()]->position()[2] != temp_hits[muontraj.nmeas()-1]->position()[2]) kount = muontraj.nmeas();
	
    }
    }
    else {
      
    _RadIt = _vRadCount.begin();
    temp_hitIt = temp_hits.begin();
      
    while ( temp_hitIt != temp_hits.end() && (*_RadIt) == 1 ){
    muontraj.add_measurement( *(*temp_hitIt) );
    (*temp_hitIt)->set_name(candHit, hit_in);
    temp_hitIt++;
    _RadIt++;
	
    ///
    if (muontraj.nmeas()>=1) m.message(" 2:temp_hits[muontraj.nmeas()]->position()[2] "<<temp_hits[muontraj.nmeas()]->position()[2]<<"  temp_hits[muontraj.nmeas()-1]->position()[2]="<<temp_hits[muontraj.nmeas()-1]->position()[2], bhep::DETAILED);

    if (muontraj.nmeas()>=1) cout<<" 2: temp_hits[muontraj.nmeas()]->position()[2] "<<temp_hits[muontraj.nmeas()]->position()[2]<<"  temp_hits[muontraj.nmeas()-1]->position()[2]="<<temp_hits[muontraj.nmeas()-1]->position()[2]<<endl;
    if (muontraj.nmeas()>=1 && temp_hits[muontraj.nmeas()]->position()[2] != temp_hits[muontraj.nmeas()-1]->position()[2]) kount = muontraj.nmeas();
	
    }
    }
    
    if(_vRadCount.size() == 0) return false;
    _meanOcc = temp_hits.size()/(double)_vRadCount.size(); 
    //cout<<temp_hits.size()<<"  "<<muontraj.size()<<endl;
    
    cout<<"kount="<<kount<<"  temp_hits end Z="<<temp_hits.back()->position()[2]<<" and begin ="<<temp_hits.front()->position()[2]<<"  size="<<temp_hits.size()<<endl;

    m.message("kount="<<kount<<"  temp_hits end Z="<<temp_hits.back()->position()[2]<<" and begin ="<<temp_hits.front()->position()[2]<<"  size="<<temp_hits.size()bhep::DETAILED);
   
    ///Rough estimation to instantiate the iterators
    _hitIt = hits.end() - muontraj.nmeas() -1;
    _planeIt = _hitsPerPlane.end() - kount -1;

        
    }*/
  
  
  
  
  //information from the trajectory 
  _lastIso = (int)muontraj.nmeas();
  m.message(" _lastIso = ",_lastIso,bhep::DETAILED);
  

  ///freeplanes inside the candidate muon trajectory
  if ( !_outLike ) _freeplanes = (int)muontraj.nmeas();
  //cout<<" _freeplanes ="<< _freeplanes<<"  _meanOcc="<< _meanOcc<<endl;
  
  
  //set to reconstruction mode.
  if ( !MINDfitman::instance().in_rec_mode() )
    MINDfitman::instance().rec_mode();


  ///for free track with atleast 5 hits  
  if ( _meanOcc == 1 ){
    _intType = 2;
    /// if ( muontraj.size() !=0 )
    if ( muontraj.size() >= _min_seed_hits )///?? how many atleast
      ok = muon_extraction( hits, muontraj, hads);
    else ok = false;
    
  } else if ((int)muontraj.nmeas() >= _min_seed_hits ){
    
    if ( _badplanes !=0 ) _intType = 3;
    ok = muon_extraction( hits, muontraj, hads);
  }

  return ok;///for memory leak ??
  
}




//***********************************************************************
bool event_classif::chargeCurrent_analysis(vector<cluster*>& hits,
						      Trajectory& muontraj, vector<cluster*>& hads){
  //***********************************************************************
  m.message("++++ Performing muon reconstruction ++++",bhep::VERBOSE);
  //cout<<"inside muon_extraction_through_PatterRec"<<endl;
  
  muontraj.reset();
  bool ok = true;
  
  _recChi = EVector(3,0);
  _recChi[1] = 100000; //SetLarge dummy value for minChiHadron bit.
  
  
  
  /*//looking into beginning of planes to guess vertex hitno, to exclude backward particle hits and if not found then  excluded_hits = 0; _exclPlanes = 0; i.e, vertGuess =0*/
  _vertGuess = exclude_backwards_particle();
  
  
  ///Zpos of vertex
  _vertGuessZ = hits[_vertGuess]->position()[2];
  
  //cout<<"_vertGuess ="<<_vertGuess<<"  "<<hits[_vertGuess]->position()[2]<<endl;
  
  const dict::Key candHit = "inMu";
  const dict::Key hit_in = "True";
  const dict::Key not_in = "False";
  
  
  
  //add measurements to muontraj 
  if ( _longestSingle >= _min_seed_hits && _endLongPlane > 0.5*(double)_nplanes ){
    while ( _planeIt >= _hitsPerPlane.begin()+_exclPlanes+_min_check && (*_planeIt) == 1 ){
      muontraj.add_measurement( *(*_hitIt) );
      (*_hitIt)->set_name(candHit, hit_in);
      _hitIt--;
      _planeIt--;
      
    }
  }
  
  //cout<<" after Z search ::traj.nmeas()="<<muontraj.nmeas()<<endl; 
  
    
  
  //information from the trajectory 
  _lastIso = (int)muontraj.nmeas();
  m.message(" _lastIso = ",_lastIso,bhep::DETAILED);
  
  
  ///freeplanes inside the candidate muon trajectory
  if ( !_outLike ) _freeplanes = (int)muontraj.nmeas();
  //cout<<" _freeplanes ="<< _freeplanes<<"  _meanOcc="<< _meanOcc<<endl;
  
  
  //set to reconstruction mode.
  if ( !MINDfitman::instance().in_rec_mode() )
    MINDfitman::instance().rec_mode();
  
  
  ///for free track with atleast 5 hits  
  if ( _meanOcc == 1){
    _intType = 2;
    /// if ( muontraj.size() !=0 )
    if((int)muontraj.nmeas() >= _min_seed_hits)
      ok = muon_extraction( hits, muontraj, hads);
    else ok = false;
    
  } else {
    if ( (int)muontraj.nmeas() < _min_seed_hits ) {
      _intType = 5;
      ok = invoke_cell_auto( hits, muontraj, hads);
      if ( !ok ) _failType = 4;
      //ok = false; _failType = 4;
       if(!ok) cout<<" invoke_cell_auto not ok"<<endl;
      
    } else {
      if ( _badplanes !=0 ) _intType = 3;
      ok = muon_extraction( hits, muontraj, hads);
      
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
  //cout<<"inside exclude_backwards_particle"<<endl;

  vector<int>::iterator measIt = _hitsPerPlane.begin();
  vector<double>::iterator zIt;
  int check_planes;


  if ( _nplanes >= 50 ) check_planes = 21;
  else check_planes = (int)( 0.4 * _nplanes ) + 1;
  int excluded_hits = 0;
  _exclPlanes = 0;
  bool found = false;
 
  //for (measIt = _hitsPerPlane.begin()+1;measIt != _hitsPerPlane.end();measIt++){
  ///here zIt starts looking from 2nd plane from beginning upto several planes (40% or 50) and 
  ///measIt for no of hits in each plane,  so to exclude backward particle hits we need to look into both plane and hits
  ///together
  for (zIt = _planeZ.begin()+1;zIt != _planeZ.begin()+check_planes;zIt++,measIt++){
   

    ///total no hits and planes
    excluded_hits += (*measIt);
    _exclPlanes++;
       
   
    ///why backwards if separation b/w two planes >104 mm ??   
    if ( abs( (*zIt) - (*(zIt-1)) ) >  2*_pieceLength ){
      //    //cout<<"found"<<endl;
      found = true;
      break;
    }
   
  }
  
  
  if ( !found ){
    excluded_hits = 0;
    _exclPlanes = 0;
  }
  
  return excluded_hits;
}





//
//***********************************************************************
void event_classif::use_mini_cellAuto(const int occ, Trajectory& muontraj){
  //***********************************************************************
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
  
  m.message("++++ event_classif::muon_extraction +++++++++++++", bhep::DETAILED);
  //  std::cout<<"++++ event_classif::muon_extraction +++++++++++++"<<std::endl;
  bool ok;


  ///if vertex is otherthan 0th hit position, then those hits will be candidate hadron 
  if (_vertGuess != 0)
    for (int i = 0;i < _vertGuess;i++){
      const dict::Key candHit = "inhad";
      const dict::Key hit_in = "True";
      hits[i]->set_name(candHit, hit_in);
      hads.push_back( hits[i] );
    }
  
  State patternSeed;
  
  ok = get_patternRec_seed( patternSeed, muontraj, hits);
  if(!ok)  m.message(" get_patternRec_seed not ok",bhep::DETAILED); 
  
  if ( ok )
    ok = perform_muon_extraction( patternSeed, hits, muontraj, hads);
  if(!ok) m.message("perform_muon_extraction not ok",bhep::DETAILED);


  ///assign the seed state
  if ( ok )
    _seedState = patternSeed;
  
  if(ok) m.message(" event_class: traj nmeas=",muontraj.nmeas()," intType =",_intType,"  && PR seed is=",_seedState,bhep::DETAILED);


  /// to get all the PR seed for reseed_traj inside fitter
  
  if(ok) _vPR_seed.push_back(_seedState);
  
  
  return ok;
}

//***********************************************************************
bool event_classif::get_patternRec_seed(State& seed, Trajectory& muontraj,
					vector<cluster*>& hits) {
  //***********************************************************************
  //std::cout<<"++++ event_classif::get_patternRec_seed +++++++++++++"<<std::endl;
  
  EVector V(6,0); EVector V2(1,0);
  EMatrix M(6,6,0); EMatrix M2(1,1,0);
  
  ///x,y,z values of the 1st node
  V[0] = muontraj.nodes()[0]->measurement().vector()[0];
  V[1] = muontraj.nodes()[0]->measurement().vector()[1];
  V[2] = muontraj.nodes()[0]->measurement().position()[2];

  //direction
  double dqtot = fit_parabola( V, muontraj);
  
  //Momentum. Estimate from empirical extent function.
  ///if cluster contains hits then _Xtent is the diff in z b/w 0th and last hit in cluster 
  ///if not, then _Xtent is the diff in z b/w 1st and one before last node (? not last node)

  if ( hits.size() != 0 ){
    _Xtent = hits[_endLongSing]->position()[2]
      - hits[_vertGuess]->position()[2];
    
  } 
  else {     
    _Xtent = muontraj.nodes()[(int)muontraj.nmeas()-1]->measurement().position()[2]
      - muontraj.nodes()[0]->measurement().position()[2];    
  }
  
  //cout<<"xtent="<<_Xtent<<endl;

  //2 pSeed = 668 + 1.06*_Xtent; //estimate in bhep::MeV, log for fit.
  double pSeed = (9180-6610*_FeWeight) + (-2.76+4.01*_FeWeight)*_Xtent; //best for 3/2 seg

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
    z[(const int)nMeas], u[(const int)nMeas];/// r[(const int)nMeas];
  int minindex = nMeas;
  double minR = 99999.999, pdR = 0.0, sumdq=0;
  double pathlength=0;/// tminR=9999.9;
  

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
  
  
  ///double wFe = geom.get_Fe_prop();
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
  int meansign = (int)(qtilde/fabs(qtilde));

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
  
  double de_dx = -( 12.37 * _FeWeight * pow( mom, 0.099) 
		    + (1 - _FeWeight) * 2.16 * pow( mom, 0.075) );
  de_dx *= bhep::MeV/bhep::cm;
  
  man().geometry_svc().setup().set_volume_property_to_sons("mother","de_dx",de_dx);

}


//***********************************************************************
bool event_classif::perform_kalman_fit(State& seed, Trajectory& track) {
  //***********************************************************************
  //cout<<" I am ******************event_classif :: before perform_kalman_fit,  man().fitting_svc().fit(seed,_track)"<<endl;    
  


  ///fit the track using the seed state                             
  bool ok = man().fitting_svc().fit(seed, track);
   
  ///print first state here, which is the 1st measurement in the trajectory
  if (ok)
    seed = track.state(track.first_fitted_node());

  //cout<<"first fitted node done"<<endl; 
  return ok;
  ///if ok then perform_muon_extraction, if not ok failtype=5
}

//***********************************************************************
bool event_classif::perform_muon_extraction(const State& seed, vector<cluster*>& hits,
					    Trajectory& muontraj, vector<cluster*>& hads) {
  //***********************************************************************


  //Loop through multiple occupancy planes finding the best match to the muon
  //traj and filtering it into the trajectory.

  //cout<<"I am  event_classif::perform_muon_extraction "<<endl; 


  bool ok;
  long ChiMin;
  int nConsecHole = 0;
 
  ///
 
  while (_hitIt >= hits.begin() + _vertGuess) {

    double Chi2[(const int)(*_planeIt)];
    
    for (int iMat = (*_planeIt)-1;iMat >= 0;iMat--, _hitIt--){
      //cout<<iMat<<"  *_planeIt="<<*_planeIt<<"  before ******************event_classif:: perform_muon_extraction,  man().fitting_svc().filter(*(*(_hitIt)), seed, muontraj)"<<endl;                                   	
      ok = man().matching_svc().match_trajectory_measurement(muontraj, (*(*_hitIt)), Chi2[iMat]);
    
      if ( !ok )
	Chi2[iMat] = 999999999;
    }
    
    ChiMin = TMath::LocMin( (const int)(*_planeIt), Chi2);
    
    for (int iFil = 0;iFil < (*_planeIt);iFil++){
    

      ///For min Chi2 the hit will be added to the Trajectory, else considered as Hadron
      if ( iFil == (int)ChiMin) {


	ok = man().fitting_svc().filter(*(*(_hitIt+iFil+1)), seed, muontraj);
	//cout<<" After***************event_classif:: perform_muon_extraction,  filter(*(*(_hitIt+iFil+1)), seed, muontraj)"<<endl;                                   	

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
    
    if ( nConsecHole > _max_consec_missed_planes ) {
      if ( _endProj ) check_forwards( seed, hits, muontraj );
      if ( muontraj.nmeas() < _min_plane_prop*(double)(_nplanes-_badplanes) ){
	_failType = 6;
	return false;
      } else {
	sort_hits( hits, muontraj, hads );
	return true;
      }
    }
   
    //cout<<" Next  *_planeIt="<<*_planeIt<<endl;
  }

  if ( _endProj ) check_forwards( seed, hits, muontraj );
  
  return true;
}







//***********************************************************************
void event_classif::check_forwards(const State& seed, vector<cluster*>& hits,
				   Trajectory& muontraj){
  //***********************************************************************

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
      m.message(" I am ******************event_classif:: check_forward,  man().fitting_svc().filter(*(*(_hitIt -2)), seed, muontraj",bhep::DETAILED);                                   	
      ok = man().fitting_svc().filter(*(*(_hitIt-2)), seed, muontraj);
      
      if ( !ok ) (*(_hitIt-2))->set_name(skipped,hit_in);
      else (*(_hitIt-2))->set_name(candHit, hit_in);
      (*(_hitIt-1))->set_name(skipped,hit_in);
    } else {
      ok = man().fitting_svc().filter(*(*(_hitIt-1)), seed, muontraj);
      m.message(" I am ******************event_classif ::check_forward,  man().fitting_svc().filter(*(*(_hitIt -1)), seed, muontraj",bhep::DETAILED);                                   	
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
  _badplanes = counter;
  
}




//***********************************************************************
bool event_classif::invoke_cell_auto(vector<cluster*>& hits,
				     Trajectory& muontraj, vector<cluster*>& hads){
  //***********************************************************************
  //uses cellular automaton to try and retrieve more complicated events.
  cout<<"+++ Performing Cellular Automaton analysis +++"<<endl;
  m.message("+++ Performing Cellular Automaton analysis +++",bhep::VERBOSE);
  
  bool ok;
  
  /// For radail search
  ///  if ( (_nplanes-_exclPlanes < __min_hits) && (_vRadCount.size() - _exclPlanes < _min_hits) ) return false;
  

  if (_nplanes-_exclPlanes < _min_hits)  return false;
  //cout<<"_nplanes ="<<_nplanes<<"   _exclPlanes="<<_exclPlanes<<"   _min_hits="<<_min_hits<<endl;
  //TEST!!
  measurement_vector hit_meas;
  get_cluster_meas( hits, hit_meas );
  
  
  //vector of trajectories from CA
  std::vector<Trajectory*> trajs;
  ok = man().matching_svc().find_trajectories( hit_meas, trajs );
  //cout<<" traj ="<<_trajs.size()<<endl;
  //cout<<" I am here 1"<<endl;


  ///if no track found
  if ( !ok || trajs.size() == 0) return false;
  //cout<<" I am here 2"<<endl;
  // output_results_tree( hits, trajs );
  

  if ( trajs.size() == 1 ) {
    //cout<<" I am here 3"<<endl;
    if ( trajs[0]->nmeas() >= _min_hits ){
      //cout<<" I am here 4"<<endl;
      muontraj.reset();
      
      muontraj = *trajs[0];
      

      ///to assign candHit & hadrons from muontraj
      sort_hits( hits, muontraj, hads);
      
      ok = true;
    } else ok = false;
    
  } else {
    
    ok = sort_trajs( muontraj, trajs);
    //cout<<" I am here 5"<<endl;
    if ( ok )
      sort_hits( hits, muontraj, hads);

  }
  
  

  if ( ok )
    delete_bad_trajs( muontraj, trajs );
  else  stc_tools::destroy( trajs );
  
  
  /// delete trajectories at the end of the event
  //  stc_tools::destroy( _trajs );
  // else if (_hitIt == hits.end() - 1 )stc_tools::destroy( _trajs );
 
  if (ok) cout<<"+++ traj infor +++  , nmeas="<<muontraj.nmeas()<<endl;
  //cout<<" end of CA"<<endl;
  m.message("+++ End of Cellular automaton +++",bhep::VERBOSE);
  
  return ok;
}

//***********************************************************************
void event_classif::get_cluster_meas(const vector<cluster*>& hits,
				     measurement_vector& meas)
  //***********************************************************************
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


  /* /// At the end of the event destroy trajectories only valid for SINGLE TRACK FINDING CODE
  vector<cluster*> hits;
  while ( !found && dIt != trajs.end() ){
    if ( muontraj.quality( traj_index ) == (*dIt)->quality( traj_index ) ){
      found = true;

      trajs.erase( dIt );
    }
    dIt++;
  }
  if( _hitIt == hits.end()-1) stc_tools::destroy( trajs ); */ 
 

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

    if ( (int)(*it1)->nmeas() >= _min_hits )
      trajs2.push_back( (*it1) );

  }

  /* if ( trajs2.size() == 0 ) return false;

  return true;*/

  /// if trajs2 exist then only true 06072012
  if ( trajs2.size() == 0 ) return false;
  else  return true;
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
    
    if ( ok1 && temp_traj.quality() < _chi2_max ){

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

    while ( it2 != temp_trajs.end() && coincidence < _max_coincedence ){

      vector<Node*>& node1 = (*it2)->nodes();

      coincidence = compare_nodes( node1, node2);

      it2++;
    }
    
    if ( coincidence < _max_coincedence )
      temp_trajs.push_back( (*it1) );

  }
  trajs.clear();
  
  if ( temp_trajs.size() == 1 ) {
    
    muontraj.reset();

    muontraj = *temp_trajs[0];

  } else {

    /// sort trajectory by hits instead of length
     sort( temp_trajs.begin(), temp_trajs.end(), sortTrajByHits() );
     muontraj = *temp_trajs[0];
    // select_trajectory( temp_trajs, muontraj);

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
  //std:://cout<<"I am in output_liklihood_info"<<std::endl;  
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
    // //cout<<"in outlke::_nhit="<<_nhit<<" _visEng="<<_visEng<< "  _freeplanes="<<_freeplanes<<"   _occ["<<counter<<"]="<< _occ[counter]<<"  _plEng[counter]="<<_plEng[counter]<<endl;   
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

    ///if it is a candidate muon
    if ( (*hitIt3)->names().has_key( candHit ) ){
      if ( (*hitIt3)->name( candHit ).compare("True") == 0 ){
	_trajhit++;
	//_trajEng += bhep::double_from_string( (*hitIt3)->name( Edep ) ) * bhep::GeV;
	_trajEng += (*hitIt3)->get_eng() * bhep::MeV;
	_trajEngPlan[counter] = (*hitIt3)->get_eng() * bhep::MeV;
	// if ( (*hitIt3)->name("MotherParticle").compare("mu+") == 0
	// 	     || (*hitIt3)->name("MotherParticle").compare("mu-") == 0 )

	///if the muon candidate is true muon
	if ( (*hitIt3)->get_mu_prop() > 0.8 )//still to be decided.
	  _trajpur++;
	_trclusthits[counter] = (*hitIt3)->get_nVox();
	////cout<<"_trajhit="<<_trajhit<<"  _trajEng="<<_trajEng<<" _trajEngPlan["<<counter<<"]="<<_trajEngPlan[counter]<<"    _trclusthits"<<"["<<counter<<"]= "<<_trclusthits[counter] <<endl;
	counter--;
      }
    }
  }
  _trajpur /= _trajhit;

 
 
}

//***********************************************************************
void event_classif::out_like(){
  //***********************************************************************
 
  int fillcatch = _likeTree->Fill();
  
}




//***********************************************************************
//Temp function for likelihoods with tru interaction in tree.
void event_classif::set_int_type(const string name){
  //***********************************************************************

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


//***********************************************************************
double event_classif::correctEdep(double edep, double X, double Y)
  //***********************************************************************
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






/***************************************************************************************/
double event_classif::RangeMomentum(double length){
  /***************************************************************************************/
  // Computing the momentum of a track based on the range
  // Parameters calculated from "Atomic Data and Nuclear Data Tables 78, 183-356 (2001)"
  double p = (_FeWeight*(0.011844*bhep::GeV * pow((length/bhep::cm),1.03235))
	      + (1- _FeWeight)*(0.0023705067*bhep::GeV* pow(length/bhep::cm,1.00711577)));
  return p;
  // Should be compared to G4 simulation directly.
}




//***********************************************************************
void event_classif::fill_traj_info( Trajectory& muontraj) {
  //***********************************************************************
  m.message("++++event_classif::fill_traj_info ", bhep::VERBOSE);

  muontraj.set_quality("intType",_intType);
  muontraj.set_quality("failType",_failType);
  // int intType = (int)(muontraj.quality("intType"));
  muontraj.set_quality("nplanes",_nplanes);
  muontraj.set_quality("freeplanes",_freeplanes);
  muontraj.set_quality("badplanes",_badplanes);
  muontraj.set_quality("longestSingle",_longestSingle);
  muontraj.set_quality("radialLongest",_radialLongest);
  muontraj.set_quality("xtent",_Xtent);
  muontraj.set_quality("vertZ",_vertGuessZ);

  /// for CA tracks it is not sure the free section belongs to the same traj or not
  if((int)muontraj.quality("intType")!=5) muontraj.set_quality("lastIso", _lastIso);
  else muontraj.set_quality("lastIso", 0);
  
  cout<<" event_clss failType ="<<_failType<<endl;;
    
  //reset the containers
  _nplanes = 0;
  _meanOcc = 0;
  _badplanes = 0;
  _longestSingle = 0;
  _endProj = false;
  _hitsPerPlane.clear();
  _energyPerPlane.clear();
  _planeZ.clear();
  _vRadCount.clear();
  _radialLongest = 0;
  _freeplanes = 0;
  _endLongSing = 0;
  _endLongPlane = 0;
  _lastIso = 0;
  ///25062012
  _intType = 0;
  _Xtent = 0;
  _vertGuess = 0;
  _vertGuessZ = 0;
  _failType = 0;

}
