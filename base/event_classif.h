/* -*- mode: c++ -*- */
#ifndef _event_classif___
#define _event_classif___

#include <mind/MINDsetup.h>
#include <mind/SetupSk.h>

#include <recpack/RecpackManager.h>
#include <recpack/RayTool.h>
#include <recpack/KalmanFitter.h>
#include <recpack/ParticleState.h>

#include <bhep/gstore.h>

#include <TGraph.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>

using namespace bhep;

class event_classif{
  
public:
  
  event_classif();
  
  virtual ~event_classif();
  
  //-------------- main functions --------------//
  bool initialize(const bhep::gstore& pstore, bhep::prlevel vlevel, Setup& det, double wFe); //needs to take store with all info for pat rec/more?
  bool execute(measurement_vector& hits,
	       Trajectory& muontraj, measurement_vector& hads); //more arguments needed?
  bool finalize();
  //-------------------------------------------//
  
  //Grabbers for monitoring etc.
  int get_int_type(){ return _intType; }
  int get_vertex(){ return _vertGuess; }
  int get_fail_type(){ return _failType; }
  EVector& get_PatRec_Chis(){ return _recChi; }
  State& get_patRec_seed(){ return _seedState; }
  int get_last_iso(){ return _lastIso; }
  //
  
protected:
  
  void readParam();
  void set_extract_properties(Setup& det);
  void reset();
  
  bool get_plane_occupancy(measurement_vector& hits);
  
  //Functions to be performed on CC mu candidates.  
  bool chargeCurrent_analysis(measurement_vector& hits,
			      Trajectory& muontraj, measurement_vector& hads);
  int exclude_backwards_particle();
  bool muon_extraction(measurement_vector& hits,
		       Trajectory& muontraj, measurement_vector& hads);
  bool get_patternRec_seed(State& seed, Trajectory& muontraj, measurement_vector& hits);
  void fit_parabola(EVector& vec, Trajectory& track);
  void set_de_dx(double mom);
  bool perform_kalman_fit(State& seed, Trajectory& track);
  bool perform_muon_extraction(const State& seed, measurement_vector& hits,
			       Trajectory& muontraj, measurement_vector& hads);
  //
  
  RecpackManager& man(){return _man;}
  
  bhep::gstore _infoStore;
  
  //bhep::prlevel level;
  
  bhep::messenger m;
  
  RecpackManager _man;
  
  //Members to store plane occupancy and mean energy.
  int _nplanes;
  double _meanOcc;
  vector<int> _hitsPerPlane;
  vector<double> _energyPerPlane;
  double _tolerance; //required 'closeness' to be considered in plane.
  
  //integer for type candidate (NC etc.)
  int _intType;
  int _lastIso;
  
  //interator for hits and container for estimated 'vertex' hit.
  measurement_vector::iterator _hitIt;
  vector<int>::iterator _planeIt;
  int _vertGuess;
  
  //Monitoring variables.
  int _failType;
  EVector _recChi;
  State _seedState;

  //Properties for muon extraction
  string model;
  string kfitter;

  double patRec_maxChi;
  double FeWeight;
  int patRec_max_outliers;
  int max_consec_missed_planes;
  int min_seed_hits;
  int min_check;

  int vfit,vnav,vmod,vmat,vsim;

  //Output Likilihood info?
  bool _outLike;
  TFile *_outFileEv;
  TTree *_likeTree;

  int _nhit;
  int _Occ[500], _freeplanes;
  double _EngP[500];

  void set_branches();
  void output_liklihood_info();
  //

};

#endif
