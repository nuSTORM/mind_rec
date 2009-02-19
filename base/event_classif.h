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

using namespace bhep;

class event_classif{

public:

  event_classif(const bhep::gstore& pstore, bhep::prlevel vlevel);

  virtual ~event_classif();

  //-------------- main functions --------------//
  bool initialize(Setup& det); //needs to take store with all info for pat rec/more?
  bool execute(measurement_vector& hits,
	       Trajectory& muontraj, measurement_vector& hads); //more arguments needed?
  bool finalize();
  //-------------------------------------------//

  int get_int_type(){ return _intType; }

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
  bool perform_kalman_fit(State& seed, Trajectory& track);
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

  int _intType;

  //interator for hits and container for estimated 'vertex' hit.
  measurement_vector::iterator _hitIt;
  int _vertGuess;

  //Properties for muon extraction
  string model;
  string kfitter;

  double patRec_maxChi;
  int patRec_max_outliers;
  int max_consec_missed_planes;

  int vfit,vnav,vmod,vmat,vsim;

};

#endif
