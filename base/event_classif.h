/* -*- mode: c++ -*- */
#ifndef _event_classif___
#define _event_classif___

#include <mind/MINDsetup.h>
#include <mind/SetupSk.h>
#include <mind/MINDfitman.h>
#include <mind/cluster.h>

#include <recpack/RecpackManager.h>
#include <recpack/RayTool.h>
#include <recpack/KalmanFitter.h>
#include <recpack/ParticleState.h>
#include <recpack/CellularAutomaton.h>

#include <bhep/gstore.h>

#include <TGraph.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>

class event_classif{
  
public:
  
  event_classif();
  
  virtual ~event_classif();
  
  //-------------- main functions --------------//
  //void initialize(const bhep::gstore& pstore, bhep::prlevel vlevel, Setup& det, double wFe);
  void initialize(const bhep::gstore& pstore, bhep::prlevel vlevel, double wFe);
  bool execute(vector<cluster*>& hits,
	       Trajectory& muontraj, vector<cluster*>& hads);
  void finalize();
  //-------------------------------------------//
  
  //Grabbers for monitoring etc.
  int get_int_type(){ return _intType; }
  int get_vertex(){ return _vertGuess; }
  int get_fail_type(){ return _failType; }
  EVector& get_PatRec_Chis(){ return _recChi; }
  State& get_patRec_seed(){ return _seedState; }
  int get_last_iso(){ return _lastIso; }
  double get_vis_eng(){ return _visEng; }
  int get_planes(){ return _nplanes; }
  int get_free_planes(){ return _freeplanes; }
  //
  
  bool get_plane_occupancy(vector<cluster*>& hits);
  void assess_event(vector<cluster*>& hits);
  
  //Tempory for likelihoods.
  void set_int_type(const string name);
  //
  
protected:
  
  void readParam();
  //void set_extract_properties(Setup& det);
  void reset();
  
  //Functions to be performed on CC mu candidates.  
  bool chargeCurrent_analysis(vector<cluster*>& hits,
			      Trajectory& muontraj, vector<cluster*>& hads);
  int exclude_backwards_particle();
  bool muon_extraction(vector<cluster*>& hits,
		       Trajectory& muontraj, vector<cluster*>& hads);
  bool get_patternRec_seed(State& seed, Trajectory& muontraj, vector<cluster*>& hits);
  void fit_parabola(EVector& vec, Trajectory& track);
  void set_de_dx(double mom);
  bool perform_kalman_fit(State& seed, Trajectory& track);
  bool perform_muon_extraction(const State& seed, vector<cluster*>& hits,
			       Trajectory& muontraj, vector<cluster*>& hads);
  void check_forwards(const State& seed, vector<cluster*>& hits,
		      Trajectory& muontraj);
  void use_mini_cellAuto(const int occ, Trajectory& muontraj);
  //specific functions using cellular automaton.
  bool invoke_cell_auto(vector<cluster*>& hits,
			Trajectory& muontraj, vector<cluster*>& hads);
  void sort_hits(vector<cluster*>& hits, Trajectory& muontraj, vector<cluster*>& hads);
  void get_cluster_meas(const vector<cluster*>& hits, measurement_vector& meas);
  void delete_bad_trajs(const Trajectory& muontraj, vector<Trajectory*>& trajs);
  bool sort_trajs(Trajectory& muontraj, vector<Trajectory*>& trajs);
  bool reject_small(vector<Trajectory*>& trajs, vector<Trajectory*>& trajs2);
  bool reject_high(vector<Trajectory*>& trajs, vector<Trajectory*>& trajs2);
  bool reject_final(vector<Trajectory*>& trajs, Trajectory& muontraj);
  double compare_nodes(const vector<Node*>& n1, const vector<Node*>& n2);
  void select_trajectory(vector<Trajectory*>& trajs, Trajectory& muontraj);
  //

  RecpackManager& man(){
    return MINDfitman::instance().manager();}
  
  bhep::gstore _infoStore;
  
  bhep::messenger m;
  
  //Members to store plane occupancy and mean energy.
  int _nplanes;
  double _meanOcc;
  vector<int> _hitsPerPlane;
  vector<double> _energyPerPlane;
  vector<double> _planeZ;
  double _tolerance; //required 'closeness' to be considered in plane.
  double _voxEdge; //Voxel edge size to check for badly reconstructed points at end.
  
  //integer for type candidate (NC etc.)
  int _intType;
  int _lastIso;
  
  //interator for hits and container for estimated 'vertex' hit.
  vector<cluster*>::iterator _hitIt;
  vector<int>::iterator _planeIt;
  int _vertGuess;
  int _exclPlanes;
  int badplanes;
  int _longestSingle;//length (in hits) of longest 'free' section.
  int _endLongSing;//End point of the above (position in hit vector).
  int _endLongPlane;//Plane position of above;
  double _maxBlobSkip;//max proportion of planes for basic skip.
  double _minBlobOcc;
  bool _endProj;//bit to tell if a forwards projection is needed
  
  //Monitoring variables.
  int _failType;
  EVector _recChi;
  State _seedState;

  double FeWeight;
  int max_consec_missed_planes;
  double min_plane_prop;
  int min_seed_hits;
  int min_check;
  int min_hits;
  double chi2_max;
  double max_coincedence;
  double _pieceLength;

  //Output Likilihood info?
  bool _outLike;
  TFile *_outFileEv;
  TTree *_likeTree;
  //

  int _nhit, _trajhit, _truInt;
  int _freeplanes, _occ[1000], _trclusthits[1000];
  double _visEng;
  double _trajpur, _trajEng;
  double _plEng[1000], _trajEngPlan[1000];

  void set_branches();
  void output_liklihood_info(const vector<cluster*>& hits);
  void traj_like(const vector<cluster*>& hits, const Trajectory& muontraj);
  void out_like();
  //


};

class chiSorter{
public:
  bool operator()(const Trajectory* T1, const Trajectory* T2){
    if ( T2->quality() > T1->quality() ) return true;
    return false;
  }
};

#endif
