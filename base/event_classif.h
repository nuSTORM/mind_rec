/* -*- mode: c++ -*- */
#ifndef _event_classif___
#define _event_classif___

#include <mind/MINDsetup.h>
#include <mind/SetupSk.h>
#include <mind/MINDfitman.h>

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

//using namespace bhep;

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
  double get_vis_eng(){ return _visEng; }
  int get_planes(){ return _nplanes; }
  int get_free_planes(){ return _freeplanes; }
  //
  
  bool get_plane_occupancy(measurement_vector& hits);
  
protected:
  
  void readParam();
  //void set_extract_properties(Setup& det);
  void reset();
  
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
  //specific functions using cellular automaton.
  bool invoke_cell_auto(measurement_vector& hits,
			Trajectory& muontraj, measurement_vector& hads);
  void sort_hits(measurement_vector& hits, Trajectory& muontraj, measurement_vector& hads);
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
  
  //RecpackManager _man;
  
  //Members to store plane occupancy and mean energy.
  int _nplanes;
  double _meanOcc;
  vector<int> _hitsPerPlane;
  vector<double> _energyPerPlane;
  vector<double> _planeZ;
  double _tolerance; //required 'closeness' to be considered in plane.
  //double _max_sep; //maximum transverse separation for cell auto neighbour.
  //int _max_traj; //maximum no. of trajectories from cell auto.
  
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
  //string model;
  //string kfitter;

  //double patRec_maxChi;
  double FeWeight;
  //int patRec_max_outliers;
  int max_consec_missed_planes;
  int min_seed_hits;
  int min_check;
  int min_hits;
  double chi2_max;
  double max_coincedence;

  //int vfit,vnav,vmod,vmat,vsim;

  //Output Likilihood info?
  bool _outLike;
  TFile *_outFileEv;
  TTree *_likeTree;
  //

  int _nhit, _trajhit;
  int _freeplanes, _occ[500];
  double _visEng;
  double _trajpur, _trajEng;
  double _plEng[500], _trajEngPlan[500];

  void set_branches();
  void output_liklihood_info(const measurement_vector& hits);
  void traj_like(const measurement_vector& hits, const Trajectory& muontraj);
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
