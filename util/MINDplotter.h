/* -*- mode: c++ -*- */
#ifndef _mind_plotter___
#define _mind_plotter___

#include <bhep/event.h>
#include <bhep/particle.h>

#include <TFile.h>
#include <TTree.h>

#include <recpack/RecpackManager.h>
//#include <recpack/Measurement.h>
#include <recpack/Ring.h>

#include <mind/fitter.h>
#include <mind/cluster.h>

using namespace std;
using namespace Recpack;
using namespace bhep;

/* A Class with functions to plot the results of fits to MIND data */

class MINDplotter{

 public:

  MINDplotter(const bhep::gstore& pstore, bhep::prlevel vlevel=bhep::NORMAL);

  virtual ~MINDplotter();

  //Main functions for output initialization.
  void initialize(string outFileName, bool patRec, bool clust);
  void execute(fitter& Fit, const bhep::event& evt, bool success);
  void finalize();
//   //
  //To calculate expected vertex position give trajectory and vertex location.
  bool extrap_to_vertex(const Trajectory& traj, 
			const bhep::Point3D& vertexLoc, fitter& fitObj, State& ste);

  //To calculate pos, charge, direction and momentum without extrapolating to vertex //tapasi
  bool fill_kinematics(const Trajectory& traj, fitter& fitObj, State& ste);

  /*Requires a vector of the node of interest, corresponding covariance Matrix
    and the event under study to calculate either position or momentum pulls */
  void position_pulls();
  void momentum_pulls();
  void direction_pulls();

  /* Function to find maximum local chi2 in a particular trajectory
     and enter that value in a histogram */
  void max_local_chi2(const Trajectory& traj);

  /* Function to plot stats about pattern recogntion, 1: no clust, 2: clust */
  void patternStats1(fitter& Fit);
  void patternStats2(fitter& Fit);
  double correctEdep(double edep, double X, double Y);

  /*Function to record quality of hadron fit*/
  void hadron_direction(fitter& fit);

  /*tru interaction type*/
  void get_tru_int_type(const bhep::event& evt);

  // void extract_State_History(std::vector<EVector> steHistory);

protected:

  bhep::gstore store;

  bhep::prlevel level;
    
  bhep::messenger m;

  TFile *outFile;

  TTree *statTree;

  bool _patR, _clu;

  double _detX, _detY, _WLSAtten;

private:

  EVector vert;
  EMatrix vertMat;
  
  bhep::particle* _truPart;
  int _evNo;
  bool _Fit;
  bool _reFit;
  int _fail;
  int _intType;
  double _vert[3];
  double _nuEng;
  int _charm[2];
  int _npi[3];
  int _nk[3];
  double _Q2;
  int _pdg[3];
  double _visEng;
  double _engTraj;
  double _engvar[4];
  //Information on hadrons.
  double _hadP[3];
  double _engTrans;
  int _nhad[2];
  double _chadP[4];
  double _nhadP[4];
  double _hadE[2];
  //double _haddot;
  double _X[6][2];
  double _Th[6][2];
  double _qP[6];
  double _leng;
  double _rangqP[3];
  int _Q[3];
  double _Chi[2];
  int _plns[2];
  int _nhits;
  int _hitType[5];
  double _XPos[5000];
  double _YPos[5000];
  double _ZPos[5000];
  double _Edep[5000];
  bool _mus[5000];
  bool _cand[5000];
  bool _node[5000];
  bool _had[5000];
  double _nchi2[5000];
  double _pChi[3];
  double _dchi2p;
  bool _isProton;
  // double hist[4];
  // bool passrec[4];
  //TString _intName;
  int _truInt;
  double IRON_z;
  double SCINT_z;
  double AIR_z;
  int nScint;
  double _pieceWidth;
  int _npieces;
  double MIND_z;
  double rel_denSI;
  double rel_denAS;
  double _wFe;
  

  void define_tree_branches();
  //1: Old method. 2: with clustering etc.
  bool extract_true_particle1(const bhep::event& evt, fitter& Fit);
  bool extract_true_particle2(const bhep::event& evt, fitter& Fit);
  void add_to_hads(const bhep::particle& part);
  //
  int _truMuHitIndex[2]; 
  int HitIndex; 
};

#endif
