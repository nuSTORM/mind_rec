/* -*- mode: c++ -*- */
#ifndef _mind_plotter___
#define _mind_plotter___

#include <bhep/event.h>
#include <bhep/particle.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h> //TH2F.h>
#include <TProfile.h>

#include <recpack/RecpackManager.h>
#include <recpack/Measurement.h>
#include <recpack/Ring.h>
#include <mind/fitter.h>

using namespace std;
using namespace Recpack;
using namespace bhep;

/* A Class with functions to plot the results of fits to MIND data */

class MINDplotter{

 public:

  MINDplotter();

  virtual ~MINDplotter();

  //Main functions for output initialization.
  bool initialize(TString outFileName, bhep::prlevel vlevel=bhep::NORMAL);
  bool execute(const EVector& V, 
	       const EMatrix& M, const bhep::event& evt);
  bool finalize();
  //
  //To calculate expected vertex position give trajectory and vertex location.
  bool extrap_to_vertex(const Trajectory& traj, 
			const bhep::Point3D& vertexLoc, fitter* fitObj, State& ste);

  /*Requires a vector of the node of interest, corresponding covariance Matrix
    and the event under study to calculate either position or momentum pulls */
  void position_pulls(const EVector& V, const EMatrix& M, const bhep::event& evt);
  void momentum_pulls(const EVector& V, const EMatrix& M);
  void direction_pulls(const EVector& V, const EMatrix& M);

  /* Functions to calculate charge misid efficiency as function of momentum
     and recorded hits */
  void momentum_efficiency(const EVector& V, const EMatrix& M);
  void hit_efficiency(const EVector& V, const EMatrix& M);

  /* Function to find maximum local chi2 in a particular trajectory
     and enter that value in a histogram */
  void max_local_chi2(const Trajectory& traj, double maxChi, const EVector& V);

protected:

  bhep::prlevel level;
    
  bhep::messenger m;

  TFile *outFile;

  //Histograms.
  TH1F* FitX;
  TH1F* truX;
  TH1F* xPull;
  TH1F* FitY;
  TH1F* truY;
  TH1F* yPull;
  TH1F* FitTx;
  TH1F* truTx;
  TH1F* tXpull;
  TH1F* FitTy;
  TH1F* truTy;
  TH1F* tYpull;
  TH1F* FitP;
  TH1F* truMom;
  TH1F* pPull;
  TH1F* pSpec;
  TH1F* misIDp;
  TH1F* hitSpec;
  TH1F* misIDhit;
  TH1F* locChi;
  TProfile* locVp;
  TProfile* trajVp;

private:

  //Particle momentum with corresponding error and charge.
  double p_;
  double d_p_;
  int tru_q_;
  bool truP;
  bhep::particle* part;

  void lowChi_MisID_track(const Trajectory& traj, double trajMax);
  void highChi_ID_track(const Trajectory& traj, double trajMax);

  bool extract_true_particle(const EVector& V, 
			     const EMatrix& M, const bhep::event& evt);

};

#endif
