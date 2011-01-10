/* -*- mode: c++ -*- */
#ifndef _mind_setup___
#define _mind_setup___

//#include <recpack/RecpackManager.h>

#include <mind/SetupSk.h>


//#include <string.h>

using namespace Recpack;

class MINDsetup: public NSetupSk{

public:
    
  MINDsetup();
  
  virtual ~MINDsetup();
    
  Setup& setup();
   
  // void init(bhep::gstore store,bhep::sstore,
// 	    bhep::prlevel level=bhep::NORMAL);
  void init(bhep::gstore store, bhep::prlevel level=bhep::NORMAL);
  
  //info to build virtual planes
  
  double getPlaneX(){return MIND_x;};
  double getPlaneY(){return MIND_y;};
  double getPlaneZ(){return MIND_z;};
  EVector getXaxis(){return xaxis;};
  EVector getYaxis(){return yaxis;};
  EVector getZaxis(){return zaxis;};
  string getMeasType(){return meastype;};
  int getMeasDim(){return meas_dim;} 
  EVector getResolution(){return resolution;} 
  EMatrix getCov(){return cov;};
  double get_Fe_prop(){return _wFe;}
  double& getDeDx(){return de_dx;}
  void setDeDx(double d){de_dx = d;}

protected:
    
  void readParam();
  void createGeom();
  void add_slab(int plane, const dict::Key det_vol);
  // void add_slab(int plane, Volume& det_vol);
  void setResolution();
  void addProperties();
  

protected:
  
  // parameter store

  bhep::gstore _pstore;
  
  // detector axis

  EVector xaxis,yaxis,zaxis;
    
  // -------------------------------------------------------------//
  //                       |  DIMENSIONS |                        //
  // -------------------------------------------------------------//
    
  //--------------------------- VOLUMES -------------------------//
     
  double MOTHER_x, MIND_x;
  double MOTHER_y, MIND_y;
  double MOTHER_z, MIND_z;
  double IRON_z, SCINT_z, AIR_z;
  double rel_denAS, rel_denSI;//AIR/Scint, Scint/Fe.
  int nScint;
  int _npieces;
  double _pieceWidth;

  // -------------------------------------------------------------//
  //                         |  PHYSICS |                         //
  // -------------------------------------------------------------//

  //------------------------ MAGNETIC FIELD ----------------------//
    
  //double B_int;
  EVector BField;
   
  //-------------------------------------------------------------//
  
  //------------------- PROPERTIES OF MATERIALS -----------------//
    
  double X0Fe, X0Sc, X0AIR, X0Eff;//members for if/when geom more strict.
  double _wFe;
  double de_dx;
  EVector _zaxis;
  
  //-------------------------------------------------------------//
  

  //------------------------- MEASUREMENTS ----------------------//
     
  int meas_dim;
  string meastype;
  EVector resolution;
  EMatrix cov;
  double resx,resy,resz;



};



#endif 
