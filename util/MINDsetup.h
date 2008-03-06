/* -*- mode: c++ -*- */
#ifndef _mind_setup___
#define _mind_setup___

#include <recpack/RecpackManager.h>
#include <nemora/NSetupSk.h>
//#include <string.h>

using namespace Recpack;

class MINDsetup: public NSetupSk{

public:
    
  MINDsetup();
  
  virtual ~MINDsetup();
    
  Setup& setup();
   
  void init(bhep::gstore store,bhep::sstore,
	    bhep::prlevel level=bhep::NORMAL);
  
  //info to build virtual planes
  
  double getPlaneX(){return MIND_x;};
  double getPlaneY(){return MIND_y;};
  EVector getXaxis(){return xaxis;};
  EVector getYaxis(){return yaxis;};
  EVector getZaxis(){return zaxis;};
  string getMeasType(){return meastype;};
  int getMeasDim(){return meas_dim;} 
  EVector getResolution(){return resolution;} 
  EMatrix getCov(){return cov;};
  
protected:
    
  void readParam();
  void createGeom();
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
    

  // -------------------------------------------------------------//
  //                         |  PHYSICS |                         //
  // -------------------------------------------------------------//

  //------------------------ MAGNETIC FIELD ----------------------//
    
  double B_int;
  EVector BField;
   
  //-------------------------------------------------------------//
  
  //------------------- PROPERTIES OF MATERIALS -----------------//
    
  double X0;
  
  //-------------------------------------------------------------//
  

  //------------------------- MEASUREMENTS ----------------------//
     
  int meas_dim;
  string meastype;
  EVector resolution;
  EMatrix cov;
  double resx,resy,resz;



};



#endif 
