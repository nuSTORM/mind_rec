#include <MINDsetup.h>
#include <recpack/string_tools.h>


//*************************************************************
MINDsetup::MINDsetup() {
//*************************************************************
  
}


//*************************************************************
MINDsetup::~MINDsetup() {
//*************************************************************

}


//*************************************************************
void MINDsetup::init(bhep::gstore pstore,bhep::sstore gstore,
		       bhep::prlevel level) {
//*************************************************************
    
    _msetup=bhep::messenger(level);
    
    _msetup.message("++MINDsetup Messenger generated++",bhep::VERBOSE);

    _pstore=pstore; _store=gstore;
    
    readParam();

    //--------------- generate recpack setup -----------//
    
    _gsetup.set_name("main");
    
    // create volumes and surfaces

    createGeom();
    
    // define detector resolutions

    setResolution();
    
    // add properties to volumes and surfaces

    addProperties();
    
    _msetup.message("++ Setup has been generated !! ++",bhep::VERBOSE);

    _msetup.message("MIND Setup:", _gsetup,bhep::VERBOSE);
}


//*************************************************************
Setup& MINDsetup::setup() {
//*************************************************************
     
    return _gsetup;

}


//*************************************************************
void MINDsetup::createGeom(){
//*************************************************************    
    
  _msetup.message("+++ CreatGeom function +++",bhep::VERBOSE);

  //----- axes for definition of volumes and surfaces ----//

  xaxis=EVector(3,0); xaxis[0] = 1.; 
  yaxis=EVector(3,0); yaxis[1] = 1.; 
  zaxis=EVector(3,0); zaxis[2] = 1.;
    
  //----------- Create a mandatory mother volume ---------//

  EVector pos(3,0); pos[2]=MOTHER_z/2;
  
  Volume* mother = new Box(pos,xaxis,yaxis,
			   MOTHER_x/2,MOTHER_y/2,MOTHER_z/2);
    
  _msetup.message("Mother volume generated",bhep::VERBOSE);

  // add mother volume
  
  _gsetup.set_mother(mother);
   
  _msetup.message("Mother added to setup",bhep::VERBOSE);

  // Create detector volume
  
  pos[2]=MIND_z/2;

  const dict::Key vol_name = "Detector";
  
  Volume* det = new Box(pos,xaxis,yaxis,MIND_x/2,MIND_y/2,MIND_z/2);
    
  _msetup.message("MIND volume generated",bhep::VERBOSE);

  // add volume

  _gsetup.add_volume("mother",vol_name,det);
   

}



//*************************************************************
void MINDsetup::setResolution(){
//*************************************************************    
        
    /*
      
    */
  
  _msetup.message("+++ setResolution function +++",bhep::VERBOSE);

  resolution = EVector(meas_dim,0);
  resolution[0] = resx*mm;  
  resolution[1] = resy*mm;  
  
    
  // cov of measurements

  cov = EMatrix(meas_dim,meas_dim,0);
  for (size_t i = 0; i < (size_t)meas_dim; i++) 
    cov[i][i] = resolution[i]*resolution[i];
    
}



//*************************************************************
void MINDsetup::addProperties(){
//*************************************************************    
  
  _msetup.message("+++ addProperties function +++",bhep::VERBOSE);

  //-------------------- magnetic field ------------------//
    
  BField = EVector(3,0);
  BField[1] = B_int;

  _gsetup.set_volume_property("mother","BField",BField);
  _msetup.message("+++B Field added to MOTHER:",BField,bhep::VERBOSE);
    
  const dict::Key vol_name = "Detector";
  _gsetup.set_volume_property(vol_name,"BField",BField);
  _msetup.message("+++B Field added to MIND:",BField,bhep::VERBOSE);
  
  _gsetup.set_volume_property(vol_name,"X0",X0);
  _msetup.message("+++X0 added to MIND:",X0,bhep::VERBOSE);
  
 
}

void MINDsetup::readParam(){

    
    // -------------------------------------------------------------//
    //                       |  DIMENSIONS |                        //
    // -------------------------------------------------------------//
    
    bhep::prlevel c = bhep::VERBOSE;
    

    MIND_x = 5000; 
      
    _msetup.message("MIND height:",MIND_x/cm,"cm",c);
    
    MIND_y = 5000; 
      
    _msetup.message("MIND width:",MIND_y/cm,"cm",c);
    
    MIND_z = 10000; 
      
    _msetup.message("MIND length:",MIND_z/cm,"cm",c);
    
  
    //-------------------------------------------------------------//

    //--------------------------- VOLUMES ------------------------//
    
    MOTHER_x = MIND_x + 10 * cm;
    MOTHER_y = MIND_y + 10 * cm;  
    MOTHER_z = MIND_z + 10 * cm;  

    // -------------------------------------------------------------//
    //                       |  MAGNETIC FIELD |                    //
    // -------------------------------------------------------------//
    

    //B_int=_pstore.fetch_dstore("B") * mm;
    //B_int=bhep::double_from_string(_store.fetch("GEOM_GG_CELL_diam"))*tesla;
    
    B_int = 0.2 * tesla;
      
    _msetup.message("Magnetic field intensity:",B_int/tesla,"tesla",c);
    
    // -------------------------------------------------------------//
    //                       |  RADIATION LENGTH |                    //
    // -------------------------------------------------------------//
    
    //X0 = _pstore.fetch_dstore("x0") * mm;
    X0 = 1e9 *mm;

    _msetup.message("Radiation length:",X0/cm,"cm",c);

    // -------------------------------------------------------------//
    //                       |  MEASUREMENTS |                    //
    // -------------------------------------------------------------//

    meas_dim = 2;
    meastype = "xy";
    
    resx = 1.;
    resy = 1.;
    
}
