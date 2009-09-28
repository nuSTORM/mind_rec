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
// void MINDsetup::init(bhep::gstore pstore,bhep::sstore gstore,
// 		       bhep::prlevel level) {
void MINDsetup::init(bhep::gstore pstore, bhep::prlevel level) {
//*************************************************************
    
  _msetup=bhep::messenger(level);
  
  _msetup.message("++MINDsetup Messenger generated++",bhep::VERBOSE);
  
  _pstore=pstore; //_store=gstore;
  
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

    //_msetup.message("MIND Setup:", _gsetup,bhep::VERBOSE);
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

  EVector pos(3,0); pos[2]=0;
  
  Volume* mother = new Box(pos,xaxis,yaxis,
			   MOTHER_x/2,MOTHER_y/2,MOTHER_z/2);
    
  _msetup.message("Mother volume generated",bhep::VERBOSE);

  // add mother volume
  
  _gsetup.set_mother(mother);
   
  _msetup.message("Mother added to setup",bhep::VERBOSE);

  // Create detector volume
  
  pos[2]=0;

  const dict::Key vol_name = "Detector";
  
  Volume* det = new Box(pos,xaxis,yaxis,MIND_x/2,MIND_y/2,MIND_z/2);
    
  _msetup.message("MIND volume generated",bhep::VERBOSE);

  // add volume

  _gsetup.add_volume("mother",vol_name,det);
   
  //Introduce IRON scintillator sandwiches.
  // int nplanes = (int)( MIND_z / (IRON_z + nScint * SCINT_z) );

//   for (int iplane = 0;iplane < nplanes;iplane++) {

//     add_slab(iplane, "mother");

//   }

  std::cout << _gsetup << std::endl;

}

//*************************************************************
void MINDsetup::add_slab(int plane, const dict::Key det_vol){
//*************************************************************

  EVector plane_pos(3,0);
  //Names for particular sections
  const dict::Key iron_name = "IRON_plane"+bhep::to_string(plane);
  
  //Define positions
  double slab_width = IRON_z + nScint * SCINT_z;
  double mind_front = -MIND_z/2;

  //Scintillator.
  plane_pos[2] = mind_front + slab_width * plane + IRON_z + SCINT_z/2;
  const dict::Key scint_name = "SCINT_plane"+bhep::to_string(plane_pos[2]);

  Volume* Sc_slab = new Box(plane_pos,xaxis,yaxis,MIND_x/2,MIND_y/2,SCINT_z/2);

  _gsetup.add_volume(det_vol,scint_name,Sc_slab);

  _gsetup.set_volume_property(scint_name,"X0",X0Sc);
  
  //IRON
  plane_pos[2] = mind_front + slab_width * plane + IRON_z/2;

  Volume* Fe_slab = new Box(plane_pos,xaxis,yaxis,MIND_x/2,MIND_y/2,IRON_z/2);

  _gsetup.add_volume(det_vol,iron_name,Fe_slab);

  _gsetup.set_volume_property(iron_name,"X0",X0Fe);
  
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

  _zaxis = EVector(3,0);
  _zaxis[2]=1;

  // _gsetup.set_volume_property("mother","BField",BField);
  _msetup.message("+++B Field added to MOTHER:",BField,bhep::VERBOSE);

  _gsetup.set_volume_property_to_sons("mother","BField",BField);
  _gsetup.set_volume_property_to_sons("mother","de_dx",de_dx);
  _gsetup.set_volume_property_to_sons("mother",RP::SurfNormal,_zaxis);
    
  const dict::Key vol_name = "Detector";
//   _gsetup.set_volume_property(vol_name,"BField",BField);
  _msetup.message("+++B Field added to MIND:",BField,bhep::VERBOSE);
  
  _gsetup.set_volume_property(vol_name,"X0",X0Eff);
  _msetup.message("+++X0 added to MIND:",X0,bhep::VERBOSE);

  // _gsetup.set_volume_property(vol_name,"de_dx",de_dx);
  _msetup.message("+++de/dx added to MIND:",de_dx,bhep::VERBOSE);
  
 
}

void MINDsetup::readParam(){

    
    // -------------------------------------------------------------//
    //                       |  DIMENSIONS |                        //
    // -------------------------------------------------------------//
    
    bhep::prlevel c = bhep::VERBOSE;
    

    MIND_x = _pstore.fetch_dstore("MIND_x") * m; 
      
    _msetup.message("MIND height:",MIND_x/cm,"cm",c);
    
    MIND_y =  _pstore.fetch_dstore("MIND_y") * m;
      
    _msetup.message("MIND width:",MIND_y/cm,"cm",c);
    
    MIND_z =  _pstore.fetch_dstore("MIND_z") * m;
      
    _msetup.message("MIND length:",MIND_z/cm,"cm",c);
    
    //-------------------------------------------------------------//
    //                      | INNER DIMENSIONS |                   //
    //-------------------------------------------------------------//

    IRON_z = _pstore.fetch_dstore("widthI") * cm;
    SCINT_z = _pstore.fetch_dstore("widthS") * cm;
    nScint = _pstore.fetch_istore("nplane");
    rel_den = _pstore.fetch_dstore("rel_den");

    //--------------------------- VOLUMES ------------------------//
    
    MOTHER_x = MIND_x + 10 * cm;
    MOTHER_y = MIND_y + 10 * cm;  
    MOTHER_z = MIND_z + 10 * cm;  

    // -------------------------------------------------------------//
    //                       |  MAGNETIC FIELD |                    //
    // -------------------------------------------------------------//
    

    //B_int=_pstore.fetch_dstore("B") * mm;
    //B_int=bhep::double_from_string(_store.fetch("GEOM_GG_CELL_diam"))*tesla;
    
    B_int = 1.0 * tesla;
      
    _msetup.message("Magnetic field intensity:",B_int/tesla,"tesla",c);
    
    // -------------------------------------------------------------//
    //            |  RADIATION LENGTH AND ENERGY LOSS |             //
    // -------------------------------------------------------------//
    
    X0Fe = _pstore.fetch_dstore("x0Fe") * mm;
    X0Sc = _pstore.fetch_dstore("x0Sc") * mm;

    _wFe = IRON_z/(IRON_z + SCINT_z*nScint*rel_den);
    double wSc = 1-_wFe;
    X0Eff = 1./(_wFe/X0Fe + wSc/X0Sc);

    de_dx = _pstore.fetch_dstore("de_dx") * MeV/cm;
    //X0 = 1e9 *mm;

    _msetup.message("Radiation length:",X0Fe/cm,"cm",c);

    // -------------------------------------------------------------//
    //                       |  MEASUREMENTS |                    //
    // -------------------------------------------------------------//

    meas_dim = 2;
    meastype = "xy";
    
    resx = _pstore.fetch_dstore("pos_res") * cm;
    resy = _pstore.fetch_dstore("pos_res") * cm;
    
}
