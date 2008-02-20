/* -*- mode: c++ -*- */
#ifndef _root_dst___
#define _root_dst___


//! root2dst Class
/*!
  Converts ROOT tree into bhep DST file
*/

class root2dst{

 public:
  
  root2dst(){};

  ~root2dst(){};
  
  bool initialize();
  bool execute();
  bool finalize();

protected:

  //verbosity level
  bhep::prlevel level;
  
  //messenger
  bhep::messenger m;
  
  //event counter
  size_t nevt;

};





#endif
