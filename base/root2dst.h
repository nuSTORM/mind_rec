/* -*- mode: c++ -*- */
#ifndef _root_dst___
#define _root_dst___

#include <bhep/event.h>
#include <bhep/messenger.h>

//! root2dst Class
/*!
  Converts ROOT tree into bhep DST file
*/

class root2dst{

 public:
  
  root2dst(bhep::prlevel);

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
