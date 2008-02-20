
#include <root2dst.h>


//*************************************************************
root2dst::root2dst(bhep::prlevel vlevel){
//*************************************************************
  
  level = vlevel;
  
  m = bhep::messenger(level);
  
  m.message("+++ root2dst Constructor +++",bhep::VERBOSE); 
  
}


//*************************************************************
bool root2dst::initialize(const bhep::sstore& run_store) {
//*************************************************************
    
    m.message("+++ root2dst init  function ++++",bhep::NORMAL);
        
    nevt=0;

    return true;
}


//*************************************************************
bool root2dst::execute(bhep::event& event){
//*************************************************************
    
  /*
    Take a ROOT tree event and convert it into bhep DST event
   */
  
  m.message("+++ root2dst execute function ++++",bhep::VERBOSE);

  nevt++;
 
  return true;  

}



//*************************************************************
bool root2dst::finalize() {
//*************************************************************
  
  m.message("+++ root2dst finalize function ++++",bhep::NORMAL);
    
  m.message("++ Number of analyzed events: ",nevt,bhep::NORMAL);

  return true;
}
