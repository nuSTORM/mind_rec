
#include <mind/root2dst.h>

int main(int argc, char* argv[]){
    
    int nevents = 10;

    if (argc==2) nevents = atoi(argv[1]);
  
    bhep::prlevel c = bhep::VERBOSE;
    
    root2dst* cvt = new root2dst(c);
    
    cvt->initialize();

    for(int i=1; i < nevents+1; i++) {
          
	if (i%100==0) cout<< "Number of events read "<<i<<endl;
        
	cvt->execute();
	
   }
   
   cvt->finalize();
   
   return 0;

}
