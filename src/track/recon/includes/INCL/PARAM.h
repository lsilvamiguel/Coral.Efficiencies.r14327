C---- all parameters - arrays dimensions etc.
C----  warning! - always as a first include
C                                             maximal numbers of:
 
            parameter (mxpl=300)         ! nb of planes in detectors
            parameter (mxplan=100)       ! nb of planes in detectors for rec. tracks check
            parameter (mpbbin=25)        ! nb of bins for MCS effect        
            parameter (mxlin=200)        ! nb of lines (pattern recognition)
            parameter (mxlin1=20)        ! effective nb of line found is mxlin-mxlin1 
            parameter (nmhit=500)        ! nb of clusters in planes 
            parameter (nmghit=5000)      ! nb of hits in all planes in zone
            parameter (ntrack=300)       ! nb of tracks 
            parameter (ncth1=100)        ! multiplicity in plane (for segments finding)
            parameter (maxiclus=10)      ! max nb of tracks in one cluster
            parameter (montecarlo=9000)  ! emergency parameter - must be the same as in C++
        
            

           
    
          
