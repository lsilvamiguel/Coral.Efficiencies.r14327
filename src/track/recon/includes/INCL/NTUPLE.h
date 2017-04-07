*--- ntuple stuff
*
*     mc parameter defines kind of output
*     mc=2 & 4 - MC info is not used at all, run for data
*     mc=3 & 5 - makes correspodence reconstructed - generated track


       parameter (ncolum_mc1=37)        ! tracks ntuple (mc=3)
       parameter (ncolum_mc2=5)         ! events ntuple (mc=3)
       parameter (ncolum_mc3=240)       ! tracks ntuple (mc=5)
       parameter (ncolum_mc4=6)         ! events ntuple (mc=5)
  
       parameter (ncolum1=22)           ! data tracks ntuple (mc=2)
       parameter (ncolum2=3)            ! data events ntuple (mc=2) 
       parameter (ncolum3=222)          ! data tracks ntuple (mc=4)
       parameter (ncolum4=95)           ! data events ntuple (mc=4)

       parameter (nhbmem=150000)

       common /pawc/ hmem(nhbmem)
