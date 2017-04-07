      subroutine reconend

**************************************************
*                                                *
*       recon_2001/2002 closing module           *
*                                                *
**************************************************

#include "INI_STAT.h"

       logical fringe,extrapol,druk 

        
       if(mmc.gt.1) then
*------- ntuple staff - end
          call hrout(0,icycle,' ')
          call hrend('ntuple')
       endif

********************************************
*                                          *
*          OUTPUT & STATISTICS             *
*                                          *
********************************************

      print *
      print *
      print *,'=== STATISTICS of RECON_2001/2002, version 10.0  ==='
      print *
      print *

      if(mmc.eq.6)then
       print*,' no reconstruction - hit profiles only ==> ended'
       return
      endif
    
      print *,' END.DET.1..Nb.of events: ',nevnta
      print *,' END.DET.2..Nb.of events: ',nevntb

        
*     print *, 'total number of input events:', max_ev 

      print *, 'number of "bad" events (reconstruction suspicious): ',
     &           l_bad

      if(druk)then
      If(l_bad.le.100)then
      print *,'list of numbers of "bad" events:',(ib(i),i=1,l_bad) 
      else
      print *,'list of numbers of first 100 "bad" events:',
     &         (ib(i),i=1,100)
      endif
      endif

      print *, 'events processed:', iv_count

       if(mmc.eq.0.or.mmc.eq.2.or.mmc.eq.4)go to 100

      print *, ' acceptance: n_tel1,n_tel2: ',nn_tel1,nn_tel2
      print *,'total number of generated tracks',ntr_licz
      print *,'total number of accepted generated tracks:',ntr_ac_licz
      print *,'total number of accepted generated tracks in new proc.:',
     &    ntr_ac_licz_new
      print *,'number of non-accepted generated tracks in new proc.:',
     &    nos(66) 

 100   continue
     
      print *, 'total number of associated tracks:', maxtr
      print *, 'total number of associated tracks with extrapolation:', 
     & maxt_extrap
      print *, 'total number of associated tracks in fringe field:',
     & maxt_fringe
      print *, 'total number of associated tracks in ff with bridging:',
     & maxt_fringe1
*      print *,'total number of assoc. tracks in f.f. before cleaning:',
*     & mlassc1_fringe


      print *
      print *, ' number of tracks/per event:   ',
     &           float(maxtr)/float(iv_count) 
      print *, ' number of tracks/per event for different triggers:   '
      print *, '   IT        ', float(maxtrIT)/float(iv_countIT)
      print *, '   MT        ', float(maxtrMT)/float(iv_countMT)
      print *, '   LT        ', float(maxtrLT)/float(iv_countLT)
      print *, '   IT+MT     ', float(maxtrITMT)/float(iv_countITMT)
      print *, '   LT+MT     ', float(maxtrLTMT)/float(iv_countLTMT)
      print *, '   LT+IT     ', float(maxtrLTIT)/float(iv_countLTIT)
      print *, '   IT+LT+MT  ', float(maxtrLTITMT)/float(iv_countLTITMT)
      print *, '   OT        ', float(maxtrOT)/float(iv_countOT)
      print *, '   OT+IT     ', float(maxtrOTIT)/float(iv_countOTIT)
      print *, '   OT+MT     ', float(maxtrOTMT)/float(iv_countOTMT)
      print *, '   incl.MT   ', float(maxtrINCMT)/float(iv_countINCMT)
      print *, '   others:   ', float(maxtrAT)/float(iv_countAT)

      print *
      print *, ' Nb.of golden events: ',m_track(7), ' with ',m_track(6),
     &            ' tracks/event'
         if(m_track(7).gt.100)then
      print *, ' First 100 Golden events:  '
      print *,   (lgolden(i),i=1,100)
         else
      print *, ' Golden events: ', (lgolden(i),i=1,m_track(7))
         endif
      print *
      print *, 'nb. of events with 1 track/event: ',m_track(1)
      print *, 'nb. of events with 2 tracks/event: ',m_track(2)
      print *, 'nb. of events with 3 tracks/event: ',m_track(3)
      print *, 'nb. of events with 4 tracks/event: ',m_track(4)
      print *, 'nb. of events with more than 4 tracks/event: ',
     &             m_track(5)



        if(mmc.eq.0.or.mmc.eq.2.or.mmc.eq.4)then
          print *
          print *,' no MC comparison - real data  '
          print *
        elseif(mmc.eq.1.or.mmc.eq.3.or.mmc.eq.5)then
          print *
          print *,' now Monte Carlo comparison'
          print *


      print *,'total number of glued (missing in ntuple) tracks: '
     &     ,nos(99),'  before cleaning : ',nos(73)
    
          do i=1,51
          if(nos(i).ne.0)then
      print *,'total number of ',i,' glued tracks: ',nos(i)
          endif
          enddo

      print *
      print *, 'tracks accepted due to used criteria:', itr_count
      print *, 'tracks accepted due to used criteria in ff :', nos(87)
      print *, 'tracks accepted in ff but continued behind SM1 :', 
     &  nos(88)

      print *, 'tracks accepted and extrapolated to SM2:', nos(83)
      print *, 'reconstructed pilup tracks: ', nos(81)
      print *
        endif

 
       do i=1,9

       if(i.le.inpar2) 
     & print *,' nb of tracks: proj: ',i ,maxt_proj(i),
     &      ' before cleaning: ',mlassc1_proj(i)

       enddo




       print *
       print *
         if(mmc.eq.0.or.mmc.eq.2.or.mmc.eq.4)go to 200
       print *, '    "ghosts":'
       print *
       
      
       print*, 'tracks bad reconstructed (due to hits criteria):'
     &        ,nos(100),' removed by cleaning: ',nos(72)
       print*, '         ghosts with p < 0.5   ',nos(96)
       print*, '         ghosts with p < 1.0   ',nos(95)
       print*, '         ghosts with p < 1.5   ',nos(94)
       print*, '         ghosts with p < 2.0   ',nos(93)
       print*, '         ghosts with p < 2.5   ',nos(92)
       print*
       print*, '         ghosts with p >100.0   ',nos(91)
       print*, '         ghosts with p>80.0 ', nos(90)
       print*
       print*,'          ghosts from fringe field rec.:',nos(89)
       print*,'          ghosts with extrapolation to SM2:',nos(82)
       print*



       print*, 'tracks bad rec.- improper nb of hits in tel.1:'
     &        ,nos(98)
       print*, 'tracks bad rec.- improper nb of hits in tel.2:'
     &        ,nos(97)

       print*, 'tracks with nb of hits in tel.1 > nb of planes:'
     &        ,nos(86)
       print*, 'tracks with nb of hits in tel.2 > nb of planes:'
     &        ,nos(85)
       print*, 'tracks with nb of supp. hits  > nb of planes:'
     &        ,nos(84)
       print*,'total nb of clusters with more than allowed tracks '
     &        ,' associated with before SM1: ',nos(67)
       print*,'total nb of clusters with more than allowed tracks '
     &        ,' associated with behind SM1: ',nos(68)
       print*,'total nb of clusters with more than allowed tracks '
     &        ,' associated with behind Rich: ',nos(69)

        g1 = float(nos(100))/float(itr_count)*100.
        g2 = float(nos(100))/float(maxtr)*100.
        eb = float(maxtr)/float(ntr_ac_licz)*100.
        er = float(itr_count)/float(ntr_ac_licz)*100.
        er_new = float(itr_count)/float(ntr_ac_licz_new)*100.

        eb1 = float(maxtra(1))/float(ntra_ac(1))*100.
        er1 = float(nos(80))/float(ntra_ac(1))*100.
        er1_new = float(nos(65))/float(ntra_ac(8))*100.

        eb2 = float(maxtra(2))/float(ntra_ac(2))*100.
        er2 = float(nos(79))/float(ntra_ac(2))*100.
        er2_new = float(nos(64))/float(ntra_ac(9))*100.


        eb3 = float(maxtra(3))/float(ntra_ac(3))*100.
        er3 = float(nos(78))/float(ntra_ac(3))*100.
        er3_new = float(nos(63))/float(ntra_ac(10))*100.

* -- now pilup_out calculations

       
        er4 = float(nos(77))/float(ntra_ac(4))*100.       
        er5 = float(nos(76))/float(ntra_ac(5))*100.       
        er6 = float(nos(75))/float(ntra_ac(6))*100.       
        er7 = float(nos(74))/float(ntra_ac(7))*100.    

        er4_new = float(nos(62))/float(ntra_ac(11))*100.       
        er5_new = float(nos(61))/float(ntra_ac(12))*100.       
        er6_new = float(nos(60))/float(ntra_ac(13))*100.       
        er7_new = float(nos(59))/float(ntra_ac(14))*100.    


      print *
      print *,'              summary: '
      print *
      print'("               ghosts in %:          ",2f7.2)', g1,g2
      print *
      print'("  range (in GeV)      efficiency (average in %)     ",
     &     " bridging (average in %)")' 
      print'("----------------------------------------------------",
     &     "-------------------------")'
      print'("   full range               ",f6.2,"          ",f6.2,
     &       "                    ",f6.2)',er,er_new,eb
      print'("      p > 2                 ",f6.2,"          ",f6.2,
     &       "                    ",f6.2)',er1,er1_new,eb1
      print'("   2 < p < 50               ",f6.2,"          ",f6.2,
     &       "                    ",f6.2)',er2,er2_new,eb2
      print'("   2 < p < 5                ",f6.2,"          ",f6.2,
     &       "                    ",f6.2)',er3,er3_new,eb3
      print'("----------------------------------------------------",
     &     "-------------------------")'
      print *,'              calculations without pileup           '
      print *
      print'("   full range         ",f6.2,"           ",f6.2)',
     &        er4,er4_new
      print'("      p > 2           ",f6.2,"           ",f6.2)',
     &        er5,er5_new
      print'("   2 < p < 50         ",f6.2,"           ",f6.2)',
     &        er6,er6_new
      print'("   2 < p < 5          ",f6.2,"           ",f6.2)',
     &        er7,er7_new
      print'("----------------------------------------------------",
     &     "-------------------------")'

 200   continue

      print *
      print *
      print *,'== END of OUTPUT from RECON_2001/2002, version 10.0 =='
      print *
      print * 



*----------------------------------------------------- end of basic output


            return
            end
         





