      subroutine mag_field(pos,h_vec)

***********************************************************************
*                                                                     *
* ----- computes magnetic field for RECON ---------                   *
*           ---input:  [m] pos(1)....Jura                             *
*                          pos(2)....up                               *
*                          pos(3)....beam                             *
*        (double precision)                                           *
*                                                                     *
*        uses  GUFLD procedure which uses TAlgo::Field function       *
*                                            from Traffic (S.G.)      *
*                                                                     *
*           ---output:  [T] h_vec(1)....Jura                          *
*                           h_vec(2)....up                            *
*                           h_vec(3)....beam                          *
*                                                                     *
*   WARNING: double precision ABSOLUTELY needed for this routine      *
*                as well as for routine where it is called!!!         *
*                                                                     *
***********************************************************************

      double precision pos(3),h_vec(3),hvec_temp(3),pos_temp(3)


       pos_temp(1)=pos(3)*100.d0   ! meters -> cm, recalculation cm->mm is performed in TAlgo.h
       pos_temp(2)=pos(1)*100.d0
       pos_temp(3)=pos(2)*100.d0  
 

      CALL GUFLD(pos_temp,hvec_temp)

*      print*,pos(1),pos(2),pos(3),' cm ',
*     1        (hvec_temp(i),i=1,3),' T'

      h_vec(1)=hvec_temp(2)/10.d0
      h_vec(2)=hvec_temp(3)/10.d0
      h_vec(3)=hvec_temp(1)/10.d0   

   
      end
