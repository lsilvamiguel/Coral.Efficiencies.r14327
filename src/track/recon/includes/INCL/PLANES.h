*   plane description.......first detector.(in front of M1).

      common/planes/ ntpl,xplane(mxpl),cosn(mxpl),sinn(mxpl),
     1 jyzt(mxpl),sgpln(mxpl),fd(mxpl),ntpl_tmp,ntpl_cor 

      common/plantyp/nypl,nzpl,nthpl,ityp1(mxpl),
     1    jypl(mxpl),jzpl(mxpl),jthpl(mxpl)

      common/logic_subset/sub,mp1,mp2

      common/planes_set/ ntpl_s1,nypl_s1,nzpl_s1,nthpl_s1,
     1                   ntpl_s2,nypl_s2,nzpl_s2,nthpl_s2,
     2                   nypl_tmp,nzpl_tmp,nthpl_tmp
 
      common/plans2_set/ ntpl2_s1,nypl2_s1,nzpl2_s1,nthpl2_s1,
     1                   ntpl2_s2,nypl2_s2,nzpl2_s2,nthpl2_s2,
     2                   nypl2_tmp,nzpl2_tmp,nthpl2_tmp

*   plane description.......second detector.(behind M1).

      common/plans2/ ntpl2,xplan2(mxpl),cosn2(mxpl),sinn2(mxpl),
     1 jyzt2(mxpl),sgpln2(mxpl),fd2(mxpl),ntpl2_tmp,ntpl2_cor

      common/plan2typ/nypl2,nzpl2,nthpl2,ityp2(mxpl),
     1    jypl2(mxpl),jzpl2(mxpl),jthpl2(mxpl),modul(mxpl)


*   plane description ...... third detector (behid Rich)
      common/plan_sup/ntpl_sup,xplan_sup(mxpl),cosn_sup(mxpl),
     1  sinn_sup(mxpl),jyzt3(mxpl),sgpln_sup(mxpl),npla33,
     2  size_y,size_z,size_yg,size_zg,size_ysf,size_zsf 

      common/plan_sup_typ/ nypl3,nzpl3,nthpl3,ityp3(mxpl),
     1    jypl3(mxpl),jzpl3(mxpl),jthpl3(mxpl)

*-----------------------
      common/efficiency/eff_tel1(mxpl),eff_tel2(mxpl),
     1                  eff_tel3(mxpl)
      common/detectors/i_s1(mxpl),i_s2(mxpl),
     1           iuid(mxpl),iuid2(mxpl),iuid3(mxpl)

      
      common/logic_planes_off/n_black1,n_black2,pat_hit
  

*    multiple coulomb scattering parameters

     

      common/mcs/tetdet(2*mxpl),xxx(2*mxpl),sss(2*mxpl),ccc(2*mxpl),
     1 vmat(2*mxpl,2*mxpl),fudge,mib(2*mxpl),mib_lab(2*mxpl),
     2 plow_mcs,vrscal(5)
      
   








