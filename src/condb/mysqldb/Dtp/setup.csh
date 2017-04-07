setenv QTDIR /afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/Qt/3.0.1
setenv PATH ${QTDIR}/bin:$PATH

if(! $?LD_LIBRARY_PATH) then 
    setenv LD_LIBRARY_PATH
endif
setenv LD_LIBRARY_PATH ${QTDIR}/lib:${LD_LIBRARY_PATH}
