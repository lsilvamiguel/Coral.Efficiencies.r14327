# To be sourced running (t)csh

set RootVersion = 3.03.09
set sysversion = `cat /etc/redhat-release | awk '{print $5}'`
set gccversion = `gcc --version`

if( $sysversion == "6.1" && $gccversion == "2.95.2" ) then
  set syspath = rh61_gcc2952
else if( $sysversion == "7.2" && $gccversion == "2.95.2" ) then
  set syspath = rh72_gcc2953
else if( $sysversion == "7.2" && $gccversion == "2.95.3" ) then
  set syspath = rh72_gcc2953
else if( $sysversion == "7.3" && $gccversion == "2.95.2" ) then
  set syspath = rh72_gcc2953
else if( $sysversion == "7.3" && $gccversion == "2.95.3" ) then
  set syspath = rh72_gcc2953
else 
  echo "RedHat $sysversion, gcc $gccversion not supported."
  exit 1
endif

set newbinpath = `printenv PATH | sed "s/compass\/delivery\/tools\/root\/bin/sw\/root\/v$RootVersion\/$syspath\/root\/bin/g"`

set newlibpath = `printenv LD_LIBRARY_PATH | sed "s/compass\/delivery\/tools\/root\/lib/sw\/root\/v$RootVersion\/$syspath\/root\/lib/g"`

setenv PATH ${newbinpath}
setenv LD_LIBRARY_PATH ${newlibpath}
setenv ROOTSYS /afs/cern.ch/sw/root/v$RootVersion/$syspath/root

unset RootVersion
unset newbinpath
unset newlibpath
