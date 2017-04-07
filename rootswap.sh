# To be . running (b)ash

RootVersion=3.03.09
sysversion=`cat /etc/redhat-release | awk '{print $5}'`
gccversion=`gcc --version`

if [ "$sysversion" = "6.1" -a "$gccversion" = "2.95.2" ]; then
  syspath=rh61_gcc2952
elif [ "$sysversion" = "7.2" -a "$gccversion" = "2.95.2" ]; then
  syspath=rh72_gcc2953
elif [ "$sysversion" = "7.2" -a "$gccversion" = "2.95.3" ]; then
  syspath=rh72_gcc2953
elif [ "$sysversion" = "7.3" -a "$gccversion" = "2.95.2" ]; then
  syspath=rh72_gcc2953
elif [ "$sysversion" = "7.3" -a "$gccversion" = "2.95.3" ]; then
  syspath=rh72_gcc2953
else
  echo "RedHat $sysversion, gcc $gccversion not supported."
  exit 1
fi

newbinpath=`printenv PATH | sed "s/compass\/delivery\/tools\/root\/bin/sw\/root\/v$RootVersion\/$syspath\/root\/bin/g"`

newlibpath=`printenv LD_LIBRARY_PATH | sed "s/compass\/delivery\/tools\/root\/lib/sw\/root\/v$RootVersion\/$syspath\/root\/lib/g"`

export PATH=${newbinpath}
export LD_LIBRARY_PATH=${newlibpath}
export ROOTSYS=/afs/cern.ch/sw/root/v$RootVersion/$syspath/root

unset RootVersion
unset newbinpath
unset newlibpath