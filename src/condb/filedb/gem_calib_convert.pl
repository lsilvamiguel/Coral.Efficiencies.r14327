#!/usr/local/bin/perl -w
require 5.004;
use strict;
use diagnostics;
use Time::Local;

#
#Parameters
my $det = "GM";
#my $oldDB = "/afs/cern.ch/compass/scratch/d01/DB";
my $oldDB = "oldDB";
#my $newDB = "/afs/cern.ch/compass/scratch/d01/newDB";
my $newDB = "newDB";
#
#
my @flist = `ls $oldDB | grep $det`;

for(@flist){
  print $_;
  /(\w+)__~~start-(\w+)-(\w+)-(\w+)-(\w+):(\w+):(\w+)~~finish-(\w+)-(\w+)-(\w+)-(\w+):(\w+):(\w+)/;
  my $time = $2.".".$3.".".$4.":".$5.".".$6.".".$7;
  print $time,"\n";
  my $ltime = timelocal($7,$6,$5,$4,$3 - 1, $2 - 1900);
  print $ltime,"\n";

  my @log = `Linux/gem_calib_convert $1__ $ltime $oldDB`;
  open(ODAT,">$newDB/$_");
  
  for(@log){
    chomp($_);
    my @line = split / /,$_;
    print ODAT $line[0]," ",$line[2]," ",$line[3]," ",$line[1],"\n";
  }
  #last;
   close ODAT;
}

