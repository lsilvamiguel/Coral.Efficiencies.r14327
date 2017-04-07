#!/usr/local/bin/perl -w

require 5.004;

use strict;
use diagnostics;

if($#ARGV != 0) {
    die "USAGE: ./Backup_logs.pl <run_number>\n";
}

unless(defined $ENV{CUR_DIR}) {
    $ENV{CUR_DIR} = $ENV{PWD};
}

my $run_number = $ARGV[0];
`mv $ENV{CUR_DIR}/logFiles/optdst*-$run_number.log $ENV{CUR_DIR}/Run_$run_number/`;
my $ret_value = system("cd $ENV{CUR_DIR}; tar cf - ./Run_$run_number/* | gzip -c > $ENV{CUR_DIR}/Run_$run_number.backup.tgz");

if($ret_value == 0) {
    `rm $ENV{CUR_DIR}/Run_$run_number/*`;
    `rm $ENV{CUR_DIR}/Run_$run_number/.info`;
    `mv $ENV{CUR_DIR}/Run_$run_number.backup.tgz $ENV{CUR_DIR}/Run_$run_number/`;
    print "Backup of log files has been done successfully.\n";
}


