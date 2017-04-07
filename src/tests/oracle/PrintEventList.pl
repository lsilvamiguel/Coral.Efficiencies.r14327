#!/usr/local/bin/perl -w

require 5.004;

use strict;
use diagnostics;

if($#ARGV != 0) {
    die("Usage: PrintEventList.pl <directory_name>\n");
}

my $list_file_name = "J_Psi.list";
my $exec_file = "./Linux/PrintEventList";
my $directory = $ARGV[0];

my @file_list = `ls -1 $directory | grep raw | grep Run_`;

foreach my $file_name (@file_list) {
    chomp($file_name);
    print("$file_name\n");
    unless(-z "$directory/$file_name") {
	`$exec_file $directory/$file_name $list_file_name`;
    }
}
