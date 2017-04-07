#!/usr/bin/perl
#
# check which local files have changed without having to run
# 'cvs up'

open(CVS, "cvs status 2>&1 |") or die("cannot start 'cvs status'\n");
while(<CVS>) {
    if(/^cvs status: Examining (.*)/) {
	$lastdir=$1;
    } elsif(/^cvs/) {
	printf("In %s:\n", $lastdir);
	print;
    } elsif(/^File:\s+(no file )?((\S+|\s+\S+)+?)\s+Status:\s+(?!Up-to-date)(.+)/) {
	printf("%-20s %s/%s\n", $4, $lastdir, $2);
    }
}
