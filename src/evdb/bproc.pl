#!/usr/local/bin/perl -w
# Author: Vinicio duic duic@cern.ch

require 5.004;

use strict;
#use diagnostics;

#*****************************************************************#
#
#  HARD CODINGS
#
my $pgmName = "batchread";
my $BOOTFile = "noperiod.2001.BOOT";
my $pgmParallelize = "Linux/parallelCheck";
my $pgmExecute = "Linux/dstreadTest";
my $optSuff = ".opt";
my $outSuff = ".out";
my $outDir = "out";
my $optDir = "opt";
my $lsfDir = "lsf";
my $bDumpExe = ""; # "export CSOBJSTORE_DUMP=1;"
my $syncToPgm = "synchroTo.pl";
my $LSFsleep = 10;
#
#*################################################################
#
#  OVERRIDABLES (not implemented yet 
#

my $pgmArguments = "";
my $showLimits = "ulimit -a; ";
my $q_1 = "-q vy_1nh";
#my $q_1 = "-q vy_8nh";
my $q_2 = "-R compass";
my $queue;

#
##################################################################

my $sumRootHistograms = 1;
my $sync = "";
my $nkd0 = getNakedName($0);
my $synopsis = "$nkd0: synopsis:\n $nkd0 {program_name coral_options_file [/ [--sync syncTo] -queue queueName -Resource resourceName] | ".
               "--summarize absWrkDir collDir optDir lsfDir outDir CardFile}\n\t".
               "program_name: is the file that processes data\n\n\t".
               "if syncTo = 1, automatic synchronization well be set.\n\t".
               "Set the environment variable PARALLEL_BATCH_DUMP to dump\n".
               "Every --option can be shortened to -o, where \"o\" is the ".
               "corresponding initial letter.\n";

unless (defined $ARGV[1]) {die "$synopsis\n";}

if ($ARGV[0] eq "--summarize" || $ARGV[0] eq "-s") {
  shift @ARGV;
  summarize(@ARGV);
} else {
  $pgmExecute = shift @ARGV;
  die "$nkd0: $pgmExecute does not exist.\n" unless (-e $pgmExecute);
  my @pgmOpts = ();
  foreach (@ARGV) {
    if ($_ eq "\/") {
      shift @ARGV;
      last;
    }
    push @pgmOpts, shift @ARGV;
  }

  setOptionsValues(@ARGV);

  $queue = "$q_1 $q_2";
  launchBatches(@pgmOpts);
}

#***************************************************************************#
sub launchBatches {
  my $thisDir = $ENV{PWD} || 
    die "$nkd0: cannot determine current directory.\n";

  my $bDump = mustDump();
  my $runN;
  my @nakedCardsFiles=();

  my $CORALDir = $ENV{CORAL} || undef;
  die "$nkd0: cannot find CORAL environment variable.\n" 
    unless (defined $CORALDir);

  my @now = localtime(time());

  my $optionsFile = shift @_;
  my $pgmArguments = join(" ", @_) if defined @_;

  my $rootHistFileName = getRootHistFrom($optionsFile);
  open (OPTFILE, $optionsFile) || die "$nkd0: cannot open $optionsFile. $!\n";
  while(<OPTFILE>) {
    if (/\s?Data\s+run/) {
      if (/select\s+(\d+)/) {
	$runN = int($1) || undef;
	last;
      }
    }
  }
  close OPTFILE;
  defined $runN or die "Cannot determine run number\n";

  my $rootHFileName = getRootHistFrom($optionsFile);
  my $collDir = sprintf("prll_%5.5d-%2.2d%2.2d%2.2d%2.2d%2.2d", 
			$runN, $now[5]-100, $now[4]+1, $now[3], 
			$now[2], $now[1]);

  my $absCollDir = "$thisDir/$collDir";
  die "$nkd0: $absCollDir already exists; no overriding is allowed.\n" 
    if -e "$absCollDir";

  mkdir ("$absCollDir", 0600) || 
    die "$nkd0: cannot create directory $absCollDir. $!\n";
  mkdir ("$absCollDir/$optDir", 0600) ||
    die "$nkd0: cannot create directory $absCollDir/$optDir. $!\n";
  mkdir ("$absCollDir/$outDir", 0600) ||
    die "$nkd0: cannot create directory $absCollDir/$outDir. $!\n";
  mkdir ("$absCollDir/$lsfDir", 0600) ||
    die "$nkd0: cannot create directory $absCollDir/$lsfDir. $!\n";

  `echo \"\" > $absCollDir/$lsfDir/chk-ended 2>&1` &&
    die "$nkd0: cannot write to log file $absCollDir/$lsfDir/chk-ended. $!\n";

  my @cmd = `$pgmParallelize $optionsFile $BOOTFile $absCollDir/$optDir 2>&1`;
  print "@cmd" if defined $bDump;
  foreach (@cmd) {
    chomp;
    if (/Generated \"$absCollDir\/$optDir\/parallel-(cdr\d+-)(\d+)$optSuff/) {
      #$runN = int($2) unless (defined $runN);
      push @nakedCardsFiles, "$1$2";
    }
  }

  my @jobIDS = ();

  # check if sync is actually available
  unless (-e "$thisDir/$syncToPgm") {
    print "$nkd0: WARNING! $thisDir/$syncToPgm was not found. ",
          "Synchronization will not be performed\n";
    $sync = "";
  } else {
    print "$nkd0: jobs will be synchronized by $sync\n" if $sync;
  }

  if ($sync eq "1") {
    my $runAtTime = time()+ $LSFsleep * @nakedCardsFiles+100;
    my @now = localtime($runAtTime);
    $sync = sprintf("%2.2d:%2.2d:%2.2d:%2.2d", $now[4]+1, $now[3], 
		    $now[2], $now[1]); 
    print "$nkd0: synchronization is automatically set to $sync\n";
  }

  foreach (@nakedCardsFiles) {   
    my $jobName = $_;
    $jobName =~ s/\-\d+$//;
    $jobName = "$collDir-" . $jobName;
    my $cpHistos = "";
    $cpHistos = "; sync; cp $rootHFileName $absCollDir/$outDir/$_-$rootHFileName"
      if defined $rootHFileName;

    my @bsubArgs = (
		    "<<EOF",
		    "#BSUB $queue",
		    "#BSUB -r",
		    "#BSUB -o $absCollDir/$outDir/$_$outSuff",
		    "#BSUB -J $jobName",
                    #"#BSUB -b $sync",

		    "(ulimit -c 0; $showLimits cd $CORALDir; . setup.sh; cd $thisDir; $sync $bDumpExe $thisDir/$pgmExecute $absCollDir/$optDir/parallel-$_$optSuff $pgmArguments $cpHistos)",

		    #"(ulimit -c 0; cd $CORALDir; . setup.sh; cd $thisDir; $bDumpExe $thisDir/$pgmExecute $absCollDir/$optDir/parallel-$_$optSuff)",

		    #"(cd $CORALDir; . setup.sh; cd $thisDir; echo \"run $absCollDir/$optDir/parallel-$_$optSuff\nbt\nquit\n\" | gdb -q $thisDir/$pgmExecute)",

		    "echo \"$_ is done\" >> $absCollDir/$lsfDir/chk-ended",
		    "EOF\n"
		    );

    my $bsubArgs = join("\n", @bsubArgs);
    @cmd = `bsub $bsubArgs 2>&1`;
    sleep($LSFsleep);

    foreach my $job (@cmd) {
      if ($job =~ /Job \<(\d+)\> is submitted to queue /) {
	push @jobIDS, " && ended($1)";
	open (LSFLOG, ">$absCollDir/$lsfDir/$_") || 
	  die "$nkd0: cannot write lsf log for job $1. $!\n";
	print LSFLOG "JOBID: $1\n";
	print LSFLOG "$bsubArgs";
	close LSFLOG;
	last;
      }
    }
    print "@cmd";
  }

  #prepare & run script that gathers information from batched jobs output
  die "No job was launched (non existing run?)\n" if (@jobIDS < 1);
  $jobIDS[0] =~ s/ \&\&\s+//;
  my $jobCondition = join("", @jobIDS);
  my @bsubArgs = (
		  "<<EOF",
		  "#BSUB $queue",
		  "#BSUB -r",
		  "#BSUB -o $absCollDir/$lsfDir/summarize.out",
		  "#BSUB -J summarize",
		  "#BSUB -w \"$jobCondition\"",
		  "$0 --summarize $thisDir $collDir $optDir $lsfDir $outDir parallel-$nakedCardsFiles[0]",
		  "EOF\n"
		  );
  my $bsubArgs = join("\n", @bsubArgs);
  @cmd = `bsub $bsubArgs 2>&1`;
  print "@cmd\n";

}  
#***************************************************************************# 
sub summarize {
  my $thisDir = shift @_ ||
    die "summarize: cannot determine absolute working directory path.\n";
  my $collDir = shift @_ || 
    die "summarize: cannot determine data collecting directory.\n";
  my $absCollDir = "$thisDir/$collDir";
  -e $absCollDir or die "summarize: $absCollDir does not exist.\n";

  my $optDir = shift @_ 
    or die "summarize: cannot determine batch option directory.\n"; 

  -e "$absCollDir/$optDir"
    or die "summarize: $absCollDir/$optDir does not exist.\n";

  my $lsfDir = shift @_ 
    or die "summarize: cannot determine batch lsf job log directory.\n"; 

  -e "$absCollDir/$lsfDir"
    or die "summarize: $absCollDir/$lsfDir does not exist.\n";

  my $outDir = shift @_ 
    or die "summarize: cannot determine batch output directory.\n"; 

  -e "$absCollDir/$outDir"
    or die "summarize: $absCollDir/$outDir does not exist.\n";

  my $thisContainer = $_[0];
  my $unitarityOptFileName = "$absCollDir/$optDir/$_[0]$optSuff";
  my @allContainers = summarize_getContainers($unitarityOptFileName);
  defined @allContainers 
    or die "summarize: cannot retrieve unitarity container directive.\n"; 

  sleep(60);                   # wait that all jobs are likely to be finished
  my $cnt = 0;
  my $absOutDir = "$absCollDir/$outDir";
  my %missingJobs;
  foreach my $cntnr (@allContainers) {
    if (-e "$absOutDir/$cntnr$outSuff") {
      #print "summarize: $cntnr found at first chance\n";
      $cnt++;
    } else {
      print "\tsummarize: look for $cntnr bookmark... ";
      if (summarize_hasBookMarked("$absCollDir/$lsfDir/chk-ended", $cntnr)) { 
	print "found.\n";
	$cnt++;            #process has finished; wait...
      } else {
	print "NOT found.\n";
	$missingJobs{$cntnr} = 1; 
      }     
    }
  }

  sleep(60);
  if ($cnt == @allContainers) {
    my $nevents = 0;
    my $nchunks = 0;
    my $totsize = 0;
    foreach (keys %missingJobs) {
      my $maxIter = 0;
      while ($maxIter < 10) {
	if (-e "$absOutDir/$_$outSuff") {
	  $missingJobs{$_} = 0;
	  $maxIter = 10;
	  $cnt++;
	} else {
	  sleep(30);
	  $maxIter++;
	}
      }
    }
    foreach (keys %missingJobs) {
      print "summarize: WARNING! missing job: $_\n" if $missingJobs{$_};
    }

    if ($sumRootHistograms) {
      my $rootHistFileName = getRootHistFrom($unitarityOptFileName);
      $rootHistFileName || 
	die "$nkd0: cannot determine root file from $unitarityOptFileName\n";
      my @histos=();
      for my $key (@allContainers) {
	my $histName = "$key-$rootHistFileName";
	if (-e "$absOutDir/$histName") {
	  my $sz = (-s "$absOutDir/$histName");
	  if ($sz > 1000) {
	    push @histos, $histName;
	  } else {
	    print "summarize: file $histName has size $sz; skipped\n";
	  }
	} else {
	  print "summarize: cannot find $histName; it will not be summed.\n";
	}
      }
      if ($#histos > 0) {
	SumHistograms($thisDir, $absOutDir, $rootHistFileName, @histos);
	if (-e "$absOutDir/$rootHistFileName") {
	  my $cmdmv = `mv $absOutDir/$rootHistFileName $absOutDir/../ 2>&1`;
	}
      } else {
	print "summarize: only ", $#histos+1, " files found.\n";
      }
    } else {
      for my $key (@allContainers) {
	open(IN,"$absOutDir/$key$outSuff") ||
	  die "summarize: ($thisContainer) file ",
	  "$absOutDir/$key$outSuff could not be opened. \n";
	$nchunks++;
	while(<IN>) {
	  chomp;
	  
###########################################################################
#                                                                         #
#                  CUSTOMIZE THIS PART OF THE SCRIPT                      #
#                         TO ANALYZE PROCESSED DATA                       #
#                                                                         #
#   

	  if (/egmentation /) {
	    print "summarize: WARNING! Segmentation fault found in $key\n";
	    next;
	  }
	  /^\s*Total events scanned\s+(\d+)/ && ($nevents+=$1);
	  if (/^\s*Total raw events size:\s+(\d+) read at an average of: (\d+[\.]*)/) {
	    $totsize += $1/1E9;
	  }
	}
	close (IN);
      }
      open(OUT, ">$absCollDir/report$outSuff") || 
	die "summarize: cannot write results.\n";
      
      print OUT "Total number of chunks $nchunks\n",
                "Total events scanned in all chunks $nevents\n",
                "Total size read $totsize GB\n";
      close OUT;
      print "Results summarized in $absCollDir/report$outSuff\n";
    }
#                                                                        #
#                                                                        #
##########################################################################
  } else {

    print "summarize: ERROR! Not all jobs ended (only $cnt found).",
          " Cannot analyse results\n";
  }
}
#*****************************************************************************#
sub getNakedName {
 my @pathAndName = split /\//, $_[0];
 return $pathAndName[$#pathAndName];
}
#*****************************************************************************#
sub mustDump {
  return  $ENV{PARALLEL_BATCH_DUMP} || undef
}

sub summarize_getContainers {
  my @cntnrs = ();
  unless (open (REFCARD, "$_[0]")) {
    print STDERR "summarize: cannot open $_[0].\n";
    return ();
  } else {
    while (<REFCARD>) {
      chomp;
      next if /\s+\/\//;
      if (s/\s*Data container_unitarity\s+//i) {
	push @cntnrs, $_;
      }
    }
    my $tmpCntnrs = join (" ", @cntnrs);
    @cntnrs = ( split " ", $tmpCntnrs);
    return @cntnrs;
  }
}

sub summarize_submitAgain {
  print "submitAgain Called!!!!\n";
  #open (REPFILE, $_[0]) || die "summarize: cannot open $_[0]. $!\n";
#foreach 
}

sub summarize_hasBookMarked {
  my $found = 0;
  open (INBKMRK, $_[0]) or 
    die "summarize: cannot open $_[0]. $!\n";
  #print "\n\t\tLooking for $_[1]:\n";
  while(<INBKMRK>) {
    chomp;
    #print "\t\t>$_\n";
    if (/$_[1]/) {
      $found = 1;
      last;
    }
  }
  close INBKMRK;
  return $found;
}
#--------------------------------------------------------------------------#
sub getRootHistFrom {
  my $found = 0;
  my $name;
  open (INOPTFILE, $_[0]) or 
    die "summarize\:\:getRootHistFrom: cannot open $_[0]. $!\n";
  while(<INOPTFILE>) {
    chomp;
    if (/\s*histograms\s+package\s+ROOT/) {
      $found = 1;  # could be after "histograms home" so only aknowledge...
    } elsif (/\s*histograms\s+home/) {
      my @tokens = split /\s+/, $_;
      $name = $tokens[2];
    }
  }
  close INOPTFILE;
  return $name if $found;
  return undef;
}
#--------------------------------------------------------------------------#
sub SumHistograms {
  my $thisDir = shift @_;
  my $absOutDir = shift @_;
  #my $filesToProcess = join ("\",\"", @_);
  my @cmdRoot = ();
  my $upTo = $#_+2;
  print "summing to $_[0]\n";
  push @cmdRoot, (
		  "root -b -l -q -n $absOutDir << END_ROOT_EXECUTE", 
		  ".L $thisDir/histoadd.C",
		  "char \*cf[$upTo];",
		  "cf[0] = new char[50]; strcpy(cf[0], \"\");",
		  "cf[1] = new char[50]; strcpy(cf[1], \"$_[0]\")"
		  );
  
  for (my $i = 2; $i <= $#_+1; $i++) {
    push @cmdRoot, "cf[$i] = new char[50]; strcpy(cf[$i], \"$_[$i-1]\");";
  }
    
	#	  "const char \*cfiles[] = {\"\", \"$filesToProcess\"};",
  push @cmdRoot, (
		  "main($upTo, cf);",
		  "END_ROOT_EXECUTE",
		  "\n"
		  );
  my $cmdRoot = join ("\n", @cmdRoot);
  print "SUMMING HISTOGRAMS...\n";
  my @sumCmd = `$cmdRoot 2>&1`;
  printOut(@sumCmd);
}

#-----------------------------------------------------------------------------
#
sub printOut {
  foreach (@_) {
    next if /^\s+$/;
    next if /^\*\*\*/;
    next if /^  * /;
    next if /^Compiled with/;
    next if /^CINT\/ROOT/;
    next if /^Type \? for/;
    next if /^Enclose mul/;
    print "$_";
  }
}
#----------------------------------------------------------------------------#

sub setOptionsValues {
  my $i = 0;
  while ($i < @_) {
    if ($_[$i] =~/-s|--sync/) {
      if ($#_ >= $i) {
        unless ($_[$i+1] =~ /^\d+$/) {
	  die "$nkd0: \"$_[$i+1]\" is illegal option for \"sync\" switch\n";
	}
	my $syncTime = $_[$i+1];
	if ($syncTime == 1) { 
	  $sync = "1";
	} else {
	  my $now = 0 + time();
	  if ($syncTime < $now) {
	    print "$nkd0: WARNING! synchronization will not be performed, as ",
	    "given time $syncTime is less than $now\n";
	  } else {
	    $sync = "$syncToPgm $syncTime;"; # Watch out: ";" must be there!!!
	  }
	}
      }
    } elsif ($_[$i] =~/-q|--queue/) {
      if ($#_ >= $i) {
	$q_1 = "-q $_[$i+1]";
      } else {
	die "$nkd0: \"queue\" switch requires an LSF queue name.\n";
      }
    } elsif ($_[$i] =~/-R|--Resource/) {
      if ($#_ >= $i) {
	$q_1 = "-q $_[$i+1]";
      } else {
	die "$nkd0: \"Resource\" switch requires an LSF Resource name.\n";
      }     
    } else {
      die "$nkd0: \"$_[$i]\" is not a valid switch.\n$synopsis\n";
    }
    $i++;
  }
}


