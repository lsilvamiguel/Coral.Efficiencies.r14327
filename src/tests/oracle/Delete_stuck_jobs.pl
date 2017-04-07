#!/usr/local/bin/perl -w

my $username    = `whoami`;
my $slot_number = 1;
my $year        = 2003;
my $period      = "P1J";
my $log_dir     = "/afs/cern.ch/compass/scratch/mdst2/production/logs/$year/$period/slot$slot_number";

chomp $username;

print "USERNAME: $username\n";
my @job_list = `bjobs -w -r`;

my $n_kiled_jobs = 0;
my @list_of_killed_jobs = ();

foreach my $job (@job_list) {
    if($job =~ /(\d+)\s+$username\s+RUN\s+(\S+)\s+\S+\s+(\S+)\s+(\d+)\-(\d+)/) {
	my $job_id       = $1;
	my $queue        = $2;
	my $machine      = $3;
	my $chunk_number = $4;
	my $run_number   = $5;
	unless(-e "./Run_$run_number/$chunk_number-$run_number.lsf") { #this is another period's run
	    next;
	}
	my $log_age = (-M "$log_dir/$run_number/$run_number-$chunk_number-$slot_number.log");
	unless(defined $log_age) {
	    print "Could not find file: $log_dir/$run_number/$run_number-$chunk_number-$slot_number.log\n";
	    $log_age = 1;
	}
	if($log_age > 0.1) {
	    my $is_finished = 0;
	    if($log_age != 1) {
		open LOG_FILE, "$log_dir/$run_number/$run_number-$chunk_number-$slot_number.log";
		while(<LOG_FILE>) {
		    if($_ =~ /CORAL END\:/) {
			$is_finished = 1;
			last;
		    }
		}
		close LOG_FILE;
	    }
	    unless($is_finished) {
		$n_kiled_jobs++;
		print "Job $chunk_number-$run_number looks to be stuck. Its last changing was $log_age day ago. It will be killed now.\n";
		push(@list_of_killed_jobs,sprintf("%7d %7s %10s %05d-%05d %.3f",$job_id, $queue, $machine, $chunk_number, $run_number, $log_age));
		if(system("bkill $job_id") == 0) {
		    unlink "$log_dir/$run_number/$run_number-$chunk_number-$slot_number.log";
		}
	    }
	}
    }
}

print "Total number of killed jobs: $n_kiled_jobs\n";
foreach my $info (@list_of_killed_jobs) {
    print "$info\n";
}
