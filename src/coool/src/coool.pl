#! /usr/bin/perl -w

#$ENV{ROOTSYS} = "/online/soft/root/root_v5.10.00";
$ENV{ROOTSYS} = "/online/soft/root/root_v5.28.00d";
$ENV{LD_LIBRARY_PATH} = "$ENV{ROOTSYS}/lib:$ENV{ROOTSYS}";
$ENV{PATH} = "$ENV{ROOTSYS}/bin:$ENV{PATH}";

# my $geomfile = "/online/detector/geometry/2008/detectors.69549.hadron.dat";
my $geomfile = "/online/detector/geometry/2011/detectors.91413.plus.dat";

my $EVB = $ENV{COOOL_DEFAULT_EVB};

if($#ARGV == 0) {
    $EVB = $ARGV[0];
}
else {
    my $mysql_cmd = "select EBelems FROM tb_actEB WHERE actEBid=(SELECT actEBid FROM tb_run WHERE runnb=(SELECT max(runnb) FROM tb_run WHERE 1))";
    my $evb_id_list = `echo \"$mysql_cmd\" | mysql -h pccodb00 -u onl -pna58db runlb | tail -1`;

    chop $evb_id_list;

    #print "$evb_id_list\n";

    my $isEnd = 0;
    my $str = $evb_id_list;
    my @evb_list;
    do {
	my $evb = 0;
	if($str =~ s/(\d+),(\S+)//) {
	    $evb = $1;
	    $str = $2;
	}
	else {
	    $evb = $str;
	    $isEnd = 1;
	}
	push(@evb_list,"pccoeb$evb");
    } until($isEnd);

    my $nEvbs = 0;

    my $min_num_clients = 100;

    foreach (@evb_list) {
	my $evb = $_;
	$nEvbs++;
#check monitoring clients
	my $n_clients = `ssh daq\@$evb \"/date/monitoring/Linux/monitorClients | grep yes | wc -l\"`;
	chop $n_clients;
	#print "$evb: number of clients:  $n_clients\n";
	if($n_clients == 0) {
	    $EVB = $evb;
	    last;
	}
	else {
	    if($n_clients < $min_num_clients) {
		$min_num_clients = $n_clients;
		$EVB = $evb;
	    }
	}
    }
    $EVB = "\@$EVB:";
}

my $CMD = "./coool -map /online/detector/maps -geom $geomfile -group /online/detector/monitor/groups.xml -cfg /online/detector/monitor/default_params -workenv DAQ $EVB";

print("$CMD \n");

`$CMD`;
