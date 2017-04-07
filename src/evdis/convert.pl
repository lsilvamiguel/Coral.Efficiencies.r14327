#!/usr/bin/perl -w

# small script to modify a detectors.C create by COMGeant and g2root so that
# it is usable with the EVE event display
# * transform axis:
#    CG | new
#     X |  Z
#     Y |  X
#     Z |  Y

use strict;

sub usage {
    print "$0 [original detectors.C] [modified detectors.C]\n";
    print "If the second file is not specified, then 'detectors.C' is used.\n";
}

if (@ARGV < 1 || @ARGV > 2) {
    usage;
    exit;
}

my $orgfile = $ARGV[0];
my $newfile = "detectors.C";
if (@ARGV==2) {
    $newfile = $ARGV[1];
}

my @rotmats;
# largest index of a rotation matrix to know which numbers to use for the matrices added
my $rotmatLastIdx = 0;
# index of rotation matrix 'unity' in new coordinate system
my $rotmatNewUnity = 0;
my $rotmatNewUnityStr = "0,0,90,0,90,90";
# rotation matrices used to add something to the hall (top volume)
my @rotmatsToMod;
# ... and their replacement
my @rotmatsToModRepl;

# first pass
# get rotation matrices
open(IN,  "<$orgfile");
while(<IN>) {
    my $line = $_;
    if ($line =~ /new TGeoRotation/) {
        # extract rotation matrix
        $line =~ /new TGeoRotation\("rot(\d+)",([-+\d\.E]+,[-+\d\.E]+,[-+\d\.E]+,[-+\d\.E]+,[-+\d\.E]+,[-+\d\.E]+)\)/;
        if ($1 > $rotmatLastIdx) {
            $rotmatLastIdx = $1;
        }
        if ($2 eq $rotmatNewUnityStr) {
            $rotmatNewUnity = $1;
        }
        $rotmats[$1] = $2;
    } elsif ($line =~ /HALL->AddNode\(\w+,\d+,new TGeoCombiTrans\([-+\d\.E]+,[-+\d\.E]+,[-+\d\.E]+,rot(\d+)\)\);/) {
        my $found=0;
        for (my $i=0; $i<@rotmatsToMod; $i++) {
            if ($rotmatsToMod[$i] == $1) {
                $found=1;
            }
        }
        if ($found==0) {
            $rotmatsToMod[@rotmatsToMod] = $1;
        }
    }
}
close(IN);

# add new required rotation matrices
if ($rotmatNewUnity == 0) {
    $rotmatNewUnity=@rotmats;
    $rotmats[$rotmatNewUnity] = $rotmatNewUnityStr;
}
for (my $i=0; $i<@rotmatsToMod; $i++) {
    my $newrotmat = $rotmats[$rotmatsToMod[$i]];

    if (-e "convert.C") {
        print "File 'convert.C' does exist, a temporary script cannot be written. Cancelling!\n";
        exit(1);
    }

    open(TMP, ">convert.C");
    print TMP <<FINISHED;
void convert() {
    TGeoRotation A("rot1",0,0,90,0,90,90);
    TGeoRotation AT=A.Inverse();
    TGeoRotation M("rotM",$newrotmat);

    std::cout << A.Determinant() << " " << AT.Determinant() << " " << M.Determinant() << std::endl;

    M.Print();

    A.MultiplyBy(&M);
    //A.MultiplyBy(&AT);

    A.Print();

    double t1,p1,t2,p2,t3,p3;
    A.GetAngles(t1,p1,t2,p2,t3,p3);
    std::cout << t1 << "," << p1 << "," << t2 << "," << p2 << "," << t3 << "," << p3 << std::endl;
}
FINISHED
    close(TMP);

    #system "root -q -b -l convert.C";
    $newrotmat = `root -q -b -l convert.C | tail -n1`;
    chomp $newrotmat;
    unlink "convert.C";

    my $found=0;
    for (my $j=0; $j<@rotmats; $j++) {
        if (defined($rotmats[$j]) && $rotmats[$j] eq $newrotmat) {
            $rotmatsToModRepl[$rotmatsToMod[$i]] = $j;
            $found=1;
        }
    }
    if ($found==0) {
        $rotmatsToModRepl[$rotmatsToMod[$i]] = @rotmats;
        $rotmats[@rotmats] = $newrotmat;
    }
}

# second pass
# actually modify the file
open(IN,  "<$orgfile");
open(OUT, ">$newfile");

while(<IN>) {
    my $line = $_;
    
    if ($line =~ /HALL/) {
        if ($line =~ /^gGeoManager->SetTopVolume\(HALL\);$/) {
            # do nothing
        } elsif ($line =~ /TGeoVolume \*HALL = gGeoManager->MakeBox\("HALL"/) {
            my $oldline = $line;
            $line =~ s/MakeBox\("HALL",med(\d+),([-+\d\.E]+),([-+\d\.E]+),([-+\d\.E]+)\)/MakeBox\("HALL",med$1,$3,$4,$2\)/g;
            if ($oldline eq $line) {
                print "Unchanged line:\n$line";
            }
        } elsif ($line =~ /^ HALL->AddNode\(\w+,\d+,new TGeoTranslation\(([-+\d\.E]+),([-+\d\.E]+),([-+\d\.E]+)\)\);$/) {
            my $oldline = $line;
            $line =~ s/new TGeoTranslation\(([-+\d\.E]+),([-+\d\.E]+),([-+\d\.E]+)\)/new TGeoCombiTrans\($2,$3,$1,rot$rotmatNewUnity\)/g;
            if ($oldline eq $line) {
                print "Unchanged line:\n$line";
            }
        } elsif ($line =~ /^ HALL->AddNode\(\w+,\d+,new TGeoCombiTrans\(([-+\d\.E]+),([-+\d\.E]+),([-+\d\.E]+),rot\d+\)\);$/) {
            my $oldline = $line;
            $line =~ s/new TGeoCombiTrans\(([-+\d\.E]+),([-+\d\.E]+),([-+\d\.E]+),rot(\d+)\)/new TGeoCombiTrans\($2,$3,$1,rot$rotmatsToModRepl[$4]\)/g;
            if ($oldline eq $line) {
                print "Unchanged line:\n$line";
            }
        } else {
            print "Unknown line:\n$line";
        }
    }

    print OUT $line;

    # append new rotation matrices to 
    if ($line =~ /new TGeoRotation\("rot$rotmatLastIdx",/) {
        for (my $i=$rotmatLastIdx+1; $i<@rotmats; $i++) {
            print OUT "TGeoRotation *rot$i = new TGeoRotation(\"rot$i\",$rotmats[$i]);\n";
        }
    }
}

close(IN);
close(OUT);
