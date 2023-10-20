
#!/usr/bin/perl -w
use strict;
use threads;
#require "baseFun.pl";

sub c_primerLoc{
        my ($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand)=@_;
        my ($e_up_a,$e_up_b,$e_down_a,$e_down_b,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len);
        if($t_type eq "ctag"){
                if($strand eq "+"){
                        $e_up_a=$cdsEnd-$t_upGeneFar;
                        $e_up_b=$cdsEnd-$t_upGeneNear;
                        $v_up_a=$cdsEnd-$t_upVfar;
                        $v_up_b=$cdsEnd-$t_upVnear;
                        
                        if($e_up_b<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (upGeneNear,upGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($e_up_a<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (upGeneNear,upGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_up_b<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (upVnear,upVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_up_a<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (upVnear,upVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }
                        #
                        if($v_up_b-$v_up_a<17){
                                print "Error: The interval was too narrow.\n";
                                die "Try to modify the config file (upVerifyNear,upVerifyFar).\n";
                        }
                        if($e_up_b-$e_up_a<17){
                                print "Error: The interval was too narrow.\n";
                                die "Try to modify the config file (upGeneNear,upGeneFar).\n";
                        }
                #
                        $e_down_a=$cdsEnd+$t_downGeneFar;
                        $e_down_b=$cdsEnd+$t_downGeneNear;
                        
                        $v_down_a=$cdsEnd+$t_downVfar;
                        $v_down_b=$cdsEnd+$t_downVnear;
                        $v_len=$t_upVfar+$t_downVfar+1;
                        
                        if($e_down_b>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (downGeneNear,downGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($e_down_a>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (downGeneNear,downGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_down_b>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (downVnear,downVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_down_a>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (downVnear,downVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }
                        
                        if($e_down_a-$e_down_b<17){
                                print "Error: The interval was too narrow.\n";
                                print "Try to modify the config file (downGeneNear,downGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }
                        
                        if($v_down_a-$v_down_b<17){
                                print "Error: The interval was too narrow.\n";
                                print "Try to modify the config file (downVnear,downVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }
                        
                }else{
                        $e_up_a=$cdsStart-$t_downGeneFar;
                        $e_up_b=$cdsStart-$t_downGeneNear;
                        $v_up_a=$cdsStart-$t_downVfar;
                        $v_up_b=$cdsStart-$t_downVnear;
                        #
                        if($e_up_b<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (downGeneNear,downGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($e_up_a<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (downGeneNear,downGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_up_b<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (downVnear,downVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_up_a<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (downVnear,downVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }
                        #
                        if($v_up_b-$v_up_a<17){
                                print "Error: The interval was too narrow.\n";
                                die "Try to modify the config file (downVerifyNear,downVerifyFar).\n";
                        }
                        if($e_up_b-$e_up_a<17){
                                print "Error: The interval was too narrow.\n";
                                die "Try to modify the config file (downGeneNear,downGeneFar).\n";
                        }
                #
                        $e_down_a=$cdsStart+$t_upGeneFar;
                        $e_down_b=$cdsStart+$t_upGeneNear;
                        
                        $v_down_a=$cdsStart+$t_upVfar;
                        $v_down_b=$cdsStart+$t_upVnear;
                        $v_len=$t_upVfar+$t_downVfar+1;
                        
                        if($e_down_b>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (upGeneNear,upGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($e_down_a>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (upGeneNear,upGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_down_b>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (upVnear,upVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_down_a>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr cds_start: $cdsStart cds_end: $cdsEnd strand: $strand).\n";
                                print "Try to modify the config file (upVnear,upVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }
                        
                        if($e_down_a-$e_down_b<17){
                                print "Error: The interval was too narrow.\n";
                                die "Try to modify the config file (upGeneNear,upGeneFar).\n";
                        }
                        
                        if($v_down_a-$v_down_b<17){
                                print "Error: The interval was too narrow.\n";
                                die "Try to modify the config file (upVnear,upVfar).\n";
                        }
                }
                return($e_up_a,$e_up_b,$e_down_a,$e_down_b,$cdsStart,$cdsEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len);
        }
}

sub c_geneNameInterval{
        my ($t_gff,$t_geneName,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,$refLen)=@_;
		my @out;
        print "Reading GFF file. Search gene\n";
		if($t_gff=~/\.gz$/){
		        (open IN,"gzip -dc $t_gff|") || die "$t_gff open failed. $!\n";
		}else{
		        (open IN,"$t_gff") || die "$t_gff open failed. $!\n";
		}
		
		my ($chr,$c_end,$geneName,$a_geneName);
        my $id;
        my $dbxref;
        
        my ($gStart,$gEnd,$strand);
        my ($cdsStart,$cdsEnd);
        my @cds;
        my $indicate=0;
        my $c_mRNA=0;
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                        if($indicate==1){
                                if(exists ${$refLen}{$chr}){
                                        $c_end=$refLen->{$chr};
                                        if(uc($t_geneName) eq uc($geneName) || $t_geneName eq $id || $t_geneName eq $dbxref){
                                                print "Gene information: $a_geneName\n";
                                                my $n=scalar(@cds);
                                                if($n<1){
                                                        print "Warning: no CDS.\n";
                                                }else{
                                                        if($strand eq "+"){   
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[$n-1][1];
                                                        }else{
                                                                if($cds[0][0]>$cds[$n-1][0]){
                                                                        $cdsStart=$cds[$n-1][0];
                                                                        $cdsEnd=$cds[$n-1][1];
                                                                }else{
                                                                        $cdsStart=$cds[0][0];
                                                                        $cdsEnd=$cds[0][1];
                                                                }
                                                        }
                                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                                if($tout[0] == 0){
                                                                        die "Gene locate in the terminal of chromosome.\n";
                                                                }
                                                                @out=($a_geneName,$chr,@tout,$strand);
                                                        }else{
                                                                print "Warning: region of CDS out of gene's.\n";
                                                        }
                                                }
                                                #
                                                if($c_mRNA>1){
                                                        print "$a_geneName mRNA : $c_mRNA\n";
                                                }
                                                $indicate=0;
                                                last;
                                        }
                                }else{
                                        print "Warning: chromosome length unknown (chr:$chr).\n";
                                }
                        }
                        ############################
                        $chr=$line[0];
                        $strand=$line[6];
                        
                        if($line[8]=~/Name=(.+?);/){
						        $geneName=$1;
						}elsif($line[8]=~/Name=(.+)/){
						        $geneName=$1;
						}else{
                                $geneName="unknown";
                        }
                        
                        if($line[8]=~/ID=gene[:,-](.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=gene[:,-](.+)/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+)/){
						        $id=$1;
						}else{
                                $id="unknown";
                        }
                        
                        if($line[8]=~/Dbxref=GeneID:(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=GeneID:(.+)/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+)/){
						        $dbxref=$1;
						}else{
                                $dbxref="unknown";
                        }
                        $a_geneName="${id}|${geneName}|$dbxref";
                        $a_geneName=~tr/\,/_/;
						# a b ===== b a
						# =====a b=== b a
						if($line[3]>$line[4]){
						        $gStart=$line[4];
								$gEnd=$line[3];
						}else{
						        $gStart=$line[3];
								$gEnd=$line[4];
						}
                        ###########
                        @cds=();
                        $indicate=0;
                        $c_mRNA=0;
                }elsif($line[2] eq "mRNA"){
                        $c_mRNA++;
                }elsif($line[2] eq "CDS"){
                        if($c_mRNA>1){
                                next;
                        }
                        push(@cds,[$line[3],$line[4]]);
                        $indicate=1;
                }				
		}
        #
        if($indicate==1){
                if(exists ${$refLen}{$chr}){
                        $c_end=$refLen->{$chr};
                        if(uc($t_geneName) eq uc($geneName) || $t_geneName eq $id || $t_geneName eq $dbxref){
                                print "Gene information: $a_geneName\n";
                                my $n=scalar(@cds);
                                if($n<1){
                                        print "Warning: no CDS.\n";
                                }else{
                                        if($strand eq "+"){   
                                                $cdsStart=$cds[$n-1][0];
                                                $cdsEnd=$cds[$n-1][1];
                                        }else{
                                                if($cds[0][0]>$cds[$n-1][0]){
                                                        $cdsStart=$cds[$n-1][0];
                                                        $cdsEnd=$cds[$n-1][1];
                                                }else{
                                                        $cdsStart=$cds[0][0];
                                                        $cdsEnd=$cds[0][1];
                                                }
                                        }
                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                if($tout[0] == 0){
                                                        die "Gene locate in the terminal of chromosome.\n";
                                                }
                                                @out=($a_geneName,$chr,@tout,$strand);
                                        }else{
                                                print "Warning: region of CDS out of gene's.\n";
                                        }
                                }
                                #
                                if($c_mRNA>1){
                                        print "$a_geneName mRNA : $c_mRNA\n";
                                }
                        }
                }else{
                        print "Warning: chromosome length unknown (chr:$chr).\n";
                }
        }
		close IN;
		if(scalar(@out)>0){
		        return(@out);
        }else{
                die "GeneName cat't be found in the gff file.\n";
		}
}

sub c_nameListInterval{
        my ($t_gff,$t_geneList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,$refLen)=@_;
		# read name list
		my %nameHash;
        my %r_hash;
        print "Reading gene list\n";
		(open IN,"$t_geneList") || die "$t_geneList $!\n";
		while(<IN>){
		        chomp;
				next if(/^\s*$/);
				my @line=split;
				$nameHash{uc($line[0])}=1;
                #push(@$r_ordGene,$line[0]);
		}
		close IN;
		if(scalar(keys %nameHash)<1){
                die "Error: $t_geneList has no valid line.\n";
        }
		# read gff
        print "Reading GFF file. Search gene\n";
		if($t_gff=~/\.gz$/){
		        (open IN,"gzip -dc $t_gff|") || die "$t_gff open failed. $!\n";
		}else{
		        (open IN,"$t_gff") || die "$t_gff open failed. $!\n";
		}
		
        my ($chr,$c_end,$geneName,$a_geneName);
        my $id;
        my $dbxref;
        my @out;
        
        my ($gStart,$gEnd,$strand);
        my ($cdsStart,$cdsEnd);
        my @cds;
        my $indicate=0;
        my $c_mRNA=0;
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                        if($indicate==1){
                                if(exists ${$refLen}{$chr}){
                                        $c_end=$refLen->{$chr};
                                        if(exists $nameHash{uc($geneName)} || $nameHash{$id} || $nameHash{$dbxref}){
                                                print "Gene information: $a_geneName\n";
                                                my $n=scalar(@cds);
                                                if($n<1){
                                                        print "Warning: no CDS.\n";
                                                }else{
                                                        if(exists $r_hash{$a_geneName}){
                                                                print "Warning: gene name duplicate in the GFF file.\n";
                                                                print "The gene on chromosome $chr, location $gStart -- $gEnd will be skipped.\n"; 
                                                        }else{
                                                                #print "Gene information: $a_geneName\n";
                                                                if($strand eq "+"){
                                                                        $cdsStart=$cds[$n-1][0];
                                                                        $cdsEnd=$cds[$n-1][1];
                                                                }else{
                                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                                $cdsStart=$cds[$n-1][0];
                                                                                $cdsEnd=$cds[$n-1][1];
                                                                        }else{
                                                                                $cdsStart=$cds[0][0];
                                                                                $cdsEnd=$cds[0][1];
                                                                        }
                                                                }
                                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                        my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                                        if($tout[0] > 0){
                                                                                my $str=join ",",$a_geneName,$chr,@tout,$strand;
                                                                                $r_hash{$a_geneName}=1;
                                                                                push(@out,"$str");
                                                                        }
                                                                }else{
                                                                        print "Warning: region of CDS was out of gene's.\n";
                                                                }
                                                        }
                                                }
                                                #
                                                if($c_mRNA>1){
                                                        print "$a_geneName mRNA : $c_mRNA\n";
                                                }
                                        }
                                }else{
                                        print "Warning: chromosome length unknown (chr:$chr).\n";
                                }
                        }
                        ############################
                        $chr=$line[0];
                        $strand=$line[6];
                        if($line[8]=~/Name=(.+?);/){
						        $geneName=$1;
						}elsif($line[8]=~/Name=(.+)/){
						        $geneName=$1;
						}else{
                                $geneName="unknown";
                        }
                        
                        if($line[8]=~/ID=gene[:,-](.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=gene[:,-](.+)/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+)/){
						        $id=$1;
						}else{
                                $id="unknown";
                        }
                        
                        if($line[8]=~/Dbxref=GeneID:(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=GeneID:(.+)/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+)/){
						        $dbxref=$1;
						}else{
                                $dbxref="unknown";
                        }
                        $a_geneName="${id}|${geneName}|$dbxref";
                        $a_geneName=~tr/\,/_/;
						# a b ===== b a
						# =====a b=== b a
						if($line[3]>$line[4]){
						        $gStart=$line[4];
								$gEnd=$line[3];
						}else{
						        $gStart=$line[3];
								$gEnd=$line[4];
						}
						###########
                        @cds=();
                        $indicate=0;
                        $c_mRNA=0;
                }elsif($line[2] eq "mRNA"){
                        $c_mRNA++;
                }elsif($line[2] eq "CDS"){
                        if($c_mRNA>1){
                                next;
                        }
                        push(@cds,[$line[3],$line[4]]);
                        $indicate=1;
                }	
		}
        #
        if($indicate==1){
                if(exists ${$refLen}{$chr}){
                        $c_end=$refLen->{$chr};
                        if(exists $nameHash{uc($geneName)} || $nameHash{$id} || $nameHash{$dbxref}){
                                print "Gene information: $a_geneName\n";
                                my $n=scalar(@cds);
                                if($n<1){
                                        print "Warning: no CDS.\n";
                                }else{
                                        if(exists $r_hash{$a_geneName}){
                                                print "Warning: gene name duplicate in the GFF file.\n";
                                                print "The gene on chromosome $chr, location $gStart -- $gEnd will be skipped.\n"; 
                                        }else{
                                                #print "Gene information: $a_geneName\n";
                                                if($strand eq "+"){
                                                        $cdsStart=$cds[$n-1][0];
                                                        $cdsEnd=$cds[$n-1][1];
                                                }else{
                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[$n-1][1];
                                                        }else{
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }
                                                }
                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                        my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                        if($tout[0] > 0){
                                                                my $str=join ",",$a_geneName,$chr,@tout,$strand;
                                                                $r_hash{$a_geneName}=1;
                                                                push(@out,"$str");
                                                        }
                                                }else{
                                                        print "Warning: region of CDS was out of gene's.\n";
                                                }
                                        }
                                }
                                #
                                if($c_mRNA>1){
                                        print "$a_geneName mRNA : $c_mRNA\n";
                                }
                        }
                }else{
                        print "Warning: chromosome length unknown (chr:$chr).\n";
                }
        }
		close IN;
		if(scalar (@out)<1){
                die "GeneName cat't be found in the GFF file.\n";
		}else{
                return(@out);
        }
}

sub c_coordInterval{
        my ($t_gff,$t_coord,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,$refLen)=@_;
		my @out;
        print "Reading GFF file. Search gene\n";
		if($t_gff=~/\.gz$/){
		        (open IN,"gzip -dc $t_gff|") || die "$t_gff open failed. $!\n";
		}else{
		        (open IN,"$t_gff") || die "$t_gff open failed. $!\n";
		}
		
        my @s=split /:/,$t_coord;
        if(scalar(@s) != 2){
                die "Error: '--coordinate $t_coord' format error.\n";
        }
        my ($s_chr,$s_pos)=@s;
		my ($chr,$c_end,$geneName,$a_geneName);
        my $id;
        my $dbxref;
        
        my ($gStart,$gEnd,$strand);
        my ($cdsStart,$cdsEnd);
        my @cds;
        my $indicate=0;
        my $c_mRNA=0;
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                        if($indicate==1){
                                if(exists ${$refLen}{$chr}){
                                        $c_end=$refLen->{$chr};
                                        if($chr eq $s_chr){
                                                if($s_pos >= $gStart && $s_pos <= $gEnd){
                                                        print "Gene information: $a_geneName\n";
                                                        if($a_geneName ne "unknown|unknown|unknown"){
                                                                my $n=scalar(@cds);
                                                                if($n<1){
                                                                        print "Warning: no CDS.\n";
                                                                }else{
                                                                        if($strand eq "+"){   
                                                                                $cdsStart=$cds[$n-1][0];
                                                                                $cdsEnd=$cds[$n-1][1];
                                                                        }else{
                                                                                if($cds[0][0]>$cds[$n-1][0]){
                                                                                        $cdsStart=$cds[$n-1][0];
                                                                                        $cdsEnd=$cds[$n-1][1];
                                                                                }else{
                                                                                        $cdsStart=$cds[0][0];
                                                                                        $cdsEnd=$cds[0][1];
                                                                                }
                                                                        }
                                                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                                my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                                                if($tout[0] == 0){
                                                                                        die "Gene locate in the terminal of chromosome.\n";
                                                                                }
                                                                                @out=($a_geneName,$chr,@tout,$strand);
                                                                        }else{
                                                                                print "Warning: region of CDS out of gene's.\n";
                                                                        }
                                                                }
                                                        }
                                                        #
                                                        if($c_mRNA>1){
                                                                print "$a_geneName mRNA : $c_mRNA\n";
                                                        }
                                                        $indicate=0;
                                                        last;
                                                }
                                        }
                                }else{
                                        print "Warning: chromosome length unknown (chr:$chr).\n";
                                }
                        }
                        ############################
                        $chr=$line[0];
                        $strand=$line[6];
                        
                        if($line[8]=~/Name=(.+?);/){
						        $geneName=$1;
						}elsif($line[8]=~/Name=(.+)/){
						        $geneName=$1;
						}else{
                                $geneName="unknown";
                        }
                        
                        if($line[8]=~/ID=gene[:,-](.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=gene[:,-](.+)/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+)/){
						        $id=$1;
						}else{
                                $id="unknown";
                        }
                        
                        if($line[8]=~/Dbxref=GeneID:(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=GeneID:(.+)/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+)/){
						        $dbxref=$1;
						}else{
                                $dbxref="unknown";
                        }
                        $a_geneName="${id}|${geneName}|$dbxref";
                        $a_geneName=~tr/\,/_/;
						# a b ===== b a
						# =====a b=== b a
						if($line[3]>$line[4]){
						        $gStart=$line[4];
								$gEnd=$line[3];
						}else{
						        $gStart=$line[3];
								$gEnd=$line[4];
						}
						###########
                        @cds=();
                        $indicate=0;
                        $c_mRNA=0;
                }elsif($line[2] eq "mRNA"){
                        $c_mRNA++;
                }elsif($line[2] eq "CDS"){
                        if($c_mRNA>1){
                                next;
                        }
                        push(@cds,[$line[3],$line[4]]);
                        $indicate=1;
                }				
		}
        #
		if($indicate==1){
                if(exists ${$refLen}{$chr}){
                        $c_end=$refLen->{$chr};
                        if($chr eq $s_chr){
                                if($s_pos >= $gStart && $s_pos <= $gEnd){
                                        print "Gene information: $a_geneName\n";
                                        if($a_geneName ne "unknown|unknown|unknown"){
                                                my $n=scalar(@cds);
                                                if($n<1){
                                                        print "Warning: no CDS.\n";
                                                }else{
                                                        if($strand eq "+"){   
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[$n-1][1];
                                                        }else{
                                                                if($cds[0][0]>$cds[$n-1][0]){
                                                                        $cdsStart=$cds[$n-1][0];
                                                                        $cdsEnd=$cds[$n-1][1];
                                                                }else{
                                                                        $cdsStart=$cds[0][0];
                                                                        $cdsEnd=$cds[0][1];
                                                                }
                                                        }
                                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                                if($tout[0] == 0){
                                                                        die "Gene locate in the terminal of chromosome.\n";
                                                                }
                                                                @out=($a_geneName,$chr,@tout,$strand);
                                                        }else{
                                                                print "Warning: region of CDS out of gene's.\n";
                                                        }
                                                }
                                        }
                                        #
                                        if($c_mRNA>1){
                                                print "$a_geneName mRNA : $c_mRNA\n";
                                        }

                                }
                        }
                }else{
                        print "Warning: chromosome length unknown (chr:$chr).\n";
                }
        }
		close IN;
		if(scalar(@out)>0){
		        return(@out);
        }else{
                die "GeneName cat't be found in the gff file.\n";
		}
}

sub c_coordListInterval{ 
        my ($t_gff,$t_coordList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,$refLen)=@_;
		# read name list
		my %posHash;
        my %r_hash;
        print "Reading coordinate list\n";
		(open IN,"$t_coordList") || die "$t_coordList open failed. $!\n";
		while(<IN>){
		        chomp;
				next if(/^\s*$/);
				my @line=split /:/,$_;
                if(scalar(@line) != 2){
                        print "Warning: bad format. $_\n";
                }else{
                        $posHash{$line[0]}{$line[1]}=1;
                }
		}
		close IN;
        if(scalar(keys %posHash)<1){
                die "Error: $t_coordList has no valid line.\n";
        }		
		# read gff
        print "Reading GFF file. Search gene\n";
		if($t_gff=~/\.gz$/){
		        (open IN,"gzip -dc $t_gff|") || die "$t_gff open failed. $!\n";
		}else{
		        (open IN,"$t_gff") || die "$t_gff open failed. $!\n";
		}
		
		my ($chr,$c_end,$geneName,$a_geneName);
        my ($gStart,$gEnd,$strand);
        my ($cdsStart,$cdsEnd);
        my @cds;
        my @out;
        my $id;
        my $dbxref;
        my $indicate=0;
        my $c_mRNA=0;
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                        if($indicate==1){
                                if(exists ${$refLen}{$chr}){
                                        $c_end=$refLen->{$chr};
                                        my $n=scalar(@cds);
                                        #if($n<1){print "Warning: no CDS.\n";}
                                        if($n>=1){
                                                if($strand eq "+"){       
                                                        $cdsStart=$cds[$n-1][0];
                                                        $cdsEnd=$cds[$n-1][1];
                                                }else{
                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[$n-1][1];
                                                        }else{
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }
                                                }
                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                        foreach my $pos(keys %{$posHash{$chr}}){
                                                                if($pos>=$gStart && $pos<=$gEnd){
                                                                        print "Gene information: $chr $pos $a_geneName\n";
                                                                        if($a_geneName eq "unknown|unknown|unknown"){
                                                                                next;
                                                                        }
                                                                        
                                                                        if(exists $r_hash{$a_geneName}){
                                                                                print "Warning: gene name duplicate in the GFF file.\n";
                                                                                print "The gene on chromosome $chr, location $gStart -- $gEnd will be skipped.\n"; 
                                                                        }else{
                                                                                my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                                                if($tout[0] > 0){
                                                                                        my $str=join ",",$a_geneName,$chr,@tout,$strand;
                                                                                        $r_hash{$a_geneName}=1;
                                                                                        push(@out,"$str");
                                                                                }
                                                                        }
                                                                        
                                                                }
                                                                #
                                                                if($c_mRNA>1){
                                                                        print "$a_geneName mRNA : $c_mRNA\n";
                                                                }
                                                        }
                                                        
                                                }else{
                                                        print "Gene information: $chr ${gStart}_$gEnd $a_geneName\n";
                                                        print "Warning: region of CDS out of gene's.\n";
                                                }
                                        }
                                }else{
                                        print "Warning: chromosome length unknown (chr:$chr).\n";
                                }
                        }
                        ############################
                        $chr=$line[0];
                        $strand=$line[6];
                        if(exists ${$refLen}{$chr}){
                                $c_end=$refLen->{$chr};
                        }else{
                                print "Warning: chromosome length unknown (chr:$chr).\n";
                                next;
                        }
                        if($line[8]=~/Name=(.+?);/){
						        $geneName=$1;
						}elsif($line[8]=~/Name=(.+)/){
						        $geneName=$1;
						}else{
                                $geneName="unknown";
                        }
                        
                        if($line[8]=~/ID=gene[:,-](.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=gene[:,-](.+)/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+)/){
						        $id=$1;
						}else{
                                $id="unknown";
                        }
                        
                        if($line[8]=~/Dbxref=GeneID:(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=GeneID:(.+)/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+)/){
						        $dbxref=$1;
						}else{
                                $dbxref="unknown";
                        }
                        $a_geneName="${id}|${geneName}|$dbxref";
                        $a_geneName=~tr/\,/_/;
						# a b ===== b a
						# =====a b=== b a
						if($line[3]>$line[4]){
						        $gStart=$line[4];
								$gEnd=$line[3];
						}else{
						        $gStart=$line[3];
								$gEnd=$line[4];
						}
						
                        ###########
                        @cds=();
                        $indicate=0;
                        $c_mRNA=0;                        
                }elsif($line[2] eq "mRNA"){
                        $c_mRNA++;
                }elsif($line[2] eq "CDS"){
                        if($c_mRNA>1){
                                next;
                        }
                        push(@cds,[$line[3],$line[4]]);
                        $indicate=1;
                }				
		}
        #
        if($indicate==1){
                if(exists ${$refLen}{$chr}){
                        $c_end=$refLen->{$chr};
                        my $n=scalar(@cds);
                        #if($n<1){print "Warning: no CDS.\n";}
                        if($n>=1){
                                if($strand eq "+"){
                                                
                                        $cdsStart=$cds[$n-1][0];
                                        $cdsEnd=$cds[$n-1][1];
                                }else{
                                        if($cds[0][0]>$cds[$n-1][0]){
                                                $cdsStart=$cds[$n-1][0];
                                                $cdsEnd=$cds[$n-1][1];
                                        }else{
                                                $cdsStart=$cds[0][0];
                                                $cdsEnd=$cds[0][1];
                                        }
                                }
                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                        foreach my $pos(keys %{$posHash{$chr}}){
                                                if($pos>=$gStart && $pos<=$gEnd){
                                                        print "Gene information: $chr $pos $a_geneName\n";
                                                        if($a_geneName eq "unknown|unknown|unknown"){
                                                                next;
                                                        }
                                                        
                                                        if(exists $r_hash{$a_geneName}){
                                                                print "Warning: gene name duplicate in the GFF file.\n";
                                                                print "The gene on chromosome $chr, location $gStart -- $gEnd will be skipped.\n"; 
                                                        }else{
                                                                my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                                if($tout[0] > 0){
                                                                        my $str=join ",",$a_geneName,$chr,@tout,$strand;
                                                                        $r_hash{$a_geneName}=1;
                                                                        push(@out,"$str");
                                                                }
                                                        }
                                                        
                                                }
                                                #
                                                if($c_mRNA>1){
                                                        print "$a_geneName mRNA : $c_mRNA\n";
                                                }
                                        }
                                        
                                }else{
                                        print "Gene information: $chr ${gStart}_$gEnd $a_geneName\n";
                                        print "Warning: region of CDS out of gene's.\n";
                                }
                        }
                }else{
                        print "Warning: chromosome length unknown (chr:$chr).\n";
                }
        }
        #
		close IN;
		if(scalar (@out)<1){
                die "GeneName cat't be found in the GFF file.\n";
		}else{
                return(@out);
        }
}

sub c_allInterval{
        my ($t_gff,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,$refLen)=@_;
        # read gff
        print "Reading GFF file\n";
		if($t_gff=~/\.gz$/){
		        (open IN,"gzip -dc $t_gff|") || die "$t_gff open failed. $!\n";
		}else{
		        (open IN,"$t_gff") || die "$t_gff open failed. $!\n";
		}
		
		my ($chr,$c_end,$geneName,$a_geneName);
		my ($e_up_a,$e_up_b,$e_down_a,$e_down_b);
        my ($v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len);
        my ($gStart,$gEnd,$strand);
        my ($cdsStart,$cdsEnd);
        my @cds;
        my @out;        
        my %r_hash;
        my $id;
        my $dbxref;
        my $indicate=0;
        my $c_mRNA=0;
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				#if($line[2] eq "gene" || $line[2] eq "ncRNA_gene" || $line[2] eq "pseudogene"){
                if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                        if($indicate==1){
                                if(exists $r_hash{$a_geneName}){
                                        print "Warning: gene name duplicate in the GFF file.\n";
                                        print "The gene on chromosome $chr, location $gStart -- $gEnd will be skipped.\n"; 
                                }else{
                                        my $n=scalar(@cds);
                                        if($n<1){
                                                print "Gene information: $chr ${gStart}_$gEnd $a_geneName\n";
                                                print "Warning: no CDS.\n";
                                        }else{
                                                if($strand eq "+"){
                                                        $cdsStart=$cds[$n-1][0];
                                                        $cdsEnd=$cds[$n-1][1];
                                                }else{
                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[$n-1][1];
                                                        }else{
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }
                                                }
                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                        my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                        if($tout[0] > 0){
                                                                my $str=join ",",$a_geneName,$chr,@tout,$strand;
                                                                $r_hash{$a_geneName}=1;
                                                                push(@out,"$str");
                                                        }
                                                        
                                                }else{
                                                        print "Gene information: $chr ${gStart}_$gEnd $a_geneName\n";
                                                        print "Warning: region of CDS out of gene's.\n";
                                                }
                                        }
                                        #
                                        if($c_mRNA>1){
                                                print "$a_geneName mRNA : $c_mRNA\n";
                                        }
                                }
                        }
                        ############################
                        $chr=$line[0];
                        $strand=$line[6];
                        if(exists ${$refLen}{$chr}){
                                $c_end=$refLen->{$chr};
                        }else{
                                print "Warning: chromosome length unknown (chr:$chr).\n";
                                next;
                        }
                        if($line[8]=~/Name=(.+?);/){
						        $geneName=$1;
						}elsif($line[8]=~/Name=(.+)/){
						        $geneName=$1;
						}else{
                                $geneName="unknown";
                        }
                        
                        if($line[8]=~/ID=gene[:,-](.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=gene[:,-](.+)/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+?);/){
						        $id=$1;
						}elsif($line[8]=~/ID=(.+)/){
						        $id=$1;
						}else{
                                $id="unknown";
                        }
                        
                        if($line[8]=~/Dbxref=GeneID:(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=GeneID:(.+)/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+?);/){
						        $dbxref=$1;
						}elsif($line[8]=~/Dbxref=(.+)/){
						        $dbxref=$1;
						}else{
                                $dbxref="unknown";
                        }
                        $a_geneName="${id}|${geneName}|$dbxref";
                        $a_geneName=~tr/\,/_/;
                        if($a_geneName eq "unknown|unknown|unknown"){
                                next;
                        }
						#
						if($line[3]>$line[4]){
						        $gStart=$line[4];
								$gEnd=$line[3];
						}else{
						        $gStart=$line[3];
								$gEnd=$line[4];
						}
                        ###########
                        @cds=();
                        $indicate=0;
                        $c_mRNA=0;
                }elsif($line[2] eq "mRNA"){
                        $c_mRNA++;
                }elsif($line[2] eq "CDS"){
                        if($c_mRNA>1){
                                next;
                        }
                        push(@cds,[$line[3],$line[4]]);
                        $indicate=1;
                }
		}
        #
        if($indicate==1){
                if(exists $r_hash{$a_geneName}){
                        print "Warning: gene name duplicate in the GFF file.\n";
                        print "The gene on chromosome $chr, location $gStart -- $gEnd will be skipped.\n"; 
                }else{
                        my $n=scalar(@cds);
                        if($n<1){
                                print "Gene information: $chr ${gStart}_$gEnd $a_geneName\n";
                                print "Warning: no CDS.\n";
                        }else{
                                if($strand eq "+"){
                                        $cdsStart=$cds[$n-1][0];
                                        $cdsEnd=$cds[$n-1][1];
                                }else{
                                        if($cds[0][0]>$cds[$n-1][0]){
                                                $cdsStart=$cds[$n-1][0];
                                                $cdsEnd=$cds[$n-1][1];
                                        }else{
                                                $cdsStart=$cds[0][0];
                                                $cdsEnd=$cds[0][1];
                                        }
                                }
                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                        my @tout=&c_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                        if($tout[0] > 0){
                                                my $str=join ",",$a_geneName,$chr,@tout,$strand;
                                                $r_hash{$a_geneName}=1;
                                                push(@out,"$str");
                                        }
                                        
                                }else{
                                        print "Gene information: $chr ${gStart}_$gEnd $a_geneName\n";
                                        print "Warning: region of CDS out of gene's.\n";
                                }
                                #
                                if($c_mRNA>1){
                                        print "$a_geneName mRNA : $c_mRNA\n";
                                }
                        }
                }
        }
        #
		close IN;
        if(scalar (@out)<1){
                die "GeneName cat't be found in the GFF file.\n";
		}else{
                return(@out);
        }
}

sub c_oneSeq{
        my ($t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$t_ref,$t_dir,$t_info,$salt,$dsalt,$dntp,$pCon)=@_;
        if(! -d "$t_dir/primer_para"){
                my $r=system("mkdir $t_dir/primer_para");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_para\n";
                }
        }
        if(! -d "$t_dir/primer_blast"){
                my $r=system("mkdir $t_dir/primer_blast");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_blast\n";
                }
        }
        #
        my ($a_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$cdsStart,$cdsEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=@$t_info;
        my ($t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len);
        if($strand eq "+"){
                $t_upR_start=$cdsEnd-$v_up_a-$t_PrimerMax-1;
                $t_upR_len=$t_PrimerMax;
                #$t_downF_start=$cdsEnd-$v_up_a-1;
                $t_downF_start=$cdsEnd-$v_up_a+2;
                $t_downF_len=$t_PrimerMax;
                
                $t_upF_start=$t_up_a-$v_up_a+1;
                $t_upF_len=$t_up_b-$t_up_a+1;
                $t_downR_start=$t_down_b-$v_up_a+1;
                $t_downR_len=$t_down_a-$t_down_b+1;
        }else{
                $t_upF_start=$cdsStart-$v_up_a+4;
                $t_upF_len=$t_PrimerMax;
                #$t_downR_start=$cdsStart-$v_up_a-$t_PrimerMax+4;
                $t_downR_start=$cdsStart-$v_up_a-$t_PrimerMax+1;
                $t_downR_len=$t_PrimerMax;
                        
                $t_upR_start=$t_down_b-$v_up_a+1;
                $t_upR_len=$t_down_a-$t_down_b+1;
                $t_downF_start=$t_up_a-$v_up_a+1;
                $t_downF_len=$t_up_b-$t_up_a+1;        
        }
        
        #
        my $v_up_start=1;
        my $v_up_len=$v_up_b-$v_up_a+1;
        #
        my $v_down_start=$v_down_b-$v_up_a+1;
        my $v_down_len=$v_down_a-$v_down_b+1;
        #
        #############
        my $seq=&chrSeq($t_ref,$chr);
        my $vseq=substr($seq,$v_up_a-1,$v_len);
        $vseq=~s/[^ATCGNatcgn]/N/g;
        #############
        my $para;
        (open OUT,">$t_dir/primer_para/$t_geneName.target.up.para") || die "$!\n";
        $para=&tUpPara($vseq,$a_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUT "$para";
        close OUT;
        
        (open OUG,">$t_dir/primer_para/$t_geneName.target.down.para") || die "$!\n";
        $para=&tDownPara($vseq,$a_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUG "$para";
        close OUG;
        #
        (open OUD,">$t_dir/primer_para/$t_geneName.verify.para") || die "$!\n";
        $para=&vpPara($vseq,$a_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_start,$v_up_len,$v_down_start,$v_down_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUD "$para";
        close OUD;
        #
}

sub c_listSeq{
        my ($t_dir,$t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_bund,$t_pNum,$t_mask,$t_ref,$geneInfo,$salt,$dsalt,$dntp,$pCon)=@_;       
        if(! -d "$t_dir/primer_para"){
                my $r=system("mkdir $t_dir/primer_para");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_para\n";
                }
        }
        if(! -d "$t_dir/primer_blast"){
                my $r=system("mkdir $t_dir/primer_blast");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_blast\n";
                }
        }
        #
        #my $fir=shift(@$geneInfo);
        my $fir=$geneInfo->[0];
        my ($t_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$cdsStart,$cdsEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=split /\,/,$fir;
        my ($t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len);
        if($strand eq "+"){
                $t_upR_start=$cdsEnd-1-$v_up_a-$t_PrimerMax;
                $t_upR_len=$t_PrimerMax;
                #$t_downF_start=$cdsEnd-$v_up_a-1;
                $t_downF_start=$cdsEnd-$v_up_a+2;
                $t_downF_len=$t_PrimerMax;
                
                $t_upF_start=$t_up_a-$v_up_a+1;
                $t_upF_len=$t_up_b-$t_up_a+1;
                $t_downR_start=$t_down_b-$v_up_a+1;
                $t_downR_len=$t_down_a-$t_down_b+1;
        }else{
                $t_upF_start=$cdsStart-$v_up_a+4;
                $t_upF_len=$t_PrimerMax;
                #$t_downR_start=$cdsStart-$v_up_a-$t_PrimerMax+4;
                $t_downR_start=$cdsStart-$v_up_a-$t_PrimerMax+1;
                $t_downR_len=$t_PrimerMax;
                        
                $t_upR_start=$t_down_b-$v_up_a+1;
                $t_upR_len=$t_down_a-$t_down_b+1;
                $t_downF_start=$t_up_a-$v_up_a+1;
                $t_downF_len=$t_up_b-$t_up_a+1;        
        }
        #
        my $v_up_start=1;
        my $v_up_len=$v_up_b-$v_up_a+1;
        #
        my $v_down_start=$v_down_b-$v_up_a+1;
        my $v_down_len=$v_down_a-$v_down_b+1;
        my $preChr=$chr;
        ##################################
        #############
        my $seq=&chrSeq($t_ref,$chr);
        my $vseq=substr($seq,$v_up_a-1,$v_len);
        $vseq=~s/[^ATCGNatcgn]/N/g;
        #############
        my $para;
        
        my $bundle=1;
        my $count=1;
        (open OUT,">$t_dir/primer_para/bundle.$bundle.target.up.para") || die "$!\n";
        (open OUG,">$t_dir/primer_para/bundle.$bundle.target.down.para") || die "$!\n";
        
        (open OUA,">$t_dir/primer_para/bundle.$bundle.verify.para") || die "$!\n";

        $para=&tUpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUT "$para";
        
        $para=&tDownPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUG "$para";  
        
        $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_start,$v_up_len,$v_down_start,$v_down_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUA "$para";
        #################################
        #foreach my $k(@$geneInfo){
        my $nInfo=scalar(@$geneInfo);
        for(my $k=1;$k<$nInfo;$k++){
                ($t_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$cdsStart,$cdsEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=split /\,/,$geneInfo->[$k];
                if($strand eq "+"){
                        $t_upR_start=$cdsEnd-1-$v_up_a-$t_PrimerMax;
                        $t_upR_len=$t_PrimerMax;
                        #$t_downF_start=$cdsEnd-$v_up_a-1;
                        $t_downF_start=$cdsEnd-$v_up_a+2;
                        $t_downF_len=$t_PrimerMax;
                        
                        $t_upF_start=$t_up_a-$v_up_a+1;
                        $t_upF_len=$t_up_b-$t_up_a+1;
                        $t_downR_start=$t_down_b-$v_up_a+1;
                        $t_downR_len=$t_down_a-$t_down_b+1;
                }else{
                        $t_upF_start=$cdsStart-$v_up_a+4;
                        $t_upF_len=$t_PrimerMax;
                        #$t_downR_start=$cdsStart-$v_up_a-$t_PrimerMax+4;
                        $t_downR_start=$cdsStart-$v_up_a-$t_PrimerMax+1;
                        $t_downR_len=$t_PrimerMax;
                                
                        $t_upR_start=$t_down_b-$v_up_a+1;
                        $t_upR_len=$t_down_a-$t_down_b+1;
                        $t_downF_start=$t_up_a-$v_up_a+1;
                        $t_downF_len=$t_up_b-$t_up_a+1;        
                }
                #
                $v_up_start=1;
                $v_up_len=$v_up_b-$v_up_a+1;
                #
                $v_down_start=$v_down_b-$v_up_a+1;
                $v_down_len=$v_down_a-$v_down_b+1;
                #
                if($preChr ne $chr){
                        $seq=&chrSeq($t_ref,$chr);
                }
                $vseq=substr($seq,$v_up_a-1,$v_len);
                $vseq=~s/[^ATCGNatcgn]/N/g;
                $count++;
                if($count>$t_bund){
                        close OUT;
                        close OUG;
                        close OUA;
                        $bundle++;
                        (open OUT,">$t_dir/primer_para/bundle.$bundle.target.up.para") || die "$!\n";
                        (open OUG,">$t_dir/primer_para/bundle.$bundle.target.down.para") || die "$!\n";
                        (open OUA,">$t_dir/primer_para/bundle.$bundle.verify.para") || die "$!\n";
                        $count=1;
                }
                #################################                
                $para=&tUpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUT "$para";
                
                $para=&tDownPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUG "$para";
                
                $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_start,$v_up_len,$v_down_start,$v_down_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUA "$para";
                #
                $preChr=$chr;
        }
        close OUT;
        close OUG;
        close OUA;
}

sub c_multiRun1{
        my ($t_dir,$s_geneName,$s_prog,$s_setFile,$t_thread)=@_;
        if(! -d "$t_dir/primer_raw"){
                my $r=system("mkdir $t_dir/primer_raw");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_raw\n";
                }
        }        
        #
        my $t_up_inFile="$t_dir/primer_para/${s_geneName}.target.up.para";
        my $t_down_inFile="$t_dir/primer_para/${s_geneName}.target.down.para";
        my $v_inFile="$t_dir/primer_para/${s_geneName}.verify.para";
        
        my $t_up_outFile="$t_dir/primer_raw/${s_geneName}.target.up.out";
        my $t_down_outFile="$t_dir/primer_raw/${s_geneName}.target.down.out";
        my $v_outFile="$t_dir/primer_raw/${s_geneName}.verify.out";
        
        my $t_up_errFile="$t_dir/primer_raw/${s_geneName}.target.up.err";
        my $t_down_errFile="$t_dir/primer_raw/${s_geneName}.target.up.err";
        my $v_errFile="$t_dir/primer_raw/${s_geneName}.verify.err";      
        
        my $task=0;
        my @allTask=([$t_up_outFile,$t_up_errFile,$t_up_inFile],[$t_down_outFile,$t_down_errFile,$t_down_inFile],[$v_outFile,$v_errFile,$v_inFile]);
        if($t_thread<2){
                &runPrimer($s_prog,$s_setFile,$t_up_outFile,$t_up_errFile,$t_up_inFile);
                &runPrimer($s_prog,$s_setFile,$t_down_outFile,$t_down_errFile,$t_down_inFile);
                &runPrimer($s_prog,$s_setFile,$v_outFile,$v_errFile,$v_inFile);
        }else{
                while(1){
                        last if($task==3);
                        while(scalar(threads->list(threads::all)) < $t_thread){
                                threads->create(\&runPrimer,$s_prog,$s_setFile,$allTask[$task][0],$allTask[$task][1],$allTask[$task][2]);
                                $task++;
                                last if($task==3); 
                        }
                        foreach my $k(threads->list(threads::all)){
                                if($k->is_joinable()){
                                        $k->join();
                                }
                        }
                }
                foreach my $k(threads->list(threads::all)){
                        $k->join();
                }
        }
}


sub c_multiRun2{
        my ($t_dir,$s_prog,$s_setFile,$t_thread,$n_bundle)=@_;
        if(! -d "$t_dir/primer_raw"){
                my $r=system("mkdir $t_dir/primer_raw");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_raw\n";
                }
        } 
        # bundle.1.target.para; bundle.1.verify.para        
        my $n=$n_bundle*3;
        if($t_thread<2){
                for(my $i=1;$i<=$n_bundle;$i++){
                        my $t_up_inFile="$t_dir/primer_para/bundle.$i.target.up.para";
                        my $t_down_inFile="$t_dir/primer_para/bundle.$i.target.down.para";
                        my $v_inFile="$t_dir/primer_para/bundle.$i.verify.para";
                        
                        my $t_up_outFile="$t_dir/primer_raw/bundle.$i.target.up.out";
                        my $t_down_outFile="$t_dir/primer_raw/bundle.$i.target.down.out";
                        my $v_outFile="$t_dir/primer_raw/bundle.$i.verify.out";
                        
                        my $t_up_errFile="$t_dir/primer_raw/bundle.$i.target.up.err";
                        my $t_down_errFile="$t_dir/primer_raw/bundle.$i.target.down.err";
                        my $v_errFile="$t_dir/primer_raw/bundle.$i.verify.err";
                        
                        &runPrimer($s_prog,$s_setFile,$t_up_outFile,$t_up_errFile,$t_up_inFile);
                        &runPrimer($s_prog,$s_setFile,$t_down_outFile,$t_down_errFile,$t_down_inFile);
                        &runPrimer($s_prog,$s_setFile,$v_outFile,$v_errFile,$v_inFile);
                }
        }else{
                my $j=0;
                while(1){
                        last if($j>=$n);
                        while(scalar(threads->list(threads::all))<$t_thread){
                                $j++;
                                my $i=int(($j+2)/3);
                                if($j % 3 == 1){
                                        my $t_inFile="$t_dir/primer_para/bundle.$i.target.up.para";
                                        my $t_outFile="$t_dir/primer_raw/bundle.$i.target.up.out";
                                        my $t_errFile="$t_dir/primer_raw/bundle.$i.target.up.err";
                                        threads->create(\&runPrimer,$s_prog,$s_setFile,$t_outFile,$t_errFile,$t_inFile);
                                        print "Bundle $i target primers (upstream) running\n";
                                }elsif($j % 3 == 2){
                                        my $t_inFile="$t_dir/primer_para/bundle.$i.target.down.para";
                                        my $t_outFile="$t_dir/primer_raw/bundle.$i.target.down.out";
                                        my $t_errFile="$t_dir/primer_raw/bundle.$i.target.down.err";
                                        threads->create(\&runPrimer,$s_prog,$s_setFile,$t_outFile,$t_errFile,$t_inFile);
                                        print "Bundle $i target primers (downstream) running\n";
                                }else{
                                        my $v_inFile="$t_dir/primer_para/bundle.$i.verify.para";
                                        my $v_outFile="$t_dir/primer_raw/bundle.$i.verify.out";
                                        my $v_errFile="$t_dir/primer_raw/bundle.$i.verify.err";
                                        threads->create(\&runPrimer,$s_prog,$s_setFile,$v_outFile,$v_errFile,$v_inFile);
                                        print "Bundle $i verification primers running\n";
                                        if($i==$n_bundle){last};
                                }
                        }
                        #sleep(5);
                        foreach my $k(threads->list(threads::all)){
                                if($k->is_joinable()){
                                        $k->join();
                                }
                        }
                }
                foreach my $k(threads->list(threads::all)){
                        $k->join();
                }
        }   
}

sub c_vpGeneFormat{
        my ($t_dir,$gName,$refine)=@_;
        my @out=("$gName.verify.out");
        my $left_fa="$gName.verify.left.fa";
        my $right_fa="$gName.verify.right.fa";
        my $formatOut="$gName.verify.format.out";
        my $sta="$gName.verify.stat";
        my $vforce=0;
        &formatPrimerMulti($t_dir,\@out,$left_fa,$right_fa,$formatOut,$sta,$vforce,$refine);
}

sub c_vpFormat{
        my ($t_dir,$bundle,$refine)=@_;
        my @out=map {"bundle.${_}.verify.out"} (1..$bundle);       
        my $left_fa="bundle.verify.left.fa";
        my $right_fa="bundle.verify.right.fa";
        my $formatOut="bundle.verify.format.out";
        my $sta="bundle.verify.stat";
        my $vforce=0;
        &formatPrimerMulti($t_dir,\@out,$left_fa,$right_fa,$formatOut,$sta,$vforce,$refine);
}

sub c_primer2blast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lpFa,$rpFa,$lvFa,$rvFa,$thread)=@_;
        &checkBlastdb($t_ref,$makeblastdb);
        my $r=0;              
        $r=system("$blastn  -task blastn-short -query $lmixFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lmix.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lmixFa\n";
        }
        
        $r=system("$blastn  -task blastn-short -query $rmixFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rmix.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rmixFa\n";
        }
               
        $r=system("$blastn  -task blastn-short -query $lpFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lp.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lpFa\n";
        }
        
        $r=system("$blastn  -task blastn-short -query $rpFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rp.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rpFa\n";
        }
        #
        $r=system("$blastn  -task blastn-short -query $lvFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lv.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lvFa\n";
        }
        
        $r=system("$blastn  -task blastn-short -query $rvFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rv.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rvFa\n";
        }
}

sub c_oneBlast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$thread,$geneName)=@_;
        my $dir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $x_ref="$dir/$refName";
        
        my $lmixFa="$t_dir/primer_blast/$geneName.target.up.reverse.fa";
        my $rmixFa="$t_dir/primer_blast/$geneName.target.down.forward.fa";
        my $lpFa="$t_dir/primer_blast/$geneName.target.up.forward.fa";
        my $rpFa="$t_dir/primer_blast/$geneName.target.down.reverse.fa";
        my $lvFa="$t_dir/primer_blast/$geneName.verify.left.fa";
        my $rvFa="$t_dir/primer_blast/$geneName.verify.right.fa";
        &c_primer2blast($t_dir,$x_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lpFa,$rpFa,$lvFa,$rvFa,$thread);
}

sub c_bundleBlast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$thread)=@_;
        my $dir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $x_ref="$dir/$refName";
        
        my $lpFa="$t_dir/primer_blast/bundle.target.up.forward.fa";
        my $lmixFa="$t_dir/primer_blast/bundle.target.up.reverse.fa";
        my $rmixFa="$t_dir/primer_blast/bundle.target.down.forward.fa";
        my $rpFa="$t_dir/primer_blast/bundle.target.down.reverse.fa";
        my $lvFa="$t_dir/primer_blast/bundle.verify.left.fa";
        my $rvFa="$t_dir/primer_blast/bundle.verify.right.fa";
        &c_primer2blast($t_dir,$x_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lpFa,$rpFa,$lvFa,$rvFa,$thread);
}

sub c_vpFilt{
        my ($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen)=@_;
        my $fblast="$t_dir/primer_blast/blastn.lv.out";
        my $rblast="$t_dir/primer_blast/blastn.rv.out";
        my $outFile="$t_dir/primer_blast/verify.length";
        &vPairFilt($t_dir,$fblast,$rblast,$outFile,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
}

sub c_vpMergeInfo{
        my ($t_dir,$vForm,$refine)=@_;        
        my $r=0;
        if(! -d "$t_dir/primer_result"){
                $r=system("mkdir $t_dir/primer_result");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_result\n";
                }
        }
        
        my %tg_hash;
        (open IN,"$t_dir/primer_result/target.primer.txt") || die "$!\n";
        <IN>;
        while(<IN>){
                chomp;
                my ($id,$geneName,$gid,$chr,$gStart,$gEnd,$v_up_a)=(split /\s+/,$_)[0,1,2,4,5,6,7];
                my $tag="${id}|${geneName}|$gid";
                if(! exists $tg_hash{$tag}){
                        $tg_hash{$tag}=[$chr,$gStart,$gEnd,$v_up_a];
                }
        }        
        close IN;
        #########################################
        (open OUA,">$t_dir/primer_result/verify.primer.txt") || die "$!\n";
        print OUA "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tPrimer_ID\tV1\tV2\tV1_Loc\tV2_Loc\tV1_Tm\tV2_Tm\tV1_V2_Products_Counts\tV1_V2_Blast_Len\tV1_V2_Expect_Len\tV1_V2_Penalty\tTier\n";
        (open AN,"$vForm") || die "$!\n";
        (open BN,"$t_dir/primer_blast/verify.length") || die "$!\n";
        <AN>;<BN>;

        my %vpStat;
        my $tier = 0;
        my $wstat;
        my $num = 0;
        my $pre_gene = "";
        my @upout = ();
        while(<AN>){
                chomp;
                my @line=split;
                my ($geneName,$strand,$penalty,$verify_f,$start_f,$len_f,$tm_f,$verify_r,$start_r,$len_r,$tm_r)=@line[0,1,3,4,5,6,7,13,14,15,16];
                if($geneName ne $pre_gene){
                        if(scalar(@upout)>0){
                                my @upSort;
                                if($refine){
                                       @upSort = sort {$a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] or $a->[2] <=> $b->[2]} @upout;
                                }else{
                                       @upSort = sort {$a->[3] <=> $b->[3] or $a->[2] <=> $b->[2]} @upout;
                                }
                                
                                $num = 0;
                                foreach my $k(@upSort){
                                        $num++;
                                        print OUA "$k->[0]\t$num\t$k->[1]\t$k->[2]\t$k->[3]\n";
                                }
                        }
                        @upout = ();
                        
                }
                $pre_gene = $geneName;
                
                my $m_line=<BN>;
                chomp($m_line);
                my @mLine=split /\s+/,$m_line;
                my $mistake_start=$mLine[4];
                my $pcr=$mLine[3];
                if($pcr eq "NA"){
                        next;
                }
                
                if(exists $tg_hash{$geneName}){
                        my $x_geneName=$geneName;
                        $x_geneName=~s/\|/\t/g;
                        
                        my $p_len=$start_r-$start_f+1;
                        $wstat = &pFlutter($pcr,$p_len);
                        
                        my ($chr,$gStart,$gEnd,$v_up_a)=@{$tg_hash{$geneName}};
                        my ($v1_loc,$v2_loc)=&vloc($v_up_a,$start_f,$start_r,$len_f,$len_r);
                        $tier = 0;
                        my $n_mis = ($pcr =~ tr/;//);
                        $n_mis += 1;
                        if($n_mis == 1){
                                $tier = 1;
                        }else{
                                if($wstat == 1){
                                        $tier = 2;
                                }
                        }
                        if($tier > 0){
                                my $taBase = 0;
                                if($refine){
                                        $taBase = &vEndBase($verify_f,$verify_r);
                                }
                                if($strand eq "plus"){
                                        #print OUA "$x_geneName\t+\t$chr\t$gStart\t$gEnd\t$num\t$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$mistake_start\t$pcr\t$p_len\t$penalty\t$tier\n";
                                        push(@upout,["$x_geneName\t+\t$chr\t$gStart\t$gEnd","$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$mistake_start\t$pcr\t$p_len",$penalty,$tier,$taBase]);
                                }else{
                                        #print OUA "$x_geneName\t-\t$chr\t$gStart\t$gEnd\t$num\t$verify_r\t$verify_f\t$v2_loc\t$v1_loc\t$tm_r\t$tm_f\t$mistake_start\t$pcr\t$p_len\t$penalty\t$tier\n";
                                        push(@upout,["$x_geneName\t-\t$chr\t$gStart\t$gEnd","$verify_r\t$verify_f\t$v2_loc\t$v1_loc\t$tm_r\t$tm_f\t$mistake_start\t$pcr\t$p_len",$penalty,$tier,$taBase]);
                                }
                                $vpStat{$geneName} = 1;
                        }   
                }
        }
        close AN;
        close BN;
        if(scalar(@upout)>0){
                my @upSort;
                if($refine){
                       @upSort = sort {$a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] or $a->[2] <=> $b->[2]} @upout;
                }else{
                       @upSort = sort {$a->[3] <=> $b->[3] or $a->[2] <=> $b->[2]} @upout;
                }
                
                $num = 0;
                foreach my $k(@upSort){
                        $num++;
                        print OUA "$k->[0]\t$num\t$k->[1]\t$k->[2]\t$k->[3]\n";
                }
        }
        @upout = ();
        #######################
        print "###################################\n";
        print "Verification primers not exist: \n"; 
        foreach my $k(keys %tg_hash){
                if(! exists($vpStat{$k})){
                        print "$k\n";
                }
        }
        print "###################################\n";
        #######################
        close OUA;
}

sub c_vpStrandTrans{
        my ($t_dir,$vForm,$refine)=@_;
        my $r=0;
        if(! -d "$t_dir/primer_result"){
                $r=system("mkdir $t_dir/primer_result");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_result\n";
                }
        }
        
        ####################################
        my %tg_hash;
        (open IN,"$t_dir/primer_result/target.primer.unblast.txt") || die "$!\n";
        <IN>;
        while(<IN>){
                chomp;
                my ($id,$geneName,$gid,$chr,$gStart,$gEnd,$v_up_a)=(split /\s+/,$_)[0,1,2,4,5,6,7];
                my $tag="${id}|${geneName}|$gid";
                if(! exists $tg_hash{$tag}){
                        $tg_hash{$tag}=[$chr,$gStart,$gEnd,$v_up_a];
                }
        }        
        close IN;
        ####################################
        (open OUA,">$t_dir/primer_result/verify.primer.unblast.txt") || die "$!\n";
        print OUA "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tPrimer_ID\tV1\tV2\tV1_Loc\tV2_Loc\tV1_Tm\tV2_Tm\tV1_V2_Expect_Len\tPenalty\n";
        
        (open AN,"$vForm") || die "$!\n";
        <AN>;
        
        my %vpStat;
        my $num=0;
        my $pre_gene="";
        my @upout = ();
        while(<AN>){
                chomp;
                my @line=split;
                my ($geneName,$strand,$penalty,$verify_f,$start_f,$len_f,$tm_f,$verify_r,$start_r,$len_r,$tm_r)=@line[0,1,3,4,5,6,7,13,14,15,16];
                if($geneName ne $pre_gene){
                        if(scalar(@upout)>0){
                                my @upSort;
                                if($refine){
                                       @upSort = sort {$a->[3] <=> $b->[3] or $a->[2] <=> $b->[2]} @upout;
                                }else{
                                       @upSort = sort {$a->[2] <=> $b->[2]} @upout;
                                }
                                
                                $num = 0;
                                foreach my $k(@upSort){
                                        $num++;
                                        print OUA "$k->[0]\t$num\t$k->[1]\t$k->[2]\n";
                                }
                        }
                        @upout = ();
                }
                $pre_gene = $geneName;
                
                if(exists $tg_hash{$geneName}){
                        my $x_geneName=$geneName;
                        $x_geneName=~s/\|/\t/g;
                        
                        my $p_len=$start_r-$start_f+1;
                        my ($chr,$gStart,$gEnd,$v_up_a)=@{$tg_hash{$geneName}};
                        my ($v1_loc,$v2_loc)=&vloc($v_up_a,$start_f,$start_r,$len_f,$len_r);

                        my $taBase = 0;
                        if($refine){
                                $taBase = &vEndBase($verify_f,$verify_r);
                        }                        
                        if($strand eq "plus"){
                                #print OUA "$x_geneName\t+\t$chr\t$gStart\t$gEnd\t$num\t$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$p_len\t$penalty\n";
                                push(@upout,["$x_geneName\t+\t$chr\t$gStart\t$gEnd","$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$p_len",$penalty,$taBase]);
                        }else{
                                #print OUA "$x_geneName\t-\t$chr\t$gStart\t$gEnd\t$num\t$verify_r\t$verify_f\t$v2_loc\t$v1_loc\t$tm_r\t$tm_f\t$p_len\t$penalty\n";
                                push(@upout,["$x_geneName\t-\t$chr\t$gStart\t$gEnd","$verify_r\t$verify_f\t$v2_loc\t$v1_loc\t$tm_r\t$tm_f\t$p_len",$penalty,$taBase]);
                        }
                        $vpStat{$geneName} = 1;
                }
        }
        
        if(scalar(@upout)>0){
                my @upSort;
                if($refine){
                       @upSort = sort {$a->[3] <=> $b->[3] or $a->[2] <=> $b->[2]} @upout;
                }else{
                       @upSort = sort {$a->[2] <=> $b->[2]} @upout;
                }
                
                $num = 0;
                foreach my $k(@upSort){
                        $num++;
                        print OUA "$k->[0]\t$num\t$k->[1]\t$k->[2]\n";
                }
        }
        @upout = ();
                                                
        print "###################################\n";
        print "Verification primers not exist: \n"; 
        foreach my $k(keys %tg_hash){
                if(! exists($vpStat{$k})){
                        print "$k\n";
                }
        }
        print "###################################\n";
        #######################
        close AN;
        close OUA;
}

sub c_vpGeneTrans{        
        my ($t_dir,$geneName,$refine)=@_;
        my $vForm="$t_dir/primer_blast/$geneName.verify.format.out";
        &c_vpStrandTrans($t_dir,$vForm,$refine);
}

sub c_vpGeneMerge{        
        my ($t_dir,$geneName,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$refine)=@_;
        my $vForm="$t_dir/primer_blast/$geneName.verify.format.out";
        &c_vpFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        &c_vpMergeInfo($t_dir,$vForm,$refine);
}

sub c_vpBundleTrans{
        my ($t_dir,$refine)=@_;
        my $vForm="$t_dir/primer_blast/bundle.verify.format.out";
        &c_vpStrandTrans($t_dir,$vForm,$refine);
}

sub c_vpBundleMerge{
        my ($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$dERegion,$dEMatch,$refine)=@_;
        my $vForm="$t_dir/primer_blast/bundle.verify.format.out";
        &c_vpFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        &c_vpMergeInfo($t_dir,$vForm,$refine);
}

sub ctag_main{
        # common
        my ($binDir,$config,$t_ref,$t_gff,$plasmidFa,$t_dir,$t_geneName,$t_coord,$t_geneList,$t_coordList,$all,$t_mask,$t_pNum,$t_bundle_size,$t_thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$force,$rec,$refine,$thermo)=@_;
        my ($t_primer,$t_blastn,$t_makeblastdb,$t_windowmasker,$primerPara,$ntthal)=&getProg($binDir,$thermo);
        
        my ($r_primer,$r_ntthal,$r_blastn,$r_makeblastdb,$r_windowmasker,$t_upVfar,$t_upVnear,
        $t_upGeneFar,$t_upGeneNear,$t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_upMixPlasmidSize,
        $t_downVfar,$t_downVnear,$t_downGeneFar,$t_downGeneNear,$t_downMixGenomeSize,$t_downMixPlasmidSize,
        $t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,
        $d_GminMinus,$d_endRegion,$d_endMatch,$settings)=&read_config($config);
        
        if($settings ne "*"){
                $primerPara=$settings;
        }
        my $d_max_tm=&dMaxTm($primerPara,$d_GminMinus);
        
        if($r_primer ne "*"){
                $t_primer=$r_primer;
        }
        if(! -r $t_primer){
                print "Program: primer3_core\n";
                die "Error: invalid program. $t_primer\n";
        }
        if($r_ntthal ne "*"){
                $ntthal=$r_ntthal;
        }
        if(! -r $ntthal){
                print "Program: ntthal\n";
                die "Error: invalid program. $ntthal\n";
        }
        if($r_blastn ne "*"){
                $t_blastn=$r_blastn;
        }
        if(! -r $t_blastn){
                print "Program: blastn\n";
                die "Error: invalid program. $t_blastn\n";
        }
        if($r_makeblastdb ne "*"){
                $t_makeblastdb=$r_makeblastdb;
        }
        if(! -r $t_makeblastdb){
                print "Program: makeblastdb\n";
                die "Error: invalid program. $t_makeblastdb\n";
        }
        if($r_windowmasker ne "*"){
                $t_windowmasker=$r_windowmasker;
        }
        if(! -r $t_windowmasker){
                print "Program: windowmasker\n";
                die "Error: invalid program. $t_windowmasker\n";
        }
        
        my $rename_ref;
        $rename_ref=&unzRef($t_ref);
        my %refHash=&getRefLen($t_gff,$rename_ref);
        if($t_mask){
                $rename_ref=&maskGenome($t_windowmasker,$t_ref);
        }
        my @plasid_m_seq;        
        my ($upSeq,$downSeq,$seqLen);
        my ($p1_tm,$p2_tm);
        if($rec){
                @plasid_m_seq=&pRedesign($t_dir,$t_primer,$primerPara,$plasmidFa,$t_PrimerMin,$t_PrimerMax,$salt,$dsalt,$dntp,$pCon);
                ($p1_tm,$p2_tm)=&getPtm(\@plasid_m_seq,$pCon,$salt,$dsalt,$dntp,$method);
                
                my $p_mask=0;
                &plasmidPrimer($binDir,$config,$plasmidFa,$t_dir,$p_mask,$t_pNum,$t_thread,$salt,$dsalt,$dntp,$pCon,$thermo,$refine);
        }else{
                ($upSeq,$downSeq,$seqLen)=&plasmid2primer($t_dir,$plasmidFa,$t_upMixPlasmidSize,$t_downMixPlasmidSize);
                @plasid_m_seq=($upSeq,$downSeq);
                ($p1_tm,$p2_tm)=&getPtm(\@plasid_m_seq,$pCon,$salt,$dsalt,$dntp,$method);
                &plasmidCheck($t_dir,$upSeq,$downSeq,$seqLen,$t_upMixPlasmidSize,$t_downMixPlasmidSize,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$p1_tm,$p2_tm);
        }
        #my @plasid_m_plus=&plasmidPlus($plasmidFa,$t_upMixPlasmidSize,$t_downMixPlasmidSize);
        #my @plasid_m_minus=&plasmidMinus($plasmidFa,$t_upMixPlasmidSize,$t_downMixPlasmidSize);
        #&plasmid_primer($t_dir,$plasmidFa,$t_upMixPlasmidSize,$t_downMixPlasmidSize);
        #################
        
        my $vforce=0;
        
        # geneName
        if($t_geneName || $t_coord){
                my @geneInfo;
                if($t_geneName){
                        @geneInfo=&c_geneNameInterval($t_gff,$t_geneName,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coord){
                        my $loc;
                        ($loc=$t_coord)=~s/:/_/g;
                        $t_geneName="LOC_$loc";
                        @geneInfo=&c_coordInterval($t_gff,$t_coord,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }

                &c_oneSeq($t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$rename_ref,$t_dir,\@geneInfo,$salt,$dsalt,$dntp,$pCon);

                &c_multiRun1($t_dir,$t_geneName,$t_primer,$primerPara,$t_thread);
                
                &tpUpGeneFormat($t_dir,$t_geneName,$force,$refine);
                &tpDownGeneFormat($t_dir,$t_geneName,$force,$refine);

                &c_vpGeneFormat($t_dir,$t_geneName,$refine);
                
                my @info;
                push(@info,join ",",@geneInfo);
                if($run_blast){
                        &c_oneBlast($t_dir,$rename_ref,$t_blastn,$t_makeblastdb,$t_thread,$t_geneName);
                        #
                        &tpGeneMerge($t_dir,$t_geneName,\@plasid_m_seq,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$d_endRegion,$d_endMatch,\@info,$t_pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                        &c_vpGeneMerge($t_dir,$t_geneName,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$refine);
                }else{
                        &tpGeneTrans($t_dir,$t_geneName,\@plasid_m_seq,\@info,$d_endRegion,$d_endMatch,$t_pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                        &c_vpGeneTrans($t_dir,$t_geneName,$refine);
                }
                return(0);
        }
        ####################
        # bundle
        if($t_geneList || $t_coordList || $all){
                my @geneInfo;
                my $bundle_counts;
                
                if($t_geneList){
                        @geneInfo=&c_nameListInterval($t_gff,$t_geneList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coordList){
                        @geneInfo=&c_coordListInterval($t_gff,$t_coordList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($all){
                        @geneInfo=&c_allInterval($t_gff,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }
                &c_listSeq($t_dir,$t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_bundle_size,$t_pNum,$t_mask,$rename_ref,\@geneInfo,$salt,$dsalt,$dntp,$pCon);
                $bundle_counts=&bundleCounts($t_dir);
                &c_multiRun2($t_dir,$t_primer,$primerPara,$t_thread,$bundle_counts);
                &tpUpFormat($t_dir,$bundle_counts,$force,$refine);
                &tpDownFormat($t_dir,$bundle_counts,$force,$refine);
                &c_vpFormat($t_dir,$bundle_counts,$refine);
                if($run_blast){
                        &c_bundleBlast($t_dir,$rename_ref,$t_blastn,$t_makeblastdb,$t_thread);
                        #
                        &tpBundleMerge($t_dir,\@plasid_m_seq,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$d_endRegion,$d_endMatch,\@geneInfo,$t_pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                        &c_vpBundleMerge($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$d_endRegion,$d_endMatch,$refine);
                }else{
                        &tpBundleTrans($t_dir,\@plasid_m_seq,\@geneInfo,$d_endRegion,$d_endMatch,$t_pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                        &c_vpBundleTrans($t_dir,$refine);
                }
                return(0);
        }
        
        die "Error: lack of parameter (--geneName or --geneList or coordinate or coordList or --all)\n";
}

########################
#       short arm
########################

sub c_geneSeq{
        my ($seq,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand)=@_;
        my $vseq=substr($seq,$v_up_a-1,$v_len);
        $vseq=~s/[^ATCGNatcgn]/N/g;
        #my $t_gene_upseq=substr($$seq,$geneStart-$t_upMixGenomeSize-1,$t_upMixGenomeSize);
        my ($t_gene_upseq,$gene_mix_upseq,$gene_mix_downseq);
        if($strand eq "+"){
                $t_gene_upseq=substr($vseq,$cdsEnd-2-$v_up_a-$t_upMixGenomeSize,$t_upMixGenomeSize);
                #$gene_mix_downseq=substr($vseq,$cdsEnd-$v_up_a-2,$t_downMixGenomeSize);
                $gene_mix_downseq=substr($vseq,$cdsEnd-$v_up_a+1,$t_downMixGenomeSize);
        }else{
                #$t_gene_upseq=substr($vseq,$cdsStart-$v_up_a-$t_downMixGenomeSize+3,$t_downMixGenomeSize);
                $t_gene_upseq=substr($vseq,$cdsStart-$v_up_a-$t_downMixGenomeSize,$t_downMixGenomeSize);
                $gene_mix_downseq=substr($vseq,$cdsStart-$v_up_a+3,$t_upMixGenomeSize);
        }
        $gene_mix_upseq=reverse($t_gene_upseq);
        #$gene_mix_upseq=~s/[^ATCGNatcgn]/N/g;
        $gene_mix_upseq=~tr/ATCGatcg/TAGCtagc/;
        #
        #$gene_mix_downseq=~s/[^ATCGNatcgn]/N/g;
        return($vseq,$gene_mix_upseq,$gene_mix_downseq);
}

sub c_geneSeq2{
        my ($t_ref,$chr,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand)=@_;
        my $seq=&chrSeq($t_ref,$chr);
        &c_geneSeq($seq,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand);
}

sub sc_oneSeq{
        my ($t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$t_ref,$t_dir,$t_info,$salt,$dsalt,$dntp,$pCon)=@_;
        if(! -d "$t_dir/primer_para"){
                my $r=system("mkdir $t_dir/primer_para");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_para\n";
                }
        }
        if(! -d "$t_dir/primer_blast"){
                my $r=system("mkdir $t_dir/primer_blast");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_blast\n";
                }
        }
        #
        my ($a_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$cdsStart,$cdsEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=@$t_info; 
        #
        my $v_up_start=1;
        my $v_up_len=$v_up_b-$v_up_a+1;
        #
        my $v_down_start=$v_down_b-$v_up_a+1;
        my $v_down_len=$v_down_a-$v_down_b+1;
        #
        my ($vseq,$mix_upseq,$mix_downseq)=&c_geneSeq2($t_ref,$chr,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand);
        my $upN=($mix_upseq=~tr/Nn//);
        my $downN=($mix_downseq=~tr/Nn//);
        if($upN>5 || $downN >5){
                die "flanking sequence of gene has too many Ns\n";
        }
        my $para;
        #
        (open OUD,">$t_dir/primer_para/$t_geneName.verify.para") || die "$!\n";
        $para=&vpPara($vseq,$a_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_start,$v_up_len,$v_down_start,$v_down_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUD "$para";
        close OUD;
        #
        (open OUB,">$t_dir/primer_blast/$t_geneName.lmix.fa") || die "$!\n";
        (open OUC,">$t_dir/primer_blast/$t_geneName.rmix.fa") || die "$!\n";
        print OUB ">${a_geneName}_lmix\n";
        print OUB "$mix_upseq\n";
        print OUC ">${a_geneName}_rmix\n";
        print OUC "$mix_downseq\n";
        close OUB;
        close OUC;
}

sub sc_listSeq{
        my ($t_dir,$t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_bund,$t_pNum,$t_mask,$t_ref,$geneInfo,$salt,$dsalt,$dntp,$pCon)=@_;       
        if(! -d "$t_dir/primer_para"){
                my $r=system("mkdir $t_dir/primer_para");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_para\n";
                }
        }
        if(! -d "$t_dir/primer_blast"){
                my $r=system("mkdir $t_dir/primer_blast");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_blast\n";
                }
        }
        (open OUB,">$t_dir/primer_blast/bundle.lmix.fa") || die "$!\n";
        (open OUC,">$t_dir/primer_blast/bundle.rmix.fa") || die "$!\n";
        #
        #my $fir=shift(@$geneInfo);
        my $fir=$geneInfo->[0];
        my ($t_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$cdsStart,$cdsEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=split /\,/,$fir;
        #
        my $v_up_start=1;
        my $v_up_len=$v_up_b-$v_up_a+1;
        #
        my $v_down_start=$v_down_b-$v_up_a+1;
        my $v_down_len=$v_down_a-$v_down_b+1;
        my $preChr=$chr;
        ##################################
        my $seq=&chrSeq($t_ref,$chr);
        my ($vseq,$mix_upseq,$mix_downseq)=&c_geneSeq($seq,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand);
        my $para;
        
        my $bundle=1;
        my $count=1;
        
        (open OUA,">$t_dir/primer_para/bundle.$bundle.verify.para") || die "$!\n";
        $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_start,$v_up_len,$v_down_start,$v_down_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUA "$para";
        
        #
        print OUB ">${t_geneName}_lmix\n";
        print OUB "$mix_upseq\n";
        print OUC ">${t_geneName}_rmix\n";
        print OUC "$mix_downseq\n";
        #################################
        #foreach my $k(@$geneInfo){
        my $nInfo=scalar(@$geneInfo);
        for(my $k=1;$k<$nInfo;$k++){
                ($t_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$cdsStart,$cdsEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=split /\,/,$geneInfo->[$k];
                #
                $v_up_start=1;
                $v_up_len=$v_up_b-$v_up_a+1;
                #
                $v_down_start=$v_down_b-$v_up_a+1;
                $v_down_len=$v_down_a-$v_down_b+1;
                #
                if($preChr ne $chr){
                        $seq=&chrSeq($t_ref,$chr);
                }
                ($vseq,$mix_upseq,$mix_downseq)=&c_geneSeq($seq,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand);
                my $upN=($mix_upseq=~tr/Nn//);
                my $downN=($mix_downseq=~tr/Nn//);
                if($upN>5 || $downN >5){
                        next;
                }
                $count++;
                if($count>$t_bund){
                        close OUA;
                        $bundle++;
                        (open OUA,">$t_dir/primer_para/bundle.$bundle.verify.para") || die "$!\n";
                        $count=1;
                }
                #################################                
                
                $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_start,$v_up_len,$v_down_start,$v_down_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUA "$para";
                #
                print OUB ">${t_geneName}_lmix\n";
                print OUB "$mix_upseq\n";
                print OUC ">${t_geneName}_rmix\n";
                print OUC "$mix_downseq\n";
                #
                $preChr=$chr;
        }
        close OUA;
        close OUB;
        close OUC;
}

sub sc_multiRun1{
        my ($t_dir,$s_geneName,$s_prog,$s_setFile,$t_thread)=@_;
        if(! -d "$t_dir/primer_raw"){
                my $r=system("mkdir $t_dir/primer_raw");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_raw\n";
                }
        }        
        #
        my $v_inFile="$t_dir/primer_para/${s_geneName}.verify.para";
        
        my $v_outFile="$t_dir/primer_raw/${s_geneName}.verify.out";
        
        my $v_errFile="$t_dir/primer_raw/${s_geneName}.verify.err";      
        &runPrimer($s_prog,$s_setFile,$v_outFile,$v_errFile,$v_inFile);
}


sub sc_oneBlast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$thread,$geneName)=@_;
        my $dir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $x_ref="$dir/$refName";
        
        my $lmixFa="$t_dir/primer_blast/$geneName.lmix.fa";
        my $rmixFa="$t_dir/primer_blast/$geneName.rmix.fa";
        my $lvFa="$t_dir/primer_blast/$geneName.verify.left.fa";
        my $rvFa="$t_dir/primer_blast/$geneName.verify.right.fa";
        &sc_primer2blast($t_dir,$x_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lvFa,$rvFa,$thread);
}

sub sc_multiRun2{
        my ($t_dir,$s_prog,$s_setFile,$t_thread,$n_bundle)=@_;
        if(! -d "$t_dir/primer_raw"){
                my $r=system("mkdir $t_dir/primer_raw");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_raw\n";
                }
        } 
        # bundle.1.target.para; bundle.1.verify.para        
        my $n=$n_bundle;
        if($t_thread<2){
                for(my $i=1;$i<=$n_bundle;$i++){
                        my $v_inFile="$t_dir/primer_para/bundle.$i.verify.para";
                        
                        my $v_outFile="$t_dir/primer_raw/bundle.$i.verify.out";
                        
                        my $v_errFile="$t_dir/primer_raw/bundle.$i.verify.err";
                        
                        &runPrimer($s_prog,$s_setFile,$v_outFile,$v_errFile,$v_inFile);
                }
        }else{
                my $j=0;
                while(1){
                        last if($j>=$n);
                        while(scalar(threads->list(threads::all))<$t_thread){
                                $j++;
                                my $v_inFile="$t_dir/primer_para/bundle.$j.verify.para";
                                my $v_outFile="$t_dir/primer_raw/bundle.$j.verify.out";
                                my $v_errFile="$t_dir/primer_raw/bundle.$j.verify.err";
                                threads->create(\&runPrimer,$s_prog,$s_setFile,$v_outFile,$v_errFile,$v_inFile);
                                print "Bundle $j verification primers running\n";
                                if($j==$n_bundle){last};
                        }
                        #sleep(5);
                        foreach my $k(threads->list(threads::all)){
                                if($k->is_joinable()){
                                        $k->join();
                                }
                        }
                }
                foreach my $k(threads->list(threads::all)){
                        $k->join();
                }
        }   
}

sub sc_primer2blast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lvFa,$rvFa,$thread)=@_;
        &checkBlastdb($t_ref,$makeblastdb);
        my $r=0; 
        if(length($lmixFa)<50){
                $r=system("$blastn  -task blastn-short -query $lmixFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lmix.out");
        }else{
                $r=system("$blastn  -task blastn -query $lmixFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lmix.out");
        }
        if($r){
                print "Error: blastn error\n";
                die "query -- $lmixFa\n";
        }
        
        if(length($rmixFa)<50){
                $r=system("$blastn  -task blastn-short -query $rmixFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rmix.out");
        }else{
                $r=system("$blastn  -task blastn -query $rmixFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rmix.out");
        }
        if($r){
                print "Error: blastn error\n";
                die "query -- $rmixFa\n";
        }
        #
        $r=system("$blastn  -task blastn-short -query $lvFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lv.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lvFa\n";
        }
        
        $r=system("$blastn  -task blastn-short -query $rvFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rv.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rvFa\n";
        }
}

sub sc_bundleBlast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$thread)=@_;
        my $dir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $x_ref="$dir/$refName";
        
        my $lmixFa="$t_dir/primer_blast/bundle.lmix.fa";
        my $rmixFa="$t_dir/primer_blast/bundle.rmix.fa";
        my $lvFa="$t_dir/primer_blast/bundle.verify.left.fa";
        my $rvFa="$t_dir/primer_blast/bundle.verify.right.fa";
        &sc_primer2blast($t_dir,$x_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lvFa,$rvFa,$thread);
}

sub short_ctag_main{
        # common
        my ($binDir,$config,$t_ref,$t_gff,$plasmidFa,$t_dir,$t_geneName,$t_coord,$t_geneList,$t_coordList,$all,$t_mask,$t_pNum,$t_bundle_size,$t_thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$rec,$refine,$thermo)=@_;
        my ($t_primer,$t_blastn,$t_makeblastdb,$t_windowmasker,$primerPara,$ntthal)=&getProg($binDir,$thermo);
        my ($r_primer,$r_ntthal,$r_blastn,$r_makeblastdb,$r_windowmasker,$t_upVfar,$t_upVnear,
        $t_upGeneFar,$t_upGeneNear,$t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_upMixPlasmidSize,
        $t_downVfar,$t_downVnear,$t_downGeneFar,$t_downGeneNear,$t_downMixGenomeSize,$t_downMixPlasmidSize,
        $t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,
        $d_GminMinus,$d_endRegion,$d_endMatch,$settings)=&read_config($config);
        if($t_upMixGenomeSize<40){
                $t_upMixGenomeSize=40;
        }
        if($t_downMixGenomeSize<40){
                $t_downMixGenomeSize=40;
        }
        
        if($settings ne "*"){
                $primerPara=$settings;
        }        
        if($r_primer ne "*"){
                $t_primer=$r_primer;
        }
        if(! -r $t_primer){
                print "Program: primer3_core\n";
                die "Error: invalid program. $t_primer\n";
        }
        if($r_ntthal ne "*"){
                $ntthal=$r_ntthal;
        }
        if(! -r $ntthal){
                print "Program: ntthal\n";
                die "Error: invalid program. $ntthal\n";
        }
        if($r_blastn ne "*"){
                $t_blastn=$r_blastn;
        }
        if(! -r $t_blastn){
                print "Program: blastn\n";
                die "Error: invalid program. $t_blastn\n";
        }
        if($r_makeblastdb ne "*"){
                $t_makeblastdb=$r_makeblastdb;
        }
        if(! -r $t_makeblastdb){
                print "Program: makeblastdb\n";
                die "Error: invalid program. $t_makeblastdb\n";
        }
        if($r_windowmasker ne "*"){
                $t_windowmasker=$r_windowmasker;
        }
        if(! -r $t_windowmasker){
                print "Program: windowmasker\n";
                die "Error: invalid program. $t_windowmasker\n";
        }
        
        my $rename_ref;
        $rename_ref=&unzRef($t_ref);
        my %refHash=&getRefLen($t_gff,$rename_ref);
        if($t_mask){
                $rename_ref=&maskGenome($t_windowmasker,$t_ref);
        }
        
        my @plasid_m_seq;        
        my ($upSeq,$downSeq,$seqLen);
        my ($p1_tm,$p2_tm);
        if($rec){
                @plasid_m_seq=&pRedesign($t_dir,$t_primer,$primerPara,$plasmidFa,$t_PrimerMin,$t_PrimerMax,$salt,$dsalt,$dntp,$pCon);
                ($p1_tm,$p2_tm)=&getPtm(\@plasid_m_seq,$pCon,$salt,$dsalt,$dntp,$method);
                
                my $p_mask=0;
                &plasmidPrimer($binDir,$config,$plasmidFa,$t_dir,$p_mask,$t_pNum,$t_thread,$salt,$dsalt,$dntp,$pCon,$thermo,$refine);
        }else{
                ($upSeq,$downSeq,$seqLen)=&plasmid2primer($t_dir,$plasmidFa,$t_upMixPlasmidSize,$t_downMixPlasmidSize);
                @plasid_m_seq=($upSeq,$downSeq);
                ($p1_tm,$p2_tm)=&getPtm(\@plasid_m_seq,$pCon,$salt,$dsalt,$dntp,$method);
                &plasmidCheck($t_dir,$upSeq,$downSeq,$seqLen,$t_upMixPlasmidSize,$t_downMixPlasmidSize,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$p1_tm,$p2_tm);
        }
        #my @plasid_m_plus=&plasmidPlus($plasmidFa,$t_upMixPlasmidSize,$t_downMixPlasmidSize);
        #my @plasid_m_minus=&plasmidMinus($plasmidFa,$t_upMixPlasmidSize,$t_downMixPlasmidSize);
        #&plasmid_primer($t_dir,$plasmidFa,$t_upMixPlasmidSize,$t_downMixPlasmidSize);
        #################
        # geneName
        if($t_geneName || $t_coord){
                my @geneInfo;
                if($t_geneName){
                        @geneInfo=&c_geneNameInterval($t_gff,$t_geneName,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coord){
                        my $loc;
                        ($loc=$t_coord)=~s/:/_/g;
                        $t_geneName="LOC_$loc";
                        @geneInfo=&c_coordInterval($t_gff,$t_coord,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }

                &sc_oneSeq($t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$rename_ref,$t_dir,\@geneInfo,$salt,$dsalt,$dntp,$pCon);

                &sc_multiRun1($t_dir,$t_geneName,$t_primer,$primerPara,$t_thread);

                &c_vpGeneFormat($t_dir,$t_geneName,$refine);
                
                my @info;
                push(@info,join ",",@geneInfo);
                if($run_blast){
                        &sc_oneBlast($t_dir,$rename_ref,$t_blastn,$t_makeblastdb,$t_thread,$t_geneName);
                        my $upAlign=sprintf("%d",$t_upMixGenomeSize*0.8);
                        my $downAlign=sprintf("%d",$t_downMixGenomeSize*0.8);
                        my $identity=80.0;
                        &s_tpGeneMerge($t_dir,$t_geneName,$upAlign,$downAlign,$identity,\@plasid_m_seq,\@info,$p1_tm,$p2_tm);
                        &c_vpGeneMerge($t_dir,$t_geneName,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$refine);
                }else{
                        &s_tpGeneTrans($t_dir,$t_geneName,\@plasid_m_seq,\@info,$p1_tm,$p2_tm);
                        &c_vpGeneTrans($t_dir,$t_geneName,$refine);
                }
                return(0);
        }
        ####################
        # bundle
        if($t_geneList || $t_coordList || $all){
                my @geneInfo;
                my $bundle_counts;
                
                if($t_geneList){
                        @geneInfo=&c_nameListInterval($t_gff,$t_geneList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coordList){
                        @geneInfo=&c_coordListInterval($t_gff,$t_coordList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($all){
                        @geneInfo=&c_allInterval($t_gff,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }
                &sc_listSeq($t_dir,$t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_bundle_size,$t_pNum,$t_mask,$rename_ref,\@geneInfo,$salt,$dsalt,$dntp,$pCon);
                $bundle_counts=&bundleCounts($t_dir);
                &sc_multiRun2($t_dir,$t_primer,$primerPara,$t_thread,$bundle_counts);
                &c_vpFormat($t_dir,$bundle_counts,$refine);
                if($run_blast){
                        &sc_bundleBlast($t_dir,$rename_ref,$t_blastn,$t_makeblastdb,$t_thread);
                        my $upAlign=sprintf("%d",$t_upMixGenomeSize*0.8);
                        my $downAlign=sprintf("%d",$t_downMixGenomeSize*0.8);
                        my $identity=80.0;
                        &s_tpBundleMerge($t_dir,$upAlign,$downAlign,$identity,\@plasid_m_seq,\@geneInfo,$p1_tm,$p2_tm);
                        &c_vpBundleMerge($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$d_endRegion,$d_endMatch,$refine);
                }else{
                        &s_tpBundleTrans($t_dir,\@plasid_m_seq,\@geneInfo,$p1_tm,$p2_tm);
                        &c_vpBundleTrans($t_dir,$refine);
                }
                return(0);
        }
        
        die "Error: lack of parameter (--geneName or --geneList or coordinate or coordList or --all)\n";
}


1;


