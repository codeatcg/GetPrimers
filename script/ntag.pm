

#!/usr/bin/perl -w
use strict;
use threads;
#require "ctag.pl";

sub n_primerLoc{
        my ($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand)=@_;
        my ($e_up_a,$e_up_b,$e_down_a,$e_down_b,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len);
        if($t_type eq "ntag"){
                if($strand eq "+"){
                        $e_up_a=$cdsStart-$t_upGeneFar;
                        $e_up_b=$cdsStart-$t_upGeneNear;
                        $v_up_a=$cdsStart-$t_upVfar;
                        $v_up_b=$cdsStart-$t_upVnear;
                        
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
                        $e_down_a=$cdsStart+$t_downGeneFar;
                        $e_down_b=$cdsStart+$t_downGeneNear;
                        
                        $v_down_a=$cdsStart+$t_downVfar;
                        $v_down_b=$cdsStart+$t_downVnear;
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
                        $e_up_a=$cdsEnd-$t_downGeneFar;
                        $e_up_b=$cdsEnd-$t_downGeneNear;
                        $v_up_a=$cdsEnd-$t_downVfar;
                        $v_up_b=$cdsEnd-$t_downVnear;
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
                        $e_down_a=$cdsEnd+$t_upGeneFar;
                        $e_down_b=$cdsEnd+$t_upGeneNear;
                        
                        $v_down_a=$cdsEnd+$t_upVfar;
                        $v_down_b=$cdsEnd+$t_upVnear;
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

sub n_geneNameInterval{
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
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }else{
                                                                if($cds[0][0]>$cds[$n-1][0]){
                                                                        $cdsStart=$cds[0][0];
                                                                        $cdsEnd=$cds[0][1];
                                                                }else{
                                                                        $cdsStart=$cds[$n-1][0];
                                                                        $cdsEnd=$cds[$n-1][1];
                                                                }
                                                        }
                                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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
                        if($t_geneName eq $geneName || $t_geneName eq $id || $t_geneName eq $dbxref){
                                print "Gene information: $a_geneName\n";
                                my $n=scalar(@cds);
                                if($n<1){
                                        print "Warning: no CDS.\n";
                                }else{
                                        if($strand eq "+"){   
                                                $cdsStart=$cds[0][0];
                                                $cdsEnd=$cds[0][1];
                                        }else{
                                                if($cds[0][0]>$cds[$n-1][0]){
                                                        $cdsStart=$cds[0][0];
                                                        $cdsEnd=$cds[0][1];
                                                }else{
                                                        $cdsStart=$cds[$n-1][0];
                                                        $cdsEnd=$cds[$n-1][1];
                                                }
                                        }
                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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

sub n_nameListInterval{
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
                                        if(exists $nameHash{uc($geneName)} || (exists $nameHash{$id}) || (exists $nameHash{$dbxref})){
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
                                                                        $cdsStart=$cds[0][0];
                                                                        $cdsEnd=$cds[0][1];
                                                                }else{
                                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                                $cdsStart=$cds[0][0];
                                                                                $cdsEnd=$cds[0][1];
                                                                        }else{
                                                                                $cdsStart=$cds[$n-1][0];
                                                                                $cdsEnd=$cds[$n-1][1];
                                                                        }
                                                                }
                                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                        my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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
                                                        $cdsStart=$cds[0][0];
                                                        $cdsEnd=$cds[0][1];
                                                }else{
                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }else{
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[$n-1][1];
                                                        }
                                                }
                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                        my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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

sub n_coordInterval{
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
                                                                                $cdsStart=$cds[0][0];
                                                                                $cdsEnd=$cds[0][1];
                                                                        }else{
                                                                                if($cds[0][0]>$cds[$n-1][0]){
                                                                                        $cdsStart=$cds[0][0];
                                                                                        $cdsEnd=$cds[0][1];
                                                                                }else{
                                                                                        $cdsStart=$cds[$n-1][0];
                                                                                        $cdsEnd=$cds[$n-1][1];
                                                                                }
                                                                        }
                                                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                                my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }else{
                                                                if($cds[0][0]>$cds[$n-1][0]){
                                                                        $cdsStart=$cds[0][0];
                                                                        $cdsEnd=$cds[0][1];
                                                                }else{
                                                                        $cdsStart=$cds[$n-1][0];
                                                                        $cdsEnd=$cds[$n-1][1];
                                                                }
                                                        }
                                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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

sub n_coordListInterval{ 
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
                                                        $cdsStart=$cds[0][0];
                                                        $cdsEnd=$cds[0][1];         
                                                }else{
                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }else{
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[$n-1][1];
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
                                                                                my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                                                if($tout[0] > 0){
                                                                                        my $str=join ",",$a_geneName,$chr,@tout,$strand;
                                                                                        $r_hash{$a_geneName}=1;
                                                                                        push(@out,"$str");
                                                                                }
                                                                        }
                                                                        #
                                                                        if($c_mRNA>1){
                                                                                print "$a_geneName mRNA : $c_mRNA\n";
                                                                        } 
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
                                        $cdsStart=$cds[0][0];
                                        $cdsEnd=$cds[0][1];        
                                }else{
                                        if($cds[0][0]>$cds[$n-1][0]){
                                                $cdsStart=$cds[0][0];
                                                $cdsEnd=$cds[0][1];
                                        }else{
                                                $cdsStart=$cds[$n-1][0];
                                                $cdsEnd=$cds[$n-1][1];
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
                                                                my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                                if($tout[0] > 0){
                                                                        my $str=join ",",$a_geneName,$chr,@tout,$strand;
                                                                        $r_hash{$a_geneName}=1;
                                                                        push(@out,"$str");
                                                                }
                                                        }
                                                        #
                                                        if($c_mRNA>1){
                                                                print "$a_geneName mRNA : $c_mRNA\n";
                                                        }  
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

sub n_allInterval{
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
                                                        $cdsStart=$cds[0][0];
                                                        $cdsEnd=$cds[0][1];
                                                }else{
                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }else{
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[$n-1][1];
                                                        }
                                                }
                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                        my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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
                        #
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
                                        $cdsStart=$cds[0][0];
                                        $cdsEnd=$cds[0][1];
                                }else{
                                        if($cds[0][0]>$cds[$n-1][0]){
                                                $cdsStart=$cds[0][0];
                                                $cdsEnd=$cds[0][1];
                                        }else{
                                                $cdsStart=$cds[$n-1][0];
                                                $cdsEnd=$cds[$n-1][1];
                                        }
                                }
                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                        my @tout=&n_primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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
        #
		close IN;
        if(scalar (@out)<1){
                die "GeneName cat't be found in the GFF file.\n";
		}else{
                return(@out);
        }
}

sub n_oneSeq{
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
                #$t_upR_start=$cdsStart-$v_up_a-$t_PrimerMax+4;
                $t_upR_start=$cdsStart-$v_up_a-$t_PrimerMax+1;
                $t_upR_len=$t_PrimerMax;
                #$t_downF_start=$cdsStart-$v_up_a+4;
                $t_downF_start=$cdsStart-$v_up_a+1;
                $t_downF_len=$t_PrimerMax;
                
                $t_upF_start=$t_up_a-$v_up_a+1;
                $t_upF_len=$t_up_b-$t_up_a+1;
                $t_downR_start=$t_down_b-$v_up_a+1;
                $t_downR_len=$t_down_a-$t_down_b+1;
        }else{
                #$t_upF_start=$cdsEnd-$v_up_a-1;
                $t_upF_start=$cdsEnd-$v_up_a+2;
                $t_upF_len=$t_PrimerMax;
                #$t_downR_start=$cdsEnd-$v_up_a-$t_PrimerMax-1;
                $t_downR_start=$cdsEnd-$v_up_a-$t_PrimerMax+2;
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
}

sub n_listSeq{
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
        #my $fir=shift(@$geneInfo);
        my $fir=$geneInfo->[0];
        my ($t_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$cdsStart,$cdsEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=split /\,/,$fir;
        my ($t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len);
        if($strand eq "+"){
                #$t_upR_start=$cdsStart-$v_up_a-$t_PrimerMax+4;
                $t_upR_start=$cdsStart-$v_up_a-$t_PrimerMax+1;
                $t_upR_len=$t_PrimerMax;
                #$t_downF_start=$cdsStart-$v_up_a+4;
                $t_downF_start=$cdsStart-$v_up_a+1;
                $t_downF_len=$t_PrimerMax;

                
                $t_upF_start=$t_up_a-$v_up_a+1;
                $t_upF_len=$t_up_b-$t_up_a+1;
                $t_downR_start=$t_down_b-$v_up_a+1;
                $t_downR_len=$t_down_a-$t_down_b+1;
        }else{
                #$t_upF_start=$cdsEnd-$v_up_a-1;
                $t_upF_start=$cdsEnd-$v_up_a+2;
                $t_upF_len=$t_PrimerMax;
                #$t_downR_start=$cdsEnd-$v_up_a-$t_PrimerMax-1;
                $t_downR_start=$cdsEnd-$v_up_a-$t_PrimerMax+2;
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
                        #$t_upR_start=$cdsStart-$v_up_a-$t_PrimerMax+4;
                        $t_upR_start=$cdsStart-$v_up_a-$t_PrimerMax+1;
                        $t_upR_len=$t_PrimerMax;
                        #$t_downF_start=$cdsStart-$v_up_a+4;
                        $t_downF_start=$cdsStart-$v_up_a+1;
                        $t_downF_len=$t_PrimerMax;
                        
                        $t_upF_start=$t_up_a-$v_up_a+1;
                        $t_upF_len=$t_up_b-$t_up_a+1;
                        $t_downR_start=$t_down_b-$v_up_a+1;
                        $t_downR_len=$t_down_a-$t_down_b+1;
                }else{
                        #$t_upF_start=$cdsEnd-$v_up_a-1;
                        $t_upF_start=$cdsEnd-$v_up_a+2;
                        $t_upF_len=$t_PrimerMax;
                        #$t_downR_start=$cdsEnd-$v_up_a-$t_PrimerMax-1;
                        $t_downR_start=$cdsEnd-$v_up_a-$t_PrimerMax+2;
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

sub ntag_main{
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
        # geneName
        if($t_geneName || $t_coord){
                my @geneInfo;
                if($t_geneName){
                        @geneInfo=&n_geneNameInterval($t_gff,$t_geneName,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coord){
                        my $loc;
                        ($loc=$t_coord)=~s/:/_/g;
                        $t_geneName="LOC_$loc";
                        @geneInfo=&n_coordInterval($t_gff,$t_coord,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }

                &n_oneSeq($t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$rename_ref,$t_dir,\@geneInfo,$salt,$dsalt,$dntp,$pCon);

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
                        @geneInfo=&n_nameListInterval($t_gff,$t_geneList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coordList){
                        @geneInfo=&n_coordListInterval($t_gff,$t_coordList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($all){
                        @geneInfo=&n_allInterval($t_gff,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }
                &n_listSeq($t_dir,$t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_bundle_size,$t_pNum,$t_mask,$rename_ref,\@geneInfo,$salt,$dsalt,$dntp,$pCon);
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

#####################
# short arm
#####################
sub n_geneSeq{
        my ($seq,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand)=@_;
        my $vseq=substr($seq,$v_up_a-1,$v_len);
        $vseq=~s/[^ATCGNatcgn]/N/g;
        #my $t_gene_upseq=substr($$seq,$geneStart-$t_upMixGenomeSize-1,$t_upMixGenomeSize);
        my ($t_gene_upseq,$gene_mix_upseq,$gene_mix_downseq);
        if($strand eq "+"){
                #$t_gene_upseq=substr($vseq,$cdsStart-$v_up_a-$t_upMixGenomeSize+3,$t_upMixGenomeSize);
                $t_gene_upseq=substr($vseq,$cdsStart-$v_up_a-$t_upMixGenomeSize,$t_upMixGenomeSize);
                #$gene_mix_downseq=substr($vseq,$cdsStart-$v_up_a+3,$t_downMixGenomeSize);
                $gene_mix_downseq=substr($vseq,$cdsStart-$v_up_a,$t_downMixGenomeSize);
        }else{
                #$t_gene_upseq=substr($vseq,$cdsEnd-$v_up_a-$t_downMixGenomeSize-2,$t_downMixGenomeSize);
                #$gene_mix_downseq=substr($vseq,$cdsEnd-$v_up_a-2,$t_upMixGenomeSize);
                
                $t_gene_upseq=substr($vseq,$cdsEnd-$v_up_a-$t_downMixGenomeSize+1,$t_downMixGenomeSize);
                $gene_mix_downseq=substr($vseq,$cdsEnd-$v_up_a+1,$t_upMixGenomeSize);
        }
        $gene_mix_upseq=reverse($t_gene_upseq);
        #$gene_mix_upseq=~s/[^ATCGNatcgn]/N/g;
        $gene_mix_upseq=~tr/ATCGatcg/TAGCtagc/;
        #
        #$gene_mix_downseq=~s/[^ATCGNatcgn]/N/g;
        return($vseq,$gene_mix_upseq,$gene_mix_downseq);
}

sub n_geneSeq2{
        my ($t_ref,$chr,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand)=@_;
        my $seq=&chrSeq($t_ref,$chr);
        &n_geneSeq($seq,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand);
}

sub sn_oneSeq{
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

        my ($a_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$cdsStart,$cdsEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=@$t_info;
        #
        my $v_up_start=1;
        my $v_up_len=$v_up_b-$v_up_a+1;
        #
        my $v_down_start=$v_down_b-$v_up_a+1;
        my $v_down_len=$v_down_a-$v_down_b+1;
        #
        #############
        my ($vseq,$mix_upseq,$mix_downseq)=&n_geneSeq2($t_ref,$chr,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand);
        my $upN=($mix_upseq=~tr/Nn//);
        my $downN=($mix_downseq=~tr/Nn//);
        if($upN>5 || $downN >5){
                die "flanking sequence of gene has too many Ns\n";
        }
        #############
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

sub sn_listSeq{
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
        #############
        my $seq=&chrSeq($t_ref,$chr);
        my ($vseq,$mix_upseq,$mix_downseq)=&n_geneSeq($seq,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand);
        #############
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
                ($vseq,$mix_upseq,$mix_downseq)=&n_geneSeq($seq,$v_up_a,$v_len,$cdsStart,$cdsEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand);
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

sub short_ntag_main{
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
                        @geneInfo=&n_geneNameInterval($t_gff,$t_geneName,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coord){
                        my $loc;
                        ($loc=$t_coord)=~s/:/_/g;
                        $t_geneName="LOC_$loc";
                        @geneInfo=&n_coordInterval($t_gff,$t_coord,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }

                &sn_oneSeq($t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$rename_ref,$t_dir,\@geneInfo,$salt,$dsalt,$dntp,$pCon);

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
                        @geneInfo=&n_nameListInterval($t_gff,$t_geneList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coordList){
                        @geneInfo=&n_coordListInterval($t_gff,$t_coordList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($all){
                        @geneInfo=&n_allInterval($t_gff,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }
                &sn_listSeq($t_dir,$t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_bundle_size,$t_pNum,$t_mask,$rename_ref,\@geneInfo,$salt,$dsalt,$dntp,$pCon);
                $bundle_counts=&bundleCounts($t_dir);
                &sc_multiRun2($t_dir,$t_primer,$primerPara,$t_thread,$bundle_counts);
                &c_vpFormat($t_dir,$bundle_counts,$refine);
                if($run_blast){
                        &sc_bundleBlast($t_dir,$rename_ref,$t_blastn,$t_makeblastdb,$t_thread);
                        #
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


