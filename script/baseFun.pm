#!/usr/bin/perl -w
use strict;
use threads;
use File::Basename;

sub getProg{
        # bin dir
        my ($dir,$thermo)=@_;
        my ($r_primer,$r_blastn,$r_makeblastdb,$r_windowmasker,$r_primerPara,$ntthal);
	    #my $third_dir;
	    #($third_dir=$dir)=~s/script$/third_party/;
        #my $dir;
        #if($script=~/(.*)\//){
        #        $dir=$1;
        #}
        my $t_dir;
        ($t_dir=$dir)=~s/\/script\d*$//;
	    $r_primer="$t_dir/bin/primer3_core";
        $ntthal="$t_dir/bin/ntthal";
        
	    $r_blastn="$t_dir/bin/blastn";
        $r_makeblastdb="$t_dir/bin/makeblastdb";
        $r_windowmasker="$t_dir/bin/windowmasker";
        if($thermo){
                $r_primerPara="$t_dir/settings_files/primer3web_v4_0_0_modify_settings.txt";
        }else{
                $r_primerPara="$t_dir/settings_files/primer3web_v0_4_0_modify_settings.txt";
        }
        
        return($r_primer,$r_blastn,$r_makeblastdb,$r_windowmasker,$r_primerPara,$ntthal);
}

sub read_config{
        print "Reading config file\n";
        my $t_config=$_[0];
        (open IN,"$t_config") || die "file $t_config open failed. $!\n";
		my %hash;
		while(<IN>){
		        chomp;
				next if(/^#/ || /^\s*$/);
				$_=~s/\s+//g;
				my @line=split /=/,$_;
				if(scalar(@line)!=2 || $line[1] eq ''){    
				        print "Error: config file error.\n";
						die "Trace: $_\n";
				}
				$hash{$line[0]}=$line[1];
		}
		close IN;
        my $r_blastn=$hash{"blastn"};
        my $r_makeblastdb=$hash{"makeblastdb"};
        my $r_primer=$hash{"primer"};
        my $r_ntthal=$hash{"ntthal"};
        my $r_windowmasker=$hash{"windowmasker"};
        #
        my $t_PrimerMin=$hash{"PrimerMin"};
        my $t_PrimerMax=$hash{"PrimerMax"};
        
        my $t_upVFar=$hash{"vUpGeneFar"};
        my $t_upVnear=$hash{"vUpGeneNear"};
        my $t_downVfar=$hash{"vDownGeneFar"};
        my $t_downVnear=$hash{"vDownGeneNear"};

		my $t_upGeneFar=$hash{"upGeneFar"};
		my $t_upGeneNear=$hash{"upGeneNear"};		
		my $t_upMixGenomeSize=$hash{"upMixGenomeSize"};
		my $t_upMixPlasmidSize=$hash{"upMixPlasmidSize"};
        
		my $t_downGeneFar=$hash{"downGeneFar"};
		my $t_downGeneNear=$hash{"downGeneNear"};
		my $t_downMixGenomeSize=$hash{"downMixGenomeSize"};
		my $t_downMixPlasmidSize=$hash{"downMixPlasmidSize"};

		my $t_mismatchMore=$hash{"mismatchMore"};
        my $t_endRegion=$hash{"endRegion"};
        my $t_endMismatch=$hash{"endMismatch"};
        my $t_mismatchAtLeast=$hash{"mismatchAtLeast"};
        my $t_productMinLen=$hash{"productMinLen"};
        my $t_productMaxLen=$hash{"productMaxLen"};
        
        #my $d_matchMost=$hash{"dMatchMost"};
        my $d_GminMinus=$hash{"dGminMinus"};
        my $d_endRegion=$hash{"dEndRegion"};
        my $d_endMatch=$hash{"dEndMatch"};
        
        my $settings=$hash{"primerSettings"};
      
        return($r_primer,$r_ntthal,$r_blastn,$r_makeblastdb,$r_windowmasker,$t_upVFar,$t_upVnear,
        $t_upGeneFar,$t_upGeneNear,$t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_upMixPlasmidSize,
        $t_downVfar,$t_downVnear,$t_downGeneFar,$t_downGeneNear,$t_downMixGenomeSize,$t_downMixPlasmidSize,
        $t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,
        $d_GminMinus,$d_endRegion,$d_endMatch,$settings);
}

sub dMaxTm{
        #my $paraFile=$_[0];
        my ($paraFile,$d_GminMinus)=@_;
        my $d_max_tm=0;
        open(IN,"$paraFile") || die "$paraFile\n$!\n";
        while(<IN>){
                chomp;
                if(/PRIMER_MIN_TM=(.+)/){
                        $d_max_tm=$1-$d_GminMinus;
                        last;
                }
        }        
        close IN;
        if($d_max_tm < 20){
                die "Error: minimum allowable melt temperature was too low!. Please check the PRIMER_MIN_TM in $paraFile\n";
        }
        return($d_max_tm);
}

sub unzRef{
        my $t_ref=$_[0];
        my $refDir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $r=0;
        if($t_ref=~/\.gz$/){
                if(! -s "$refDir/$refName.fa"){
                        $r=system("gzip -dc $t_ref > $refDir/$refName.fa");
                        if($r){
                                die "Error: reference file decompress failed.\n";
                        }
                }
        }else{
                if(! -s "$refDir/$refName.fa"){
                        $r=system("ln -s $t_ref $refDir/$refName.fa");
                        if($r){
                                die "Error: reference symlink failed.\n";
                        }
                }
        }
        return("$refDir/$refName.fa");
}

sub maskGenome{
        my ($t_prog,$t_ref)=@_;
        # linux
        my $refDir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $r=0;
        if(! -s "$refDir/$refName.softmask.fa"){
                if($t_ref=~/\.gz$/){
                        if(! -s "$refDir/$refName.fa"){
                                $r=system("gzip -dc $t_ref > $refDir/$refName.fa");
                                if($r){
                                        die "Error: reference file decompress failed.\n";
                                }
                        }
                }else{
                        if(! -s "$refDir/$refName.fa"){
                                $r=system("ln -s $t_ref $refDir/$refName.fa");
                                if($r){
                                        die "Error: reference symlink failed.\n";
                                }
                        }
                }
                
                if(! -s "$refDir/$refName.fa.unit"){
                        $r=system("$t_prog -checkdup true -mk_counts -in $refDir/$refName.fa -out $refDir/$refName.fa.unit");
                        if($r){
                                die "Error: windowmasker stage 1.\n";
                        }
                        #
                        $r=system("$t_prog -ustat $refDir/$refName.fa.unit -dust true -in $refDir/$refName.fa -out $refDir/$refName.softmask.fa -outfmt fasta");
                        if($r){
                                die "Error: windowmasker stage 2.\n";
                        }
                }else{
                        $r=system("$t_prog -ustat $refDir/$refName.fa.unit -dust true -in $refDir/$refName.fa -out $refDir/$refName.softmask.fa -outfmt fasta");
                        if($r){
                                die "Error: windowmasker stage 2.\n";
                        }   
                } 
        }
        return("$refDir/$refName.softmask.fa");
}

#

sub faRefLen{
        my ($ref,$outfile)=@_;
        (open AN,"$ref") || die "$!\n";
        my $tag_chr=0;
        my $tag_seq=0;
        my $chr;
        my %hash;
        (open OUT,">$outfile") || die "output file (chromosome length) open failed. $!\n";
        my $seqLen=0;
        while(<AN>){
                chomp;
                if(/^>/){
                        if($tag_chr){
                                if($tag_seq){
                                        $hash{$chr}=$seqLen;
                                        print OUT "$chr\t$hash{$chr}\n";
                                }else{
                                        die "Error: fasta file format error. chromosome $chr\n";
                                }
                        }
                        my $tchr=(split /\s+/,$_)[0];
                        $chr=substr($tchr,1);
                        $tag_chr=1;
                        $tag_seq=0;
                        $seqLen=0;
                }else{
                        $tag_seq=1;
                        $seqLen+=length($_);
                }
        }
        close AN;
        if($tag_chr){
                if($tag_seq){
                        $hash{$chr}=$seqLen;
                        print OUT "$chr\t$hash{$chr}\n";
                        close OUT;
                }else{
                        die "Error: fasta file format error\n";
                }
        }
        return(%hash);
}

sub gffRefLen{
        my ($t_gff,$ref,$outfile)=@_;
        if($t_gff=~/\.gz$/){
		        (open IN,"gzip -dc $t_gff|") || die "$t_gff open failed. $!\n";
		}else{
		        (open IN,"$t_gff") || die "$t_gff open failed. $!\n";
		}
        my $indicate=0;
        my %hash;
        while(<IN>){
                chomp;
                if(/^##sequence-region/){
                        my @line=split;
                        if(scalar(@line) != 4){
                                print "Warning: header line format error. $_\n";
                                $indicate=0;
                                last;
                        }else{
                                $hash{$line[1]}=$line[3];
                        }
                        $indicate=1;
                }
        
        }        
        close IN;
        if($indicate){
                return(%hash);
        }else{
                return(&faRefLen($ref,$outfile));               
        }
}

sub gffRefLen2{
        my $len_file=$_[0];
        print "Reading chromosome length\n";
		(open IN,"$len_file") || die "$len_file open failed. $!\n";
        my %hash;
        while(<IN>){
                chomp;
                next if(/^\s*$/);
                my @line=split;
                $hash{$line[0]}=$line[1];
        }
        close IN;
        return(%hash);
}

sub getRefLen{
        my ($gff,$t_ref)=@_;
        my $refDir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $len_file="$refDir/$refName.length";
        if(-s "$refDir/$refName.length"){
                return(&gffRefLen2($len_file));
        }else{
                return(&gffRefLen($gff,$t_ref,$len_file));
        }  
}

# sequence of selection marker
sub plasmid2primer{
        my ($t_dir,$t_pfa,$t_upMixPlasmidSize,$t_downMixPlasmidSize)=@_;
        my ($upSeq,$downSeq);
        
        if($t_pfa =~ /\.gz$/){
                (open IN,"gzip -dc $t_pfa|") || die "$t_pfa open failed. $!\n";
        }else{
                (open IN,"$t_pfa") || die "$t_pfa open failed. $!\n";
        }
        my $seq='';
        while(<IN>){
                chomp;
                next if(/^\s*$/ || />/);
                s/[^ATCGNatcgn]/N/g;
                $seq .=$_;
        }
        
        my $seqLen = length($seq);
        if($seqLen < $t_upMixPlasmidSize+$t_downMixPlasmidSize){
                die "Error: length of insertion sequence (plasmid) is too short!\n";
        }
        my $t_up=substr($seq,0,$t_upMixPlasmidSize);
        my $s_up=reverse($t_up); 
        ($upSeq=$s_up)=~tr/ATCGatcg/TAGCtagc/;
        my $seqlen=length($seq);
        $downSeq=substr($seq,$seqlen-$t_downMixPlasmidSize,$t_downMixPlasmidSize);
        close IN;
        return($upSeq,$downSeq,$seqLen);
}

sub plasmidCheck{
        my ($t_dir,$upSeq,$downSeq,$seqLen,$t_upMixPlasmidSize,$t_downMixPlasmidSize,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$p1_tm,$p2_tm)=@_;
        
        my $v5_loc="1-$t_upMixPlasmidSize:-";
        my $v6_start=$seqLen - $t_downMixPlasmidSize+1;
        my $v6_loc="$v6_start-$seqLen:+";
                
        my $up_range_s=$t_upVnear+$t_upMixPlasmidSize;
        my $up_range_e=$t_upVfar+$t_upMixPlasmidSize;
        my $v1_v5_range="${up_range_s}-$up_range_e";
        
        my $down_range_s=$t_downVnear+$t_downMixPlasmidSize;
        my $down_range_e=$t_downVfar+$t_downMixPlasmidSize;
        my $v6_v2_range="${down_range_s}-$down_range_e";
        
        if(! -d "$t_dir/primer_result"){
                my $r=system("mkdir $t_dir/primer_result");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_result\n";
                }
        }
        
        open(OUA, ">$t_dir/primer_result/plasmid.verify.reverse.txt") || die "Error: $t_dir/primer_result/plasmid.verify.reverse.txt open failed!\n";
        print OUA "#Gene_Name\tPrimer_ID\tV5\tV5_Loc\tV5_Tm\tV1_V5_Len_Range\tV5_Penalty\n";
        print OUA "plasmid_verify\t1\t$upSeq\t$v5_loc\t$p1_tm\t$v1_v5_range\tNA\n";
        close OUA;
        
        open(OUB, ">$t_dir/primer_result/plasmid.verify.forward.txt") || die "Error: $t_dir/primer_result/plasmid.verify.forward.txt open failed!\n";
        print OUB "#Gene_Name\tPrimer_ID\tV6\tV6_Loc\tV6_Tm\tV6_V2_Len_Range\tV6_Penalty\n";
        print OUB "plasmid_verify\t1\t$downSeq\t$v6_loc\t$p2_tm\t$v6_v2_range\tNA\n";
        close OUB;
}

sub chrSeq{
        my ($t_ref,$chr)=@_;
        (open IN,"$t_ref") || die "$!\n";
        my $indicate=0;
        my $seq='';
        while(<IN>){
                chomp;
                next if(/^#/);
                if(/^>/){
                        if($indicate == 1){last;}
                        my $tchr=(split /\s+/,$_)[0];
                        my $xchr=substr($tchr,1,length($tchr)-1);
                        if($xchr eq $chr){  
                                $indicate=1;
                        }
                }else{
                        if($indicate == 1){
                                $seq="$seq" . "$_";
                        }
                
                }
        }        
        close IN;
        #
        if(! $indicate){
                print "Error: chromosome can't be found.\n";
                die "The GFF file may not match the reference file.\n";
        }
        return($seq);
}

sub tPPara{
        my ($vseq,$t_PrimerMin,$t_PrimerMax,$salt,$dsalt,$dntp,$pCon)=@_;
        my $right_start=length($vseq);
        my $t_down_start=$right_start-$t_PrimerMax+1;
        my $para="SEQUENCE_ID=plasmid_target\n" .
                "SEQUENCE_TEMPLATE=$vseq\n" .
                "PRIMER_TASK=generic\n" .
                "PRIMER_PICK_LEFT_PRIMER=1\n" .
                "PRIMER_PICK_RIGHT_PRIMER=1\n" .
                "PRIMER_PRODUCT_SIZE_RANGE=100-10000\n" .
                "PRIMER_PICK_INTERNAL_OLIGO=0\n" .
                "PRIMER_OPT_SIZE=20\n" .
                "PRIMER_MIN_SIZE=$t_PrimerMin\n" .
                "PRIMER_MAX_SIZE=$t_PrimerMax\n" .
                "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,$t_PrimerMax,$t_down_start,$t_PrimerMax\n" .
                "SEQUENCE_FORCE_LEFT_START=1\n" .
                "SEQUENCE_FORCE_RIGHT_START=$right_start\n" .
                
                "PRIMER_DNA_CONC=$pCon\n" . 
                "PRIMER_SALT_MONOVALENT=$salt\n" .
                "PRIMER_SALT_DIVALENT=$dsalt\n" .
                "PRIMER_DNTP_CONC=$dntp\n" . 
                
                "PRIMER_NUM_RETURN=1\n" .
                "=\n";
        return($para);
}

sub tUpPara{
        my ($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_up_start,$t_up_len,$t_down_start,$t_down_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon)=@_;
        my $tstrand;
        my $fix_start;
        if($strand eq "+"){
                $tstrand="plus";
                my $right_start=$t_down_start+$t_down_len-1;
                $fix_start="SEQUENCE_FORCE_RIGHT_START=$right_start";
        }else{
                $tstrand="minus";
                $fix_start="SEQUENCE_FORCE_LEFT_START=$t_up_start";
        }
        
        
        my $para="SEQUENCE_ID=${t_geneName}_target_$tstrand\n" .
                "SEQUENCE_TEMPLATE=$vseq\n" .
                "PRIMER_LOWERCASE_MASKING=$t_mask\n" .
                "PRIMER_TASK=generic\n" .
                "PRIMER_PICK_LEFT_PRIMER=1\n" .
                "PRIMER_PICK_RIGHT_PRIMER=1\n" .
                "PRIMER_PRODUCT_SIZE_RANGE=${t_productMinLen}-$t_productMaxLen\n" .
                "PRIMER_PICK_INTERNAL_OLIGO=0\n" .
                "PRIMER_OPT_SIZE=20\n" .
                "PRIMER_MIN_SIZE=$t_PrimerMin\n" .
                "PRIMER_MAX_SIZE=$t_PrimerMax\n" .
                "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=$t_up_start,$t_up_len,$t_down_start,$t_down_len\n" .
                "$fix_start\n" .
                
                "PRIMER_DNA_CONC=$pCon\n" . 
                "PRIMER_SALT_MONOVALENT=$salt\n" .
                "PRIMER_SALT_DIVALENT=$dsalt\n" .
                "PRIMER_DNTP_CONC=$dntp\n" . 
                
                "PRIMER_NUM_RETURN=$t_pNum\n" .
                "=\n";
        return($para);
}

sub tDownPara{
        my ($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_up_start,$t_up_len,$t_down_start,$t_down_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon)=@_;
        my $tstrand;
        my $fix_start;
        if($strand eq "+"){
                $tstrand="plus";
                $fix_start="SEQUENCE_FORCE_LEFT_START=$t_up_start";
        }else{
                $tstrand="minus";
                my $right_start=$t_down_start+$t_down_len-1;
                $fix_start="SEQUENCE_FORCE_RIGHT_START=$right_start";
        }
        my $para="SEQUENCE_ID=${t_geneName}_target_$tstrand\n" .
                "SEQUENCE_TEMPLATE=$vseq\n" .
                "PRIMER_LOWERCASE_MASKING=$t_mask\n" .
                "PRIMER_TASK=generic\n" .
                "PRIMER_PICK_LEFT_PRIMER=1\n" .
                "PRIMER_PICK_RIGHT_PRIMER=1\n" .
                "PRIMER_PRODUCT_SIZE_RANGE=${t_productMinLen}-$t_productMaxLen\n" .
                "PRIMER_PICK_INTERNAL_OLIGO=0\n" .
                "PRIMER_OPT_SIZE=20\n" .
                "PRIMER_MIN_SIZE=$t_PrimerMin\n" .
                "PRIMER_MAX_SIZE=$t_PrimerMax\n" .
                "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=$t_up_start,$t_up_len,$t_down_start,$t_down_len\n" .
                "$fix_start\n" .
                
                "PRIMER_DNA_CONC=$pCon\n" . 
                "PRIMER_SALT_MONOVALENT=$salt\n" .
                "PRIMER_SALT_DIVALENT=$dsalt\n" .
                "PRIMER_DNTP_CONC=$dntp\n" . 
                
                "PRIMER_NUM_RETURN=$t_pNum\n" .
                "=\n";
        return($para);
}

sub vpPara{
        my ($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_start,$v_up_len,$v_down_start,$v_down_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon)=@_;
        my $tstrand;
        if($strand eq "+"){
                $tstrand="plus";
        }else{
                $tstrand="minus";
        }
        my $para="SEQUENCE_ID=${t_geneName}_verify_$tstrand\n" .
                "SEQUENCE_TEMPLATE=$vseq\n" .
                "PRIMER_LOWERCASE_MASKING=$t_mask\n" .
                "PRIMER_TASK=generic\n" .
                "PRIMER_PICK_LEFT_PRIMER=1\n" .
                "PRIMER_PICK_RIGHT_PRIMER=1\n" .
                "PRIMER_PRODUCT_SIZE_RANGE=${t_productMinLen}-$t_productMaxLen\n" .
                "PRIMER_PICK_INTERNAL_OLIGO=0\n" .
                "PRIMER_OPT_SIZE=20\n" .
                "PRIMER_MIN_SIZE=$t_PrimerMin\n" .
                "PRIMER_MAX_SIZE=$t_PrimerMax\n" .
                "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=$v_up_start,$v_up_len,$v_down_start,$v_down_len\n" .
                
                "PRIMER_DNA_CONC=$pCon\n" . 
                "PRIMER_SALT_MONOVALENT=$salt\n" .
                "PRIMER_SALT_DIVALENT=$dsalt\n" .
                "PRIMER_DNTP_CONC=$dntp\n" . 
                
                "PRIMER_NUM_RETURN=$t_pNum\n" .
                "=\n";
        return($para);
}

sub runPrimer{
        my ($t_prog,$setFile,$outFile,$errFile,$inFile)=@_;
        system("$t_prog --p3_settings_file $setFile --output $outFile --error $errFile $inFile");
}

sub pRedesign{
        my ($t_dir,$s_prog,$s_setFile,$t_pfa,$t_PrimerMin,$t_PrimerMax,$salt,$dsalt,$dntp,$pCon)=@_;
        if(! -d "$t_dir/primer_plasmid"){
                my $r=system("mkdir $t_dir/primer_plasmid");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_plasmid\n";
                }
        }
        #
        if($t_pfa =~ /\.gz$/){
                (open IN,"gzip -dc $t_pfa|") || die "$t_pfa open failed. $!\n";
        }else{
                (open IN,"$t_pfa") || die "$t_pfa open failed. $!\n";
        }
        my $seq='';
        while(<IN>){
                chomp;
                next if(/^\s*$/);
                s/[^ATCGNatcgn]/N/g;
                $seq .=$_;
        }
        close IN;
        
        my $seqLen = length($seq);
        if($seqLen < 70){
                die "Error: length of insertion sequence (plasmid) is too short!\n";
        }
        #
        my $t_inFile="$t_dir/primer_plasmid/plasmid.target.para";
        my $t_outFile="$t_dir/primer_plasmid/plasmid.target.out";
        my $t_errFile="$t_dir/primer_plasmid/plasmid.target.err";
        my $para=&tPPara($seq,$t_PrimerMin,$t_PrimerMax,$salt,$dsalt,$dntp,$pCon);
        #
        (open OUT,">$t_inFile") || die "$!\n";
        print OUT "$para";
        close OUT;
        &runPrimer($s_prog,$s_setFile,$t_outFile,$t_errFile,$t_inFile);
        #
        my ($seq1,$seq2)=('','');
        if(-s $t_outFile){
                open(IN,"$t_outFile") || die "$!\n";
                while(<IN>){
                        chomp;
                        my @line=split /=/,$_;
                        
                        if($line[0]=~/LEFT_\d+_SEQUENCE$/){
                                $seq1=$line[1];
                        }elsif($line[0]=~/RIGHT_\d+_SEQUENCE$/){
                                $seq2=$line[1];                                
                        }

                }
                close IN;
                if($seq1 ne '' && $seq2 ne ''){
                        my $upSeq=reverse($seq1);
                        my $downSeq=reverse($seq2);
                        $upSeq=~tr/ATCGatcg/TAGCtagc/;
                        $downSeq=~tr/ATCGatcg/TAGCtagc/;
                        return($upSeq,$downSeq);
                }else{
                        die "Error: common sequence can't be redesigned! Please check the insertion sequence or try to get common sequence by length.\n";
                }
        }else{
                die "Error: file is not created or empty!\n$t_outFile\n";
        }
}

sub bundleCounts{
        my $t_dir=$_[0];
        my @allwork=glob "$t_dir/primer_para/bundle.*.verify.up.para";
        my $n_bundle=scalar(@allwork);
        if($n_bundle > 0){
                return($n_bundle);
        }else{
                @allwork=glob "$t_dir/primer_para/bundle.*.verify.para";
                $n_bundle=scalar(@allwork);
                if($n_bundle > 0){
                        return($n_bundle);
                }else{
                        die "Error: failed to create parameter files!\n";
                }
        }
}

# @$pair_list, 1..$bundle
sub formatPrimerMulti{
        my ($t_dir,$pair_list,$left_fa,$right_fa,$formatOut,$sta,$force,$refine)=@_;
        if(! -d "$t_dir/primer_blast"){
                my $r=system("mkdir $t_dir/primer_blast");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_blast\n";
                }
        }
        
        (open OUA,">$t_dir/primer_blast/$left_fa") || die "$!\n";
        (open OUB,">$t_dir/primer_blast/$right_fa") || die "$!\n";
        (open OUT,">$t_dir/primer_blast/$formatOut") || die "$!\n";
        (open OUN,">$t_dir/primer_blast/$sta") || die "$!\n";
        
        if($force){
                print OUT "#GeneName\tStrand\tPrimerID\tPenalty_pair\tSeq_F\tStart_F\tLen_F\tTm_F\tGC%_F\tSelf_any_F\tSelf_end_F\tEnd_stab_F\tPenalty_F\tSeq_R\tStart_R\tLen_R\tTm_R\tGC%_R\tSelf_any_R\tSelf_end_R\tEnd_stab_R\tPenalty_R\tWarning\n";
        }else{
                print OUT "#GeneName\tStrand\tPrimerID\tPenalty_pair\tSeq_F\tStart_F\tLen_F\tTm_F\tGC%_F\tSelf_any_F\tSelf_end_F\tEnd_stab_F\tPenalty_F\tSeq_R\tStart_R\tLen_R\tTm_R\tGC%_R\tSelf_any_R\tSelf_end_R\tEnd_stab_R\tPenalty_R\n";
        }
        print OUN "#GeneName\tStat\n";
        #        
        my ($pairPenalty,$num)=(0,0);
        my ($penalty1,$seq1,$start1,$len1,$tm1,$gc1,$selAny1,$selEnd1,$endSta1)=(0,0,0,0,0,0,0,0,0);
        my ($penalty2,$seq2,$start2,$len2,$tm2,$gc2,$selAny2,$selEnd2,$endSta2)=(0,0,0,0,0,0,0,0,0);
        my $geneName;
        my $strand;
        my $tSta=0;
        my ($flag,$ok)=(0,-1);
        my $warnInfo="";
        foreach my $k(@$pair_list){
                (open AN,"$t_dir/primer_raw/$k") || die "$!\n";
                while(<AN>){
                        chomp;
                        my @line=split /=/,$_;
                        
                        if($line[0] eq "SEQUENCE_ID"){
                                if($ok>-1){
                                        if($force){
                                                print OUN "$geneName\t$tSta\n";
                                        }else{
                                                print OUN "$geneName\t$ok\n";
                                        } 
                                }
                                $ok=0;
                                #($geneName,$strand)=(split /_/,$line[1])[0,2];
                                my @tag=split /_/,$line[1];
                                $strand=pop(@tag);
                                pop(@tag);
                                $geneName=join "_",@tag;
                        }elsif($line[0] eq "PRIMER_PAIR_EXPLAIN"){
                                if($line[1]=~/ok\s+(\d+)/){
                                        if($refine){
                                                $tSta=0;
                                                next;
                                        }
                                        
                                        if($1>0){
                                                $tSta=$1;
                                        }else{
                                                $tSta=0;
                                                $pairPenalty="NA";
                                        }
                                }
                        }elsif($line[0]=~/PAIR_(\d+)_PENALTY$/){
                                $pairPenalty=sprintf("%.3f",$line[1]);
                                #$num=$1;
                                $flag=0;
                                $warnInfo="";
                        }elsif($line[0]=~/LEFT_(\d+)_PENALTY$/){
                                $penalty1=sprintf("%.3f",$line[1]);
                                $num=$1;
                        }elsif($line[0]=~/RIGHT_\d+_PENALTY$/){
                                $penalty2=sprintf("%.3f",$line[1]);
                                #$num=$1;
                        }elsif($line[0]=~/_\d+_PROBLEMS$/){
                                $flag=1;
                                #
                                $warnInfo .= "$line[1]";
                        }elsif($line[0]=~/LEFT_\d+_SEQUENCE$/){
                                $seq1=$line[1];
                        }elsif($line[0]=~/RIGHT_\d+_SEQUENCE$/){
                                $seq2=$line[1];                                
                        }elsif($line[0]=~/LEFT_\d+$/){
                                ($start1,$len1)=split /\,/,$line[1];
                        }elsif($line[0]=~/RIGHT_\d+$/){
                                ($start2,$len2)=split /\,/,$line[1];
                        }elsif($line[0]=~/LEFT_\d+_TM$/){
                                $tm1=$line[1];
                        }elsif($line[0]=~/RIGHT_\d+_TM$/){
                                $tm2=$line[1];
                        }elsif($line[0]=~/LEFT_\d+_GC_PERCENT$/){
                                $gc1=$line[1];
                        }elsif($line[0]=~/RIGHT_\d+_GC_PERCENT$/){
                                $gc2=$line[1];
                        }elsif($line[0]=~/LEFT_\d+_SELF_ANY(_TH)?$/){
                                $selAny1=$line[1];
                        }elsif($line[0]=~/RIGHT_\d+_SELF_ANY(_TH)?$/){
                                $selAny2=$line[1];
                        }elsif($line[0]=~/LEFT_\d+_SELF_END(_TH)?$/){
                                $selEnd1=$line[1];
                        }elsif($line[0]=~/RIGHT_\d+_SELF_END(_TH)?$/){
                                $selEnd2=$line[1];
                        }elsif($line[0]=~/LEFT_\d+_END_STABILITY$/){
                                $endSta1=$line[1];
                        }elsif($line[0]=~/RIGHT_\d+_END_STABILITY$/){
                                $endSta2=$line[1];
                                if($force){
                                        if($refine){
                                                my $four_f=substr($seq1,-4);
                                                $four_f=uc($four_f);
                                                if($four_f eq "AAAA" || $four_f eq "TTTT" || $four_f eq "GGGG" || $four_f eq "CCCC"){
                                                        next;
                                                }
                                                
                                                my $four_r=substr($seq2,-4);
                                                $four_r=uc($four_r);
                                                if($four_r eq "AAAA" || $four_r eq "TTTT" || $four_r eq "GGGG" || $four_r eq "CCCC"){
                                                        next;
                                                }
                                                
                                                $tSta++;
                                        }
                                        print OUA ">${geneName}_${num}_left\n";
                                        print OUA "$seq1\n";
                                        print OUB ">${geneName}_${num}_right\n";
                                        print OUB "$seq2\n";
                                        print OUT "$geneName\t$strand\t$num\t$pairPenalty\t$seq1\t$start1\t$len1\t$tm1\t$gc1\t$selAny1\t$selEnd1\t$endSta1\t$penalty1\t$seq2\t$start2\t$len2\t$tm2\t$gc2\t$selAny2\t$selEnd2\t$endSta2\t$penalty2\t$warnInfo\n";
                                }else{
                                        if(! $flag){
                                                if($refine){
                                                        my $four_f=substr($seq1,-4);
                                                        $four_f=uc($four_f);
                                                        if($four_f eq "AAAA" || $four_f eq "TTTT" || $four_f eq "GGGG" || $four_f eq "CCCC"){
                                                                next;
                                                        }
                                                        
                                                        my $four_r=substr($seq2,-4);
                                                        $four_r=uc($four_r);
                                                        if($four_r eq "AAAA" || $four_r eq "TTTT" || $four_r eq "GGGG" || $four_r eq "CCCC"){
                                                                next;
                                                        }
                                                        
                                                }
                                                #
                                                $ok++;
                                                print OUA ">${geneName}_${num}_left\n";
                                                print OUA "$seq1\n";
                                                print OUB ">${geneName}_${num}_right\n";
                                                print OUB "$seq2\n";
                                                print OUT "$geneName\t$strand\t$num\t$pairPenalty\t$seq1\t$start1\t$len1\t$tm1\t$gc1\t$selAny1\t$selEnd1\t$endSta1\t$penalty1\t$seq2\t$start2\t$len2\t$tm2\t$gc2\t$selAny2\t$selEnd2\t$endSta2\t$penalty2\n";
                                        }
                                }                                
                        }

                }
                if($ok>-1){
                        if($force){
                                print OUN "$geneName\t$tSta\n";
                        }else{
                                print OUN "$geneName\t$ok\n";
                        } 
                }
                close AN;
        }
        
        close OUA;
        close OUB;
        close OUT;
        close OUN;
}

sub tpUpFormat{
        my ($t_dir,$bundle,$force,$refine)=@_;
        my @out=map {"bundle.${_}.target.up.out"} (1..$bundle);       
        my $left_fa="bundle.target.up.forward.fa";
        my $right_fa="bundle.target.up.reverse.fa";
        my $formatOut="bundle.target.up.format.out";
        my $sta="bundle.target.up.stat";
        &formatPrimerMulti($t_dir,\@out,$left_fa,$right_fa,$formatOut,$sta,$force,$refine);
}

sub tpDownFormat{
        my ($t_dir,$bundle,$force,$refine)=@_;
        my @out=map {"bundle.${_}.target.down.out"} (1..$bundle);       
        my $left_fa="bundle.target.down.forward.fa";
        my $right_fa="bundle.target.down.reverse.fa";
        my $formatOut="bundle.target.down.format.out";
        my $sta="bundle.target.down.stat";
        &formatPrimerMulti($t_dir,\@out,$left_fa,$right_fa,$formatOut,$sta,$force,$refine);
}

sub tpUpGeneFormat{
        my ($t_dir,$gName,$force,$refine)=@_;
        my @out=("$gName.target.up.out");
        my $left_fa="$gName.target.up.forward.fa"; 
        my $right_fa="$gName.target.up.reverse.fa";
        my $formatOut="$gName.target.up.format.out";
        my $sta="$gName.target.up.stat";
        &formatPrimerMulti($t_dir,\@out,$left_fa,$right_fa,$formatOut,$sta,$force,$refine);
}

sub tpDownGeneFormat{
        my ($t_dir,$gName,$force,$refine)=@_;
        my @out=("$gName.target.down.out");
        my $left_fa="$gName.target.down.forward.fa"; 
        my $right_fa="$gName.target.down.reverse.fa";
        my $formatOut="$gName.target.down.format.out";
        my $sta="$gName.target.down.stat";
        &formatPrimerMulti($t_dir,\@out,$left_fa,$right_fa,$formatOut,$sta,$force,$refine);
}

sub checkBlastdb{
        my ($t_ref,$makeblastdb)=@_;
        if((-s "$t_ref.fa.ndb") && (-s "$t_ref.fa.nhd") && (-s "$t_ref.fa.nhi") && (-s "$t_ref.fa.nhr") && (-s "$t_ref.fa.nin")  && (-s "$t_ref.fa.nog") && (-s "$t_ref.fa.not") && (-s "$t_ref.fa.nsq") && (-s "$t_ref.fa.ntf") && (-s "$t_ref.fa.nto")){
                print "balstdb files have existed.\n";
        }else{
                system("$makeblastdb  -in $t_ref.fa  -dbtype nucl -hash_index");
        }
}

sub tpUpFilt{
        my ($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen)=@_;
        my $fblast="$t_dir/primer_blast/blastn.lp.out";
        my $rblast="$t_dir/primer_blast/blastn.lmix.out";
        my $outFile="$t_dir/primer_blast/target.up.length";
        &vPairFilt($t_dir,$fblast,$rblast,$outFile,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
}

sub tpDownFilt{
        my ($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen)=@_;
        my $fblast="$t_dir/primer_blast/blastn.rmix.out";
        my $rblast="$t_dir/primer_blast/blastn.rp.out";
        my $outFile="$t_dir/primer_blast/target.down.length";
        &vPairFilt($t_dir,$fblast,$rblast,$outFile,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
}

sub vEndBase{
        my ($seq1,$seq2) = @_;
        my $taBase = 0;
        my $base1 = substr($seq1,-1);
        my $base2 = substr($seq2,-1);
        if($base1 eq 'A' || $base1 eq 'T' || $base1 eq 'a' || $base1 eq 't'){
                $taBase++;
        }
        if($base2 eq 'A' || $base2 eq 'T' || $base2 eq 'a' || $base2 eq 't'){
                $taBase++;
        }
        return($taBase);
}

# Tm, dimer
# my ($seq1,$seq2,$dMMost,$regionLen,$matchLen)
# info -> 
sub tpTrans{
        my ($r_up_primer,$r_down_primer,$gene_info,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo) = @_;
        my @group1;
        my @group2;
        my @out;
        
        my ($geneName,$strand,$g1_g3_penalty,$g1_seq,$g1_start,$g1_len,$g1_tm,$g3_seq,$g3_start,$g3_len,$g3_tm,$g3_warn);
        my ($g4_g2_penalty,$g4_seq,$g4_start,$g4_len,$g4_tm,$g2_seq,$g2_start,$g2_len,$g2_tm,$g4_warn);
        my ($g1_loc_s,$g1_loc_e,$g2_loc_s,$g2_loc_e,$g3_loc_s,$g3_loc_e,$g4_loc_s,$g4_loc_e);
        my ($g1_loc,$g2_loc,$g3_loc,$g4_loc);
        my ($left_d_len,$right_d_len);
        my $up_n = scalar(@$r_up_primer);
        my $down_n = scalar(@$r_down_primer);
        my $outName;
        my ($chr,$gStart,$gEnd,$v_up_a) = @$gene_info;
        my ($up_rmix,$down_fmix)=@$plasSeq;
        #
        for(my $i = 0; $i < $up_n; $i++){
                my @up_arr = split /\t+/,$r_up_primer->[$i];
                my ($geneName,$strand) = @up_arr[0,1];
                ($outName=$geneName)=~s/\|/\t/g;
                if($strand eq "plus"){
                        $strand = "+";
                        ($g1_g3_penalty,$g1_seq,$g1_start,$g1_len,$g1_tm,$g3_seq,$g3_start,$g3_len,$g3_tm) = @up_arr[3,4,5,6,7,13,14,15,16];
                        $g1_loc_s=$v_up_a+$g1_start-1;
                        $g1_loc_e=$g1_loc_s+$g1_len-1;
                        
                        $g3_loc_e=$v_up_a+$g3_start-1;
                        $g3_loc_s=$g3_loc_e-$g3_len+1;
                        
                        $g1_loc="${g1_loc_s}-${g1_loc_e}:+";
                        $g3_loc="${g3_loc_s}-${g3_loc_e}:-";
                        
                        $left_d_len=$g3_loc_e-$g1_loc_s+1;
                }else{
                        $strand = "-";
                        ($g1_g3_penalty,$g3_seq,$g3_start,$g3_len,$g3_tm,$g1_seq,$g1_start,$g1_len,$g1_tm) = @up_arr[3,4,5,6,7,13,14,15,16];
                        $g3_loc_s=$v_up_a+$g3_start-1;
                        $g3_loc_e=$g3_loc_s+$g3_len-1;
                        
                        $g1_loc_e=$v_up_a+$g1_start-1;
                        $g1_loc_s=$g1_loc_e-$g1_len+1;
                        
                        $g1_loc="${g1_loc_s}-${g1_loc_e}:-";
                        $g3_loc="${g3_loc_s}-${g3_loc_e}:+";
                        
                        $left_d_len=$g1_loc_e-$g3_loc_s+1;
                }
                my $left_r="$up_rmix" . "$g3_seq";
                my $fplasmid=reverse($left_r);
                $fplasmid=~tr/ATCGatcg/TAGCtagc/;
                
                $left_d_len += length($up_rmix);
                
                if(scalar(@up_arr) == 23){
                        $g3_warn=$up_arr[22];
                }else{
                        $g3_warn = "";
                }
                #
                for(my $j = 0; $j < $down_n; $j++){
                        my @down_arr = split /\t+/,$r_down_primer->[$j];
                        if($strand eq "+"){
                                ($g4_g2_penalty,$g4_seq,$g4_start,$g4_len,$g4_tm,$g2_seq,$g2_start,$g2_len,$g2_tm) = @down_arr[3,4,5,6,7,13,14,15,16];
                                $g4_loc_s=$v_up_a+$g4_start-1;
                                $g4_loc_e=$g4_loc_s+$g4_len-1;
                                
                                $g2_loc_e=$v_up_a+$g2_start-1;
                                $g2_loc_s=$g2_loc_e-$g2_len+1;
                                
                                $g4_loc="${g4_loc_s}-${g4_loc_e}:+";
                                $g2_loc="${g2_loc_s}-${g2_loc_e}:-";
                                
                                $right_d_len=$g2_loc_e-$g4_loc_s+1;
                                
                        }else{
                                ($g4_g2_penalty,$g2_seq,$g2_start,$g2_len,$g2_tm,$g4_seq,$g4_start,$g4_len,$g4_tm) = @down_arr[3,4,5,6,7,13,14,15,16];
                                $g2_loc_s=$v_up_a+$g2_start-1;
                                $g2_loc_e=$g2_loc_s+$g2_len-1;
                                
                                $g4_loc_e=$v_up_a+$g4_start-1;
                                $g4_loc_s=$g4_loc_e-$g4_len+1;
                                
                                $g4_loc="${g4_loc_s}-${g4_loc_e}:-";
                                $g2_loc="${g2_loc_s}-${g2_loc_e}:+";
                                
                                $right_d_len=$g4_loc_e-$g2_loc_s+1;
                                
                        }
                        my $right_f="$down_fmix" . "$g4_seq";
                        my $rplasmid=reverse($right_f);
                        $rplasmid=~tr/ATCGatcg/TAGCtagc/;
                        
                        $right_d_len += length($down_fmix);
                        #
                        if($g1_tm - $g2_tm <= 6 && $g1_tm - $g2_tm >= -6){
                                my $g1_g2_dimer;
                                if($thermo){
                                        $g1_g2_dimer = &th_dimer($ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$g1_seq,$g2_seq);
                                }else{
                                        $g1_g2_dimer = &dimer($g1_seq,$g2_seq,$dERegion,$dEMatch);
                                }
                                my $part1 = "$outName\t$strand\t$chr\t$gStart\t$gEnd\t$v_up_a";
                                #my $part1 = "$outName\t$strand\t$chr\t$gStart\t$gEnd";
                                my $part2 = "$g1_seq\t$left_r\t$right_f\t$g2_seq\t$fplasmid\t$rplasmid\t$g1_loc\t$g3_loc\t$g4_loc\t$g2_loc\t$g1_tm\t$g3_tm\t$g4_tm\t$g2_tm\t$p1_tm\t$p2_tm\t$left_d_len\t$right_d_len\t$g1_g2_dimer\t$g1_g3_penalty\t$g4_g2_penalty";
                                
                                my $warn = "";
                                if($force){
                                        if(scalar(@down_arr)==23){
                                                $g4_warn = $down_arr[22];
                                        }else{
                                                $g4_warn = "";
                                        }
                                        if($g3_warn eq ""){
                                                if($g4_warn eq ""){
                                                        $warn = "None";
                                                }else{
                                                        $warn = "G4:$g4_warn";
                                                }
                                        }else{
                                                if($g4_warn eq ""){
                                                        $warn = "G3:$g3_warn";
                                                }else{
                                                        $warn = "G3:$g3_warn G4:$g4_warn";
                                                }
                                        }
                                }
                                my $m_penalty = ($g1_g3_penalty + $g4_g2_penalty) / 2;        
                                if($g1_g2_dimer){
                                        push(@group2,[$part1,$part2,2,$warn,$m_penalty]);
                                }else{
                                        push(@group1,[$part1,$part2,1,$warn,$m_penalty]);
                                }
                                
                        }
                }
        }
        #
        my $t = 0;
        my @allprimer = (\@group1,\@group2);
        my $flag = 0;
        my @pos = (0,3,4,5);
        foreach my $k(@allprimer){
                if($flag == 1){
                        last;
                }
                my @sp;
                if($refine){
                        foreach my $x(@$k){
                                my $seqinfo = $x->[1];
                                my @xinfo = split /\t/,$seqinfo;
                                
                                my $taBase = 0;
                                foreach my $p(@pos){
                                        my $seq1 = $xinfo[$p];
                                        my $base1 = substr($seq1,-1);
                                        if($base1 eq 'A' || $base1 eq 'T' || $base1 eq 'a' || $base1 eq 't'){
                                                $taBase++;
                                        }
                                }
                                
                                push(@$x,$taBase);
                        }
                        #
                        @sp = sort {$a->[5]<=>$b->[5] or $a->[4]<=>$b->[4]} @$k;
                }else{
                        @sp = sort {$a->[4]<=>$b->[4]} @$k;
                }
                foreach my $f(@sp){
                        $t++;
                        if($force){
                                push(@out,"$f->[0]\t$t\t$f->[1]\t$f->[2]\t$f->[3]");
                        }else{
                                push(@out,"$f->[0]\t$t\t$f->[1]\t$f->[2]");
                        }
                        if($t == $pNum){
                                $flag = 1;
                                last;
                        }
                }
                
        }
        return(@out);
}


sub tcommon{
        my ($up_sta,$down_sta) = @_;
        (open BN,"$up_sta") || die "$!\n";
        (open DN,"$down_sta") || die "$!\n";
        <BN>;
        my %hash;
        my (@upno,@downno);
        while(<BN>){
                chomp;
                my ($geneName,$sta) = split /\s+/,$_;
                if($sta > 0){
                        $hash{$geneName} = 1;
                }else{
                        push(@upno,$geneName);
                }
        }
        close BN;
        #
        <DN>;
        my %common;
        #        
        while(<DN>){
                chomp;
                my ($geneName,$sta) = split /\s+/,$_;
                if($sta > 0){
                        if(exists $hash{$geneName}){
                                $common{$geneName} = 1;
                                delete($hash{$geneName});
                        }else{
                                push(@downno,$geneName);
                        }
                }else{
                        push(@downno,$geneName);
                }
        }
        close DN;
        if(scalar(keys %common) < 1){
                die "\nTargetting primers that meet the criteria do not exist!\n";
        }
        #####################
        print "Targetting primers that meet the criteria from upstream of genes don't exist:\n";
        foreach my $k(@upno){
                print "$k\n";
        }
        
        print "\nTargetting primers that meet the criteria from downstream of genes don't exist:\n";
        foreach my $s(@downno){
                print "$s\n";
        }
        return(%common);
}

sub tpStrandTrans{
        my ($t_dir,$f_format,$r_format,$up_sta,$down_sta,$info,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo)=@_;        
        my $r=0;
        my %common = &tcommon($up_sta,$down_sta);
        if(! -d "$t_dir/primer_result"){
                $r=system("mkdir $t_dir/primer_result");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_result\n";
                }
        }
        
        (open OUT,">$t_dir/primer_result/target.primer.unblast.txt") || die "$!\n";
        
        (open EN,"$f_format") || die "$!\n";
        (open FN,"$r_format") || die "$!\n";
        my $f_head = <EN>;
        chomp($f_head);
        my $force = 0;
        if($f_head =~ /Warning$/){
                print OUT "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tSeq_Start\tPrimer_ID\tG1\tG3\tG4\tG2\tP1\tP2\tG1_Loc\tG3_Loc\tG4_Loc\tG2_Loc\tG1_Tm\tG3_Tm\tG4_Tm\tG2_Tm\tP1_Tm\tP2_Tm\tG1_G3_Expect_Len\tG4_G2_Expect_Len\tG1_G2_Dimer\tG1_G3_Penalty\tG4_G2_penalty\tTier\tPrimer_Warning\n";
                $force = 1;
        }else{
                print OUT "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tSeq_Start\tPrimer_ID\tG1\tG3\tG4\tG2\tP1\tP2\tG1_Loc\tG3_Loc\tG4_Loc\tG2_Loc\tG1_Tm\tG3_Tm\tG4_Tm\tG2_Tm\tP1_Tm\tP2_Tm\tG1_G3_Expect_Len\tG4_G2_Expect_Len\tG1_G2_Dimer\tG1_G3_Penalty\tG4_G2_penalty\tTier\n";
        }
        <FN>;
        my $pre_up = "";
        my @up_p;
        my @down_p;
        my ($flag,$flag2) = (0,0);
        my @geneInfo;
        my $start = 0;
        my $end = scalar(@$info);
        my ($up_gene,$down_gene,$up_line,$down_line);
        my $fileEnd = 1;
        print "\nGenes filtered:\n";
        while(<EN>){
                chomp;
                $up_line = $_;
                $up_gene = (split /\s+/,$_)[0];
                if($up_gene ne $pre_up){
                        if($flag == 1){
                                @down_p = ();
                                $flag2 = 0;
                                if($down_gene eq $pre_up){
                                        push(@down_p,$down_line);
                                        $flag2 = 1;
                                }
                                #
                                ##
                                $fileEnd = 0;
                                while(1){
                                        $down_line = <FN>;
                                        if(! $down_line){
                                                $fileEnd = 1;
                                                last;
                                        }
                                        chomp($down_line);
                                        $down_gene = (split /\s+/,$down_line)[0];
                                        #
                                        if($down_gene eq $pre_up){
                                                push(@down_p,$down_line);
                                                $flag2 = 1;
                                        }else{
                                                if($flag2){
                                                        for(my $k = $start; $k < $end; $k++){
                                                                my $tline = $info->[$k];
                                                                my @allInfo = split /\,/,$tline;
                                                                my $geneName = $allInfo[0];
                                                                if($geneName eq $pre_up){
                                                                        @geneInfo = @allInfo[1,6,7,8];
                                                                        $start = $k + 1;
                                                                        last;
                                                                }
                                                        }
                                                        my @tout = &tpTrans(\@up_p,\@down_p,\@geneInfo,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                                                        if(scalar(@tout) > 0){
                                                                foreach my $s(@tout){
                                                                        print OUT "$s\n";
                                                                }
                                                        }else{
                                                                print "$pre_up\n";
                                                        }
                                                        last;
                                                }
                                        }
                                }
                                if($fileEnd){
                                        if($flag2){
                                                for(my $k = $start; $k < $end; $k++){
                                                        my $tline = $info->[$k];
                                                        my @allInfo = split /\,/,$tline;
                                                        my $geneName = $allInfo[0];
                                                        if($geneName eq $pre_up){
                                                                @geneInfo = @allInfo[1,6,7,8];
                                                                $start = $k + 1;
                                                                last;
                                                        }
                                                }
                                                my @tout = &tpTrans(\@up_p,\@down_p,\@geneInfo,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                                                if(scalar(@tout) > 0){
                                                        foreach my $s(@tout){
                                                                print OUT "$s\n";
                                                        }
                                                }else{
                                                        print "$pre_up\n";
                                                }
                                                #
                                                $flag = 0;
                                        }
                                        last;
                                }
                        }
                        #
                        if(exists $common{$up_gene}){
                                $flag = 1;
                                @up_p = ();
                                push(@up_p,$up_line);
                        }else{
                                $flag = 0;
                        }
                }else{
                        if($flag == 1){
                                push(@up_p,$up_line);
                        }
                
                }
                #
                $pre_up = $up_gene;
        }
        if($flag == 1){
                @down_p = ();
                $flag2 = 0;
                if($down_gene eq $pre_up){
                        push(@down_p,$down_line);
                }
                ##
                $fileEnd = 0;
                while(1){
                        $down_line = <FN>;
                        if(! $down_line){
                                $fileEnd = 1;
                                last;
                        }
                        chomp($down_line);
                        $down_gene = (split /\s+/,$down_line)[0];
                        #
                        if($down_gene eq $pre_up){
                                push(@down_p,$down_line);
                                $flag2 = 1;
                        }else{
                                if($flag2){
                                        for(my $k = $start; $k < $end; $k++){
                                                my $tline = $info->[$k];
                                                my @allInfo = split /\,/,$tline;
                                                my $geneName = $allInfo[0];
                                                if($geneName eq $pre_up){
                                                        @geneInfo = @allInfo[1,6,7,8];
                                                        $start = $k + 1;
                                                        last;
                                                }
                                        }
                                        my @tout = &tpTrans(\@up_p,\@down_p,\@geneInfo,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                                        if(scalar(@tout) > 0){
                                                foreach my $s(@tout){
                                                        print OUT "$s\n";
                                                }
                                        }else{
                                                print "$pre_up\n";
                                        }
                                        last;
                                }
                                #
                                
                        }
                }
                if($fileEnd){
                        if($flag2){
                                for(my $k = $start; $k < $end; $k++){
                                        my $tline = $info->[$k];
                                        my @allInfo = split /\,/,$tline;
                                        my $geneName = $allInfo[0];
                                        if($geneName eq $pre_up){
                                                @geneInfo = @allInfo[1,6,7,8];
                                                $start = $k + 1;
                                                last;
                                        }
                                }
                                my @tout = &tpTrans(\@up_p,\@down_p,\@geneInfo,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                                if(scalar(@tout) > 0){
                                        foreach my $s(@tout){
                                                print OUT "$s\n";
                                        }
                                }else{
                                        print "$pre_up\n";
                                }
                        }
                }
        }
        close OUT;
}

sub tpGeneTrans{
        my ($t_dir,$geneName,$plasSeq,$info,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo)=@_;
        my $up_format = "$t_dir/primer_blast/$geneName.target.up.format.out";
        my $down_format = "$t_dir/primer_blast/$geneName.target.down.format.out";
        my $up_sta = "$t_dir/primer_blast/$geneName.target.up.stat";
        my $down_sta = "$t_dir/primer_blast/$geneName.target.down.stat";
        &tpStrandTrans($t_dir,$up_format,$down_format,$up_sta,$down_sta,$info,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
}

sub tpBundleTrans{
        my ($t_dir,$plasSeq,$info,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo)=@_;
        my $up_format="$t_dir/primer_blast/bundle.target.up.format.out";
        my $down_format="$t_dir/primer_blast/bundle.target.down.format.out";     
        my $up_sta = "$t_dir/primer_blast/bundle.target.up.stat";
        my $down_sta = "$t_dir/primer_blast/bundle.target.down.stat";
        &tpStrandTrans($t_dir,$up_format,$down_format,$up_sta,$down_sta,$info,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
}

sub chPcrLen{
        my ($pcr_str,$add)=@_;
        my @line1=split /\;/,$pcr_str;
        my @out1;
        foreach my $k(@line1){
                push(@out1,$k+$add);
        }
        my $outpcr=join ";",@out1;
}

# primer concentration 50nM, monovalent cation concentration 50mM
sub atcgCount{
        my $seq=$_[0];
        my ($a,$t,$c,$g)=(0,0,0,0);
        $a=($seq=~tr/Aa//);
        $t=($seq=~tr/Tt//);
        $c=($seq=~tr/Cc//);
        $g=($seq=~tr/Gg//);
        my $len=length($seq);
        return($len,$a,$t,$c,$g);
}

sub basicTm{
        my $seq=$_[0];
        my ($len,$a,$t,$c,$g)=&atcgCount($seq);
        my $tm;
        if($len<14){
                $tm=($a+$t)*2+($g+$c)*4;
        }else{
                $tm=64.9+41*($g+$c-16.4)/($a+$t+$g+$c);
        }
        return(sprintf("%.3f",$tm));
}

sub longTm{
        my ($seq,$salt)=@_;
        $salt=$salt*1e-3;
        my ($len,$a,$t,$c,$g)=&atcgCount($seq);
        my $tm;
        my $baseLen=$a+$t+$g+$c;
        my $gc=100*($g+$c)/$baseLen;
        $tm=81.5+16.6*log($salt/(1+0.7*$salt))/log(10)+0.41*$gc-500/$baseLen;
        return(sprintf("%.3f",$tm));
}

sub saltAdj{
        my ($seq,$salt)=@_;
        $salt=$salt*1e-3;
        my ($len,$a,$t,$c,$g)=&atcgCount($seq);
        my $tm;
        if($len<14){
                $tm=($a+$t)*2+($g+$c)*4-16.6*(log(0.05)/log(10))+16.6*(log($salt)/log(10));
        }else{
                $tm=100.5+41*($g+$c)/($a+$t+$g+$c)-820/($a+$t+$g+$c)+16.6*(log($salt)/log(10));
        }
        return(sprintf("%.3f",$tm));
}

sub nnParaH{
        my %pdeltaH;
        $pdeltaH{"AA"}=-7.9; $pdeltaH{"AT"}=-7.2; $pdeltaH{"AG"}=-7.8; $pdeltaH{"AC"}=-8.4; $pdeltaH{"AN"}=-7.2;
        $pdeltaH{"TA"}=-7.2; $pdeltaH{"TT"}=-7.9; $pdeltaH{"TG"}=-8.5; $pdeltaH{"TC"}=-8.2; $pdeltaH{"TN"}=-7.2;        
        $pdeltaH{"CA"}=-8.5; $pdeltaH{"CT"}=-7.8; $pdeltaH{"CG"}=-10.6; $pdeltaH{"CC"}=-8.0; $pdeltaH{"CN"}=-7.8;        
        $pdeltaH{"GA"}=-8.2; $pdeltaH{"GT"}=-8.4; $pdeltaH{"GG"}=-8.0; $pdeltaH{"GC"}=-9.8; $pdeltaH{"GN"}=-8.0;
        $pdeltaH{"NA"}=-7.2; $pdeltaH{"NT"}=-7.2; $pdeltaH{"NG"}=-7.8; $pdeltaH{"NC"}=-8.0; $pdeltaH{"NN"}=-7.2;
        return(%pdeltaH);
}

sub nnParaS{
        my %pdeltaS;
        $pdeltaS{"AA"}=-22.2; $pdeltaS{"AT"}=-20.4; $pdeltaS{"AG"}=-21.0; $pdeltaS{"AC"}=-22.4; $pdeltaS{"AN"}=-22.4;
        $pdeltaS{"TA"}=-21.3; $pdeltaS{"TT"}=-22.2; $pdeltaS{"TG"}=-22.7; $pdeltaS{"TC"}=-22.2; $pdeltaS{"TN"}=-22.7;     
        $pdeltaS{"CA"}=-22.7; $pdeltaS{"CT"}=-21.0; $pdeltaS{"CG"}=-27.2; $pdeltaS{"CC"}=-19.9; $pdeltaS{"CN"}=-27.2;       
        $pdeltaS{"GA"}=-22.2; $pdeltaS{"GT"}=-22.4; $pdeltaS{"GG"}=-19.9; $pdeltaS{"GC"}=-24.4; $pdeltaS{"GN"}=-24.4;
        $pdeltaS{"NA"}=-16.8; $pdeltaS{"NT"}=-21.5; $pdeltaS{"NG"}=-22.0; $pdeltaS{"NC"}=-21.0; $pdeltaS{"NN"}=-22.0;
        return(%pdeltaS);
} 

sub nearNeighbor{
        my ($seq,$hH,$hS,$pCon,$salt,$dsalt,$dntp)=@_;
        $pCon=$pCon*1e-9;
        $salt=($salt+120*sqrt(($dsalt - $dntp)))*1e-3;
        
        my $initH=0.6;
        my $initS=-9.0;
        my $symS=-1.4;

        my $useq=uc($seq);
        my @line=split //,$useq;
        
        my ($initH_1,$initH_2,$initS_1,$initS_2);
        if($line[0] eq 'G' || $line[0] eq 'C'){
                $initH_1=0.1;
                $initS_1=-2.8;
        }elsif($line[0] eq 'A' || $line[0] eq 'T'){
                $initH_1=2.3;
                $initS_1=4.1;
        }else{
                $initH_1=1.2;
                $initS_1=0.65;
        }
        
        if($line[-1] eq 'G' || $line[-1] eq 'C'){
                $initH_2=0.1;
                $initS_2=-2.8;
        }elsif($line[-1] eq 'A' || $line[-1] eq 'T'){
                $initH_2=2.3;
                $initS_2=4.1;
        }else{
                $initH_2=1.2;
                $initS_2=0.65;
        }
        
        my ($deltaH,$deltaS,$iDeltaH,$iDeltaS)=(0,0,0,0);
        
        for(my $i=1;$i<=$#line;$i++){
                if($line[$i-1] eq 'N' || $line[$i] eq 'N'){
                        next;
                }
                $deltaH+=$hH->{"${line[$i-1]}$line[$i]"};
                $deltaS+=$hS->{"${line[$i-1]}$line[$i]"};
        }
        my $tm;
        
        $deltaH=$deltaH+$initH_1+$initH_2;
        $deltaS=$deltaS+$initS_1+$initS_2;
        $deltaS=$deltaS+0.368*(length($seq)-1)*log($salt);
        
        $tm=$deltaH/($deltaS+1.987*log($pCon/4))*1000-273.15;
        return(sprintf("%.3f",$tm));
}

# --tm base basesalt nn
sub tmValue{
        my ($seq,$pCon,$salt,$dsalt,$dntp,$method)=@_;
        if($method eq "longsalt"){        
                &longTm($seq,$salt);
        }elsif($method eq "base"){
                &basicTm($seq);
        }elsif($method eq "basesalt"){
                &saltAdj($seq,$salt);
        }elsif($method eq "nn"){
                my %pdeltaH=&nnParaH();
                my %pdeltaS=&nnParaS();
                my $tm=&nearNeighbor($seq,\%pdeltaH,\%pdeltaS,$pCon,$salt,$dsalt,$dntp);
        }else{
                die "Error: unsupported method.\n";
        }

}

sub getPtm{
        my ($mixSeq,$pCon,$salt,$dsalt,$dntp,$method)=@_;
        #my $up_plas=reverse($mixSeq->[0]);
        #$up_plas=~tr/ATCGatcg/TAGCtagc/;
        #my $down_plas=reverse($mixSeq->[1]);
        #$down_plas=~tr/ATCGatcg/TAGCtagc/;
        my $tm_p1=&tmValue($mixSeq->[0],$pCon,$salt,$dsalt,$dntp,$method);
        my $tm_p2=&tmValue($mixSeq->[1],$pCon,$salt,$dsalt,$dntp,$method);
        return($tm_p1,$tm_p2);
}

sub pFlutter{
        my ($pcr,$len)=@_;
        my @tag=split /\;/,$pcr;
        my $pcr_out=0;
        my $pcr_in=0;
        foreach my $m(@tag){
                my $n=$m-$len;
                if($n<100 && $n>-100){
                        $pcr_in+=1;
                }
        }
        if($pcr_in == 1){
                $pcr_out=scalar(@tag);
        }else{
                $pcr_out=0;
        }
        return($pcr_out);
}

# %upPcr, %downPcr;
# $upPcr{$num} = [$pcr_product,$pcr_m_start] ,
sub tpBlastMerge{
        my ($r_up_primer,$r_down_primer,$gene_info,$plasSeq,$upPcr,$downPcr,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo) = @_;
        my (@group1,@group2,@group3,@group4);
        my @out;
        #my $tier;
        
        my ($geneName,$strand,$g1_g3_penalty,$g1_seq,$g1_start,$g1_len,$g1_tm,$g3_seq,$g3_start,$g3_len,$g3_tm,$g3_warn);
        my ($g4_g2_penalty,$g4_seq,$g4_start,$g4_len,$g4_tm,$g2_seq,$g2_start,$g2_len,$g2_tm,$g4_warn);
        my ($g1_loc_s,$g1_loc_e,$g2_loc_s,$g2_loc_e,$g3_loc_s,$g3_loc_e,$g4_loc_s,$g4_loc_e);
        my ($g1_loc,$g2_loc,$g3_loc,$g4_loc);
        my ($left_d_len,$right_d_len);
        my $up_n = scalar(@$r_up_primer);
        my $down_n = scalar(@$r_down_primer);
        my $outName;
        my ($chr,$gStart,$gEnd,$v_up_a) = @$gene_info;
        my ($up_rmix,$down_fmix)=@$plasSeq;
        my $num = 0;
        my $m_penalty;
        for(my $i = 0; $i < $up_n; $i++){
                if($upPcr->[$i]=~ /NA/){
                        next;
                }
                my @up_arr = split /\t+/,$r_up_primer->[$i];
                my ($geneName,$strand) = @up_arr[0,1];
                ($outName=$geneName)=~s/\|/\t/g;
                if($strand eq "plus"){
                        $strand = "+";
                        ($g1_g3_penalty,$g1_seq,$g1_start,$g1_len,$g1_tm,$g3_seq,$g3_start,$g3_len,$g3_tm) = @up_arr[3,4,5,6,7,13,14,15,16];
                        $g1_loc_s=$v_up_a+$g1_start-1;
                        $g1_loc_e=$g1_loc_s+$g1_len-1;
                        
                        $g3_loc_e=$v_up_a+$g3_start-1;
                        $g3_loc_s=$g3_loc_e-$g3_len+1;
                        
                        $g1_loc="${g1_loc_s}-${g1_loc_e}:+";
                        $g3_loc="${g3_loc_s}-${g3_loc_e}:-";
                        
                        $left_d_len=$g3_loc_e-$g1_loc_s+1;
                }else{
                        $strand = "-";
                        ($g1_g3_penalty,$g3_seq,$g3_start,$g3_len,$g3_tm,$g1_seq,$g1_start,$g1_len,$g1_tm) = @up_arr[3,4,5,6,7,13,14,15,16];
                        $g3_loc_s=$v_up_a+$g3_start-1;
                        $g3_loc_e=$g3_loc_s+$g3_len-1;
                        
                        $g1_loc_e=$v_up_a+$g1_start-1;
                        $g1_loc_s=$g1_loc_e-$g1_len+1;
                        
                        $g1_loc="${g1_loc_s}-${g1_loc_e}:-";
                        $g3_loc="${g3_loc_s}-${g3_loc_e}:+";
                        
                        $left_d_len=$g1_loc_e-$g3_loc_s+1;
                }
                my $left_r="$up_rmix" . "$g3_seq";
                my $fplasmid=reverse($left_r);
                $fplasmid=~tr/ATCGatcg/TAGCtagc/;
                
                my $up_wstat = &pFlutter($upPcr->[$i],$left_d_len);
                next if($up_wstat==0);
                if(scalar(@up_arr)==23){
                        $g3_warn = $up_arr[22];
                }else{
                        $g3_warn = "";
                }
                #
                $left_d_len += length($up_rmix);
                #
                for(my $j = 0; $j < $down_n; $j++){
                        if($downPcr->[$j] =~ /NA/){
                                next;
                        }
                        my @down_arr = split /\t+/,$r_down_primer->[$j];
                        if($strand eq "+"){
                                ($g4_g2_penalty,$g4_seq,$g4_start,$g4_len,$g4_tm,$g2_seq,$g2_start,$g2_len,$g2_tm) = @down_arr[3,4,5,6,7,13,14,15,16];
                                $g4_loc_s=$v_up_a+$g4_start-1;
                                $g4_loc_e=$g4_loc_s+$g4_len-1;
                                
                                $g2_loc_e=$v_up_a+$g2_start-1;
                                $g2_loc_s=$g2_loc_e-$g2_len+1;
                                
                                $g4_loc="${g4_loc_s}-${g4_loc_e}:+";
                                $g2_loc="${g2_loc_s}-${g2_loc_e}:-";
                                
                                $right_d_len=$g2_loc_e-$g4_loc_s+1;
                                
                        }else{
                                ($g4_g2_penalty,$g2_seq,$g2_start,$g2_len,$g2_tm,$g4_seq,$g4_start,$g4_len,$g4_tm) = @down_arr[3,4,5,6,7,13,14,15,16];
                                $g2_loc_s=$v_up_a+$g2_start-1;
                                $g2_loc_e=$g2_loc_s+$g2_len-1;
                                
                                $g4_loc_e=$v_up_a+$g4_start-1;
                                $g4_loc_s=$g4_loc_e-$g4_len+1;
                                
                                $g4_loc="${g4_loc_s}-${g4_loc_e}:-";
                                $g2_loc="${g2_loc_s}-${g2_loc_e}:+";
                                
                                $right_d_len=$g4_loc_e-$g2_loc_s+1;
                                
                        }
                        my $down_wstat = &pFlutter($downPcr->[$j],$right_d_len);
                        next if($down_wstat == 0);
                        
                        my $right_f="$down_fmix" . "$g4_seq";
                        my $rplasmid=reverse($right_f);
                        $rplasmid=~tr/ATCGatcg/TAGCtagc/;
                        #
                        $right_d_len += length($down_fmix);
                        #
                        if($g1_tm - $g2_tm <= 6 && $g1_tm - $g2_tm >= -6){
                                my $g1_g2_dimer;
                                if($thermo){
                                        $g1_g2_dimer = &th_dimer($ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$g1_seq,$g2_seq);
                                }else{
                                        $g1_g2_dimer = &dimer($g1_seq,$g2_seq,$dERegion,$dEMatch);
                                }
                                my $part1 = "$outName\t$strand\t$chr\t$gStart\t$gEnd\t$v_up_a";
                                #my $part1 = "$outName\t$strand\t$chr\t$gStart\t$gEnd";
                                my $part2 = "$g1_seq\t$left_r\t$right_f\t$g2_seq\t$fplasmid\t$rplasmid\t$g1_loc\t$g3_loc\t$g4_loc\t$g2_loc\t$g1_tm\t$g3_tm\t$g4_tm\t$g2_tm\t$p1_tm\t$p2_tm\t$up_wstat\t$upPcr->[$i]\t$left_d_len\t$down_wstat\t$downPcr->[$j]\t$right_d_len\t$g1_g2_dimer\t$g1_g3_penalty\t$g4_g2_penalty";
                                my $warn = "";
                                if($force){
                                        if(scalar(@down_arr)==23){
                                                $g4_warn = $down_arr[22];
                                        }else{
                                                $g4_warn = "";
                                        }
                                        if($g3_warn eq ""){
                                                if($g4_warn eq ""){
                                                        $warn = "None";
                                                }else{
                                                        $warn = "G4:$g4_warn";
                                                }
                                        }else{
                                                if($g4_warn eq ""){
                                                        $warn = "G3:$g3_warn";
                                                }else{
                                                        $warn = "G3:$g3_warn G4:$g4_warn";
                                                }
                                        }
                                }
                                
                                my $m_penalty = ($g1_g3_penalty + $g4_g2_penalty) / 2;                                 
                                if($g1_g2_dimer){
                                        # group2, group4
                                        # 
                                        if($up_wstat > 1){
                                                push(@group4,[$part1,$part2,4,$warn,$m_penalty]);
                                        }elsif($down_wstat > 1){
                                                push(@group4,[$part1,$part2,4,$warn,$m_penalty]);
                                        }else{
                                                push(@group2,[$part1,$part2,2,$warn,$m_penalty]);
                                        }
                                        
                                }else{
                                        # group1, group3
                                        if($up_wstat > 1){
                                                push(@group3,[$part1,$part2,3,$warn,$m_penalty]);
                                        }elsif($down_wstat > 1){
                                                push(@group3,[$part1,$part2,3,$warn,$m_penalty]);
                                        }else{
                                                push(@group1,[$part1,$part2,1,$warn,$m_penalty]);
                                        }
                                }
                        }
                }
        }
        #
        my $t = 0;
        my @allprimer = (\@group1,\@group2,\@group3,\@group4);
        my $flag = 0;
        #my @pos = (0,3,4,5);
        my @pos = (0,1,2,3);
        foreach my $k(@allprimer){
                if($flag == 1){
                        last;
                }
                
                my @sp;
                if($refine){
                        foreach my $x(@$k){
                                my $seqinfo = $x->[1];
                                my @xinfo = split /\t/,$seqinfo;
                                
                                my $taBase = 0;
                                foreach my $p(@pos){
                                        my $seq1 = $xinfo[$p];
                                        my $base1 = substr($seq1,-1);
                                        if($base1 eq 'A' || $base1 eq 'T' || $base1 eq 'a' || $base1 eq 't'){
                                                $taBase++;
                                        }
                                }
                                
                                push(@$x,$taBase);
                        }
                        #
                        @sp = sort {$a->[5]<=>$b->[5] or $a->[4]<=>$b->[4]} @$k;
                }else{
                        @sp = sort {$a->[4]<=>$b->[4]} @$k;
                }
                
                foreach my $f(@sp){
                        $t++;
                        if($force){
                                push(@out,"$f->[0]\t$t\t$f->[1]\t$f->[2]\t$f->[3]");
                        }else{
                                push(@out,"$f->[0]\t$t\t$f->[1]\t$f->[2]");
                        }
                        if($t == $pNum){
                                $flag = 1;
                                last;
                        }
                }
                
        }
        return(@out);
}

sub tpMergeInfo{
        my ($t_dir,$f_format,$r_format,$up_sta,$down_sta,$info,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo)=@_;        
        my $r=0;
        if(! -d "$t_dir/primer_result"){
                $r=system("mkdir $t_dir/primer_result");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_result\n";
                }
        }
        my %common = &tcommon($up_sta,$down_sta);
        my $up_mstart="$t_dir/primer_blast/target.up.length";
        my $down_mstart="$t_dir/primer_blast/target.down.length";
        my ($up_rmix,$down_fmix)=@$plasSeq;
        my $up_len=length($up_rmix);
        my $down_len=length($down_fmix);
        
        (open OUT,">$t_dir/primer_result/target.primer.txt") || die "$!\n";
        (open IN,"$up_mstart") || die "$!\n";
        (open BN,"$down_mstart") || die "$!\n";
        
        (open EN,"$f_format") || die "$!\n";
        (open FN,"$r_format") || die "$!\n";
        <IN>;<BN>;
        my $f_head = <EN>;
        chomp($f_head);
        my $force = 0;
        if($f_head =~ /Warning$/){
                print OUT "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tSeq_Start\tPrimer_ID\tG1\tG3\tG4\tG2\tP1\tP2\tG1_Loc\tG3_Loc\tG4_Loc\tG2_Loc\tG1_Tm\tG3_Tm\tG4_Tm\tG2_Tm\tP1_Tm\tP2_Tm\tG1_G3_Products_Counts\tG1_G3_Blast_Len\tG1_G3_Expect_Len\tG4_G2_Products_Counts\tG4_G2_Blast_Len\tG4_G2_Expect_Len\tG1_G2_Dimer\tG1_G3_Penalty\tG4_G2_Penalty\tTier\tPrimer_Warning\n";
                $force = 1;
        }else{
                print OUT "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tSeq_Start\tPrimer_ID\tG1\tG3\tG4\tG2\tP1\tP2\tG1_Loc\tG3_Loc\tG4_Loc\tG2_Loc\tG1_Tm\tG3_Tm\tG4_Tm\tG2_Tm\tP1_Tm\tP2_Tm\tG1_G3_Products_Counts\tG1_G3_Blast_Len\tG1_G3_Expect_Len\tG4_G2_Products_Counts\tG4_G2_Blast_Len\tG4_G2_Expect_Len\tG1_G2_Dimer\tG1_G3_Penalty\tG4_G2_Penalty\tTier\n";
        }
        <FN>;
        
        my $pre_up = "";
        my @up_p;
        my @down_p;
        my ($flag,$flag2) = (0,0);
        my @geneInfo;
        my $start = 0;
        my $end = scalar(@$info);
        my (@upPcr,@downPcr);
        my ($up_str,$down_str);
        my ($up_gene,$down_gene,$up_line,$down_line);
        my $fileEnd = 0;
        my (@upLen,@downLen);
        print "\nGenes filtered:\n";
        while(<EN>){
                chomp;
                $up_line = $_;
                $up_gene = (split /\s+/,$_)[0];
                #
                my $up_len_line=<IN>;
                chomp($up_len_line);
                @upLen=split /\s+/,$up_len_line;
                if($up_gene ne $pre_up){
                        if($flag == 1){
                                @down_p = ();
                                #
                                @downPcr = ();
                                $flag2 = 0;
                                if($down_gene eq $pre_up){
                                        push(@down_p,$down_line);
                                        #
                                        $down_str = &chPcrLen($downLen[3],$down_len);
                                        push(@downPcr,$down_str);
                                        $flag2 = 1;
                                }
                                ##
                                $fileEnd = 0;
                                while(1){
                                        $down_line = <FN>;
                                        if(! $down_line){
                                                $fileEnd = 1;
                                                last;
                                        }
                                        chomp($down_line);
                                        $down_gene = (split /\s+/,$down_line)[0];
                                        
                                        #
                                        my $down_len_line=<BN>;
                                        chomp($down_len_line);
                                        @downLen=split /\s+/,$down_len_line;
                                        #
                                        if($down_gene eq $pre_up){
                                                push(@down_p,$down_line);
                                                #
                                                $down_str = &chPcrLen($downLen[3],$down_len);
                                                push(@downPcr,$down_str);
                                                $flag2 = 1;
                                        }else{
                                                if($flag2){
                                                        for(my $k = $start; $k < $end; $k++){
                                                                my $tline = $info->[$k];
                                                                my @allInfo = split /\,/,$tline;
                                                                my $geneName = $allInfo[0];
                                                                if($geneName eq $pre_up){
                                                                        @geneInfo = @allInfo[1,6,7,8];
                                                                        $start = $k + 1;
                                                                        last;
                                                                }
                                                        }
                                                        
                                                        my @tout = &tpBlastMerge(\@up_p,\@down_p,\@geneInfo,$plasSeq,\@upPcr,\@downPcr,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                                                        if(scalar(@tout) > 0){
                                                                foreach my $s(@tout){
                                                                        print OUT "$s\n";
                                                                }
                                                        }else{
                                                                print "$pre_up\n";
                                                        }
                                                        last;
                                                }
                                        }
                                        
                                }
                                if($fileEnd){
                                        if($flag2){
                                                for(my $k = $start; $k < $end; $k++){
                                                        my $tline = $info->[$k];
                                                        my @allInfo = split /\,/,$tline;
                                                        my $geneName = $allInfo[0];
                                                        if($geneName eq $pre_up){
                                                                @geneInfo = @allInfo[1,6,7,8];
                                                                $start = $k + 1;
                                                                last;
                                                        }
                                                }
                                                
                                                my @tout = &tpBlastMerge(\@up_p,\@down_p,\@geneInfo,$plasSeq,\@upPcr,\@downPcr,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                                                if(scalar(@tout) > 0){
                                                        foreach my $s(@tout){
                                                                print OUT "$s\n";
                                                        }
                                                }else{
                                                        print "$pre_up\n";
                                                }
                                                #
                                                $flag = 0;
                                        }
                                        last;
                                }
                        }
                        #
                        if(exists $common{$up_gene}){
                                $flag = 1;
                                @up_p = ();
                                push(@up_p,$up_line);
                                #########
                                @upPcr = ();
                                $up_str = &chPcrLen($upLen[3],$up_len);
                                push(@upPcr,$up_str);
                        }else{
                                $flag = 0;
                        }
                }else{
                        if($flag == 1){
                                push(@up_p,$up_line);
                                #########
                                $up_str = &chPcrLen($upLen[3],$up_len);
                                push(@upPcr,$up_str);
                        }
                
                }
                ######
                $pre_up = $up_gene;
        }
        
        if($flag == 1){
                @down_p = ();
                #
                @downPcr = ();
                $flag2 = 0;
                if($down_gene eq $pre_up){
                        push(@down_p,$down_line);
                        #
                        $down_str = &chPcrLen($downLen[3],$down_len);
                        push(@downPcr,$down_str);

                        $flag2 = 1;
                }
                ##
                $fileEnd = 0;
                while(1){
                        $down_line = <FN>;
                        if(! $down_line){
                                $fileEnd = 1;
                                last;
                        }
                        chomp($down_line);
                        $down_gene = (split /\s+/,$down_line)[0];
                        #
                        my $down_len_line=<BN>;
                        chomp($down_len_line);
                        @downLen=split /\s+/,$down_len_line;
                        #
                        if($down_gene eq $pre_up){
                                push(@down_p,$down_line);
                                #
                                $down_str = &chPcrLen($downLen[3],$down_len);
                                push(@downPcr,$down_str);

                                $flag2 = 1;
                        }else{
                                if($flag2){
                                        for(my $k = $start; $k < $end; $k++){
                                                my $tline = $info->[$k];
                                                my @allInfo = split /\,/,$tline;
                                                my $geneName = $allInfo[0];
                                                if($geneName eq $pre_up){
                                                        @geneInfo = @allInfo[1,6,7,8];
                                                        $start = $k + 1;
                                                        last;
                                                }
                                        }
                                        
                                        my @tout = &tpBlastMerge(\@up_p,\@down_p,\@geneInfo,$plasSeq,\@upPcr,\@downPcr,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                                        if(scalar(@tout) > 0){
                                                foreach my $s(@tout){
                                                        print OUT "$s\n";
                                                }
                                        }else{
                                                print "$pre_up\n";
                                        }
                                        last;
                                }
                        }
                }
                if($fileEnd){
                        if($flag2){
                                for(my $k = $start; $k < $end; $k++){
                                        my $tline = $info->[$k];
                                        my @allInfo = split /\,/,$tline;
                                        my $geneName = $allInfo[0];
                                        if($geneName eq $pre_up){
                                                @geneInfo = @allInfo[1,6,7,8];
                                                $start = $k + 1;
                                                last;
                                        }
                                }
                                
                                my @tout = &tpBlastMerge(\@up_p,\@down_p,\@geneInfo,$plasSeq,\@upPcr,\@downPcr,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$force,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                                if(scalar(@tout) > 0){
                                        foreach my $s(@tout){
                                                print OUT "$s\n";
                                        }
                                }else{
                                        print "$pre_up\n";
                                }
                        }
                }
        }
        close IN;
        close BN;
        close EN;
        close FN;
        close OUT;
}

sub tpGeneMerge{
        my ($t_dir,$geneName,$plasSeq,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$dERegion,$dEMatch,$info,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo)=@_;
        my $f_format="$t_dir/primer_blast/$geneName.target.up.format.out";
        my $r_format="$t_dir/primer_blast/$geneName.target.down.format.out";
        my $up_sta = "$t_dir/primer_blast/$geneName.target.up.stat";
        my $down_sta = "$t_dir/primer_blast/$geneName.target.down.stat";
        
        #my $upMix="$t_dir/primer_blast/$geneName.lmix.fa";
        #my $downMix="$t_dir/primer_blast/$geneName.rmix.fa";
        &tpUpFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        &tpDownFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        ######
        &tpMergeInfo($t_dir,$f_format,$r_format,$up_sta,$down_sta,$info,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
}

sub tpBundleMerge{
        my ($t_dir,$plasSeq,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$dERegion,$dEMatch,$info,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo)=@_;
        my $f_format="$t_dir/primer_blast/bundle.target.up.format.out";
        my $r_format="$t_dir/primer_blast/bundle.target.down.format.out";
        my $up_sta = "$t_dir/primer_blast/bundle.target.up.stat";
        my $down_sta = "$t_dir/primer_blast/bundle.target.down.stat";
        
        #my $upMix="$t_dir/primer_blast/bundle.lmix.fa";
        #my $downMix="$t_dir/primer_blast/bundle.rmix.fa";
        &tpUpFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        &tpDownFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        &tpMergeInfo($t_dir,$f_format,$r_format,$up_sta,$down_sta,$info,$plasSeq,$dERegion,$dEMatch,$pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
}

sub boolstart{
        my ($blastLine,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast)=@_;
        my $misStart=1;
        my $schr=$blastLine->[1];
        my $allignLen=$blastLine->[3];
        my $mismatch=$blastLine->[4];
        my $gap=$blastLine->[5];
        my $qstart=$blastLine->[6];
        my $qend=$blastLine->[7];
        my $sstart=$blastLine->[8];
        my $queryLen=$blastLine->[12];
        my $mapInfo=$blastLine->[14];
        
        my $unAlignLen;
        #$gap=($mapInfo=~tr/-//);
        my $etrunc=$queryLen-$qend;
        $unAlignLen=$mismatch+$gap+$etrunc+$qstart-1;
        if($unAlignLen>=$t_mismatchMore){
                $misStart=0;
        }else{
                if($etrunc >= $t_endMismatch){
                        $misStart=0;
                }else{
                        my @info=split //,$mapInfo;
                        #my $preStr=shift(@info);
                        my @mismatchSite;
                        my @gapSite;
                        my @insertSite;
                        my $tpos=$qstart;
                        my $outStr="";
                        my $nchar;
                        my $n=scalar(@info);
                        my $k=0;
                        my ($fir,$sec);
                        while(1){
                                if($info[$k]=~/\d/){
                                        $outStr="$outStr" . "$info[$k]";
                                }else{
                                        if($outStr ne ''){
                                                my $value=int($outStr);
                                                $tpos=$tpos+$value-1;
                                        }
                                        $fir=$info[$k];
                                        $k++;
                                        $sec=$info[$k];
                                        $outStr='';
                                        if($fir eq '-'){
                                                $tpos++;
                                                #-A
                                                push(@gapSite,$tpos);
                                        }else{
                                                if($sec eq '-'){
                                                        #A-
                                                        $tpos++;
                                                        push(@insertSite,$tpos);
                                                }else{
                                                        #AG
                                                        $tpos++;
                                                        push(@mismatchSite,$tpos);
                                                }
                                        }
                                }
                                $k++;
                                last if($k==$n);
                        }
                        #if($outStr ne ''){
                        #        my $value=int($outStr);
                        #        $tpos=$tpos+$value-1;
                        #}    
                        
                        my $pre_len=$queryLen-$t_endRegion;
                        my $endMis=$etrunc;
                        foreach my $j(@mismatchSite){
                                if($j>$pre_len){
                                        $endMis++;
                                }
                        }
                        
                        foreach my $j(@gapSite){
                                if($j>$pre_len){
                                        $endMis++;
                                }
                        }
                        
                        foreach my $j(@insertSite){
                                if($j>$pre_len){
                                        $endMis++;
                                }
                        }
                        #
                        if($endMis>=$t_endMismatch){ 
                                if($unAlignLen>=$t_mismatchAtLeast){
                                        $misStart=0;
                                }   
                        }
                }
        }
        return($misStart);
}

# $t_productMinLen $t_productMaxLen
sub hitcompare{
        my ($fhit_arr,$rhit_arr,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen)=@_;
        my ($fchr,$fstart,$rchr,$rend);
        my $productLen;
        my @out;
        #my ($x,$y)=(-1,-1);
        #my $min_size=$t_productMinLen-100;
        my $max_size=$t_productMaxLen+3000;
        foreach my $i(@$fhit_arr){
                #$x++;
                if(&boolstart($i,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast)){
                        $fchr=$i->[1];
                        $fstart=$i->[8];
                        
                }else{
                        next;
                }
                #
                #$y=-1;
                foreach my $j(@$rhit_arr){
                        #$y++;
                        if(&boolstart($j,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast)){
                                $rchr=$j->[1];
                                $rend=$j->[8];
                                if($fchr eq $rchr){
                                        $productLen=$rend-$fstart+1;
                                        #if($productLen>=$min_size && $productLen<=$max_size){
                                        if($productLen>0 && $productLen<=$max_size){
                                                #push(@out,"$x" . "_" . "$y" . "_" . "$productLen");
                                                push(@out,$productLen);
                                        }
                                }
                        }
                }  
        }
        return(@out);
}

# verify primer length
sub vPairFilt{
        # mix primer exist for all gene
        my ($t_dir,$fblast,$rblast,$outFile,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen)=@_;
        (open IN,"$rblast") || die "$!\n";
        (open AN,"$fblast") || die "$!\n";
        # $t_dir/primer_blast/up.target.length
        (open OUT,">$outFile") || die "$!\n";
        print OUT "#GeneName\tF_name\tR_name\tLength\tCounts\n";
           
        my ($fgeneName,$rgeneName);
        my (@fhit,@rhit);
        my $indicate=0;
        my $indicateB=0;
        my ($fpName,$rpName);
        while(<AN>){
                chomp;
                my $xtag=$_;
                if($xtag=~/^# BLASTN/){
                        $indicate++;
                        if($indicate == 1){
                                @fhit=();
                                next;
                        }else{
                                @rhit=();
                                while(<IN>){
                                        my $rtag=$_;
                                        chomp($rtag);
                                        if($rtag=~/^# BLASTN/){
                                                $indicateB++;
                                                if($indicateB==1){
                                                        @rhit=();
                                                        next;
                                                }else{ 
                                                        
                                                        if(scalar(@rhit)<1 || scalar(@fhit)<1){
                                                                $indicateB=1;
                                                                print OUT "$fgeneName\t$fpName\t$rpName\tNA\tNA\n";
                                                                last;
                                                        }
                                                        #
                                                        my @productPCR=&hitcompare(\@fhit,\@rhit,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
                                                        my $nproduct=scalar(@productPCR);
                                                        if($nproduct>0){
                                                                my $tp=join ';',@productPCR;
                                                                print OUT "$fgeneName\t$fhit[0][0]\t$rhit[0][0]\t$tp\t$nproduct\n";
                                                        }else{
                                                                print OUT "$fgeneName\t$fpName\t$rpName\tNA\tNA\n";
                                                        }
                                                        $indicateB=1;
                                                        last;
                                                }
                                        }
                                        #
                                        if($indicateB == 1){
                                                if($rtag=~/^# Query: (.+\|.+\|.+?)_/){
                                                        $rgeneName=$1;
                                                        $rpName=substr($rtag,9);
                                                }else{
                                                        next if($rtag=~/^#/);
                                                        my @rline=split /\s+/,$rtag;
                                                        push(@rhit,\@rline);
                                                }        
                                        }
                                }
                                $indicate=1;
                                @fhit=();
                        }
                }
                if($indicate == 1){
                        if($xtag=~/^# Query: (.+\|.+\|.+?)_/){
                                $fgeneName=$1;
                                $fpName=substr($xtag,9);
                        }else{
                                next if($xtag=~/^#/);
                                my @fline=split /\s+/,$xtag;
                                push(@fhit,\@fline);
                        }
                        
                }
        }
        #
        @rhit=();
        while(<IN>){
                my $rtag=$_;
                chomp($rtag);
                if($rtag=~/^# Query: (.+\|.+\|.+?)_/){
                        $rgeneName=$1;
                        $rpName=substr($rtag,9);
                }else{
                        next if($rtag=~/^#/);
                        my @rline=split /\s+/,$rtag;
                        push(@rhit,\@rline);
                }        
        }
        if(scalar(@rhit)<1 || scalar(@fhit)<1){
                print OUT "$fgeneName\t$fpName\t$rpName\tNA\tNA\n";
        }else{
                my @productPCR=&hitcompare(\@fhit,\@rhit,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
                my $nproduct=scalar(@productPCR);
                if($nproduct>0){
                        my $tp=join ';',@productPCR;
                        print OUT "$fgeneName\t$fhit[0][0]\t$rhit[0][0]\t$tp\t$nproduct\n";
                }else{
                        print OUT "$fgeneName\t$fpName\t$rpName\tNA\tNA\n";
                }
        }
        
        close IN;
        close AN;
        close OUT;
}

sub alignScore{
        # left 1, up 2, diag 3; a 1, b 2, c 3; d 12, e 13, f 23; g 123
        my ($query,$subject,$direction,$match,$mismatch,$gap)=@_;
        if($direction eq "-"){
                $query=reverse($query);
                $query=~tr/ATCGatcg/TAGCtagc/;
        }

        my @scoreMat;
        #my @ori;
        $query=uc($query);
        $subject=uc($subject);
        my $queryLen=length($query);
        my $subjectLen=length($subject);
        my @qArray=split //,$query;
        my @sArray=split //,$subject;
        my @fir=map {$_ * $gap * (-1)} (0..$subjectLen);
        push(@scoreMat,\@fir);
        my ($fromLeft,$fromUp,$fromDiag);
        for(my $i=1;$i<=$queryLen;$i++){
                my @out;
                push(@out,$i*$gap*(-1));
                push(@scoreMat,\@out);
                
                #my @dline;
                #push(@dline,'b');
                #push(@ori,\@dline);
                for(my $j=1;$j<=$subjectLen;$j++){
                        $fromLeft=$scoreMat[$i][$j-1]-$gap;
                        $fromUp=$scoreMat[$i-1][$j]-$gap;
                        if($qArray[$i-1] eq "N"){
                                if($sArray[$i-1] eq "N"){
                                        $fromDiag=$scoreMat[$i-1][$j-1];
                                }else{
                                        $fromDiag=$scoreMat[$i-1][$j-1]-$gap;
                                }
                        }else{
                                if($sArray[$i-1] eq "N"){
                                        $fromDiag=$scoreMat[$i][$j]-$gap;
                                }else{
                                        if($qArray[$i-1] eq $sArray[$j-1]){
                                                $fromDiag=$scoreMat[$i-1][$j-1]+$match;
                                        }else{
                                                $fromDiag=$scoreMat[$i-1][$j-1]-$mismatch;
                                        }
                                }
                        }
                        my $max=($fromLeft>$fromUp)?($fromLeft):($fromUp);
                        $max=($max>$fromDiag)?($max):($fromDiag);  
                        push(@{$scoreMat[$i]},$max);
                }
        }
        #
        return($scoreMat[$queryLen][$subjectLen]);
}

sub dimer{
        my ($seq1,$seq2,$regionLen,$matchLen)=@_;
        $seq1=uc($seq1);
        $seq2=uc($seq2);
        
        # long $seq1 5'-3'; short $seq2 3'-5'
        my $len1=length($seq1);
        my $len2=length($seq2);
        if($len1<$len2){
                my $tseq=$seq1;
                $seq1=$seq2;
                $seq2=$tseq;
                my $tlen=$len1;
                $len1=$len2;
                $len2=$tlen;
        }
        
        if($len2<$matchLen){
                return(0);
        }
        $seq2=reverse($seq2);
        $seq2=~tr/ATCG/TAGC/;
        
        my $cirLen=$len2;
        if($len2<$regionLen){
                $regionLen=$len2;
        }
        my $stop1=$len1-$regionLen;
        my $start1=$stop1;
        my $sStop=$len1-$matchLen;
        my $start2=$len2-$regionLen-1;
        my $start3=$regionLen-1;
        my @arr1=split //,$seq1;
        my @arr2=split //,$seq2;
        my ($match1,$match2)=(0,0);
        my ($mBase1,$mBase2)=(0,0);
        my $base1=$arr1[$len1-1];
        my $base2=$arr2[0];
        
        for(my $i=1;$i<=$sStop;$i++){
                $match1=0;
                $match2=0;
                $mBase1=0;
                $mBase2=0;
                my $k=$i;
                for(my $j=0;$j<$regionLen;$j++){
                        if($arr1[$k] eq $arr2[$j]){
                                $match2++;
                        }
                        $k++;
                        if($k == $len1){last}
                }
                if($arr1[$i] eq $base2){
                        $mBase2=1;
                }
                $cirLen+=$i;
                my $q=$start2;
                if($i<$stop1){
                        if($cirLen>$len1){
                                for(my $p=$start1;$p<$len1;$p++){
                                        if($arr1[$p] eq $arr2[$q]){
                                                $match1++;
                                        }
                                        $q++;
                                }
                                $start2--;
                                if($arr2[$q-1] eq $base1){
                                        $mBase1=1;
                                }
                        }
                        
                        
                }else{
                        $match1=$match2;
                        if($arr2[$start3] eq $base1){
                                $mBase1=1;
                        }
                        $start3--;
                }
                
                if(($match1>=$matchLen && $mBase1>0) || ($match2>=$matchLen && $mBase2>0)){
                        #print "@arr1[$stop1..($len1-1)] # @arr2[0..($regionLen-1)]\n";
                        #print "@arr1[($i)..($k-1)] ## @arr2[0..($regionLen-1)] ||| @arr1[$stop1..($len1-1)] # @arr2[($start2+1)..($q-1)]\n";
                        #print join "",@arr1,"\n";
                        #print " " x $i;
                        #print join "",@arr2,"\n";
                        return(1);
                }
        }
        
        return(0);
}

sub th_dimer{
        my ($ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$seq1,$seq2)=@_;
        my $command="$ntthal -s1 $seq1 -s2 $seq2 -r -mv $salt -dv $dsalt -n $dntp -d $pCon";
        my $line=`$command`;
        if($line <= $d_max_tm){
                return(0);
        }else{
                return(1);
        }
}

sub vloc{
        my ($v_up_a,$start_f,$start_r,$len_f,$len_r)=@_;
        
        my $v1_loc_s=$v_up_a+$start_f-1;
        my $v1_loc_e=$v1_loc_s+$len_f-1;
        
        my $v2_loc_e=$v_up_a+$start_r-1;
        my $v2_loc_s=$v2_loc_e-$len_r+1;
        
        my $v1_loc="${v1_loc_s}-${v1_loc_e}:+";
        my $v2_loc="${v2_loc_s}-${v2_loc_e}:-";
        return($v1_loc,$v2_loc);

}

#######################################
#  short arm
#######################################


sub s_tpStrandTrans{
        my ($t_dir,$upMix,$downMix,$plasSeq,$info,$p1_tm,$p2_tm)=@_;        
        my $r=0;
        if(! -d "$t_dir/primer_result"){
                $r=system("mkdir $t_dir/primer_result");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_result\n";
                }
        }
        
        my ($up_rmix,$down_fmix)=@$plasSeq;
        
        (open GN,"$upMix") || die "$!\n";
        (open HN,"$downMix") || die "$!\n";
        (open OUT,">$t_dir/primer_result/target.primer.unblast.txt") || die "$!\n";
        print OUT "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tSeq_Start\tP1\tP2\tP1_Tm\tP2_Tm\n";
        my ($geneName,$tag);
        my ($up_mix_fir,$up_mix_sec,$down_mix_fir,$down_mix_sec);
        my ($left_r,$right_f);
        my ($fplasmid,$rplasmid);
        my $infoName="";
        my ($chr,$gStart,$gEnd,$v_up_a,$strand);
        my ($tm_p1,$tm_p2);
        
        my $flag = 1;
        while($up_mix_fir=<GN>){
                $up_mix_sec=<GN>;
                $down_mix_fir=<HN>;
                $down_mix_sec=<HN>;
                chomp($up_mix_fir);
                chomp($up_mix_sec);
                chomp($down_mix_fir);
                chomp($down_mix_sec);
                if($up_mix_fir=~/>(.+)_lmix/){
                        $geneName=$1;
                        ($tag=$geneName)=~s/\|/\t/g;
                }
                #
                $flag = 1;
                while($geneName ne $infoName){
                        if(scalar(@{$info})>0){
                                my $geneTag=shift(@{$info});
                                my @geneInfo=split /\,/,$geneTag;
                                ($infoName,$chr,$gStart,$gEnd,$v_up_a,$strand)=@geneInfo[0,1,6,7,8,13];
                        }else{
                                $flag = 0;
                                last;
                        }
                }
                if($flag){
                    if($strand eq "+"){
                            $left_r="$up_rmix" . "$up_mix_sec";
                            $right_f="$down_fmix" . "$down_mix_sec";
                            
                            $fplasmid=reverse($left_r);
                            $fplasmid=~tr/ATCGatcg/TAGCtagc/;
                            $rplasmid=reverse($right_f);
                            $rplasmid=~tr/ATCGatcg/TAGCtagc/;
                            
                            #$v_up_a\t
                            print OUT "$tag\t$strand\t$chr\t$gStart\t$gEnd\t$v_up_a\t$fplasmid\t$rplasmid\t$p1_tm\t$p2_tm\n";
                    }else{
                            $left_r="$up_rmix" . "$down_mix_sec";
                            $right_f="$down_fmix" . "$up_mix_sec";
                            
                            $fplasmid=reverse($left_r);
                            $fplasmid=~tr/ATCGatcg/TAGCtagc/;
                            $rplasmid=reverse($right_f);
                            $rplasmid=~tr/ATCGatcg/TAGCtagc/;
                            
                            print OUT "$tag\t$strand\t$chr\t$gStart\t$gEnd\t$v_up_a\t$fplasmid\t$rplasmid\t$p1_tm\t$p2_tm\n";
                    }
                }
        }
        close GN;
        close HN;
        close OUT;
}

sub s_tpGeneTrans{
        my ($t_dir,$geneName,$plasSeq,$info,$p1_tm,$p2_tm)=@_;
        my $upMix="$t_dir/primer_blast/$geneName.lmix.fa";
        my $downMix="$t_dir/primer_blast/$geneName.rmix.fa";
        &s_tpStrandTrans($t_dir,$upMix,$downMix,$plasSeq,$info,$p1_tm,$p2_tm);
}

sub s_tpBundleTrans{
        my ($t_dir,$plasSeq,$info,$p1_tm,$p2_tm)=@_;
        my $upMix="$t_dir/primer_blast/bundle.lmix.fa";
        my $downMix="$t_dir/primer_blast/bundle.rmix.fa";        
        &s_tpStrandTrans($t_dir,$upMix,$downMix,$plasSeq,$info,$p1_tm,$p2_tm);
}

sub s_tpGeneMerge{
        my ($t_dir,$geneName,$upAlign,$downAlign,$identity,$plasSeq,$info,$p1_tm,$p2_tm)=@_;
        my $r=0;
        if(! -d "$t_dir/primer_result"){
                $r=system("mkdir $t_dir/primer_result");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_result\n";
                }
        }
        my $fblast="$t_dir/primer_blast/blastn.rmix.out";
        my $rblast="$t_dir/primer_blast/blastn.lmix.out";
        my $fFa="$t_dir/primer_blast/$geneName.rmix.fa";
        my $rFa="$t_dir/primer_blast/$geneName.lmix.fa";
        my $outFile="$t_dir/primer_result/target.primer.txt";
        &s_tpMergeInfo($fblast,$rblast,$fFa,$rFa,$outFile,$upAlign,$downAlign,$identity,$plasSeq,$info,$p1_tm,$p2_tm);
}

sub s_tpBundleMerge{
        my ($t_dir,$upAlign,$downAlign,$identity,$plasSeq,$info,$p1_tm,$p2_tm)=@_;
        my $r=0;
        if(! -d "$t_dir/primer_result"){
                $r=system("mkdir $t_dir/primer_result");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_result\n";
                }
        }
        my $fblast="$t_dir/primer_blast/blastn.rmix.out";
        my $rblast="$t_dir/primer_blast/blastn.lmix.out";
        my $fFa="$t_dir/primer_blast/bundle.rmix.fa";
        my $rFa="$t_dir/primer_blast/bundle.lmix.fa";
        my $outFile="$t_dir/primer_result/target.primer.txt";
        &s_tpMergeInfo($fblast,$rblast,$fFa,$rFa,$outFile,$upAlign,$downAlign,$identity,$plasSeq,$info,$p1_tm,$p2_tm);
}

sub s_tpMergeInfo{
        # mix primer exist for all gene
        my ($fblast,$rblast,$fFa,$rFa,$outFile,$upAlign,$downAlign,$identity,$plasSeq,$info,$p1_tm,$p2_tm)=@_;
        (open IN,"$rblast") || die "file contain P1 (reverse complement) blast result open failed. $!\n";
        (open AN,"$fblast") || die "file contain P2 (reverse complement) blast result open failed. $!\n";
        # $t_dir/primer_blast/up.target.length
        (open BN,"$rFa") || die "file contain P1 (reverse complement) sequence open failed. $!\n";
        (open DN,"$fFa") || die "file contain P2 (reverse complement) sequence open failed. $!\n";
        (open OUT,">$outFile") || die "$!\n";
        print OUT "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tSeq_Start\tP1\tP2\tP1_Tm\tP2_Tm\tP1_Blast_Match\tP2_Blast_Match\n";
        
        my ($fhit,$rhit)=(0,0);
        my $indicate=0;
        my $indicateB=0;
        my ($rgeneName,$fgeneName);
        my ($up_rmix,$down_fmix)=@$plasSeq;
        #my $up_plas=reverse($up_rmix);
        #$up_plas=~tr/ATCGatcg/TAGCtagc/;
        #my $down_plas=reverse($down_fmix);
        #$down_plas=~tr/ATCGatcg/TAGCtagc/;
        
        my ($up_mix_fir,$up_mix_sec,$down_mix_fir,$down_mix_sec);
        my ($left_r,$right_f);
        my ($fplasmid,$rplasmid);
        my $infoName="";
        my $geneName;
        my ($chr,$gStart,$gEnd,$v_up_a,$strand);
        while(<AN>){
                chomp;
                my $xtag=$_;
                if($xtag=~/^# BLASTN/){
                        $indicate++;
                        if($indicate == 1){
                                $fhit=0;
                                next;
                        }else{
                                while(<IN>){
                                        chomp;
                                        my $rtag=$_;
                                        if($rtag=~/^# BLASTN/){
                                                $indicateB++;
                                                if($indicateB==1){
                                                        $rhit=0;
                                                        next;
                                                }else{
                                                        $geneName=~s/\|/\t/g;
                                                        #$tm_p1=&tmValue($fplasmid,$pCon,$salt,$method);
                                                        #$tm_p2=&tmValue($rplasmid,$pCon,$salt,$method);
                                                        if($strand eq "+"){
                                                                #$v_up_a\t
                                                                print OUT "$geneName\t$strand\t$chr\t$gStart\t$gEnd\t$v_up_a\t$fplasmid\t$rplasmid\t$p1_tm\t$p2_tm\t$rhit\t$fhit\n";
                                                        }else{
                                                                print OUT "$geneName\t$strand\t$chr\t$gStart\t$gEnd\t$v_up_a\t$fplasmid\t$rplasmid\t$p1_tm\t$p2_tm\t$fhit\t$rhit\n";
                                                        }
                                                        last;
                                                }
                                        }else{
                                                if($indicateB == 1){
                                                        #if($rtag=~/^# Query:/){
                                                        #        $fgeneName=$rtag;
                                                        #        print "$fgeneName $rgeneName\n";
                                                        #}
                                                        next if($rtag=~/^#/);
                                                        my @rline=split /\s+/,$rtag;
                                                        if($strand eq "+"){
                                                                if($rline[3]>$upAlign){
                                                                        if($rline[2]>$identity){
                                                                                $rhit++;
                                                                        }
                                                                }
                                                        }else{
                                                                if($rline[3]>$downAlign){
                                                                        if($rline[2]>$identity){
                                                                                $rhit++;
                                                                        }
                                                                }
                                                        }         
                                                }
                                        }           
                                }

                                $indicate=1;
                                $indicateB=1;
                                $fhit=0;
                                $rhit=0;
                        }
                }else{
                        if($indicate == 1){
                                if($xtag=~/^# Query:/){
                                        $rgeneName=$xtag;
                                        $up_mix_fir=<BN>;
                                        $up_mix_sec=<BN>;
                                        $down_mix_fir=<DN>;
                                        $down_mix_sec=<DN>;
                                        chomp($up_mix_fir);
                                        chomp($up_mix_sec);
                                        chomp($down_mix_fir);
                                        chomp($down_mix_sec);
                                        if($up_mix_fir=~/>(.+)_lmix/){
                                                $geneName=$1;
                                        }
                                        while($geneName ne $infoName){
                                                if(scalar(@{$info})>0){
                                                        my $geneTag=shift(@{$info});
                                                        my @geneInfo=split /\,/,$geneTag;
                                                        ($infoName,$chr,$gStart,$gEnd,$v_up_a,$strand)=@geneInfo[0,1,6,7,8,13];
                                                }else{
                                                        last;
                                                }
                                        }
                        
                                        if($strand eq "+"){
                                                $left_r="$up_rmix" . "$up_mix_sec";
                                                $right_f="$down_fmix" . "$down_mix_sec";
                                                
                                                $fplasmid=reverse($left_r);
                                                $fplasmid=~tr/ATCGatcg/TAGCtagc/;
                                                $rplasmid=reverse($right_f);
                                                $rplasmid=~tr/ATCGatcg/TAGCtagc/;
                                                
                                        }else{
                                                $left_r="$up_rmix" . "$down_mix_sec";
                                                $right_f="$down_fmix" . "$up_mix_sec";
                                                
                                                $fplasmid=reverse($left_r);
                                                $fplasmid=~tr/ATCGatcg/TAGCtagc/;
                                                $rplasmid=reverse($right_f);
                                                $rplasmid=~tr/ATCGatcg/TAGCtagc/;    
                                        }
                                }else{
                                        next if($xtag=~/^#/);
                                        my @fline=split /\s+/,$xtag;
                                        if($strand eq "+"){
                                                if($fline[3]>$downAlign){
                                                        if($fline[2]>$identity){
                                                                $fhit++;
                                                        }
                                                }
                                        }else{
                                                if($fline[3]>$upAlign){
                                                        if($fline[2]>$identity){
                                                                $fhit++;
                                                        }
                                                }
                                        }
                                }
                                
                        }
                }
        }
        
        while(<IN>){
                my $rtag=$_;
                chomp($rtag);
                # $indicateB=1;
                #if($rtag=~/^# Query:/){
                #        $fgeneName=$rtag;
                #        print "$fgeneName $rgeneName\n";
                #}
                next if($rtag=~/^#/);
                my @rline=split /\s+/,$rtag;
                if($strand eq "+"){
                        if($rline[3]>$upAlign){
                                if($rline[2]>$identity){
                                        $rhit++;
                                }
                        }
                }else{
                        if($rline[3]>$downAlign){
                                if($rline[2]>$identity){
                                        $rhit++;
                                }
                        }
                }
        }
        $geneName=~s/\|/\t/g;
        #$tm_p1=&tmValue($fplasmid,$pCon,$salt,$method);
        #$tm_p2=&tmValue($rplasmid,$pCon,$salt,$method);
        if($strand eq "+"){
                #$v_up_a
                print OUT "$geneName\t$strand\t$chr\t$gStart\t$gEnd\t$v_up_a\t$fplasmid\t$rplasmid\t$p1_tm\t$p2_tm\t$rhit\t$fhit\n";
        }else{
                print OUT "$geneName\t$strand\t$chr\t$gStart\t$gEnd\t$v_up_a\t$fplasmid\t$rplasmid\t$p1_tm\t$p2_tm\t$fhit\t$rhit\n";
        }                               
        close IN;
        close AN;
        close BN;
        close DN;
        close OUT;
}

1;
















