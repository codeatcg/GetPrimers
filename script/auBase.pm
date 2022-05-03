

#!/usr/bin/perl -w
use strict;
use threads;
use Getopt::Long;

#require "baseFun.pl";

sub selectList{
        my ($allGene,$geneList,$outFile)=@_;
        (open AN,"$allGene") || die "primer file ($allGene) that contains all genes open failed. $!\n";
        (open IN,"$geneList") || die "gene list file ($geneList) open failed. $!\n";
        (open OUT,">$outFile") || die "output file open failed. $!\n";
        my @gene;
        my %hash;
        my %out;
        my $ncol;
        my $pidcol=0;
        while(<IN>){
                chomp;
                if(/^\s*$/){
                        my @line=split;
                        next;
                }
                push(@gene,$_);
                $hash{$_}=1;
        }
        
        while(<AN>){
                chomp;
                if(/^#/){
                        print OUT "$_\n";
                        my @line=split;
                        $ncol=scalar(@line);
                        for(my $k=0;$k<=$#line;$k++){
                                if($line[$k] eq "Primer_ID" || $line[$k] eq "PrimerID"){
                                        $pidcol=$k;
                                        last;
                                }
                        }
                }else{
                        my @line=split;
                        my @tag=@line[0,1,2];
                        my $tgene="";
                        foreach my $k(@tag){
                                if(exists $hash{$k}){                               
                                        $tgene=$k;
                                        last;
                                }
                        }
                        if($tgene ne ""){
                                $out{$tgene}{$line[$pidcol]}=$_;
                        }
                        
                }
                
        }        
        close IN;
        close AN;
        foreach my $k(@gene){
                if(exists $out{$k}){
                        foreach my $m(sort {$a<=>$b} keys %{$out{$k}}){
                                print OUT "$out{$k}{$m}\n";
                        }
                }else{
                        my $unLine="${k}\tunknown\tunknown" . "\tNA" x ($ncol-3);
                        print OUT "$unLine\n";
                }
        }
        close OUT;
}

sub geneLoc{
        my ($t_gff,$outFile)=@_;
		if($t_gff=~/\.gz$/){
		        (open IN,"gzip -dc $t_gff|") || die "GFF open failed. $!\n";
		}else{
		        (open IN,"$t_gff") || die "GFF open failed. $!\n";
		}
        
        (open OUT,">$outFile") || die "$!\n";
        print OUT "#Chr\tGene_Start\tGene_End\tID\tGene_Name\tDbxref_ID\tBiotype\tStrand\n";
        
		my ($chr,$geneName,$a_geneName,$type);      
        my %r_hash;
        my $id;
        my $dbxref;
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				if($line[2] eq "gene" || $line[2] eq "ncRNA_gene" || $line[2] eq "pseudogene"){
                        $chr=$line[0];
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
                        
                        if($line[8]=~/biotype=(.+?);/){
						        $type=$1;
						}elsif($line[8]=~/biotype=(.+)/){
						        $type=$1;
						}else{
                                $type="unknown";
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
                        if($a_geneName eq "unknown|unknown|unknown"){
                                next;
                        }
						# a b ===== b a
						# =====a b=== b a
						my ($gStart,$gEnd);
						if($line[3]>$line[4]){
						        $gStart=$line[4];
								$gEnd=$line[3];
						}else{
						        $gStart=$line[3];
								$gEnd=$line[4];
						}
				
                        if(exists $r_hash{$a_geneName}){
                                print "Warning: gene name duplicate in the GFF file.\n";
                                print "The gene on chromosome $chr, location $gStart -- $gEnd will be skipped.\n"; 
                        }else{
                                my $strand=$line[6];
                                print OUT "$chr\t$gStart\t$gEnd\t$id\t$geneName\t$dbxref\t$type\t$strand\n";
                        }	
                }				
		}		
		close IN;
        close OUT;
}

############
sub getPseq{
        my $t_pfa=$_[0];
        my $seq='';
        if($t_pfa =~ /\.gz$/){
                (open IN,"gzip -dc $t_pfa|") || die "plasmid sequence open failed. $!\n";
        }else{
                (open IN,"$t_pfa") || die "plasmid sequence open failed. $!\n";
        }
        while(<IN>){
                chomp;
                next if(/^\s*$/ || />/);
                s/[^ATCGNatcgn]/N/g;
                $seq .=$_;
        }
        close IN;
        
        my $seqLen = length($seq);
        if($seqLen < 70){
                die "Error: length of insertion sequence (plasmid) is too short!\n";
        }
        return($seq);
}

sub vplasmidPara{
        my ($pseq,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_pNum,$n_left,$n_right,$seq_start,$seq_len,$salt,$dsalt,$dntp,$pCon)=@_;
        my $para="SEQUENCE_ID=plasmid_verify\n" .
                "SEQUENCE_TEMPLATE=$pseq\n" .
                "PRIMER_LOWERCASE_MASKING=$t_mask\n" .
                "PRIMER_TASK=generic\n" .
                "PRIMER_PICK_LEFT_PRIMER=$n_left\n" .
                "PRIMER_PICK_RIGHT_PRIMER=$n_right\n" .
                "PRIMER_PICK_INTERNAL_OLIGO=0\n" .
                "PRIMER_OPT_SIZE=20\n" .
                "PRIMER_MIN_SIZE=$t_PrimerMin\n" .
                "PRIMER_MAX_SIZE=$t_PrimerMax\n" .
                "SEQUENCE_INCLUDED_REGION=$seq_start,$seq_len\n" .
                
                "PRIMER_DNA_CONC=$pCon\n" . 
                "PRIMER_SALT_MONOVALENT=$salt\n" .
                "PRIMER_SALT_DIVALENT=$dsalt\n" .
                "PRIMER_DNTP_CONC=$dntp\n" . 
                
                "PRIMER_NUM_RETURN=$t_pNum\n" .
                "=\n";
        return($para);
}


sub vplasmidSeq{
        my ($t_dir,$pseqFile,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_pNum,$t_productMinLen,$t_productMaxLen,$t_upVnear,$t_upVfar,$t_downVnear,$t_downVfar,$salt,$dsalt,$dntp,$pCon)=@_;
        my $para;
        if(! -d "$t_dir/primer_para"){
                my $r=system("mkdir $t_dir/primer_para");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_para\n";
                }
        }
        my $tpseq=&getPseq($pseqFile);
        my $p_len=length($tpseq);
        
        my $r_start=$t_productMinLen-$t_upVfar;
        if($r_start<1){
                $r_start=1;
        }
        
        my $r_end=$t_productMaxLen-$t_upVnear;
        if($r_end>$p_len){
                $r_end=$p_len;
        }
        my $r_len=$r_end-$r_start+1;
        if($r_len<20){
                die "Error: un-proper PCR product size.\n";
        }
        #
        my $f_start=$p_len-$t_productMaxLen+$t_downVnear;
        if($f_start<1){
                $f_start=1;
        }
        my $f_end=$p_len-$t_productMinLen+$t_downVfar;
        if($f_end>$p_len){
                $f_end=$p_len;
        }
        my $f_len=$f_end-$f_start;
        if($f_len<20){
                die "Error: un-proper PCR product size.\n";
        }
        
        (open OUA,">$t_dir/primer_para/plasmid.verify.reverse.para") || die "$!\n";
        $para=&vplasmidPara($tpseq,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_pNum,0,1,$r_start,$r_len,$salt,$dsalt,$dntp,$pCon);
        print OUA "$para";
        close OUA;
        
        (open OUB,">$t_dir/primer_para/plasmid.verify.forward.para") || die "$!\n";
        $para=&vplasmidPara($tpseq,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_pNum,1,0,$f_start,$f_len,$salt,$dsalt,$dntp,$pCon);
        print OUB "$para";
        close OUB;
        return($p_len);
}

sub multiRun3{
        my ($t_dir,$s_prog,$s_setFile,$t_thread)=@_;
        if(! -d "$t_dir/primer_raw"){
                my $r=system("mkdir $t_dir/primer_raw");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_raw\n";
                }
        }        
        #
        my $r_inFile="$t_dir/primer_para/plasmid.verify.reverse.para";
        my $f_inFile="$t_dir/primer_para/plasmid.verify.forward.para";
        
        my $r_outFile="$t_dir/primer_raw/plasmid.verify.reverse.out";
        my $f_outFile="$t_dir/primer_raw/plasmid.verify.forward.out";
        
        my $r_errFile="$t_dir/primer_raw/plasmid.verify.reverse.err";
        my $f_errFile="$t_dir/primer_raw/plasmid.verify.forward.err";
        
        if($t_thread<2){
                &runPrimer($s_prog,$s_setFile,$r_outFile,$r_errFile,$r_inFile);
                &runPrimer($s_prog,$s_setFile,$f_outFile,$f_errFile,$f_inFile);
        }else{
                threads->create(\&runPrimer,$s_prog,$s_setFile,$r_outFile,$r_errFile,$r_inFile);
                threads->create(\&runPrimer,$s_prog,$s_setFile,$f_outFile,$f_errFile,$f_inFile);
                foreach my $k(threads->list(threads::all)){
                        $k->join();
                }
        }
}

sub formatPrimerOne{
        my ($t_dir,$rawPrimer,$outFa,$formatOut,$sta,$result,$header,$strand,$p_len,$t_upVnear,$t_upVfar,$t_downVnear,$t_downVfar,$refine)=@_;
        if(! -d "$t_dir/primer_blast"){
                my $r=system("mkdir $t_dir/primer_blast");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_blast\n";
                }
        }
        
        if(! -d "$t_dir/primer_result"){
                my $r=system("mkdir $t_dir/primer_result");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_result\n";
                }
        }
        
        (open OUA,">$t_dir/primer_blast/$outFa") || die "$!\n";
        (open OUT,">$t_dir/primer_blast/$formatOut") || die "$!\n";
        (open OUN,">$t_dir/primer_blast/$sta") || die "$!\n";
        (open OUQ,">$t_dir/primer_result/$result") || die "$!\n";
        
        print OUT "#GeneName\tPrimerID\tPenalty\tSeq\tStart\tLen\tTm\tGC%\tSelf_any\tSelf_end\tEnd_stab\n";
        print OUN "#GeneName\tStat\n";
        print OUQ "$header\n";
        #        
        my $num=0;
        my ($penalty1,$seq1,$start1,$len1,$tm1,$gc1,$selAny1,$selEnd1,$endSta1)=(0,0,0,0,0,0,0,0,0);
        my $geneName;
        my $tSta=0;
        (open AN,"$t_dir/primer_raw/$rawPrimer") || die "$!\n";
        my @out;
        while(<AN>){
                chomp;
                my @line=split /=/,$_;
                
                if($line[0] eq "SEQUENCE_ID"){
                        $geneName=$line[1];
                }elsif($line[0]=~/_EXPLAIN$/){
                        if($line[1]=~/ok\s+(\d+)/){
                                if($refine){
                                        $tSta=0;
                                        next;
                                }
                                        
                                if($1>0){
                                        $tSta=$1;
                                }else{
                                        $tSta=0;
                                }
                                print OUN "$geneName\t$tSta\n";
                        }
                }elsif($line[0]=~/_(\d+)_PENALTY$/){
                        $penalty1=sprintf("%.3f",$line[1]);
                        $num=$1;
                }elsif($line[0]=~/_SEQUENCE$/){
                        $seq1=$line[1];
                        print OUA ">${geneName}_${num}\n";
                        print OUA "$seq1\n";
                }elsif($line[0]=~/_\d+$/){
                        ($start1,$len1)=split /\,/,$line[1];
                }elsif($line[0]=~/_\d+_TM$/){
                        $tm1=$line[1];
                }elsif($line[0]=~/_GC_PERCENT$/){
                        $gc1=$line[1];
                }elsif($line[0]=~/_SELF_ANY(_TH)?$/){
                        $selAny1=$line[1];
                }elsif($line[0]=~/_\d+_SELF_END(_TH)?$/){
                        $selEnd1=$line[1];
                }elsif($line[0]=~/_\d+_END_STABILITY$/){
                        $endSta1=$line[1];
                        print OUT "$geneName\t$num\t$penalty1\t$seq1\t$start1\t$len1\t$tm1\t$gc1\t$selEnd1\t$endSta1\n";
                        #
                        if($refine){
                                my $four_f=substr($seq1,-4);
                                $four_f=uc($four_f);
                                if($four_f eq "AAAA" || $four_f eq "TTTT" || $four_f eq "GGGG" || $four_f eq "CCCC"){
                                        next;
                                }
                                $tSta++;
                        }
                        #
                        my $base1 = substr($seq1,-1);
                        my $taBase = 0;
                        if($base1 eq 'A' || $base1 eq 'T' || $base1 eq 'a' || $base1 eq 't'){
                                $taBase++;
                        }
                        if($strand eq '+'){
                                my $p_start=$start1;
                                my $p_end=$p_start+$len1-1;
                                my $r_len=$p_len-$p_start+1;
                                my $r_near=$r_len+$t_downVnear;
                                my $r_far=$r_len+$t_downVfar;
                                #print OUQ "$geneName\t$num\t$seq1\t${p_start}-${p_end}:+\t$tm1\t${r_near}-${r_far}\t$penalty1\n";
                                push(@out,[$geneName,"$seq1\t${p_start}-${p_end}:+\t$tm1\t${r_near}-${r_far}",$penalty1,$taBase]);
                        }else{
                                my $p_end=$start1;
                                my $p_start=$p_end-$len1+1;
                                my $l_near=$p_end+$t_upVnear;
                                my $l_far=$p_end+$t_upVfar;
                                #print OUQ "$geneName\t$num\t$seq1\t${p_start}-${p_end}:-\t$tm1\t${l_near}-${l_far}\t$penalty1\n";
                                push(@out,[$geneName,"$seq1\t${p_start}-${p_end}:-\t$tm1\t${l_near}-${l_far}",$penalty1,$taBase]);
                        }
                }

        }
        
        if($refine){
                my @sout=sort {$a->[3] <=> $b->[3] or $a->[2] <=> $b->[2]} @out;
                $num=0;
                foreach my $k(@sout){
                        $num++;
                        print OUQ "$k->[0]\t$num\t$k->[1]\t$k->[2]\n";
                }
                print OUN "$geneName\t$tSta\n";
        }else{
                $num=0;
                foreach my $k(@out){
                        $num++;
                        print OUQ "$k->[0]\t$num\t$k->[1]\t$k->[2]\n";
                }
        }
        
        close AN;
        close OUA;
        close OUT;
        close OUN;
        close OUQ;
}

sub plasmidFormat{
        my ($t_dir,$p_len,$t_upVnear,$t_upVfar,$t_downVnear,$t_downVfar,$refine)=@_;
        my ($rawPrimer,$outFa,$formatOut,$sta,$result,$header);
        $rawPrimer="plasmid.verify.reverse.out";
        $outFa="plasmid.verify.reverse.fa";
        $formatOut="plasmid.verify.reverse.format.out";
        $sta="plasmid.verify.reverse.stat";
        $result="plasmid.verify.reverse.txt";
        $header="#Gene_Name\tPrimer_ID\tV5\tV5_Loc\tV5_Tm\tV1_V5_Len_Range\tV5_Penalty";
        &formatPrimerOne($t_dir,$rawPrimer,$outFa,$formatOut,$sta,$result,$header,'-',$p_len,$t_upVnear,$t_upVfar,$t_downVnear,$t_downVfar,$refine);
        #
        $rawPrimer="plasmid.verify.forward.out";
        $outFa="plasmid.verify.forward.fa";
        $formatOut="plasmid.verify.forward.format.out";
        $sta="plasmid.verify.forward.stat";
        $result="plasmid.verify.forward.txt";
        $header="#Gene_Name\tPrimer_ID\tV6\tV6_Loc\tV6_Tm\tV6_V2_Len_Range\tV6_Penalty";
        &formatPrimerOne($t_dir,$rawPrimer,$outFa,$formatOut,$sta,$result,$header,'+',$p_len,$t_upVnear,$t_upVfar,$t_downVnear,$t_downVfar,$refine);
}

sub plasmidPrimer{
        my ($binDir,$config,$plasmidFa,$t_dir,$t_mask,$t_pNum,$t_thread,$salt,$dsalt,$dntp,$pCon,$thermo,$refine)=@_;
        my ($t_primer,$t_blastn,$t_makeblastdb,$t_windowmasker,$primerPara,$ntthal)=&getProg($binDir,$thermo);
        my ($r_primer,$r_blastn,$r_makeblastdb,$r_windowmasker,$t_upVfar,$t_upVnear,
        $t_upGeneFar,$t_upGeneNear,$t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_upMixPlasmidSize,
        $t_downVfar,$t_downVnear,$t_downGeneFar,$t_downGeneNear,$t_downMixGenomeSize,$t_downMixPlasmidSize,
        $t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,
        $d_GminMinus,$d_endRegion,$d_endMatch,$settings)=&read_config($config);
        
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
        my $p_len=&vplasmidSeq($t_dir,$plasmidFa,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_pNum,$t_productMinLen,$t_productMaxLen,$t_upVnear,$t_upVfar,$t_downVnear,$t_downVfar,$salt,$dsalt,$dntp,$pCon);
        &multiRun3($t_dir,$t_primer,$primerPara,$t_thread);
        &plasmidFormat($t_dir,$p_len,$t_upVnear,$t_upVfar,$t_downVnear,$t_downVfar,$refine);
}

1;














