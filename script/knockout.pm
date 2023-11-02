
#!/usr/bin/perl -w
use strict;
use threads;
#require "baseFun.pl";

# generate input file of primer3
sub oneSeq{
        my ($t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$t_ref,$t_dir,$t_info,$sscodon,$salt,$dsalt,$dntp,$pCon)=@_;
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
        my ($a_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$geneStart,$geneEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=@$t_info;
        my ($t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len);
        my ($v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len);
        if($strand eq "+"){
                $v_up_f_start=1;
                $v_up_f_len=$v_up_b-$v_up_a+1;
                $v_up_r_start=$geneStart-$v_up_a+1;
                $v_up_r_len=$geneEnd-$geneStart+1;

                #############
                if($sscodon){
                        $v_up_r_start+=3;
                        $v_up_r_len-=6;
                        #
                        $t_upR_start=$geneStart+4-$v_up_a-$t_PrimerMax;
                        $t_upR_len=$t_PrimerMax;
                        $t_downF_start=$geneEnd-1-$v_up_a;
                        $t_downF_len=$t_PrimerMax;
                }else{
                        $t_upR_start=$geneStart-$v_up_a-$t_PrimerMax+1;
                        $t_upR_len=$t_PrimerMax;
                        $t_downF_start=$geneEnd-$v_up_a+2;
                        $t_downF_len=$t_PrimerMax;
                }
                #############
                $v_down_f_start=$v_up_r_start;
                $v_down_f_len=$v_up_r_len;
                $v_down_r_start=$v_down_b-$v_up_a+1;
                $v_down_r_len=$v_down_a-$v_down_b+1;
                
                $t_upF_start=$t_up_a-$v_up_a+1;
                $t_upF_len=$t_up_b-$t_up_a+1;
                $t_downR_start=$t_down_b-$v_up_a+1;
                $t_downR_len=$t_down_a-$t_down_b+1;
        }else{
                $v_down_f_start=1;
                $v_down_f_len=$v_up_b-$v_up_a+1;
                $v_down_r_start=$geneStart-$v_up_a+1;
                $v_down_r_len=$geneEnd-$geneStart+1;
                #############
                if($sscodon){
                        $v_down_r_start+=3;
                        $v_down_r_len-=6;
                        #
                        $t_upF_start=$geneEnd-1-$v_up_a;
                        $t_upF_len=$t_PrimerMax;
                        $t_downR_start=$geneStart+4-$v_up_a-$t_PrimerMax;
                        $t_downR_len=$t_PrimerMax;
                }else{
                        $t_upF_start=$geneEnd+2-$v_up_a;
                        $t_upF_len=$t_PrimerMax;
                        $t_downR_start=$geneStart-$v_up_a-$t_PrimerMax+1;
                        $t_downR_len=$t_PrimerMax;
                }
                #############
                $v_up_f_start=$v_down_r_start;
                $v_up_f_len=$v_down_r_len;
                $v_up_r_start=$v_down_b-$v_up_a+1;
                $v_up_r_len=$v_down_a-$v_down_b+1;
                
                $t_upR_start=$t_down_b-$v_up_a+1;
                $t_upR_len=$t_down_a-$t_down_b+1;
                $t_downF_start=$t_up_a-$v_up_a+1;
                $t_downF_len=$t_up_b-$t_up_a+1;
        }
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
        (open OUD,">$t_dir/primer_para/$t_geneName.verify.up.para") || die "$!\n";
        $para=&vpPara($vseq,$a_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUD "$para";
        close OUD;
        
        (open OUE,">$t_dir/primer_para/$t_geneName.verify.down.para") || die "$!\n";
        $para=&vpPara($vseq,$a_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUE "$para";
        close OUE;
}

sub listSeq{
        my ($t_dir,$t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_bund,$t_pNum,$t_mask,$t_ref,$geneInfo,$sscodon,$salt,$dsalt,$dntp,$pCon)=@_;       
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
        my ($t_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$geneStart,$geneEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=split /\,/,$fir;
        my($t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len);
        my ($v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len);
        if($strand eq "+"){
                $v_up_f_start=1;
                $v_up_f_len=$v_up_b-$v_up_a+1;
                $v_up_r_start=$geneStart-$v_up_a+1;
                $v_up_r_len=$geneEnd-$geneStart+1;

                #############
                if($sscodon){
                        $v_up_r_start+=3;
                        $v_up_r_len-=6;
                        #
                        $t_upR_start=$geneStart+4-$v_up_a-$t_PrimerMax;
                        $t_upR_len=$t_PrimerMax;
                        $t_downF_start=$geneEnd-1-$v_up_a;
                        $t_downF_len=$t_PrimerMax;
                }else{
                        $t_upR_start=$geneStart-$v_up_a-$t_PrimerMax+1;
                        $t_upR_len=$t_PrimerMax;
                        $t_downF_start=$geneEnd-$v_up_a+2;
                        $t_downF_len=$t_PrimerMax;
                }
                #############
                $v_down_f_start=$v_up_r_start;
                $v_down_f_len=$v_up_r_len;
                $v_down_r_start=$v_down_b-$v_up_a+1;
                $v_down_r_len=$v_down_a-$v_down_b+1;
                
                $t_upF_start=$t_up_a-$v_up_a+1;
                $t_upF_len=$t_up_b-$t_up_a+1;
                $t_downR_start=$t_down_b-$v_up_a+1;
                $t_downR_len=$t_down_a-$t_down_b+1;
        }else{
                $v_down_f_start=1;
                $v_down_f_len=$v_up_b-$v_up_a+1;
                $v_down_r_start=$geneStart-$v_up_a+1;
                $v_down_r_len=$geneEnd-$geneStart+1;
                #############
                if($sscodon){
                        $v_down_r_start+=3;
                        $v_down_r_len-=6;
                        #
                        $t_upF_start=$geneEnd-1-$v_up_a;
                        $t_upF_len=$t_PrimerMax;
                        $t_downR_start=$geneStart+4-$v_up_a-$t_PrimerMax;
                        $t_downR_len=$t_PrimerMax;
                }else{
                        $t_upF_start=$geneEnd+2-$v_up_a;
                        $t_upF_len=$t_PrimerMax;
                        $t_downR_start=$geneStart-$v_up_a-$t_PrimerMax+1;
                        $t_downR_len=$t_PrimerMax;
                }
                #############
                $v_up_f_start=$v_down_r_start;
                $v_up_f_len=$v_down_r_len;
                $v_up_r_start=$v_down_b-$v_up_a+1;
                $v_up_r_len=$v_down_a-$v_down_b+1;
                
                $t_upR_start=$t_down_b-$v_up_a+1;
                $t_upR_len=$t_down_a-$t_down_b+1;
                $t_downF_start=$t_up_a-$v_up_a+1;
                $t_downF_len=$t_up_b-$t_up_a+1;
        }
        my $preChr=$chr;
        ##################################
        my $seq=&chrSeq($t_ref,$chr);
        my $vseq=substr($seq,$v_up_a-1,$v_len);
        $vseq=~s/[^ATCGNatcgn]/N/g;
        #############
        my $para;

        my $bundle=1;
        my $count=1;
        (open OUT,">$t_dir/primer_para/bundle.$bundle.target.up.para") || die "$!\n";
        (open OUG,">$t_dir/primer_para/bundle.$bundle.target.down.para") || die "$!\n";
        (open OUA,">$t_dir/primer_para/bundle.$bundle.verify.up.para") || die "$!\n";
        (open OUE,">$t_dir/primer_para/bundle.$bundle.verify.down.para") || die "$!\n";
        
        $para=&tUpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUT "$para";
        
        $para=&tDownPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUG "$para";     
        
        $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUA "$para";
        
        $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUE "$para";
        #################################
        my $nInfo=scalar(@$geneInfo);
        #foreach my $k(@$geneInfo){
        for(my $k=1;$k<$nInfo;$k++){
                ($t_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$geneStart,$geneEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=split /\,/,$geneInfo->[$k];
                if($strand eq "+"){
                        $v_up_f_start=1;
                        $v_up_f_len=$v_up_b-$v_up_a+1;
                        $v_up_r_start=$geneStart-$v_up_a+1;
                        $v_up_r_len=$geneEnd-$geneStart+1;

                        #############
                        if($sscodon){
                                $v_up_r_start+=3;
                                $v_up_r_len-=6;
                                #
                                $t_upR_start=$geneStart+4-$v_up_a-$t_PrimerMax;
                                $t_upR_len=$t_PrimerMax;
                                $t_downF_start=$geneEnd-1-$v_up_a;
                                $t_downF_len=$t_PrimerMax;
                        }else{
                                $t_upR_start=$geneStart-$v_up_a-$t_PrimerMax+1;
                                $t_upR_len=$t_PrimerMax;
                                $t_downF_start=$geneEnd-$v_up_a+2;
                                $t_downF_len=$t_PrimerMax;
                        }
                        #############
                        $v_down_f_start=$v_up_r_start;
                        $v_down_f_len=$v_up_r_len;
                        $v_down_r_start=$v_down_b-$v_up_a+1;
                        $v_down_r_len=$v_down_a-$v_down_b+1;
                        
                        $t_upF_start=$t_up_a-$v_up_a+1;
                        $t_upF_len=$t_up_b-$t_up_a+1;
                        $t_downR_start=$t_down_b-$v_up_a+1;
                        $t_downR_len=$t_down_a-$t_down_b+1;
                }else{
                        $v_down_f_start=1;
                        $v_down_f_len=$v_up_b-$v_up_a+1;
                        $v_down_r_start=$geneStart-$v_up_a+1;
                        $v_down_r_len=$geneEnd-$geneStart+1;
                        #############
                        if($sscodon){
                                $v_down_r_start+=3;
                                $v_down_r_len-=6;
                                #
                                $t_upF_start=$geneEnd-1-$v_up_a;
                                $t_upF_len=$t_PrimerMax;
                                $t_downR_start=$geneStart+4-$v_up_a-$t_PrimerMax;
                                $t_downR_len=$t_PrimerMax;
                        }else{
                                $t_upF_start=$geneEnd+2-$v_up_a;
                                $t_upF_len=$t_PrimerMax;
                                $t_downR_start=$geneStart-$v_up_a-$t_PrimerMax+1;
                                $t_downR_len=$t_PrimerMax;
                        }
                        #############
                        $v_up_f_start=$v_down_r_start;
                        $v_up_f_len=$v_down_r_len;
                        $v_up_r_start=$v_down_b-$v_up_a+1;
                        $v_up_r_len=$v_down_a-$v_down_b+1;
                        
                        $t_upR_start=$t_down_b-$v_up_a+1;
                        $t_upR_len=$t_down_a-$t_down_b+1;
                        $t_downF_start=$t_up_a-$v_up_a+1;
                        $t_downF_len=$t_up_b-$t_up_a+1;
                }
                #
                if($preChr ne $chr){
                        $seq=&chrSeq($t_ref,$chr);
                }
                #############
                $vseq=substr($seq,$v_up_a-1,$v_len);
                $vseq=~s/[^ATCGNatcgn]/N/g;
                #############
                $count++;
                if($count>$t_bund){
                        close OUT;
                        close OUG;
                        close OUA;
                        close OUE;
                        $bundle++;
                        (open OUT,">$t_dir/primer_para/bundle.$bundle.target.up.para") || die "$!\n";
                        (open OUG,">$t_dir/primer_para/bundle.$bundle.target.down.para") || die "$!\n";
                        (open OUA,">$t_dir/primer_para/bundle.$bundle.verify.up.para") || die "$!\n";
                        (open OUE,">$t_dir/primer_para/bundle.$bundle.verify.down.para") || die "$!\n";
                        $count=1;
                }
                #################################                
                $para=&tUpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUT "$para";
                
                $para=&tDownPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUG "$para";
                
                $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUA "$para";
                
                $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUE "$para";
                #
                $preChr=$chr;
        }
        close OUT;
        close OUG;
        close OUA;
        close OUE;
        #return($bundle);
}


sub multiRun1{
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
        my $v_up_inFile="$t_dir/primer_para/${s_geneName}.verify.up.para";
        my $v_down_inFile="$t_dir/primer_para/${s_geneName}.verify.down.para";
        
        my $t_up_outFile="$t_dir/primer_raw/${s_geneName}.target.up.out";
        my $t_down_outFile="$t_dir/primer_raw/${s_geneName}.target.down.out";
        my $v_up_outFile="$t_dir/primer_raw/${s_geneName}.verify.up.out";
        my $v_down_outFile="$t_dir/primer_raw/${s_geneName}.verify.down.out";
        
        my $t_up_errFile="$t_dir/primer_raw/${s_geneName}.target.up.err";
        my $t_down_errFile="$t_dir/primer_raw/${s_geneName}.target.down.err";
        my $v_up_errFile="$t_dir/primer_raw/${s_geneName}.verify.up.err";
        my $v_down_errFile="$t_dir/primer_raw/${s_geneName}.verify.down.err";       
        
        my @allTask=([$t_up_outFile,$t_up_errFile,$t_up_inFile],[$t_down_outFile,$t_down_errFile,$t_down_inFile],[$v_up_outFile,$v_up_errFile,$v_up_inFile],[$v_down_outFile,$v_down_errFile,$v_down_inFile]);
        my $task=0;
        if($t_thread<2){
                &runPrimer($s_prog,$s_setFile,$t_up_outFile,$t_up_errFile,$t_up_inFile);
                &runPrimer($s_prog,$s_setFile,$t_down_outFile,$t_down_errFile,$t_down_inFile);
                &runPrimer($s_prog,$s_setFile,$v_up_outFile,$v_up_errFile,$v_up_inFile);
                &runPrimer($s_prog,$s_setFile,$v_down_outFile,$v_down_errFile,$v_down_inFile);
        }else{
                #my $indicate=0;
                while(1){
                        last if($task==4);
                        while(scalar(threads->list(threads::all))<$t_thread){
                                threads->create(\&runPrimer,$s_prog,$s_setFile,$allTask[$task][0],$allTask[$task][1],$allTask[$task][2]);
                                $task++;
                                last if($task==4);
                        }
                        foreach my $k(threads->list(threads::all)){
                                if($k->is_joinable()){
                                        $k->join();
                                        last;
                                }
                        }
                }
                foreach my $k(threads->list(threads::all)){
                        $k->join();
                }
        }
}

sub multiRun2{
        my ($t_dir,$s_prog,$s_setFile,$t_thread,$n_bundle)=@_;
        if(! -d "$t_dir/primer_raw"){
                my $r=system("mkdir $t_dir/primer_raw");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_raw\n";
                }
        } 
        # bundle.1.target.para; bundle.1.verify.para        
        my $n=$n_bundle*4;
        if($t_thread<2){
                for(my $i=1;$i<=$n_bundle;$i++){
                        my $t_up_inFile="$t_dir/primer_para/bundle.$i.target.up.para";
                        my $t_down_inFile="$t_dir/primer_para/bundle.$i.target.down.para";
                        my $v_up_inFile="$t_dir/primer_para/bundle.$i.verify.up.para";
                        my $v_down_inFile="$t_dir/primer_para/bundle.$i.verify.down.para";
                        
                        my $t_up_outFile="$t_dir/primer_raw/bundle.$i.target.up.out";
                        my $t_down_outFile="$t_dir/primer_raw/bundle.$i.target.down.out";
                        my $v_up_outFile="$t_dir/primer_raw/bundle.$i.verify.up.out";
                        my $v_down_outFile="$t_dir/primer_raw/bundle.$i.verify.down.out";
                        
                        my $t_up_errFile="$t_dir/primer_raw/bundle.$i.target.up.err";
                        my $t_down_errFile="$t_dir/primer_raw/bundle.$i.target.down.err";
                        my $v_up_errFile="$t_dir/primer_raw/bundle.$i.verify.up.err";
                        my $v_down_errFile="$t_dir/primer_raw/bundle.$i.verify.down.err";
                        
                        &runPrimer($s_prog,$s_setFile,$t_up_outFile,$t_up_errFile,$t_up_inFile);
                        &runPrimer($s_prog,$s_setFile,$t_down_outFile,$t_down_errFile,$t_down_inFile);
                        &runPrimer($s_prog,$s_setFile,$v_up_outFile,$v_up_errFile,$v_up_inFile);
                        &runPrimer($s_prog,$s_setFile,$v_down_outFile,$v_down_errFile,$v_down_inFile);
                }
        }else{
                my $j=0;
                while(1){
                        last if($j>=$n);
                        while(scalar(threads->list(threads::all))<$t_thread){
                                $j++;
                                my $i=int(($j+3)/4);
                                if($j % 4 == 1){
                                        my $t_inFile="$t_dir/primer_para/bundle.$i.target.up.para";
                                        my $t_outFile="$t_dir/primer_raw/bundle.$i.target.up.out";
                                        my $t_errFile="$t_dir/primer_raw/bundle.$i.target.up.err";
                                        threads->create(\&runPrimer,$s_prog,$s_setFile,$t_outFile,$t_errFile,$t_inFile);
                                        print "Bundle $i target primers (upstream) running\n";
                                }elsif($j % 4 == 2){
                                        my $t_inFile="$t_dir/primer_para/bundle.$i.target.down.para";
                                        my $t_outFile="$t_dir/primer_raw/bundle.$i.target.down.out";
                                        my $t_errFile="$t_dir/primer_raw/bundle.$i.target.down.err";
                                        threads->create(\&runPrimer,$s_prog,$s_setFile,$t_outFile,$t_errFile,$t_inFile);
                                        print "Bundle $i target primers (downstream) running\n";
                                }elsif($j % 4 == 3){
                                        my $v_up_inFile="$t_dir/primer_para/bundle.$i.verify.up.para";
                                        my $v_up_outFile="$t_dir/primer_raw/bundle.$i.verify.up.out";
                                        my $v_up_errFile="$t_dir/primer_raw/bundle.$i.verify.up.err";
                                        threads->create(\&runPrimer,$s_prog,$s_setFile,$v_up_outFile,$v_up_errFile,$v_up_inFile);
                                        print "Bundle $i verification primers (upstream) running\n";
                                }else{
                                        my $v_down_inFile="$t_dir/primer_para/bundle.$i.verify.down.para";
                                        my $v_down_outFile="$t_dir/primer_raw/bundle.$i.verify.down.out";
                                        my $v_down_errFile="$t_dir/primer_raw/bundle.$i.verify.down.err";
                                        threads->create(\&runPrimer,$s_prog,$s_setFile,$v_down_outFile,$v_down_errFile,$v_down_inFile);
                                        print "Bundle $i verification primers (downstream) running\n";
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

sub primerLoc{
        my ($chr,$c_end,$gStart,$gEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand)=@_;
        my ($e_up_a,$e_up_b,$e_down_a,$e_down_b,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len);
        if($t_type eq "knockout"){
                if($strand eq "+"){
                        $e_up_a=$gStart-$t_upGeneFar;
                        $e_up_b=$gStart-$t_upGeneNear;
                        $v_up_a=$gStart-$t_upVfar;
                        $v_up_b=$gStart-$t_upVnear;
                        
                        if($e_up_b<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (upGeneNear,upGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($e_up_a<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (upGeneNear,upGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_up_b<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (upVnear,upVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_up_a<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (upVnear,upVfar).\n";
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
                        $e_down_a=$gEnd+$t_downGeneFar;
                        $e_down_b=$gEnd+$t_downGeneNear;
                        
                        $v_down_a=$gEnd+$t_downVfar;
                        $v_down_b=$gEnd+$t_downVnear;
                        $v_len=$t_upVfar+$t_downVfar+$gEnd-$gStart+1;
                        
                        if($e_down_b>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (downGeneNear,downGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($e_down_a>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (downGeneNear,downGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_down_b>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (downVnear,downVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_down_a>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (downVnear,downVfar).\n";
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
                        $e_up_a=$gStart-$t_downGeneFar;
                        $e_up_b=$gStart-$t_downGeneNear;
                        $v_up_a=$gStart-$t_downVfar;
                        $v_up_b=$gStart-$t_downVnear;
                        #
                        if($e_up_b<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (downGeneNear,downGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($e_up_a<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (downGeneNear,downGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_up_b<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (downVnear,downVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_up_a<=0){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene (downVnear,downVfar).\n";
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
                        $e_down_a=$gEnd+$t_upGeneFar;
                        $e_down_b=$gEnd+$t_upGeneNear;
                        
                        $v_down_a=$gEnd+$t_upVfar;
                        $v_down_b=$gEnd+$t_upVnear;
                        $v_len=$t_upVfar+$t_downVfar+$gEnd-$gStart+1;
                        
                        if($e_down_b>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene(upGeneNear,upGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($e_down_a>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene(upGeneNear,upGeneFar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_down_b>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene(upVnear,upVfar).\n";
                                return(0,0,0,0,0,0,0,0,0,0,0);
                        }elsif($v_down_a>$c_end){
                                print "Warning: The target gene may be close to the terminal of the chromosome (chromosome: $chr gene_start: $gStart gene_end: $gEnd strand: $strand).\n";
                                print "Try to modify the config file for this gene(upVnear,upVfar).\n";
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
                return($e_up_a,$e_up_b,$e_down_a,$e_down_b,$gStart,$gEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len);
        }
}

sub geneNameInterval{
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
        my $p_cds_p="";
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				#if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                if($line[2] =~ /.*gene$/){
                        if($indicate==1){
                                if(exists ${$refLen}{$chr}){
                                        $c_end=$refLen->{$chr};
                                        if(uc($t_geneName) eq uc($geneName) || $t_geneName eq $id || $t_geneName eq $dbxref){
                                                print "Gene information: $a_geneName\n";
                                                my $n=scalar(@cds);
                                                if($n<1){
                                                        print "Warning: no CDS.\n";
                                                }else{
                                                        # my @s_cds=sort {$a->[0] <=> $b->[0]} @cds;
                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }else{
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[$n-1][1];
                                                        }
                                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
                                                                if($tout[0] == 0){
                                                                        die "Gene locates in the terminal of chromosome.\n";
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
                        #
                        if($line[8] =~ /Parent=(.+?)(;|$)/){
                            my $cds_p=$1;
                            if($indicate){
                                if($cds_p ne $p_cds_p){
                                    next;
                                }
                            }else{
                                $p_cds_p=$cds_p;
                            }
                        }
                        #
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
                                        if($cds[0][0]>$cds[$n-1][0]){
                                                $cdsStart=$cds[$n-1][0];
                                                $cdsEnd=$cds[0][1];
                                        }else{
                                                $cdsStart=$cds[0][0];
                                                $cdsEnd=$cds[$n-1][1];
                                        }
                                                        
                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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

sub nameListInterval{
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
        my $p_cds_p="";
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				#if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                if($line[2] =~ /.*gene$/){
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
                                                                if($cds[0][0]>$cds[$n-1][0]){
                                                                        $cdsStart=$cds[$n-1][0];
                                                                        $cdsEnd=$cds[0][1];
                                                                }else{
                                                                        $cdsStart=$cds[0][0];
                                                                        $cdsEnd=$cds[$n-1][1];
                                                                }
                                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                        my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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
                        #
                        if($line[8] =~ /Parent=(.+?)(;|$)/){
                            my $cds_p=$1;
                            if($indicate){
                                if($cds_p ne $p_cds_p){
                                    next;
                                }
                            }else{
                                $p_cds_p=$cds_p;
                            }
                        }
                        #
                        push(@cds,[$line[3],$line[4]]);
                        $indicate=1;
                }	
		}
        #
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
                                                if($cds[0][0]>$cds[$n-1][0]){
                                                        $cdsStart=$cds[$n-1][0];
                                                        $cdsEnd=$cds[0][1];
                                                }else{
                                                        $cdsStart=$cds[0][0];
                                                        $cdsEnd=$cds[$n-1][1];
                                                }
                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                        my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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

sub coordInterval{
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
        my $p_cds_p="";
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				#if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                if($line[2] =~ /.*gene$/){
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
                                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                                $cdsStart=$cds[$n-1][0];
                                                                                $cdsEnd=$cds[0][1];
                                                                        }else{
                                                                                $cdsStart=$cds[0][0];
                                                                                $cdsEnd=$cds[$n-1][1];
                                                                        }
                                                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                                my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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
                        #
                        if($line[8] =~ /Parent=(.+?)(;|$)/){
                            my $cds_p=$1;
                            if($indicate){
                                if($cds_p ne $p_cds_p){
                                    next;
                                }
                            }else{
                                $p_cds_p=$cds_p;
                            }
                        }
                        #
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
                                                        if($cds[0][0]>$cds[$n-1][0]){
                                                                $cdsStart=$cds[$n-1][0];
                                                                $cdsEnd=$cds[0][1];
                                                        }else{
                                                                $cdsStart=$cds[0][0];
                                                                $cdsEnd=$cds[$n-1][1];
                                                        }
                                                        if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                                my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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

sub coordListInterval{ 
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
        my $p_cds_p="";
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				#if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                if($line[2] =~ /.*gene$/){
                        if($indicate==1){
                                if(exists ${$refLen}{$chr}){
                                        $c_end=$refLen->{$chr};
                                        my $n=scalar(@cds);
                                        #if($n<1){print "Warning: no CDS.\n";}
                                        if($n>=1){
                                                if($cds[0][0]>$cds[$n-1][0]){
                                                        $cdsStart=$cds[$n-1][0];
                                                        $cdsEnd=$cds[0][1];
                                                }else{
                                                        $cdsStart=$cds[0][0];
                                                        $cdsEnd=$cds[$n-1][1];
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
                                                                                my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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
                        #
                        if($line[8] =~ /Parent=(.+?)(;|$)/){
                            my $cds_p=$1;
                            if($indicate){
                                if($cds_p ne $p_cds_p){
                                    next;
                                }
                            }else{
                                $p_cds_p=$cds_p;
                            }
                        }
                        #
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
                                if($cds[0][0]>$cds[$n-1][0]){
                                        $cdsStart=$cds[$n-1][0];
                                        $cdsEnd=$cds[0][1];
                                }else{
                                        $cdsStart=$cds[0][0];
                                        $cdsEnd=$cds[$n-1][1];
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
                                                                my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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

sub allInterval{
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
        my $p_cds_p="";
		while(<IN>){
		        chomp;
				next if(/^#/);
				my @line=split;
				#if($line[2] eq "gene" || $line[2] eq "ncRNA_gene" || $line[2] eq "pseudogene"){
                #if($line[2] eq "gene" || $line[2] eq "pseudogene"){
                if($line[2] =~ /.*gene$/){
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
                                                if($cds[0][0]>$cds[$n-1][0]){
                                                        $cdsStart=$cds[$n-1][0];
                                                        $cdsEnd=$cds[0][1];
                                                }else{
                                                        $cdsStart=$cds[0][0];
                                                        $cdsEnd=$cds[$n-1][1];
                                                }
                                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                                        my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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
                        #
                        if($line[8] =~ /Parent=(.+?)(;|$)/){
                            my $cds_p=$1;
                            if($indicate){
                                if($cds_p ne $p_cds_p){
                                    next;
                                }
                            }else{
                                $p_cds_p=$cds_p;
                            }
                        }
                        #
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
                                if($cds[0][0]>$cds[$n-1][0]){
                                        $cdsStart=$cds[$n-1][0];
                                        $cdsEnd=$cds[0][1];
                                }else{
                                        $cdsStart=$cds[0][0];
                                        $cdsEnd=$cds[$n-1][1];
                                }
                                if($cdsEnd<=$gEnd && $cdsStart>=$gStart){
                                        my @tout=&primerLoc($chr,$c_end,$cdsStart,$cdsEnd,$t_upGeneFar,$t_upGeneNear,$t_upVfar,$t_upVnear,$t_downGeneFar,$t_downGeneNear,$t_downVfar,$t_downVnear,$t_type,$strand);
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

sub vpUpFormat{
        my ($t_dir,$bundle,$refine)=@_;
        my @out=map {"bundle.${_}.verify.up.out"} (1..$bundle);       
        my $up_fa="bundle.verify.up.forward.fa";
        my $down_fa="bundle.verify.up.reverse.fa";
        my $formatOut="bundle.verify.up.format.out";
        my $sta="bundle.verify.up.stat";
        my $vforce=0;
        &formatPrimerMulti($t_dir,\@out,$up_fa,$down_fa,$formatOut,$sta,$vforce,$refine);
}

sub vpDownFormat{
        my ($t_dir,$bundle,$refine)=@_;
        my @out=map {"bundle.${_}.verify.down.out"} (1..$bundle);       
        my $up_fa="bundle.verify.down.forward.fa";
        my $down_fa="bundle.verify.down.reverse.fa";
        my $formatOut="bundle.verify.down.format.out";
        my $sta="bundle.verify.down.stat";
        my $vforce=0;
        &formatPrimerMulti($t_dir,\@out,$up_fa,$down_fa,$formatOut,$sta,$vforce,$refine);
}

sub vpUpGeneFormat{
        my ($t_dir,$gName,$refine)=@_;
        my @out=("$gName.verify.up.out");
        my $up_fa="$gName.verify.up.forward.fa";
        my $down_fa="$gName.verify.up.reverse.fa";
        my $formatOut="$gName.verify.up.format.out";
        my $sta="$gName.verify.up.stat";
        my $vforce=0;
        &formatPrimerMulti($t_dir,\@out,$up_fa,$down_fa,$formatOut,$sta,$vforce,$refine);
}

sub vpDownGeneFormat{
        my ($t_dir,$gName,$refine)=@_;
        my @out=("$gName.verify.down.out");
        my $up_fa="$gName.verify.down.forward.fa";
        my $down_fa="$gName.verify.down.reverse.fa";
        my $formatOut="$gName.verify.down.format.out";
        my $sta="$gName.verify.down.stat";
        my $vforce=0;
        &formatPrimerMulti($t_dir,\@out,$up_fa,$down_fa,$formatOut,$sta,$vforce,$refine);
}


sub primer2blast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lpFa,$rpFa,$lvfFa,$lvrFa,$rvfFa,$rvrFa,$thread)=@_;
        &checkBlastdb($t_ref,$makeblastdb);
        my $r=0;
        $r=system("$blastn -task blastn-short -query $lpFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lp.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lpFa\n";
        }
        
        $r=system("$blastn -task blastn-short -query $lmixFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lmix.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lmixFa\n";
        }
        
        $r=system("$blastn -task blastn-short -query $rmixFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rmix.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rmixFa\n";
        }
        
        $r=system("$blastn -task blastn-short -query $rpFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rp.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rpFa\n";
        }
        #
        $r=system("$blastn -task blastn-short -query $lvfFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lvf.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lvfFa\n";
        }
        
        $r=system("$blastn -task blastn-short -query $lvrFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lvr.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lvrFa\n";
        }
        #
        $r=system("$blastn -task blastn-short -query $rvfFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rvf.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rvfFa\n";
        }
        
        $r=system("$blastn -task blastn-short -query $rvrFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rvr.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rvrFa\n";
        }
}

sub oneBlast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$thread,$geneName)=@_;
        my $dir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $x_ref="$dir/$refName";
        my $lpFa="$t_dir/primer_blast/$geneName.target.up.forward.fa";
        my $lmixFa="$t_dir/primer_blast/$geneName.target.up.reverse.fa";
        my $rmixFa="$t_dir/primer_blast/$geneName.target.down.forward.fa";
        my $rpFa="$t_dir/primer_blast/$geneName.target.down.reverse.fa";
        
        my $lvfFa="$t_dir/primer_blast/$geneName.verify.up.forward.fa";
        my $lvrFa="$t_dir/primer_blast/$geneName.verify.up.reverse.fa";
        my $rvfFa="$t_dir/primer_blast/$geneName.verify.down.forward.fa";
        my $rvrFa="$t_dir/primer_blast/$geneName.verify.down.reverse.fa";
        &primer2blast($t_dir,$x_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lpFa,$rpFa,$lvfFa,$lvrFa,$rvfFa,$rvrFa,$thread);
}

sub bundleBlast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$thread)=@_;
        my $dir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $x_ref="$dir/$refName";
        my $lpFa="$t_dir/primer_blast/bundle.target.up.forward.fa";
        my $lmixFa="$t_dir/primer_blast/bundle.target.up.reverse.fa";
        my $rmixFa="$t_dir/primer_blast/bundle.target.down.forward.fa";
        my $rpFa="$t_dir/primer_blast/bundle.target.down.reverse.fa";
        
        my $lvfFa="$t_dir/primer_blast/bundle.verify.up.forward.fa";
        my $lvrFa="$t_dir/primer_blast/bundle.verify.up.reverse.fa";
        my $rvfFa="$t_dir/primer_blast/bundle.verify.down.forward.fa";
        my $rvrFa="$t_dir/primer_blast/bundle.verify.down.reverse.fa";
        &primer2blast($t_dir,$x_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lpFa,$rpFa,$lvfFa,$lvrFa,$rvfFa,$rvrFa,$thread);
}


sub vpUpFilt{
        my ($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen)=@_;
        my $fblast="$t_dir/primer_blast/blastn.lvf.out";
        my $rblast="$t_dir/primer_blast/blastn.lvr.out";
        my $outFile="$t_dir/primer_blast/verify.up.length";
        &vPairFilt($t_dir,$fblast,$rblast,$outFile,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
}

sub vpDownFilt{
        my ($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen)=@_;
        my $fblast="$t_dir/primer_blast/blastn.rvf.out";
        my $rblast="$t_dir/primer_blast/blastn.rvr.out";
        my $outFile="$t_dir/primer_blast/verify.down.length";
        &vPairFilt($t_dir,$fblast,$rblast,$outFile,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
}


sub vpStrandTrans{
        my ($t_dir,$lvForm,$rvForm,$refine)=@_;
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
                my ($id,$geneName,$gid,$strand,$chr,$gStart,$gEnd,$v_up_a)=(split /\s+/,$_)[0,1,2,3,4,5,6,7];
                my $tag="${id}|${geneName}|$gid";
                if(! exists $tg_hash{$tag}){
                        $tg_hash{$tag}=[$chr,$gStart,$gEnd,$v_up_a];
                }
        }        
        close IN;
        ####################################
        my (@up_no,@down_no,@both_no);
        (open OUA,">$t_dir/primer_result/verify.up.unblast.txt") || die "$!\n";
        print OUA "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tPrimer_ID\tV1\tV3\tV1_Loc\tV3_Loc\tV1_Tm\tV3_Tm\tV1_V3_Expect_Len\tPenalty\n";
        (open OUB,">$t_dir/primer_result/verify.down.unblast.txt") || die "$!\n";
        print OUB "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tPrimer_ID\tV4\tV2\tV4_Loc\tV2_Loc\tV4_Tm\tV2_Tm\tV4_V2_Expect_Len\tPenalty\n";
        
        (open AN,"$lvForm") || die "$!\n";
        (open DN,"$rvForm") || die "$!\n";
        <AN>;<DN>;
        
        
        my %vpStat;
        my $num=0;
        my $pre_gene="";
        my @upout = ();
        my @downout = ();
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
                                #$v_up_a\t
                                #print OUA "$x_geneName\t+\t$chr\t$gStart\t$gEnd\t$num\t$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$p_len\t$penalty\n";
                                push(@upout,["$x_geneName\t+\t$chr\t$gStart\t$gEnd","$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$p_len",$penalty,$taBase]);
                        }else{
                                #print OUA "$x_geneName\t-\t$chr\t$gStart\t$gEnd\t$num\t$verify_r\t$verify_f\t$v2_loc\t$v1_loc\t$tm_r\t$tm_r\t$p_len\t$penalty\n";
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

        $pre_gene="";
        while(<DN>){
                chomp;
                my @line=split;
                my ($geneName,$strand,$penalty,$verify_f,$start_f,$len_f,$tm_f,$verify_r,$start_r,$len_r,$tm_r)=@line[0,1,3,4,5,6,7,13,14,15,16];
                if($geneName ne $pre_gene){
                        if(scalar(@downout)>0){
                                my @downSort;
                                if($refine){
                                       @downSort = sort {$a->[3] <=> $b->[3] or $a->[2] <=> $b->[2]} @downout;
                                }else{
                                       @downSort = sort {$a->[2] <=> $b->[2]} @downout;
                                }
                                
                                $num = 0;
                                foreach my $k(@downSort){
                                        $num++;
                                        print OUB "$k->[0]\t$num\t$k->[1]\t$k->[2]\n";
                                }
                        }
                        @downout = ();
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
                                #print OUB "$x_geneName\t+\t$chr\t$gStart\t$gEnd\t$num\t$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$p_len\t$penalty\n";
                                push(@downout,["$x_geneName\t+\t$chr\t$gStart\t$gEnd","$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$p_len",$penalty,$taBase]);
                        }else{
                                #print OUB "$x_geneName\t-\t$chr\t$gStart\t$gEnd\t$num\t$verify_r\t$verify_f\t$v2_loc\t$v1_loc\t$tm_r\t$tm_f\t$p_len\t$penalty\n";
                                push(@downout,["$x_geneName\t-\t$chr\t$gStart\t$gEnd","$verify_r\t$verify_f\t$v2_loc\t$v1_loc\t$tm_r\t$tm_f\t$p_len",$penalty,$taBase]);
                                
                        }
                        if(exists $vpStat{$geneName}){
                                if($vpStat{$geneName} == 1){
                                        $vpStat{$geneName} = 2;
                                }
                        }else{
                                $vpStat{$geneName} = 3;
                        }
                }
        }
        
        if(scalar(@downout)>0){
                my @downSort;
                if($refine){
                       @downSort = sort {$a->[3] <=> $b->[3] or $a->[2] <=> $b->[2]} @downout;
                }else{
                       @downSort = sort {$a->[2] <=> $b->[2]} @downout;
                }
                
                $num = 0;
                foreach my $k(@downSort){
                        $num++;
                        print OUB "$k->[0]\t$num\t$k->[1]\t$k->[2]\n";
                }
        }
        @downout = ();
                 
        foreach my $k(keys %tg_hash){
                if(! exists($vpStat{$k})){
                        push(@both_no,$k);
                }else{
                        if($vpStat{$k} == 1){
                                push(@down_no,$k);
                        }elsif($vpStat{$k} == 3){
                                push(@up_no,$k);
                        }
                }
        }
        print "#################################\n";
        print "Verification primers not exist (both upstream and downstream):\n";
        foreach my $m(sort(@both_no)){
                print "$m\n";
        }
        
        print "\n";
        print "Verification primers not exist (upstream):\n";
        foreach my $n(sort(@up_no)){
                print "$n\n";
        }
        
        print "\n";
        print "Verification primers not exist (downstream):\n";
        foreach my $n(sort(@down_no)){
                print "$n\n";
        }
        
        close AN;
        close DN;
        close OUA;
        close OUB;
}

sub vpGeneTrans{        
        my ($t_dir,$geneName,$refine)=@_;
        my $lvForm="$t_dir/primer_blast/$geneName.verify.up.format.out";
        my $rvForm="$t_dir/primer_blast/$geneName.verify.down.format.out";
        &vpStrandTrans($t_dir,$lvForm,$rvForm,$refine);
}

sub vpBundleTrans{
        my ($t_dir,$refine)=@_;
        my $lvForm="$t_dir/primer_blast/bundle.verify.up.format.out";
        my $rvForm="$t_dir/primer_blast/bundle.verify.down.format.out";
        &vpStrandTrans($t_dir,$lvForm,$rvForm,$refine);
}

sub vpMergeInfo{
        my ($t_dir,$lvForm,$rvForm,$refine)=@_;        
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
        #my @tgName;
        <IN>;
        while(<IN>){
                chomp;
                my ($id,$geneName,$gid,$chr,$gStart,$gEnd,$v_up_a)=(split /\s+/,$_)[0,1,2,4,5,6,7];
                my $tag="${id}|${geneName}|$gid";
                if(! exists $tg_hash{$tag}){
                        $tg_hash{$tag}=[$chr,$gStart,$gEnd,$v_up_a];
                        #push(@tgName,[$tag,$strand]);
                }
        }        
        close IN;
        #
        ####################################
        my $lvLen="$t_dir/primer_blast/verify.up.length";
        my $rvLen="$t_dir/primer_blast/verify.down.length";
        
        my (@up_no,@down_no,@both_no);
        (open OUA,">$t_dir/primer_result/verify.up.txt") || die "$!\n";
        print OUA "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tPrimer_ID\tV1\tV3\tV1_Loc\tV3_Loc\tV1_Tm\tV3_Tm\tV1_V3_Products_Counts\tV1_V3_Blast_Len\tV1_V3_Expect_Len\tV1_V3_Penalty\tTier\n";
        (open OUB,">$t_dir/primer_result/verify.down.txt") || die "$!\n";
        print OUB "#ID\tGene_Name\tDbxref_ID\tStrand\tChr\tCDS_Start\tCDS_End\tPrimer_ID\tV4\tV2\tV4_Loc\tV2_Loc\tV4_Tm\tV2_Tm\tV4_V2_Products_Counts\tV4_V2_Blast_Len\tV4_V2_Expect_Len\tV4_V2_Penalty\tTier\n";
        (open AN,"$lvForm") || die "$!\n";
        (open BN,"$t_dir/primer_blast/verify.up.length") || die "$!\n";
        (open DN,"$rvForm") || die "$!\n";
        (open EN,"$t_dir/primer_blast/verify.down.length") || die "$!\n";
        <AN>;<BN>;<DN>;<EN>;
        
        my %vpStat;
        my $tier = 0;
        my $wstat;
        my $num = 0;
        my $pre_gene = "";
        my @upout = ();
        my @downout = ();
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
        #
        $pre_gene = "";
        while(<DN>){
                chomp;
                my @line=split;
                my ($geneName,$strand,$penalty,$verify_f,$start_f,$len_f,$tm_f,$verify_r,$start_r,$len_r,$tm_r)=@line[0,1,3,4,5,6,7,13,14,15,16];
                if($geneName ne $pre_gene){
                        
                        if(scalar(@downout)>0){
                                my @downSort;
                                if($refine){
                                       @downSort = sort {$a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] or $a->[2] <=> $b->[2]} @downout;
                                }else{
                                       @downSort = sort {$a->[3] <=> $b->[3] or $a->[2] <=> $b->[2]} @downout;
                                }
                                
                                $num = 0;
                                foreach my $k(@downSort){
                                        $num++;
                                        print OUB "$k->[0]\t$num\t$k->[1]\t$k->[2]\t$k->[3]\n";
                                }
                        }
                        @downout = ();
                }
                $pre_gene = $geneName;
                
                my $m_line=<EN>;
                chomp($m_line);
                my @mLine=split /\s+/,$m_line;
                my $mistake_start=$mLine[4];
                my $pcr=$mLine[3];
                if($pcr eq "NA"){
                        next;
                }
                my ($geneName,$strand,$penalty,$verify_f,$start_f,$len_f,$tm_f,$verify_r,$start_r,$len_r,$tm_r)=@line[0,1,3,4,5,6,7,13,14,15,16];
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
                                        #print OUB "$x_geneName\t+\t$chr\t$gStart\t$gEnd\t$v_up_a\t$num\t$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$mistake_start\t$pcr\t$p_len\t$penalty\t$tier\n";
                                        push(@downout,["$x_geneName\t+\t$chr\t$gStart\t$gEnd","$verify_f\t$verify_r\t$v1_loc\t$v2_loc\t$tm_f\t$tm_r\t$mistake_start\t$pcr\t$p_len",$penalty,$tier,$taBase]);
                                }else{
                                        #print OUB "$x_geneName\t-\t$chr\t$gStart\t$gEnd\t$v_up_a\t$num\t$verify_r\t$verify_f\t$v2_loc\t$v1_loc\t$tm_r\t$tm_f\t$mistake_start\t$pcr\t$p_len\t$penalty\t$tier\n";
                                        push(@downout,["$x_geneName\t-\t$chr\t$gStart\t$gEnd","$verify_r\t$verify_f\t$v2_loc\t$v1_loc\t$tm_r\t$tm_f\t$mistake_start\t$pcr\t$p_len",$penalty,$tier,$taBase]);
                                }
                                if(exists $vpStat{$geneName}){
                                        if($vpStat{$geneName} == 1){
                                                $vpStat{$geneName} = 2;
                                        }
                                }else{
                                        $vpStat{$geneName} = 3;
                                }
                        }  
                }
        }
        
        if(scalar(@downout)>0){
                my @downSort;
                if($refine){
                       @downSort = sort {$a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] or $a->[2] <=> $b->[2]} @downout;
                }else{
                       @downSort = sort {$a->[3] <=> $b->[3] or $a->[2] <=> $b->[2]} @downout;
                }
                $num = 0;
                foreach my $k(@downSort){
                        $num++;
                        print OUB "$k->[0]\t$num\t$k->[1]\t$k->[2]\t$k->[3]\n";
                }
        }
        @downout = ();            
        #
        foreach my $k(keys %tg_hash){
                if(! exists($vpStat{$k})){
                        push(@both_no,$k);
                }else{
                        if($vpStat{$k} == 1){
                                push(@down_no,$k);
                        }elsif($vpStat{$k} == 3){
                                push(@up_no,$k);
                        }
                }
        }
        print "#################################\n";
        print "Verification primers not exist (both upstream and downstream):\n";
        foreach my $m(sort(@both_no)){
                print "$m\n";
        }
        
        print "\n";
        print "Verification primers not exist (upstream):\n";
        foreach my $n(sort(@up_no)){
                print "$n\n";
        }
        
        print "\n";
        print "Verification primers not exist (downstream):\n";
        foreach my $n(sort(@down_no)){
                print "$n\n";
        }
        
        close AN;
        close DN;
        close BN;
        close EN;
        close OUA;
        close OUB;
}


sub vpGeneMerge{        
        my ($t_dir,$geneName,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$refine)=@_;
        my $lvForm="$t_dir/primer_blast/$geneName.verify.up.format.out";
        my $rvForm="$t_dir/primer_blast/$geneName.verify.down.format.out";
        &vpUpFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        &vpDownFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        &vpMergeInfo($t_dir,$lvForm,$rvForm,$refine);
}

sub vpBundleMerge{
        my ($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$refine)=@_;
        my $lvForm="$t_dir/primer_blast/bundle.verify.up.format.out";
        my $rvForm="$t_dir/primer_blast/bundle.verify.down.format.out";
        &vpUpFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        &vpDownFilt($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen);
        &vpMergeInfo($t_dir,$lvForm,$rvForm,$refine);
}


sub knock_long_main{
        # common
        my ($binDir,$config,$t_ref,$t_gff,$plasmidFa,$t_dir,$t_geneName,$t_coord,$t_geneList,$t_coordList,$all,$t_mask,$t_pNum,$t_bundle_size,$t_thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$sscodon,$force,$rec,$refine,$thermo)=@_;
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
                        @geneInfo=&geneNameInterval($t_gff,$t_geneName,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coord){
                        my $loc;
                        ($loc=$t_coord)=~s/:/_/g;
                        $t_geneName="LOC_$loc";
                        @geneInfo=&coordInterval($t_gff,$t_coord,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }

                &oneSeq($t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$rename_ref,$t_dir,\@geneInfo,$sscodon,$salt,$dsalt,$dntp,$pCon);

                &multiRun1($t_dir,$t_geneName,$t_primer,$primerPara,$t_thread);
                
                &tpUpGeneFormat($t_dir,$t_geneName,$force,$refine);
                &tpDownGeneFormat($t_dir,$t_geneName,$force,$refine);

                &vpUpGeneFormat($t_dir,$t_geneName,$refine);
                &vpDownGeneFormat($t_dir,$t_geneName,$refine);
                
                my @info;
                push(@info,join ",",@geneInfo);
                if($run_blast){
                        &oneBlast($t_dir,$rename_ref,$t_blastn,$t_makeblastdb,$t_thread,$t_geneName);
                        #
                        #my ($t_dir,$geneName,$plasSeq,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$dMMost,$dERegion,$dEMatch,$info)
                        &tpGeneMerge($t_dir,$t_geneName,\@plasid_m_seq,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$d_endRegion,$d_endMatch,\@info,$t_pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                        #($t_dir,$geneName,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$dMMost,$dERegion,$dEMatch)
                        &vpGeneMerge($t_dir,$t_geneName,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$refine);
                }else{
                        
                        # my ($t_dir,$geneName,$plasSeq,$info,$dMMost,$dERegion,$dEMatch)
                        &tpGeneTrans($t_dir,$t_geneName,\@plasid_m_seq,\@info,$d_endRegion,$d_endMatch,$t_pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                        &vpGeneTrans($t_dir,$t_geneName,$refine);
                }
                return(0);
        }
        ####################
        # bundle
        if($t_geneList || $t_coordList || $all){
                my @geneInfo;
                my $bundle_counts;
                
                if($t_geneList){
                        @geneInfo=&nameListInterval($t_gff,$t_geneList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coordList){
                        @geneInfo=&coordListInterval($t_gff,$t_coordList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($all){
                        @geneInfo=&allInterval($t_gff,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }
                &listSeq($t_dir,$t_PrimerMin,$t_PrimerMax,$t_productMinLen,$t_productMaxLen,$t_bundle_size,$t_pNum,$t_mask,$rename_ref,\@geneInfo,$sscodon,$salt,$dsalt,$dntp,$pCon);
                $bundle_counts=&bundleCounts($t_dir);
                &multiRun2($t_dir,$t_primer,$primerPara,$t_thread,$bundle_counts);
                &tpUpFormat($t_dir,$bundle_counts,$force,$refine);
                &tpDownFormat($t_dir,$bundle_counts,$force,$refine);
                &vpUpFormat($t_dir,$bundle_counts,$refine);
                &vpDownFormat($t_dir,$bundle_counts,$refine);
                if($run_blast){
                        &bundleBlast($t_dir,$rename_ref,$t_blastn,$t_makeblastdb,$t_thread);
                        #
                        #my ($t_dir,$plasSeq,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$dMMost,$dERegion,$dEMatch,$info)
                        &tpBundleMerge($t_dir,\@plasid_m_seq,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$d_endRegion,$d_endMatch,\@geneInfo,$t_pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                        #($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$dMMost,$dERegion,$dEMatch)
                        &vpBundleMerge($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$refine);
                }else{
                        #($t_dir,$plasSeq,$info,$dMMost,$dERegion,$dEMatch)
                        &tpBundleTrans($t_dir,\@plasid_m_seq,\@geneInfo,$d_endRegion,$d_endMatch,$t_pNum,$p1_tm,$p2_tm,$refine,$ntthal,$salt,$dsalt,$dntp,$pCon,$d_max_tm,$thermo);
                        &vpBundleTrans($t_dir,$refine);
                        
                }
                return(0);
        }
        
        die "Error: lack of parameter (--geneName or --geneList or coordinate or coordList or --all)\n";
}

################################
#       short arm
################################

sub geneSeq{
        my ($seq,$v_up_a,$v_len,$geneStart,$geneEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand,$sscodon)=@_;
        my $vseq=substr($seq,$v_up_a-1,$v_len);
        $vseq=~s/[^ATCGNatcgn]/N/g;
        #my $t_gene_upseq=substr($$seq,$geneStart-$t_upMixGenomeSize-1,$t_upMixGenomeSize);
        my ($t_gene_upseq,$gene_mix_upseq,$gene_mix_downseq);
        if($strand eq "+"){
                if($sscodon){
                        $t_gene_upseq=substr($vseq,$geneStart+3-$v_up_a-$t_upMixGenomeSize,$t_upMixGenomeSize);
                        $gene_mix_downseq=substr($vseq,$geneEnd-2-$v_up_a,$t_downMixGenomeSize);
                }else{
                        $t_gene_upseq=substr($vseq,$geneStart-$v_up_a-$t_upMixGenomeSize,$t_upMixGenomeSize);
                        $gene_mix_downseq=substr($vseq,$geneEnd-$v_up_a+1,$t_downMixGenomeSize);
                }
        }else{
                if($sscodon){
                        $t_gene_upseq=substr($vseq,$geneStart+3-$v_up_a-$t_downMixGenomeSize,$t_downMixGenomeSize);
                        $gene_mix_downseq=substr($vseq,$geneEnd-2-$v_up_a,$t_upMixGenomeSize);
                }else{
                        $t_gene_upseq=substr($vseq,$geneStart-$v_up_a-$t_downMixGenomeSize,$t_downMixGenomeSize);
                        $gene_mix_downseq=substr($vseq,$geneEnd-$v_up_a+1,$t_upMixGenomeSize);
                }
        }
        $gene_mix_upseq=reverse($t_gene_upseq);
        $gene_mix_upseq=~tr/ATCGatcg/TAGCtagc/;
        #
        #$gene_mix_downseq=~s/[^ATCGNatcgn]/N/g;
        return($vseq,$gene_mix_upseq,$gene_mix_downseq);
}

sub geneSeq2{
        my ($t_ref,$chr,$v_up_a,$v_len,$geneStart,$geneEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand,$sscodon)=@_;
        my $seq=&chrSeq($t_ref,$chr);
        &geneSeq($seq,$v_up_a,$v_len,$geneStart,$geneEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand,$sscodon);
}

sub sk_oneSeq{
        my ($t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$t_ref,$t_dir,$t_info,$sscodon,$salt,$dsalt,$dntp,$pCon)=@_;
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
        my ($a_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$geneStart,$geneEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=@$t_info;
        my ($t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len);
        my ($v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len);
        if($strand eq "+"){
                $v_up_f_start=1;
                $v_up_f_len=$v_up_b-$v_up_a+1;
                $v_up_r_start=$geneStart-$v_up_a+1;
                $v_up_r_len=$geneEnd-$geneStart+1;

                #############
                if($sscodon){
                        $v_up_r_start+=3;
                        $v_up_r_len-=6;
                }
                #############
                $v_down_f_start=$v_up_r_start;
                $v_down_f_len=$v_up_r_len;
                $v_down_r_start=$v_down_b-$v_up_a+1;
                $v_down_r_len=$v_down_a-$v_down_b+1;
        }else{
                $v_down_f_start=1;
                $v_down_f_len=$v_up_b-$v_up_a+1;
                $v_down_r_start=$geneStart-$v_up_a+1;
                $v_down_r_len=$geneEnd-$geneStart+1;
                #############
                if($sscodon){
                        $v_down_r_start+=3;
                        $v_down_r_len-=6;
                }
                #############
                $v_up_f_start=$v_down_r_start;
                $v_up_f_len=$v_down_r_len;
                $v_up_r_start=$v_down_b-$v_up_a+1;
                $v_up_r_len=$v_down_a-$v_down_b+1;
        }
        #
        my ($vseq,$mix_upseq,$mix_downseq)=&geneSeq2($t_ref,$chr,$v_up_a,$v_len,$geneStart,$geneEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand,$sscodon);
        my $upN=($mix_upseq=~tr/Nn//);
        my $downN=($mix_downseq=~tr/Nn//);
        if($upN>5 || $downN >5){
                die "flanking sequence of gene has too many Ns\n";
        }
        my $para;                
        #
        (open OUD,">$t_dir/primer_para/$t_geneName.verify.up.para") || die "$!\n";
        $para=&vpPara($vseq,$a_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUD "$para";
        close OUD;
        
        (open OUE,">$t_dir/primer_para/$t_geneName.verify.down.para") || die "$!\n";
        $para=&vpPara($vseq,$a_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUE "$para";
        close OUE;
        
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

sub sk_listSeq{
        my ($t_dir,$t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_bund,$t_pNum,$t_mask,$t_ref,$geneInfo,$sscodon,$salt,$dsalt,$dntp,$pCon)=@_;       
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
        my ($t_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$geneStart,$geneEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=split /\,/,$fir;
        my($t_upF_start,$t_upF_len,$t_upR_start,$t_upR_len,$t_downF_start,$t_downF_len,$t_downR_start,$t_downR_len);
        my ($v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len);
        if($strand eq "+"){
                $v_up_f_start=1;
                $v_up_f_len=$v_up_b-$v_up_a+1;
                $v_up_r_start=$geneStart-$v_up_a+1;
                $v_up_r_len=$geneEnd-$geneStart+1;

                #############
                if($sscodon){
                        $v_up_r_start+=3;
                        $v_up_r_len-=6;
                }
                #############
                $v_down_f_start=$v_up_r_start;
                $v_down_f_len=$v_up_r_len;
                $v_down_r_start=$v_down_b-$v_up_a+1;
                $v_down_r_len=$v_down_a-$v_down_b+1;
        }else{
                $v_down_f_start=1;
                $v_down_f_len=$v_up_b-$v_up_a+1;
                $v_down_r_start=$geneStart-$v_up_a+1;
                $v_down_r_len=$geneEnd-$geneStart+1;
                #############
                if($sscodon){
                        $v_down_r_start+=3;
                        $v_down_r_len-=6;
                }
                #############
                $v_up_f_start=$v_down_r_start;
                $v_up_f_len=$v_down_r_len;
                $v_up_r_start=$v_down_b-$v_up_a+1;
                $v_up_r_len=$v_down_a-$v_down_b+1;
        }
        my $preChr=$chr;
        ##################################
        my $seq=&chrSeq($t_ref,$chr);
        my ($vseq,$mix_upseq,$mix_downseq)=&geneSeq($seq,$v_up_a,$v_len,$geneStart,$geneEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand,$sscodon);
        my $para;

        my $bundle=1;
        my $count=1;
        (open OUA,">$t_dir/primer_para/bundle.$bundle.verify.up.para") || die "$!\n";
        (open OUE,">$t_dir/primer_para/bundle.$bundle.verify.down.para") || die "$!\n";   
        
        $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUA "$para";
        
        $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
        print OUE "$para";
        #
        print OUB ">${t_geneName}_lmix\n";
        print OUB "$mix_upseq\n";
        print OUC ">${t_geneName}_rmix\n";
        print OUC "$mix_downseq\n";
        #################################
        my $nInfo=scalar(@$geneInfo);
        #foreach my $k(@$geneInfo){
        for(my $k=1;$k<$nInfo;$k++){
                ($t_geneName,$chr,$t_up_a,$t_up_b,$t_down_a,$t_down_b,$geneStart,$geneEnd,$v_up_a,$v_up_b,$v_down_a,$v_down_b,$v_len,$strand)=split /\,/,$geneInfo->[$k];
                if($strand eq "+"){
                        $v_up_f_start=1;
                        $v_up_f_len=$v_up_b-$v_up_a+1;
                        $v_up_r_start=$geneStart-$v_up_a+1;
                        $v_up_r_len=$geneEnd-$geneStart+1;

                        #############
                        if($sscodon){
                                $v_up_r_start+=3;
                                $v_up_r_len-=6;
                        }
                        #############
                        $v_down_f_start=$v_up_r_start;
                        $v_down_f_len=$v_up_r_len;
                        $v_down_r_start=$v_down_b-$v_up_a+1;
                        $v_down_r_len=$v_down_a-$v_down_b+1;
                }else{
                        $v_down_f_start=1;
                        $v_down_f_len=$v_up_b-$v_up_a+1;
                        $v_down_r_start=$geneStart-$v_up_a+1;
                        $v_down_r_len=$geneEnd-$geneStart+1;
                        #############
                        if($sscodon){
                                $v_down_r_start+=3;
                                $v_down_r_len-=6;
                        }
                        #############
                        $v_up_f_start=$v_down_r_start;
                        $v_up_f_len=$v_down_r_len;
                        $v_up_r_start=$v_down_b-$v_up_a+1;
                        $v_up_r_len=$v_down_a-$v_down_b+1;
                }
                #
                if($preChr ne $chr){
                        $seq=&chrSeq($t_ref,$chr);
                }
                ($vseq,$mix_upseq,$mix_downseq)=&geneSeq($seq,$v_up_a,$v_len,$geneStart,$geneEnd,$t_upMixGenomeSize,$t_downMixGenomeSize,$strand,$sscodon);
                my $upN=($mix_upseq=~tr/Nn//);
                my $downN=($mix_downseq=~tr/Nn//);
                if($upN>5 || $downN >5){
                        next;
                }
                $count++;
                if($count>$t_bund){
                        close OUA;
                        close OUE;
                        $bundle++;
                        (open OUA,">$t_dir/primer_para/bundle.$bundle.verify.up.para") || die "$!\n";
                        (open OUE,">$t_dir/primer_para/bundle.$bundle.verify.down.para") || die "$!\n";
                        $count=1;
                }
                #################################                
                
                $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_up_f_start,$v_up_f_len,$v_up_r_start,$v_up_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUA "$para";
                
                $para=&vpPara($vseq,$t_geneName,$t_mask,$t_PrimerMin,$t_PrimerMax,$v_down_f_start,$v_down_f_len,$v_down_r_start,$v_down_r_len,$strand,$t_pNum,$t_productMinLen,$t_productMaxLen,$salt,$dsalt,$dntp,$pCon);
                print OUE "$para";
                
                #
                print OUB ">${t_geneName}_lmix\n";
                print OUB "$mix_upseq\n";
                print OUC ">${t_geneName}_rmix\n";
                print OUC "$mix_downseq\n";
                #
                $preChr=$chr;
        }
        close OUA;
        close OUE;
        #return($bundle);
}

sub sk_multiRun1{
        my ($t_dir,$s_geneName,$s_prog,$s_setFile,$t_thread)=@_;
        if(! -d "$t_dir/primer_raw"){
                my $r=system("mkdir $t_dir/primer_raw");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_raw\n";
                }
        }        
        #
        my $v_up_inFile="$t_dir/primer_para/${s_geneName}.verify.up.para";
        my $v_down_inFile="$t_dir/primer_para/${s_geneName}.verify.down.para";
        
        my $v_up_outFile="$t_dir/primer_raw/${s_geneName}.verify.up.out";
        my $v_down_outFile="$t_dir/primer_raw/${s_geneName}.verify.down.out";
        
        my $v_up_errFile="$t_dir/primer_raw/${s_geneName}.verify.up.err";
        my $v_down_errFile="$t_dir/primer_raw/${s_geneName}.verify.down.err";       
        
        if($t_thread<2){
                &runPrimer($s_prog,$s_setFile,$v_up_outFile,$v_up_errFile,$v_up_inFile);
                &runPrimer($s_prog,$s_setFile,$v_down_outFile,$v_down_errFile,$v_down_inFile);
        }else{
                threads->create(\&runPrimer,$s_prog,$s_setFile,$v_up_outFile,$v_up_errFile,$v_up_inFile);
                threads->create(\&runPrimer,$s_prog,$s_setFile,$v_down_outFile,$v_down_errFile,$v_down_inFile);
                foreach my $k(threads->list(threads::all)){
                        $k->join();
                }
        }
}

sub sk_multiRun2{
        my ($t_dir,$s_prog,$s_setFile,$t_thread,$n_bundle)=@_;
        if(! -d "$t_dir/primer_raw"){
                my $r=system("mkdir $t_dir/primer_raw");
                if($r){
                        print "Error: mkdir failed.\n";
                        die "$t_dir/primer_raw\n";
                }
        } 
        # bundle.1.target.para; bundle.1.verify.para        
        my $n=$n_bundle*2;
        if($t_thread<2){
                for(my $i=1;$i<=$n_bundle;$i++){
                        #my $t_inFile="$t_dir/primer_para/bundle.$i.target.para";
                        my $v_up_inFile="$t_dir/primer_para/bundle.$i.verify.up.para";
                        my $v_down_inFile="$t_dir/primer_para/bundle.$i.verify.down.para";
                        
                        #my $t_outFile="$t_dir/primer_raw/bundle.$i.target.out";
                        my $v_up_outFile="$t_dir/primer_raw/bundle.$i.verify.up.out";
                        my $v_down_outFile="$t_dir/primer_raw/bundle.$i.verify.down.out";
                        
                        #my $t_errFile="$t_dir/primer_raw/bundle.$i.target.err";
                        my $v_up_errFile="$t_dir/primer_raw/bundle.$i.verify.up.err";
                        my $v_down_errFile="$t_dir/primer_raw/bundle.$i.verify.down.err";
                        
                        #&runPrimer($s_prog,$s_setFile,$t_outFile,$t_errFile,$t_inFile);
                        &runPrimer($s_prog,$s_setFile,$v_up_outFile,$v_up_errFile,$v_up_inFile);
                        &runPrimer($s_prog,$s_setFile,$v_down_outFile,$v_down_errFile,$v_down_inFile);
                }
        }else{
                my $j=0;
                while(1){
                        last if($j>=$n);
                        while(scalar(threads->list(threads::all))<$t_thread){
                                $j++;
                                my $i=int(($j+1)/2);
                                if($j % 2 == 1){
                                        my $v_up_inFile="$t_dir/primer_para/bundle.$i.verify.up.para";
                                        my $v_up_outFile="$t_dir/primer_raw/bundle.$i.verify.up.out";
                                        my $v_up_errFile="$t_dir/primer_raw/bundle.$i.verify.up.err";
                                        threads->create(\&runPrimer,$s_prog,$s_setFile,$v_up_outFile,$v_up_errFile,$v_up_inFile);
                                        print "Bundle $i verification primers (upstream)\n";
                                }else{
                                        my $v_down_inFile="$t_dir/primer_para/bundle.$i.verify.down.para";
                                        my $v_down_outFile="$t_dir/primer_raw/bundle.$i.verify.down.out";
                                        my $v_down_errFile="$t_dir/primer_raw/bundle.$i.verify.down.err";
                                        threads->create(\&runPrimer,$s_prog,$s_setFile,$v_down_outFile,$v_down_errFile,$v_down_inFile);
                                        print "Bundle $i verification primers (downstream)\n";
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

sub sk_primer2blast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lvfFa,$lvrFa,$rvfFa,$rvrFa,$thread)=@_;
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
        $r=system("$blastn  -task blastn-short -query $lvfFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lvf.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lvfFa\n";
        }
        
        $r=system("$blastn  -task blastn-short -query $lvrFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.lvr.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $lvrFa\n";
        }
        #
        $r=system("$blastn  -task blastn-short -query $rvfFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand plus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rvf.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rvfFa\n";
        }
        
        $r=system("$blastn  -task blastn-short -query $rvrFa -num_threads $thread -evalue 1 -db $t_ref.fa -strand minus -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs btop' >$t_dir/primer_blast/blastn.rvr.out");
        if($r){
                print "Error: blastn error\n";
                die "query -- $rvrFa\n";
        }
}

sub sk_oneBlast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$thread,$geneName)=@_;
        my $dir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $x_ref="$dir/$refName";
        
        my $lmixFa="$t_dir/primer_blast/$geneName.lmix.fa";
        my $rmixFa="$t_dir/primer_blast/$geneName.rmix.fa";
        my $lvfFa="$t_dir/primer_blast/$geneName.verify.up.forward.fa";
        my $lvrFa="$t_dir/primer_blast/$geneName.verify.up.reverse.fa";
        my $rvfFa="$t_dir/primer_blast/$geneName.verify.down.forward.fa";
        my $rvrFa="$t_dir/primer_blast/$geneName.verify.down.reverse.fa";
        &sk_primer2blast($t_dir,$x_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lvfFa,$lvrFa,$rvfFa,$rvrFa,$thread);
}

sub sk_bundleBlast{
        my ($t_dir,$t_ref,$blastn,$makeblastdb,$thread)=@_;
        my $dir=dirname($t_ref);
        my @suffix=(".fa",".fasta",".fna",".gz");
        my $refName=basename($t_ref,@suffix);
        my $x_ref="$dir/$refName";
        
        my $lmixFa="$t_dir/primer_blast/bundle.lmix.fa";
        my $rmixFa="$t_dir/primer_blast/bundle.rmix.fa";
        my $lvfFa="$t_dir/primer_blast/bundle.verify.up.forward.fa";
        my $lvrFa="$t_dir/primer_blast/bundle.verify.up.reverse.fa";
        my $rvfFa="$t_dir/primer_blast/bundle.verify.down.forward.fa";
        my $rvrFa="$t_dir/primer_blast/bundle.verify.down.reverse.fa";
        &sk_primer2blast($t_dir,$x_ref,$blastn,$makeblastdb,$lmixFa,$rmixFa,$lvfFa,$lvrFa,$rvfFa,$rvrFa,$thread);
}

sub knock_short_main{
        # common
        my ($binDir,$config,$t_ref,$t_gff,$plasmidFa,$t_dir,$t_geneName,$t_coord,$t_geneList,$t_coordList,$all,$t_mask,$t_pNum,$t_bundle_size,$t_thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$sscodon,$rec,$refine,$thermo)=@_;
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
                        @geneInfo=&geneNameInterval($t_gff,$t_geneName,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coord){
                        my $loc;
                        ($loc=$t_coord)=~s/:/_/g;
                        $t_geneName="LOC_$loc";
                        @geneInfo=&coordInterval($t_gff,$t_coord,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }

                &sk_oneSeq($t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_geneName,$t_pNum,$t_mask,$rename_ref,$t_dir,\@geneInfo,$sscodon,$salt,$dsalt,$dntp,$pCon);

                &sk_multiRun1($t_dir,$t_geneName,$t_primer,$primerPara,$t_thread);
                
                #&tpGeneFormat($t_dir,$t_geneName);

                &vpUpGeneFormat($t_dir,$t_geneName,$refine);
                &vpDownGeneFormat($t_dir,$t_geneName,$refine);
                
                my @info;
                push(@info,join ",",@geneInfo);
                if($run_blast){
                        &sk_oneBlast($t_dir,$rename_ref,$t_blastn,$t_makeblastdb,$t_thread,$t_geneName);
                        my $upAlign=sprintf("%d",$t_upMixGenomeSize*0.8);
                        my $downAlign=sprintf("%d",$t_downMixGenomeSize*0.8);
                        my $identity=80.0;
                        &s_tpGeneMerge($t_dir,$t_geneName,$upAlign,$downAlign,$identity,\@plasid_m_seq,\@info,$p1_tm,$p2_tm);
                        &vpGeneMerge($t_dir,$t_geneName,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$refine);
                }else{
                        &s_tpGeneTrans($t_dir,$t_geneName,\@plasid_m_seq,\@info,$p1_tm,$p2_tm);
                        &vpGeneTrans($t_dir,$t_geneName,$refine);
                }
                return(0);
        }
        ####################
        # bundle
        if($t_geneList || $t_coordList || $all){
                my @geneInfo;
                my $bundle_counts;
                
                if($t_geneList){
                        @geneInfo=&nameListInterval($t_gff,$t_geneList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($t_coordList){
                        @geneInfo=&coordListInterval($t_gff,$t_coordList,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }elsif($all){
                        @geneInfo=&allInterval($t_gff,$t_upVfar,$t_upVnear,$t_downVfar,$t_downVnear,$t_upGeneFar,$t_upGeneNear,$t_downGeneFar,$t_downGeneNear,$t_type,\%refHash);
                }
                &sk_listSeq($t_dir,$t_PrimerMin,$t_PrimerMax,$t_upMixGenomeSize,$t_downMixGenomeSize,$t_productMinLen,$t_productMaxLen,$t_bundle_size,$t_pNum,$t_mask,$rename_ref,\@geneInfo,$sscodon,$salt,$dsalt,$dntp,$pCon);
                $bundle_counts=&bundleCounts($t_dir);
                &sk_multiRun2($t_dir,$t_primer,$primerPara,$t_thread,$bundle_counts);
                #&tpFormat($t_dir,$bundle_counts);
                &vpUpFormat($t_dir,$bundle_counts,$refine);
                &vpDownFormat($t_dir,$bundle_counts,$refine);
                if($run_blast){
                        &sk_bundleBlast($t_dir,$rename_ref,$t_blastn,$t_makeblastdb,$t_thread);
                        my $upAlign=sprintf("%d",$t_upMixGenomeSize*0.8);
                        my $downAlign=sprintf("%d",$t_downMixGenomeSize*0.8);
                        my $identity=80.0;
                        &s_tpBundleMerge($t_dir,$upAlign,$downAlign,$identity,\@plasid_m_seq,\@geneInfo,$p1_tm,$p2_tm);
                        &vpBundleMerge($t_dir,$t_mismatchMore,$t_endRegion,$t_endMismatch,$t_mismatchAtLeast,$t_productMinLen,$t_productMaxLen,$refine);
                }else{
                        &s_tpBundleTrans($t_dir,\@plasid_m_seq,\@geneInfo,$p1_tm,$p2_tm);
                        &vpBundleTrans($t_dir,$refine);
                        
                }
                return(0);
        }
        
        die "Error: lack of parameter (--geneName or --geneList or coordinate or coordList or --all)\n";
}

1;





