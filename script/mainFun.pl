
#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
#use threads;
use lib $Bin;
use baseFun;
use auBase;
use knockout;
use ctag;
use ntag;

# rec and method were obsolete
my ($mode,$ctag,$ntag,$knockout,$force,$refine,$sscodon,$coordinate,$coordList,$geneName,$geneList,$all,$gff,
$genome,$mask,$marker,$rec,$rek,$outDir,$config,$nump,$bundle_size,$run_blast,$thread,$pCon,$salt,$dsalt,$dntp,$thermo,$method,$help);
$mode="short";
$mask=0;
$rec=0;
$rek=0;
$thread=1;
$run_blast=0;
$nump=5;
$bundle_size=1000;
$pCon=50;
$salt=50;
$dsalt=1.5;
$dntp=0.6;
$help=0;
$thermo=0;
$config='';
$force=0;
$method='nn';
GetOptions(        
		'ctag!'=>\$ctag,
		'ntag!'=>\$ntag,
		'knockout!'=>\$knockout,
        'force!'=>\$force,
        'refine!'=>\$refine,
        'sscodon!'=>\$sscodon,
		'mode:s'=>\$mode,
		'coordinate:s'=>\$coordinate,
        'coordList:s'=>\$coordList,
		'geneName:s'=>\$geneName,
        'geneList:s'=>\$geneList,
        'all!'=>\$all,
		'gff:s'=>\$gff,
		'genome:s'=>\$genome,
        'mask!'=>\$mask,
		'marker:s'=>\$marker,
        'rec!'=>\$rec,
        'rek!'=>\$rek,
		'outDir:s'=>\$outDir,
		'config:s'=>\$config,
        'nump:i'=>\$nump,
        'bundle:i'=>\$bundle_size,
        'blast!'=>\$run_blast,
        'thread:i'=>\$thread,
        'pcon:f'=>\$pCon,
        'salt:f'=>\$salt,
        'dsalt:f'=>\$dsalt,
        'dntp:f'=>\$dntp,
        'thermo!'=>\$thermo,
        'method:s'=>\$method,
		'help!'=>\$help
);

sub usage{
        print "Version 0.1.0\n";
        print "Usage: perl $0 [parameters]\n";
		print "--ctag                  design primers for adding tags at the C terminal.\n";
		print "--ntag                  design primers for adding tags at the N terminal.\n";
		print "--knockout              design primers for gene knockout.\n";
        print "--force                 force to output G3 and G4 that don't meet the criteria to design primers\n";
        print "--refine                refine primer design.\n";
        print "--sscodon               primers for gene knockout contain start codon and stop codon\n";
        print "--mode        [string]  method for primer design, by default: short (long | short).\n";
		print "--coordinate  [string]  coordinate of target gene, any location in the gene. format: chr:coordinate.\n";
		print "--coordList   [file]    list of coordinates of target genes. one line per gene.\n";
		print "--geneName    [string]  name of target gene, name can be ID, GeneName or DbxrefID.\n";
		print "--geneList    [file]    list of names of target genes. one line per gene.\n";
		print "--all                   design primers for all genes.\n";
        print "--gff         [file]    GFF file.\n";
		print "--genome      [file]    reference genome file.\n";
        print "--mask                  mask repetitive elements of genome.\n";
		print "--marker      [file]    selection marker or insertion sequence.\n";
        print "--rec                   redesign common sequence of primers, which is homologous to the insertion cassette (section of P1 and P2, V5, V6).\n";
        #print "--rek                   redesign checking primers, which is homologous to the insertion cassette.\n";
		print "--outDir      [dir]     output directory.\n";
		print "--config      [file]    config file (in most cases there's no need to modify the file).\n";
        print "--nump        [int]     number of primers output at most, by default 5.\n";
        print "--bundle      [int]     in case of designing primers for many genes bundle a certain number of genes together, by default 1000.\n";
        print "--blast                 run blast.\n";
        print "--thread      [int]     number of threads.\n";
        print "--pcon        [float]   primer concentration, by default 50 nM (unit nM).\n";
        print "--salt        [float]   concentration of monovalent cation, by default 50 mM (unit mM).\n";
        print "--dsalt       [float]   concentration of divalent cation, by default 1.5 mM (unit mM).\n";
        print "--dntp        [float]   concentration of dNTP, by default 0.6 mM (unit mM).\n";
        print "--thermo                thermodynamic models is used for oligo-oligo interactions and hairpins.\n";
        #print "--method      [string]  method to calculate Tm ['base' | 'basesalt' | 'longsalt' | 'nn'], by default 'nn'.\n";    
		  die "--help                  print help.\n";        
}

if($help || scalar(@ARGV)==1){
        &usage();
}
if($config eq ''){
    $config="$Bin/config.txt";
}

if($mode eq "short"){
        if($knockout){
                my $t_type="knockout";
                &knock_long_main($Bin,$config,$genome,$gff,$marker,$outDir,$geneName,$coordinate,$geneList,$coordList,$all,$mask,$nump,$bundle_size,$thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$sscodon,$force,$rec,$refine,$thermo);
        }elsif($ctag){
                my $t_type="ctag";
                &ctag_main($Bin,$config,$genome,$gff,$marker,$outDir,$geneName,$coordinate,$geneList,$coordList,$all,$mask,$nump,$bundle_size,$thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$force,$rec,$refine,$thermo);
        }elsif($ntag){
                my $t_type="ntag";
                &ntag_main($Bin,$config,$genome,$gff,$marker,$outDir,$geneName,$coordinate,$geneList,$coordList,$all,$mask,$nump,$bundle_size,$thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$force,$rec,$refine,$thermo);
        }else{
                die "Error: unknown command.\n";
        }
}else{
        if($knockout){
                my $t_type="knockout";
                &knock_short_main($Bin,$config,$genome,$gff,$marker,$outDir,$geneName,$coordinate,$geneList,$coordList,$all,$mask,$nump,$bundle_size,$thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$sscodon,$rec,$refine,$thermo);
        }elsif($ctag){
                my $t_type="ctag";
                &short_ctag_main($Bin,$config,$genome,$gff,$marker,$outDir,$geneName,$coordinate,$geneList,$coordList,$all,$mask,$nump,$bundle_size,$thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$rec,$refine,$thermo);
        }elsif($ntag){
                my $t_type="ntag";
                &short_ntag_main($Bin,$config,$genome,$gff,$marker,$outDir,$geneName,$coordinate,$geneList,$coordList,$all,$mask,$nump,$bundle_size,$thread,$run_blast,$t_type,$pCon,$salt,$dsalt,$dntp,$method,$rec,$refine,$thermo);
        }else{
                die "Error: unknown command.\n";
        }
}

#if($rek){
#        $mask=0;
#        &plasmidPrimer($Bin,$config,$marker,$outDir,$mask,$nump,$thread,$salt,$dsalt,$dntp,$pCon,$thermo,$refine);
#}



