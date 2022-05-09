
prefix=`pwd`
mkdir -p $prefix/third_party
mkdir -p $prefix/bin

cd $prefix/third_party
wget -c https://github.com/primer3-org/primer3/archive/refs/tags/v2.6.1.tar.gz
tar -zxf v2.6.1.tar.gz
cd $prefix/third_party/primer3-2.6.1/src
make
#make test

cd $prefix/third_party
wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar -zxf ncbi-blast-2.11.0+-x64-linux.tar.gz

#cd $prefix/third_party
#wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/windowmasker/windowmasker
#chmod +x windowmasker

cp $prefix/third_party/primer3-2.6.1/settings_files/primer3web_v0_4_0_default_settings.txt $prefix/third_party/primer3-2.6.1/settings_files/primer3web_v0_4_0_modify_settings.txt
sed -i -e 's/PRIMER_MAX_END_GC=5/PRIMER_MAX_END_GC=4/' \
-e 's/PRIMER_TM_FORMULA=0/PRIMER_TM_FORMULA=1/' \
-e 's/PRIMER_SALT_DIVALENT=0.0/PRIMER_SALT_DIVALENT=1.5/' \
-e 's/PRIMER_DNTP_CONC=0.0/PRIMER_DNTP_CONC=0.6/' \
-e 's/PRIMER_SALT_CORRECTIONS=0/PRIMER_SALT_CORRECTIONS=1/' \
-e 's/PRIMER_MAX_POLY_X=5/PRIMER_MAX_POLY_X=4/' \
-e 's/PRIMER_MIN_THREE_PRIME_DISTANCE=-1/PRIMER_MIN_THREE_PRIME_DISTANCE=1/' $prefix/third_party/primer3-2.6.1/settings_files/primer3web_v0_4_0_modify_settings.txt

cp $prefix/third_party/primer3-2.6.1/settings_files/primer3web_v4_0_0_default_settings.txt $prefix/third_party/primer3-2.6.1/settings_files/primer3web_v4_0_0_modify_settings.txt
sed -i -e 's/PRIMER_MIN_THREE_PRIME_DISTANCE=3/PRIMER_MIN_THREE_PRIME_DISTANCE=1/' \
-e 's/PRIMER_PAIR_MAX_DIFF_TM=5.0/#PRIMER_PAIR_MAX_DIFF_TM=5.0/' \
-e 's/PRIMER_MAX_END_GC=5/PRIMER_MAX_END_GC=4/' $prefix/third_party/primer3-2.6.1/settings_files/primer3web_v4_0_0_modify_settings.txt

ln -s $prefix/third_party/primer3-2.6.1/src/primer3_core $prefix/bin/primer3_core
ln -s $prefix/third_party/primer3-2.6.1/src/ntthal $prefix/bin/ntthal
ln -s $prefix/third_party/primer3-2.6.1/settings_files $prefix/settings_files

ln -s $prefix/third_party/ncbi-blast-2.11.0+/bin/makeblastdb $prefix/bin/makeblastdb
ln -s $prefix/third_party/ncbi-blast-2.11.0+/bin/blastn $prefix/bin/blastn
ln -s $prefix/third_party/ncbi-blast-2.11.0+/bin/windowmasker $prefix/bin/windowmasker

