
prefix=`pwd`
mkdir -p $prefix/third_party
mkdir -p $prefix/bin

cd $prefix/third_party
git clone https://github.com/primer3-org/primer3.git primer3
cd $prefix/third_party/primer3/src
make
make test

cd $prefix/third_party
wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar -zxf ncbi-blast-2.11.0+-x64-linux.tar.gz

#cd $prefix/third_party
wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/windowmasker/windowmasker
chmod +x windowmasker

ln -s $prefix/third_party/primer3/src/primer3_core $prefix/bin/primer3_core
ln -s $prefix/third_party/primer3/src/ntthal $prefix/bin/ntthal
ln -s $prefix/third_party/primer3/settings_files $prefix/settings_files

ln -s $prefix/third_party/ncbi-blast-2.11.1+/bin/makeblastdb $prefix/bin/makeblastdb
ln -s $prefix/third_party/ncbi-blast-2.11.1+/bin/blastn $prefix/bin/blastn
ln -s $prefix/third_party/windowmasker $prefix/bin/windowmasker
