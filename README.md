# GetPrimers

<p align="center">
  <img src="https://github.com/codeatcg/GetPrimers/blob/main/doc/figure/getprimer_logo.png" width="50%" height="50%"/>
</p>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Description
GetPrimers was originally developed to design PCR based gene targeting primers for yeasts with no limitations of strains. It provided two strategies, short primer strategy and long primer strategy. By short primer strategy customers can construct 200-600 bp or even larger homologous arms using a set of 18-30 bp primers. This may be useful for genes that are difficult to be targeted with traditional ~100 bp homologous arms. Long primer strategy is the traditional one-step method that homologous arms are constructed with 40-100 bp primers. To obtain high quality primers a serial of criteria were taken such as the base composition at 3’ end, the difference of melt temperature, the tendency of forming dimers and the specificity. For ease of use we also developed a web interface for GetPrimers (https://www.evomicslab.org/app/getprimers/). 

# Installation
Linux 64-bit operation system, GCC, g++ compiler, gzip and perl (>=5.8) environment are pre-required. Run the following commands to install GetPrimers and third-party softwares,
```
git clone https://github.com/codeatcg/GetPrimers.git
cd GetPrimers
sh install.sh
```
**Third-party softwares**

* [WindowMasker](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/winmasker/README)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [primer3](https://github.com/primer3-org/primer3)

# User guide
The workflow of GetPrimers was well wrapped. Qualified gene targeting primers can be obtained by just running one-line command. The web interface provided a better visualization of options and results.

* [GetPrimers web application manual](https://github.com/codeatcg/GetPrimers/wiki/GetPrimers-web-application-manual) 
* [GetPrimers command line manual](https://github.com/codeatcg/GetPrimers/wiki/GetPrimers-command-line-manual) 

# Licence
MIT
