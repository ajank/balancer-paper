### Custom code used for the analysis in Ghavi-Helm et al., 2019

Custom code used for the computational analyses in Ghavi-Helm*, Jankowski*, Meiers* et al. paper "Highly rearranged chromosomes reveal uncoupling between genome topology and gene expression". The project involved whole genome sequencing, mate-pair, RNA-seq, Hi-C and Capture-C data from highly rearranged balancer chromosomes in _Drosophila melanogaster_.

The software is provided "as is", in the hope that it will be useful, but without warranty of any kind. For historical reasons, the code is organized into following subdirectories:
  * `drosophila_balancer` – authored by Sascha Meiers, involving SNV and SV calling from whole-genome and mate pair sequencing data, as well as RNA-seq analysis
  * `Capture-C` – authored by Aleksander Jankowski, involving Capture-C data analysis
  * `Hi-C` – authored by Aleksander Jankowski, involving Hi-C data analysis, as well as Hi-C, Capture-C and other data plotting
  * `Hi-C-ggbio` – authored by Sascha Meiers and Aleksander Jankowski, for plots combining Hi-C contact maps, gene expression and gene models using ggbio package.

The contents of these subdirectories were accordingly located in:
  * `/g/korbel/shared/projects/drosophila_balancer`
  * `/g/furlong/project/37_Capture-C`
  * `/g/furlong/project/33_Hi-C`
  * `/g/furlong/project/33_Hi-C/Hi-C-ggbio`.

Software dependencies are listed in the Life Sciences Reporting Summary linked to the article.

In case of any questions, please contact Aleksander Jankowski <aleksander.jankowski(at)embl.de>.
