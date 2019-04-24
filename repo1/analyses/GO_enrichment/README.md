# quick GO term analysis

I take the files gene ids (generated from `1.get_geneids.sh`) and input them online at <http://www.geneontology.org/page/go-enrichment-analysis> to compare them against the control.

Here are the results as of May 4 2017


### All ASE genes

```
Analysis Type:	PANTHER Overrepresentation Test (release 20170413)
Annotation Version and Release Date:	GO Ontology database  Released 2017-04-24
Analyzed List:	FBgn.all_ase.txt (Drosophila melanogaster)
Reference List:	FBgn.control.txt (Drosophila melanogaster)
Bonferroni correction:	true
Bonferroni count:	2333
GO biological process complete	FBgn.control.txt - REFLIST (5016)	FBgn.all_ase.txt (452)	FBgn.all_ase.txt (expected)	FBgn.all_ase.txt (over/under)	FBgn.all_ase.txt (fold Enrichment)	FBgn.all_ase.txt (P-value)
```
GO term | FBgn.control | test | expected | +/- | fold | p-value
--- |--- |--- |--- |--- |--- | ---
Unclassified (UNCLASSIFIED)|704|68|63.44|+|1.07|0.00E00
biological regulation (GO:0065007)|2092|143|188.51|-|.76|1.54E-02
regulation of biological process (GO:0050789)|1929|131|173.83|-|.75|3.61E-02
regulation of cellular metabolic process (GO:0031323)|1015|56|91.46|-|.61|1.82E-02
regulation of metabolic process (GO:0019222)|1123|60|101.20|-|.59|1.40E-03
regulation of nitrogen compound metabolic process (GO:0051171)|954|48|85.97|-|.56|1.91E-03
regulation of macromolecule metabolic process (GO:0060255)|1022|51|92.09|-|.55|4.90E-04
regulation of primary metabolic process (GO:0080090)|1002|49|90.29|-|.54|3.28E-04
regulation of gene expression (GO:0010468)|808|39|72.81|-|.54|5.36E-03

### ASE genes with FC > 1.5

```
Analysis Type:	PANTHER Overrepresentation Test (release 20170413)
Annotation Version and Release Date:	GO Ontology database  Released 2017-04-24
Analyzed List:	FBgn.ase_with_fc1.5.txt (Drosophila melanogaster)
Reference List:	FBgn.control.txt (Drosophila melanogaster)
Bonferroni correction:	true
Bonferroni count:	2333
GO biological process complete	FBgn.control.txt - REFLIST (5016)	FBgn.ase_with_fc1.5.txt (314)	FBgn.ase_with_fc1.5.txt (expected)	FBgn.ase_with_fc1.5.txt (over/under)	FBgn.ase_with_fc1.5.txt (fold Enrichment)	FBgn.ase_with_fc1.5.txt (P-value)
```
GO term | FBgn.control | test | expected | +/- | fold | p-value
--- |--- |--- |--- |--- |--- | ---
oxidation-reduction process (GO:0055114)|216|31|13.52|+|2.29|4.51E-02
Unclassified (UNCLASSIFIED)|704|56|44.07|+|1.27|0.00E00
**cellular component organization** (GO:0016043)|1392|50|87.14|-|.57|1.32E-03
-> cellular process (GO:0009987)|3151|159|197.25|-|.81|1.66E-02
-> cellular component organization or biogenesis (GO:0071840)|1467|51|91.83|-|.56|1.49E-04
**regulation of RNA metabolic process** (GO:0051252)|658|18|41.19|-|.44|3.53E-02
-> regulation of macromolecule metabolic process (GO:0060255)|1022|28|63.98|-|.44|6.56E-05
-> -> regulation of metabolic process (GO:0019222)|1123|33|70.30|-|.47|8.15E-05
-> -> -> regulation of biological process (GO:0050789)|1929|73|120.75|-|.60|1.65E-05
-> -> -> -> biological regulation (GO:0065007)|2092|82|130.96|-|.63|1.43E-05
-> regulation of nucleobase-containing compound metabolic process (GO:0019219)|686|19|42.94|-|.44|2.76E-02
-> -> regulation of cellular metabolic process (GO:0031323)|1015|31|63.54|-|.49|1.35E-03
-> -> -> regulation of cellular process (GO:0050794)|1782|69|111.55|-|.62|3.16E-04
-> -> regulation of primary metabolic process (GO:0080090)|1002|28|62.72|-|.45|1.60E-04
-> -> regulation of nitrogen compound metabolic process (GO:0051171)|954|28|59.72|-|.47|1.26E-03
**gene expression** (GO:0010467)|624|16|39.06|-|.41|2.24E-02
animal organ development (GO:0048513)|842|21|52.71|-|.40|2.10E-04
regulation of gene expression (GO:0010468)|808|22|50.58|-|.43|2.58E-03
-> anatomical structure development (GO:0048856)|1745|68|109.24|-|.62|6.43E-04
-> -> developmental process (GO:0032502)|1815|72|113.62|-|.63|6.91E-04
-> system development (GO:0048731)|1373|53|85.95|-|.62|1.99E-02
-> -> multicellular organism development (GO:0007275)|1596|61|99.91|-|.61|1.41E-03
-> -> -> single-organism developmental process (GO:0044767)|1796|71|112.43|-|.63|7.19E-04
macromolecule localization (GO:0033036)|397|7|24.85|-|.28|3.44E-02
epithelial tube morphogenesis (GO:0060562)|347|4|21.72|-| < 0.2|5.67E-03
-> tube morphogenesis (GO:0035239)|367|5|22.97|-|.22|9.52E-03



### ASE genes up-regulated on the balancer

```
Analysis Type:	PANTHER Overrepresentation Test (release 20170413)
Annotation Version and Release Date:	GO Ontology database  Released 2017-04-24
Analyzed List:	FBgn.all_up-regulated.txt (Drosophila melanogaster)
Reference List:	FBgn.control.txt (Drosophila melanogaster)
Bonferroni correction:	true
Bonferroni count:	2333
GO biological process complete	FBgn.control.txt - REFLIST (5016)	FBgn.all_up-regulated.txt (234)	FBgn.all_up-regulated.txt (expected)	FBgn.all_up-regulated.txt (over/under)	FBgn.all_up-regulated.txt (fold Enrichment)	FBgn.all_up-regulated.txt (P-value)
```
No enrichment

### ASE genes down-regulated on the balancer

```
Analysis Type:	PANTHER Overrepresentation Test (release 20170413)
Annotation Version and Release Date:	GO Ontology database  Released 2017-04-24
Analyzed List:	FBgn.all_down-regulated.txt (Drosophila melanogaster)
Reference List:	FBgn.control.txt (Drosophila melanogaster)
Bonferroni correction:	true
Bonferroni count:	2333
GO biological process complete	FBgn.control.txt - REFLIST (5016)	FBgn.all_down-regulated.txt (218)	FBgn.all_down-regulated.txt (expected)	FBgn.all_down-regulated.txt (over/under)	FBgn.all_down-regulated.txt (fold Enrichment)	FBgn.all_down-regulated.txt (P-value)
```

GO term | FBgn.control | test | expected | +/- | fold | p-value
--- |--- |--- |--- |--- |--- | ---
response to hypoxia (GO:0001666)|23|8|1.00|+|8.00|2.16E-02
response to decreased oxygen levels (GO:0036293)|24|8|1.04|+|7.67|2.92E-02
response to oxygen levels (GO:0070482)|27|9|1.17|+|7.67|8.36E-03
Unclassified (UNCLASSIFIED)|704|37|30.60|+|1.21|0.00E00