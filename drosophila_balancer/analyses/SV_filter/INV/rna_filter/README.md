# Can I use RNA-seq data to validate inversions?

I tried to do this in an automized manner by counting +/- RNA-seq 
reads inside and left and right of the inversion (see `Snakefile`).

But then I realized that the inversions are ~10kb, meaning they
span multiple genes. Hence the whole analysis is *BS* :(

