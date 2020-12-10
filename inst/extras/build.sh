#~/bin/R-patched/bin/R CMD BATCH document.R
~/bin/R-4.0.0/bin/R CMD build ../../ --resave-data #--no-examples  --no-build-vignettes 
~/bin/R-4.0.0/bin/R CMD check mia_0.0.0.9010.tar.gz #--no-build-vignettes --no-examples
~/bin/R-4.0.0/bin/R CMD BiocCheck mia_0.0.0.9010.tar.gz
~/bin/R-4.0.0/bin/R CMD INSTALL mia_0.0.0.9010.tar.gz 



