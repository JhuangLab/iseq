#! /bin/env Rscript
installed.pkgs <- .packages(all.available=T)
for(i in c('optparse', 'futile.logger', 'stringr', 'devtools', 'ngstk', "data.table")){
  if (!i %in% installed.pkgs) {
    install.packages(i)
  }
  devtools::install_github("JhuangLab/BioInstaller")
}
system('pip install . --upgrade')
