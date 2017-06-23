#! /bin/env Rscript
installed.pkgs <- .packages(all.available=T)
for(i in c('optparse', 'futile.logger', 'stringr')){
  if (!i %in% installed.pkgs) {
    install.packages(i)
  }
}
system('pip install . --upgrade')
