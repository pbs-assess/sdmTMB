PACKAGE=sdmTMB
VERSION := $(shell sed -n '/^Version: /s///p' DESCRIPTION)
DATE := $(shell sed -n '/^Date: /s///p' DESCRIPTION)
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

# Allow e.g. "make R=R-devel install"
R=R

all:
	make doc-update
	make build-package
	make install

doc-update:
	echo "roxygen2::roxygenize(\".\")" | $(R) --slave

build-package:
	$(R) CMD build --no-build-vignettes --no-manual .

install:
	$(R) CMD INSTALL --preclean --no-multiarch --with-keep.source .

cran-check:
	echo "devtools::check(\".\")" | $(R) --slave
