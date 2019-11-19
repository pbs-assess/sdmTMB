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
	$(R) CMD check --as-cran $(TARBALL)

build-parallel:
	rsync -av --exclude='build' --exclude='.git' --exclude='*.tar.gz' --exclude='*.o' --exclude='*.so' --exclude='inst/*.rds' --exclude='.Rproj.user' . build
	echo "SHLIB_OPENMP_CFLAGS=-Xpreprocessor -fopenmp" >> build/src/Makevars
	echo "SHLIB_OPENMP_CXXFLAGS=-Xpreprocessor -fopenmp" >> build/src/Makevars
	sed -i "" "s/Type jnll = 0;/parallel_accumulator<Type> jnll(this);/g" build/src/sdmTMB.cpp

install-parallel:
	make build-parallel
	$(R) CMD INSTALL --preclean --no-multiarch --with-keep.source .
	rm -rf build

install-parallel-quick:
	make build-parallel
	$(R) CMD INSTALL --install-tests --no-docs --no-multiarch --no-demo .
	rm -rf build
