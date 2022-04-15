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

build-parallel:
	rsync -qav --exclude='build' --exclude='.git' --exclude='*.tar.gz' --exclude='*.o' --exclude='*.so' --exclude='inst/*.rds' --exclude='.Rproj.user' . build
	echo 'PKG_CXXFLAGS += $$(SHLIB_OPENMP_CXXFLAGS)' >> build/src/Makevars
	echo 'PKG_LIBS = $$(SHLIB_OPENMP_CXXFLAGS)' >> build/src/Makevars
	sed -i "" "s/Type jnll = 0;/parallel_accumulator<Type> jnll(this);/g" build/src/sdmTMB.cpp

install-parallel:
	make build-parallel
	$(R) CMD INSTALL --preclean --no-multiarch --with-keep.source build
	rm -rf build
