# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================

#setenv BOOST_PATH /netopt/ncbi_tools64/boost-1.62.0-ncbi1

ifdef BOOST_PATH
BOOST_INCL := -I $(BOOST_PATH)/include
BOOST_LIB := -L $(BOOST_PATH)/lib
endif

NGS_DIR := $(CURDIR)/NGS
VDB_PATH := $(NGS_DIR)/vdb_out
NGS_PATH := $(NGS_DIR)/ngs_out
BUILD_PATH := $(NGS_DIR)/build

VDB_INCL := -I $(VDB_PATH)/include
VDB_LIB := -L $(VDB_PATH)/lib64
NGS_INCL := -I $(NGS_PATH)/include
NGS_LIB := -L $(NGS_PATH)/lib64

CC = c++ -std=c++11 -fdiagnostics-color=never
CFLAGS = -Wall -Wno-format-y2k  -pthread -fPIC -O3 -finline-functions -fstrict-aliasing \
         -fomit-frame-pointer -msse4.2 $(BOOST_INCL) $(NGS_INCL) $(VDB_INCL)

LIBS = $(VDB_LIB) -lncbi-ngs-c++-static -lncbi-vdb-static \
       $(NGS_LIB) -lngs-c++-static \
       -Wl,-Bstatic $(BOOST_LIB) \
       -lboost_program_options \
       -lboost_iostreams \
       -lboost_regex \
       -lboost_timer \
       -lboost_chrono \
       -lboost_system \
       -Wl,-Bdynamic -lrt -ldl -lm  -lpthread -lz

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

binaries=skesa gfa_connector kmercounter

all: $(binaries)

.PHONY: clean

clean: 
	rm -f $(binaries) *.o $(NGS_DIR)/ngs.done

glb_align.o: glb_align.hpp Makefile

skesa.o: common_util.hpp concurrenthash.hpp readsgetter.hpp ngs_includes.hpp counter.hpp graphdigger.hpp assembler.hpp KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp glb_align.hpp Makefile $(NGS_DIR)/ngs.done
skesa: skesa.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)

gfa_connector.o: gfa.hpp glb_align.hpp common_util.hpp concurrenthash.hpp readsgetter.hpp ngs_includes.hpp graphdigger.hpp KmerInit.hpp  DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp Makefile $(NGS_DIR)/ngs.done
gfa_connector: gfa_connector.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)

kmercounter.o: common_util.hpp concurrenthash.hpp readsgetter.hpp ngs_includes.hpp counter.hpp graphdigger.hpp assembler.hpp KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp glb_align.hpp Makefile $(NGS_DIR)/ngs.done
kmercounter: kmercounter.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)


$(NGS_DIR)/ngs.done:
	rm -fr $(NGS_DIR)
	mkdir -p $(NGS_DIR)/ngs
	mkdir $(BUILD_PATH)
	mkdir $(NGS_PATH)
	mkdir $(VDB_PATH)
	cd $(NGS_DIR)/ngs; git init; git remote add -f origin https://github.com/ncbi/ngs.git; git config core.sparseCheckout true; echo "ngs-sdk" >> .git/info/sparse-checkout; git pull origin master
	cd $(NGS_DIR)/ngs/ngs-sdk; ./configure --prefix=$(NGS_PATH) --build-prefix=$(BUILD_PATH); make; make install
	cd $(NGS_DIR); git clone https://github.com/ncbi/ncbi-vdb.git
	cd $(NGS_DIR)/ncbi-vdb; ./configure --prefix=$(VDB_PATH) --build-prefix=$(BUILD_PATH); make; make install
	touch $@
