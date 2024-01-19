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
NGS_INCL := -I $(NGS_DIR)/include
NGS_LIB := -L $(NGS_DIR)/lib64

CC = c++ -std=c++11 -fdiagnostics-color=never
CFLAGS = -Wall -Wno-format-y2k  -pthread -fPIC -O3 -finline-functions -fstrict-aliasing \
         -fomit-frame-pointer -msse4.2 $(BOOST_INCL) $(NGS_INCL)

PLATFORM=$(shell uname -s)

ifeq ($(PLATFORM),Linux)

LIBS = $(NGS_LIB) -lncbi-ngs-c++-static -lngs-c++-static -lncbi-ngs-static -lncbi-vdb-static \
       -Wl,-Bstatic $(BOOST_LIB) \
       -lboost_program_options \
       -lboost_iostreams \
       -lboost_regex \
       -lboost_timer \
       -lboost_chrono \
       -lboost_system \
       -Wl,-Bdynamic -lrt -ldl -lm  -lpthread -lz

else

LIBS = $(BOOST_LIB) \
       -lboost_program_options \
       -lboost_iostreams \
       -lboost_regex \
       -lboost_timer \
       -lboost_chrono \
       -lboost_system \
       -ldl -lm -lpthread -lz

endif

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

binaries=skesa saute saute_prot gfa_connector kmercounter

all: $(binaries)

.PHONY: clean

clean: 
	rm -f $(binaries) *.o $(NGS_DIR)/ngs.done

glb_align.o: glb_align.hpp Makefile

skesa.o: common_util.hpp concurrenthash.hpp readsgetter.hpp ngs_includes.hpp counter.hpp graphdigger.hpp assembler.hpp KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp glb_align.hpp Makefile $(NGS_DIR)/ngs.done
skesa: skesa.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)

saute.o: guidedassembler.hpp guidedgraph.hpp gfa.hpp common_util.hpp guidedpath_naa.hpp readsgetter.hpp ngs_includes.hpp concurrenthash.hpp counter.hpp graphdigger.hpp KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp Makefile $(NGS_DIR)/ngs.done
saute: saute.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)

saute_prot.o: guidedassembler.hpp guidedgraph.hpp genetic_code.hpp gfa.hpp common_util.hpp guidedpath_naa.hpp readsgetter.hpp ngs_includes.hpp concurrenthash.hpp counter.hpp graphdigger.hpp KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp Makefile $(NGS_DIR)/ngs.done
saute_prot: saute_prot.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)

gfa_connector.o: gfa.hpp glb_align.hpp common_util.hpp concurrenthash.hpp readsgetter.hpp ngs_includes.hpp graphdigger.hpp KmerInit.hpp  DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp Makefile $(NGS_DIR)/ngs.done
gfa_connector: gfa_connector.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)

kmercounter.o: common_util.hpp concurrenthash.hpp readsgetter.hpp ngs_includes.hpp counter.hpp graphdigger.hpp assembler.hpp KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp glb_align.hpp Makefile $(NGS_DIR)/ngs.done
kmercounter: kmercounter.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)

$(NGS_DIR)/ngs.done:
	rm -fr $(NGS_DIR)
	mkdir -p $(NGS_DIR)/build
	mkdir $(NGS_DIR)/out
	cd $(NGS_DIR); git clone https://github.com/ncbi/ncbi-vdb.git
	cd $(NGS_DIR); git clone https://github.com/ncbi/sra-tools.git
	cd $(NGS_DIR)/ncbi-vdb; ./configure --prefix=$(NGS_DIR)/out --build-prefix=$(NGS_DIR)/build; make; make install
	cd $(NGS_DIR)/sra-tools; ./configure --prefix=$(NGS_DIR)/out --build-prefix=$(NGS_DIR)/build
	cd $(NGS_DIR)/sra-tools/ngs/ngs-sdk; make; make install
	cd $(NGS_DIR)/out; mv lib64 include ../
	cd $(NGS_DIR); rm -rf out build sra-tools ncbi-vdb
	touch $@
