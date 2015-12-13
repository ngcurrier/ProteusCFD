include $(MAKE).opts

EXE_SOLVER = ./bin/ucs.x
EXE_DECOMP = ./bin/udecomp.x
EXE_RECOMP = ./bin/urecomp.x
EXE_FINDPOINT = ./bin/findpoint.x
EXE_PORTOPT = ./bin/opt.x
EXE_CHEMPROPS = ./bin/chemprops.x
EXE_STRUCT_SOLVER = ./bin/csd.x
EXE_STRUCT_ERROR = ./bin/esd.x

EXE_TESTS = ./bin/tests.x

# --------- BEGIN COMMON LIBRARY SECTION

SRCS_COMMON = ./ucs/timer.cpp ./ucs/endian_util.cpp ./ucs/strings_util.cpp \
	./ucs/h5layer.cpp ./ucs/parameterParser.cpp
OBJS_COMMON = $(SRCS_COMMON:.cpp=.o)

# --------- BEGIN FLUID DYNAMICS CORE SECTION

SRCS_SOLVER = ./ucs/main.cpp ./ucs/eqnset.cpp ./ucs/threaded.cpp ./ucs/parallel.cpp \
	./ucs/customics.cpp ./ucs/oddsNends.cpp ./ucs/portFileio.cpp \
	./ucs/elements.cpp ./ucs/dataInfo.cpp ./ucs/derivatives.cpp
SRCS_DECOMP = ./ucs/decomp.cpp ./ucs/mesh.cpp ./ucs/parallel.cpp \
	./ucs/oddsNends.cpp ./ucs/dataInfo.cpp
SRCS_RECOMP = ./ucs/recomp.cpp ./ucs/mesh.cpp ./ucs/parallel.cpp \
	./ucs/oddsNends.cpp ./ucs/dataInfo.cpp
SRCS_FINDPOINT = ./ucs/find_point.cpp ./ucs/mesh.cpp ./ucs/parallel.cpp \
	./ucs/oddsNends.cpp ./ucs/dataInfo.cpp

OBJS_SOLVER = $(SRCS_SOLVER:.cpp=.o)
OBJS_DECOMP = $(SRCS_DECOMP:.cpp=.o)
OBJS_RECOMP = $(SRCS_RECOMP:.cpp=.o)
OBJS_FINDPOINT = $(SRCS_FINDPOINT:.cpp=.o)

# --------- BEGIN CHEMICAL PROPERTIES UTILITY SECTION

SRCS_CHEMPROPS = ./ucs/chemprops.cpp ./ucs/elements.cpp 
OBJS_CHEMPROPS = $(SRCS_CHEMPROPS:.cpp=.o)

# --------- BEGIN OPTIMIZATION TOOLKIT SECTION

CSRCS_PORTOPT = ./ucs/portDriver.cpp ./ucs/portFunctions.cpp ./ucs/portFileio.cpp ./ucs/lineSearch.cpp
FSRCS_PORTOPT = ./ucs/portF.f
OBJS_PORTOPT = $(CSRCS_PORTOPT:.cpp=.o) $(FSRCS_PORTOPT:.cpp=.o)

# --------- BEGIN STRUCTURAL DYNAMICS SECTION

SRCS_STRUCT_SOLVER = ./structuralDynamics/main.cpp ./structuralDynamics/element_lib.cpp \
	./structuralDynamics/structparam.cpp ./structuralDynamics/utility.cpp \
	./structuralDynamics/explicit.cpp ./structuralDynamics/implicit.cpp \
	./structuralDynamics/forces.cpp ./structuralDynamics/bc.cpp ./structuralDynamics/io.cpp \
	./structuralDynamics/fluid_structure.cpp
OBJS_STRUCT_SOLVER = $(SRCS_STRUCT_SOLVER:.cpp=.o)

SRCS_ALL = $(SRCS_SOLVER) ./ucs/decomp.cpp ./ucs/recomp.cpp ./ucs/mesh.cpp ./ucs/find_point.cpp $(CSRCS_PORTOPT)\
	 $(SRCS_CHEMPROPS) $(SRCS_STRUCT_SOLVER) $(SRCS_COMMON)
OBJS_ALL = $(SRCS_ALL:.cpp=.o)

# --------- BEGIN EXECUTABLE TARGETS SECTION

$(EXE_SOLVER):  $(TINYXMLDIR)/libtinyxml.a $(HDF5_LIB)/libhdf5.a $(OBJS_SOLVER) ./ucs/libcommon.a structuralDynamics/libstructdyn.a 
	$(MPICXX) $(LINK_OPTS) -o $(EXE_SOLVER) $(LCXXFLAGS) $(OBJS_SOLVER) $(CXXLIBS) -lstructdyn -ltinyxml

$(EXE_DECOMP): $(METISINSTALLDIR)/libmetis.a $(HDF5_LIB)/libhdf5.a ./ucs/libcommon.a  $(OBJS_DECOMP) 
	$(MPICXX) $(LINK_OPTS) -o $(EXE_DECOMP) $(LCXXFLAGS) -L$(METIS_LIB) $(OBJS_DECOMP) $(CXXLIBS) -lmetis -ltinyxml

$(EXE_RECOMP): $(HDF5_LIB)/libhdf5.a $(OBJS_RECOMP) 
	$(MPICXX) $(LINK_OPTS) -o $(EXE_RECOMP) $(LCXXFLAGS) $(OBJS_RECOMP) $(CXXLIBS) -ltinyxml

$(EXE_FINDPOINT): $(OBJS_FINDPOINT)
	$(MPICXX) $(LINK_OPTS) -o $(EXE_FINDPOINT) $(LCXXFLAGS) $(OBJS_FINDPOINT) $(CXXLIBS) -ltinyxml

$(EXE_PORTOPT): $(OBJS_PORTOPT)
	$(FXX) $(FXX_LINK_OPTS) -o $(EXE_PORTOPT) $(LCXXFLAGS) $(OBJS_PORTOPT) $(FXXLIBS) -lcommon

$(EXE_CHEMPROPS): $(OBJS_CHEMPROPS)
	$(MPICXX) $(LINK_OPTS) -o $(EXE_CHEMPROPS) $(LCXXFLAGS) $(OBJS_CHEMPROPS) $(CXXLIBS)

$(EXE_STRUCT_SOLVER): $(OBJS_STRUCT_SOLVER) ./ucs/libcommon.a
	$(MPICXX) -o $(EXE_STRUCT_SOLVER) $(LCXXFLAGS) $(OBJS_STRUCT_SOLVER) $(CXXLIBS)

$(EXE_STRUCT_ERROR): ./structuralDynamics/error.o ./ucs/libcommon.a
	$(MPICXX) -o $(EXE_STRUCT_ERROR) $(LCXXFLAGS) ./structuralDynamics/error.o $(CXXLIBS)

$(EXE_TESTS): $(GTEST_LIB)/libgtest.la ./unitTest/main.o 
	$(MPICXX) -o $(EXE_TESTS) $(LCXXFLAGS) ./unitTest/main.o $(CXXLIBS) $(GTEST_LIB)/libgtest.a

SRCDIRS = ./structuralDynamics ./ucs
ROOT = $$PWD

all: $(EXE_SOLVER) $(EXE_DECOMP) $(EXE_RECOMP) $(EXE_FINDPOINT) $(EXE_TESTS)\
	$(EXE_CHEMPROPS) $(EXE_STRUCT_SOLVER) $(EXE_STRUCT_ERROR)

tools: $(EXE_DECOMP) $(EXE_RECOMP) $(EXE_CHEMPROPS) $(EXE_FINDPOINT)

./structuralDynamics/libstructdyn.a: $(OBJS_STRUCT_SOLVER) 
	ar rvs ./structuralDynamics/libstructdyn.a $(OBJS_STRUCT_SOLVER)

./ucs/libcommon.a: $(OBJS_COMMON)
	ar rvs ./ucs/libcommon.a $(OBJS_COMMON) 

$(HDF5_LIB)/libhdf5.a: 
	@orig=$$PWD;\
	cd $$PWD/$(HDF5DIR); \
	./configure --prefix $$PWD/build; \
	$(MAKE); \
	$(MAKE) install; \
	$(MAKE) lib;

$(METIS_LIB)/libmetis.a: 
	@orig=$$PWD;\
	cd $$PWD/$(METISDIR); \
	$(MAKE) config; \
	$(MAKE)

$(TINYXML_LIB)/libtinyxml.a:
	@orig=$$PWD;\
	cd $$PWD/$(TINYXMLDIR);\
	$(MAKE)

$(GTEST_LIB)/.libs/libgtest.a:
	@orig=$$PWD;\
	cd $$PWD/$(GTESTDIR); \
	./configure; \
	$(MAKE)

.cpp.o:
	$(MPICXX) -c $< $(CXX_OPTS) $(INCLUDES) -o $@

.f.o: 
	$(FXX) -c $< $(FXX_OPTS) $(INCLUDES) -o $@

clean:
	rm $(OBJS_ALL); \
	rm $(EXE_SOLVER); \
	rm $(EXE_DECOMP); \
	rm $(EXE_RECOMP); \
	rm $(EXE_FINDPOINT); \
	rm $(EXE_PORTOPT); \
	rm $(EXE_CHEMPROPS); \
	rm $(EXE_STRUCT_SOLVER); \
	rm $(EXE_STRUCT_ERROR); \
	rm ./structuralDynamics/error.o; \
	rm ./ucs/libcommon.a; \
	rm ./structuralDynamics/libstructdyn.a; 

cleanall:
	@orig=$$PWD;\
	$(MAKE) clean;\
	cd $$orig/$(HDF5DIR);\
	$(MAKE) clean;\
	$(MAKE) distclean-am;\
	cd $$orig/$(METISDIR);\
	$(MAKE) clean;\
	$(MAKE) distclean;\
	cd $$orig/$(TINYXMLDIR);\
	$(MAKE) clean;\
	rm $$orig/$(METISINSTALLDIR)/configRun;\
	rm $$orig/$(HDF5INSTALLDIR)/configRun

todo: 
	@orig=$$PWD;\
	for dir in $(SRCDIRS);  do\
	  cd $$dir;\
	  pwd;\
	  grep -iRn todo ./*.cpp ./*.h ./*.tcc;\
	  cd $$orig ;\
	done	
dist:
	tar -cvzf src.tgz ./* --exclude='src.tgz'

cleanucs:
	rm $(OBJS_ALL) $(EXE_ALL) Make.depend

depend: .depend

.depend: $(SRCS_SOLVER) $(SRCS_DECOMP) $(SRCS_RECOMP) $(SRCS_FINDPOINT) $(SRCS_CHEMPROPS)
	rm -rf .depend
	$(MPICXX) -MM $(INCLUDES) $(SRCS_SOLVER) $(SRCS_DECOMP) $(SRCS_RECOMP) \
	$(SRCS_FINDPOINT) $(SRCS_CHEMPROPS) > ./.depend

include .depend

