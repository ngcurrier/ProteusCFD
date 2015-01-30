include $(MAKE).opts
include $(MAKE).ucs
include $(MAKE).common

SRCDIRS = common structuralDynamics ucs

ROOT = $$PWD

all: ucs/ucs ucs/decomp ucs/chemprops.x ucs/urecomp


tools: ucs/udecomp ucs/urecomp ucs/chemprops.x
	$(MAKE) ucs/udecomp
	$(MAKE) ucs/chemprops.x
	$(MAKE) ucs/urecomp

structuralDynamics/libstructdyn.a: 
	@cd structuralDynamics;\
	$(MAKE)

$(HDF5INSTALLDIR)/configRun:
	cd $$PWD/$(HDF5INSTALLDIR);\
	env CC=mpicc ./configure --prefix=$(HDF5INSTALLDIR);\
	touch configRun

$(HDF5INSTALLDIR)/libhdf5.a: $(HDF5INSTALLDIR)/configRun
	@orig=$$PWD;\
	cd $$PWD/$(HDF5DIR);\
	$(MAKE);\
	$(MAKE) install;\
	$(MAKE) lib

$(METISINSTALLDIR)/configRun:
	cd $$PWD/$(METISINSTALLDIR);\
	$(MAKE) config;\
	touch configRun

$(METISINSTALLDIR)/libmetis.a: $(METISINSTALLDIR)/configRun
	@orig=$$PWD;\
	cd $$PWD/$(METISDIR);\
	$(MAKE)

$(TINYXMLDIR)/libtinyxml.a:
	@orig=$$PWD;\
	cd $$PWD/$(TINYXMLDIR);\
	$(MAKE)

.cpp.o:
	$(MPICXX) -c $< $(CXX_OPTS) $(INCLUDES)

.f.o: 
	$(FXX) -c $< $(FXX_OPTS) $(INCLUDES)


clean:
	@orig=$$PWD;\
	for dir in $(SRCDIRS);  do\
	  cd $$dir;\
	  pwd;\
	  $(MAKE) clean;\
	  cd $$orig ;\
	done;\
	cleancommon;\
	cleanucs

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

