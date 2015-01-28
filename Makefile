include $(MAKE).opts

SRCDIRS = common structuralDynamics ucs

ROOT = $$PWD

all: ucs/ucs tools

tools: ucs/udecomp ucs/urecomp ucs/chemprops.x
	$(MAKE) ucs/udecomp
	$(MAKE) ucs/chemprops.x
	$(MAKE) ucs/urecomp

ucs/ucs: common/libcommon.a structuralDynamics/libstructdyn.a
	@cd ucs;\
	$(MAKE) ucs

ucs/udecomp: common/libcommon.a $(METISINSTALLDIR)/libmetis.a $(HDF5INSTALLDIR)/libhdf5.a 
	@cd ucs;\
	$(MAKE) udecomp

ucs/urecomp: common/libcommon.a $(HDF5INSTALLDIR)/libhdf5.a
	@cd ucs;\
	$(MAKE) urecomp

ucs/chemprops.x: common/libcommon.a
	@cd ucs;\
	$(MAKE) chemprops.x

common/libcommon.a: 
	@cd common;\
	$(MAKE)

structuralDynamics/libstructdyn.a: 
	@cd structuralDynamics;\
	$(MAKE)

$(HDF5INSTALLDIR)/libhdf5.a:
	@orig=$$PWD;\
	cd $$PWD/$(HDF5DIR);\
	env CC=mpicc ./configure --prefix=$(HDF5INSTALLDIR);\
	$(MAKE);\
	$(MAKE) install;\
	$(MAKE) lib

$(METISINSTALLDIR)/libmetis.a:
	@orig=$$PWD;\
	cd $$PWD/$(METISDIR);\
	$(MAKE) config;\
	$(MAKE)

$(TINYXMLDIR)/libtinyxml.a:
	@orig=$$PWD;\
	cd $$PWD/$(TINYXMLDIR);\
	$(MAKE)

clean:
	@orig=$$PWD;\
	for dir in $(SRCDIRS);  do\
	  cd $$dir;\
	  pwd;\
	  $(MAKE) clean;\
	  cd $$orig ;\
	done	

cleanall:
	@orig=$$PWD;\
	$(MAKE) clean;\
	cd $(HDF5DIR);\
	$(MAKE) clean;\
	$(MAKE) distclean-am;\
	cd $$orig;\
	cd $(METISDIR);\
	$(MAKE) clean;\
	$(MAKE) distclean;\
	cd $(TINYXMLDIR);\
	$(MAKE) clean;\
	cd $$orig

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

