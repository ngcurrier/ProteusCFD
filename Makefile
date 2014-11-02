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

ucs/udecomp: common/libcommon.a
	@cd ucs;\
	$(MAKE) udecomp

ucs/urecomp: common/libcommon.a
	@cd ucs;\
	$(MAKE) urecomp

ucs/chemprops.x: common/libcommon.a
	@cd ucs;\
	$(MAKE) chemprops.x

config:
	@orig=$$PWD;\
	cd $$PWD/$(HDF5DIR);\
	env CC=mpicc ./configure --prefix=$(HDF5INSTALLDIR);\
	cd $$PWD/$(METISDIR);\
	$(MAKE) config;\
	cd $$orig

common/libcommon.a: force_look
	@cd common;\
	$(MAKE)

structuralDynamics/libstructdyn.a: force_look
	@cd structuralDynamics;\
	$(MAKE)

libhdf5.a:
	@orig=$$PWD;\
	cd $$PWD/$(HDF5DIR);\
	$(MAKE);\
	$(MAKE) install;\
	$(MAKE) lib

libmetis.a:
	@orig=$$PWD;\
	cd $$PWD/$(METISDIR);\
	$(MAKE)

libtinyxml.a:
	@orig=$$PWD;\
	cd $$PWD/$(TINYXMLDIR);\
	$(MAKE)

tpls:
	$(MAKE) libhdf5.a;\
	$(MAKE) libmetis.a;\
	$(MAKE) libtinyxml.a

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

force_look:
	true

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

