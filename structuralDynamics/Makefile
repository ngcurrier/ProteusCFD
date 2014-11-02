include ../make.opts

EXE_SOLVER = csd
EXE_ERROR = esd

SRCS_SOLVER = main.cpp element_lib.cpp structparam.cpp utility.cpp explicit.cpp \
	      implicit.cpp forces.cpp bc.cpp io.cpp fluid_structure.cpp

OBJS_SOLVER = $(SRCS_SOLVER:.cpp=.o)

all: $(EXE_SOLVER) $(EXE_ERROR) libstructdyn.a 

libstructdyn.a: $(OBJS_SOLVER) $(EXE_SOLVER)
		ar rvs libstructdyn.a $(OBJS_SOLVER) 


$(EXE_SOLVER): $(OBJS_SOLVER)
	$(MPICXX) -o $(EXE_SOLVER) $(LCXXFLAGS) $(OBJS_SOLVER) $(CXXLIBS)

$(EXE_ERROR): error.o
	$(MPICXX) -o $(EXE_ERROR) $(LCXXFLAGS) error.o $(CXXLIBS)

.cpp.o:
	$(MPICXX) -c $< $(CXX_OPTS) $(INCLUDES)

clean: 
	rm $(EXE_SOLVER) $(EXE_ERROR) error.o $(OBJS_SOLVER) libstructdyn.a