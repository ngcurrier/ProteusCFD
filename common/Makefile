include ../make.opts

SRCS_SOLVER = timer.cpp endian_util.cpp strings_util.cpp h5layer.cpp parameterParser.cpp
OBJS_SOLVER = $(SRCS_SOLVER:.cpp=.o)

all: libcommon.a 

libcommon.a: $(OBJS_SOLVER)
		ar rvs libcommon.a $(OBJS_SOLVER) 

.cpp.o:
	$(CXX) -c $< $(CXX_OPTS) $(INCLUDES)

clean: 
	rm $(OBJS_SOLVER) libcommon.a