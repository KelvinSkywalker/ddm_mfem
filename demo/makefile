CXXC = clang++
obj = ddmsolver.o
LIBRARY = ../lib/libddm.a
LIB_SOURCES = ddmsolver.cpp
LIB_OBJECTS = $(LIB_SOURCES:.cpp=.o)
-include ./makefile.in

all: $(LIBRARY)

$(LIBRARY): $(LIB_OBJECTS)
	ar rcs $@ $^

$(obj):%.o:%.cpp 
	$(CXXC) $(CXXFLAGS) ${incpath}  -c $^ -o $@


ex1:ex1.cpp 
	${CXXC} ${CXXFLAGS} ${incpath} $^ ${libpath} -o $@

testASM:testASM.cpp
	${CXXC} ${CXXFLAGS} ${incpath} $^ -L../lib -lddm ${libpath} -o $@
 

assemblebug:assemblebug.cpp
	${CXXC} ${CXXFLAGS} ${incpath} $^ -L../lib -lddm ${libpath} -o $@  	

twoleveltest:twoleveltest.cpp
	${CXXC} ${CXXFLAGS} ${incpath} $^ -L../lib -lddm ${libpath} -o $@

test:TestDDM.cpp
	${CXXC} ${CXXFLAGS} ${incpath} $^ -L../lib -lddm ${libpath} -o $@   

Impedancetest:Impedancetest.cpp
	${CXXC} ${CXXFLAGS} ${incpath} $^ -L../lib -lddm ${libpath} -o $@

impedance3Dtest:impedance3Dtest.cpp
	${CXXC} ${CXXFLAGS} ${incpath} $^ -L../lib -lddm ${libpath} -o $@

cavitytest:cavitytest.cpp
	${CXXC} ${CXXFLAGS} ${incpath} $^ -L../lib -lddm ${libpath} -o $@

shifted-precondition:shifted-precondition.cpp
	${CXXC} ${CXXFLAGS} ${incpath} $^ -L../lib -lddm ${libpath} -o $@

mixpoisson : mixpoisson.cpp
	${CXXC} ${CXXFLAGS} ${incpath} $^ -L../lib -lddm ${libpath} -o $@
	
.PHONY : clean
clean:
	rm -f *.o ex1 testASM test libddm.a
