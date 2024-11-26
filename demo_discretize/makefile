CXXC = g++
CXXFLAGS = -std=c++11
opt = 0

ifeq ($(opt),0)
	incpath = -I ../../../mfem/debug-install/include  -I ../../../mfem//metis/include -I /usr/include/suitesparse/
	libpath = -L ../../../mfem/debug-install/lib -lmfem -L ../../../mfem/metis/lib -lmetis -Wl,-rpath,/usr/lib/x86_64-linux-gnu/ -L/usr/lib/x86_64-linux-gnu/ -lklu -lbtf -lumfpack -lcholmod -lcolamd -lamd -lcamd -lccolamd -lsuitesparseconfig -llapack -lblas -lrt
	CXXFLAGS += -g -Wall
else
	incpath = -I ../../../mfem/install/include  -I ../../../mfem//metis/include -I /usr/include/suitesparse/
	libpath = -L ../../../mfem/install/lib -lmfem -L ../../../mfem/metis/lib -lmetis -Wl,-rpath,/usr/lib/x86_64-linux-gnu/ -L/usr/lib/x86_64-linux-gnu/ -lklu -lbtf -lumfpack -lcholmod -lcolamd -lamd -lcamd -lccolamd -lsuitesparseconfig -llapack -lblas -lrt
	CXXFLAGS += -O3
endif


cavity:cavity.cpp 
	${CXXC} ${CXXFLAGS}  ${incpath} $^ ${libpath} -lmfem -o $@  

.PHONY : clean
clean:
	rm -f *.o cavity
