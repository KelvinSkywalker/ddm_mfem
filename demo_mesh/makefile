CXXC = g++
CXXFLAGS = -std=c++11
opt = 0

ifeq ($(opt),0)
	incpath = -I ../../../mfem/debug-install/include  -I ../../../mfem/metis/include -I /Applications/anaconda/anaconda3/envs/ddm4py/include
	libpath = -L ../../../mfem/debug-install/lib -lmfem -L ../../../mfem/metis/lib -lmetis -Wl,-rpath,/Applications/anaconda/anaconda3/envs/ddm4py/lib -L//Applications/anaconda/anaconda3/envs/ddm4py/lib -lklu -lbtf -lumfpack -lcholmod -lcolamd -lamd -lcamd -lccolamd -lsuitesparseconfig -llapack -lblas
	CXXFLAGS += -g -Wall
else
	incpath = -I ../../../mfem/debug-install/include  -I ../../../mfem/metis/include -I /Applications/anaconda/anaconda3/envs/ddm4py/include
	libpath = -L ../../../mfem/debug-install/lib -lmfem -L ../../../mfem/metis/lib -lmetis -Wl,-rpath,/Applications/anaconda/anaconda3/envs/ddm4py/lib -L//Applications/anaconda/anaconda3/envs/ddm4py/lib -lklu -lbtf -lumfpack -lcholmod -lcolamd -lamd -lcamd -lccolamd -lsuitesparseconfig  -llapack -lblas
	CXXFLAGS += -O3
endif


mesh:test_mesh.cpp 
	${CXXC} ${CXXFLAGS}  ${incpath} $^ ${libpath} -o $@  

.PHONY : clean
clean:
	rm -f *.o mesh
