all: libbdsurf.a

include ../makefile.opt

LIB_SRC = 	bdsurf.cpp    ccsubmatlib.cpp  ccsurfop.cpp  gpmesh.cpp  trackball.cpp	bdsurfvr.cpp  ccsurf.cpp       ccsurfvr.cpp  IImage.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	tt.cpp

MEX_SRC = 	fmm3d_mex.cpp

DEP     = 	$(LIB_SRC:.cpp=.d)

libbdsurf.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) ../libbdsurf.a $(LIB_OBJ)
	$(RANLIB) ../libbdsurf.a

tt0: libbdsurf.a tt0.o
	${CXX} -o tt0 tt0.o ../libbdsurf.a ${LDFLAGS}

tt1: libbdsurf.a tt1.o
	${CXX} -o tt1 tt1.o ../libbdsurf.a ${LDFLAGS}

-include $(DEP)

