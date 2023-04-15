all: libbdsurf.a

include ../makefile.opt

LIB_SRC = 	vecmatop.cpp	svdrep.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

DEP     = 	$(LIB_SRC:.cpp=.d)

libbdsurf.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) ../libbdsurf.a $(LIB_OBJ)
	$(RANLIB) ../libbdsurf.a

-include $(DEP)

