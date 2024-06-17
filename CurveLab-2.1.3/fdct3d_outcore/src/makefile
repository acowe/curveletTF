include ../../makefile.opt

LIB_SRC = 	fdct3d_param.cpp	fdct3d_forward.cpp	fdct3d_inverse.cpp	cpxcrvletocr.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	test.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libfdct3d.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

test: libfdct3d.a test.o
	${CXX} -o test test.o libfdct3d.a ${LDFLAGS}
	./test options

custom: custom.cpp $(LIB_OBJ)
	cp ${FFTW_DIR}/fftw/.libs/libfftw.a .
	${CXX} custom.cpp -o custom ${INCLUDES} libfdct3d.a libfftw.a

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.out libfdct3d.a test

tags:
	etags *.hpp *.cpp
