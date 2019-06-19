# to add new objects, just put them in this list and run make depend

OBJ = bid.o BidSet.o Distribution.o arbitrary.o Legacy.o main.o matching.o normal.o Param.o paths.o regions.o scheduling.o polyModel.o featureCalc.o


########################################################
# UNCOMMENT/EDIT THESE TO COMPILE WITH CPLEX LIBRARIES
########################################################

INCLUDE = -I/usr/local/ilog/cplex80/include/ilcplex/ -I/usr/local/ilog/cplex80/include/
LIB = -lnsl -pthread -lm -lcplex
LIBDIRS = -L/usr/local/ilog/cplex80/lib/i86_linux2_glibc2.2_gcc3.0/static_pic_mt/

# NOTE: REMOVE -DLINUX IN NEXT 2 LINES FOR COMPILATION IN UNIX
# ADD: -DWIN32 AND -D_WIN32 FOR COMPILATION IN WINDOWS
RELEASEFLAGS = -O5 -DNDEBUG -DWIN32 -D_WIN32 -DUSE_CPLEX
DEBUGFLAGS = -Wall -DDEBUG -g -DWIN32 -D_WIN32 -DUSE_CPLEX

########################################################
# UNCOMMENT/EDIT THESE TO COMPILE WITH LPSOLVE
########################################################

#INCLUDE = -Ilp_solve_4.0
#LIB = -lm -llpk
#LIBDIRS = -Llp_solve_4.0

# NOTE: REMOVE -DLINUX IN NEXT 2 LINES FOR COMPILATION IN UNIX
# ADD: -DWIN32 AND -D_WIN32 FOR COMPILATION IN WINDOWS
#RELEASEFLAGS = -O5 -DNDEBUG -DLINUX
#DEBUGFLAGS = -Wall -DDEBUG -g -DLINUX
#LPMAKE = cd lp_solve_4.0; make;
#LPCLEAN = cd lp_solve_4.0; make clean;

############################################################
release:
	cd obj; make all "OBJ = ${OBJ}" "CPPFLAGS = ${RELEASEFLAGS} ${INCLUDE}";
	${LPMAKE}
	make all "VPATH = obj" "FLAGS = -o cats ${RELEASEFLAGS}"

debug:
	cd debugobj; make all "OBJ = ${OBJ}" "CPPFLAGS = ${DEBUGFLAGS} ${INCLUDE}";
	${LPMAKE}
	make all "VPATH = debugobj" "FLAGS = -o debugcats ${DEBUGFLAGS}"

all : $(OBJ)
	g++ $^ ${FLAGS} ${LIBDIRS} $(LIB)

.PHONY: clean depend again

clean:
	${LPCLEAN}
	rm *obj/*.o cats >& /dev/null;

depend:
	makedepend -Y -fobj/Makefile *.cc *.cpp # >& /dev/null
	makedepend -Y -fdebugobj/Makefile *.cc *.cpp # >& /dev/null

again:
	make clean; make depend; make
# DO NOT DELETE
