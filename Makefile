CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2 -march=native
CXXFLAGS=	$(CFLAGS) -std=c++14
CPPFLAGS=
INCLUDES=
OBJS=		kalloc.o miniwfa.o mwf-dbg.o
PROG=		test-mwf test-wfa test-wfalm
LIBS=		-lz -lpthread -lm

WFA_ROOT=WFA2-lib

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:test-mwf test-wfalm

test-mwf:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

test-wfa:main-wfa.c
		$(CC) -I$(WFA_ROOT) $(CFLAGS) $< -o $@ -L$(WFA_ROOT)/lib -lwfa $(LIBS)

test-wfalm:main-wfalm.cpp
		$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

kalloc.o: kalloc.h
main-wfa.o: ketopt.h kseq.h
main.o: ketopt.h kalloc.h miniwfa.h kseq.h
miniwfa.o: miniwfa.h kalloc.h
mwf-dbg.o: miniwfa.h
main-wfalm.o: lib/wfa_lm.hpp ketopt.h kseq.h
