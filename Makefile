CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall #-O2
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		kalloc.o blockwfa.o
PROG=		wfa-test
LIBS=		-lz -lpthread -lm
LIBS_WFA2=

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

all:$(PROG)

$(PROG):$(OBJS) main.o
		$(CXX) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

blockwfa.o: blockwfa.h kalloc.h
kalloc.o: kalloc.h
main.o: ketopt.h kalloc.h blockwfa.h kseq.h
