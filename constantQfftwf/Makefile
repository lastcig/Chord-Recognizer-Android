LOCAL_CFLAGS = -c -O3
#  -DFIXED_POINT
TARGET = test
OBJS = chroma.o constQ.o
CC = $(LOCAL_CFLAGS)
LIBS = -lm -lstdc++ -lfftw3f -lpthread

test : $(OBJS)
	g++ -o $@ $^ $(LIBS)

test_neon : chroma_neon.o constQ_neon.o
	g++ -o $@ $^ $(LIBS)

.cpp.o:
	gcc $(CC) $<

chroma.o: chroma.cpp chroma.h constQ.cpp constQ.h
chroma_neon.o: chroma_neon.cpp chroma_neon.h constQ_neon.cpp constQ_neon.h
constQ_neon.o: constQ_neon.cpp constQ_neon.h
constQ.o: constQ.cpp constQ.h
clean:
	rm *.o
	rm test


