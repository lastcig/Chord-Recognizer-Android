LOCAL_CFLAGS = -c -O2 -Wall 
#  -DFIXED_POINT
TARGET = test
OBJS = chroma.o constQ.o
CC = $(LOCAL_CFLAGS)
LIBS = -lm -lstdc++ -lfftw3 -lpthread

$(TARGET) : $(OBJS)
	g++ -o $@ $^ $(LIBS)
.cpp.o:
	gcc $(CC) $<

chroma.o: chroma.cpp chroma.h constQ.cpp constQ.h
constQ.o: constQ.cpp constQ.h
clean:
	rm *.o
	rm test


