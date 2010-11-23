LOCAL_CFLAGS = -c -O2 -Wall -ffast-math -fomit-frame-pointer -unroll-loops -dA -fverbose-asm
#  -DFIXED_POINT
TARGET = test
OBJS = constQ.o kiss_fft.o
CC = $(LOCAL_CFLAGS)

$(TARGET) : $(OBJS)
	g++ -lm -o $@ $^
.c.o:
	gcc $(CC) $<

constQ.o: constQ.cpp constQ.h
# kiss_fftr.o: kiss_fft.c kiss_fftr.c kiss_fft.h kiss_fftr.h _kiss_fft_guts.h
kiss_fft.o: kiss_fft.c kiss_fft.h _kiss_fft_guts.h
clean:
	rm *.o
	rm test


