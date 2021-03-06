/**
 * @file   constQ.h
 * @author Taemin Cho <taemin@taemin-linux>
 * @date   Fri Nov 26 06:20:19 2010
 * 
 * @brief  
 * 
 * 
 */

#include <math.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include <pthread.h>

typedef struct
{
  kiss_fft_cpx *cpx;
  int s;
  int e;
}sparKernel;

typedef struct
{
  int minNote;
  int maxNote;
  int bpo;
  int fs;
  float thresh;
}constQparm;


class ConstQ
{
 private:
  const int minNote;
  const int maxNote;
  const int bpo;
  const int fs;
  const float thresh;

  static const double pi=
    3.141592653589793238462643383279502884197169399375105820974944;

  kiss_fft_cfg fft;
  kiss_fftr_cfg fftr;

  int constQBins;
  int Q;
  double minFreq;
  double maxFreq;
  float progress;

  sparKernel *sparkernel;
  unsigned int fftLength;
  pthread_t thread;


  static void *sparseKernelThread(void *arg);
  void genSparseKernel();

  inline double noteNum2Freq(float noteNum) {
    return 440 * pow(2.0, (noteNum - 69)/12.0);
  }

  /* Finds next power of two for n. If n itself
     is a power of two then returns n*/
  inline unsigned int nextPowerOf2(unsigned int n)
  {
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
  }   

  inline double hamming(int n, int N)
  {
      return 0.54 - 0.46 * cos(2*pi*n/(N-1));
  }
  
  inline kiss_fft_cpx e(float x)
  {
    kiss_fft_cpx y;
    y.r = cos(x);
    y.i = sin(x);
    return y;
  }

  inline double mag(kiss_fft_cpx x)
  {
    return sqrt(pow(x.r,2) + pow(x.i,2));
  }

 public:
  ConstQ(int _minNote=21, int _maxNote=88, int _bpo=36, 
	 int _fs=44100, float _thresh=0.0054);

  virtual ~ConstQ();
  
  // getters
  int getConstQBins(){ return constQBins; }
  int getBPO() { return bpo; }
  int getMaxNote() { return maxNote; }
  int getMinNote() { return minNote; }
  float getKernelProgress() { return progress;};
  pthread_t getThread() { return thread; }
  unsigned int getFFTLength(){ return fftLength; }

  
  void genConstantQspec(kiss_fft_scalar *samples,
			int numSample, float *constQSpec);



};
