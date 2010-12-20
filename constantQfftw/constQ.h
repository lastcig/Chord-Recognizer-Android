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
#include <pthread.h>
#include <complex>
#include "fftw3.h"

typedef struct
{
  fftw_complex *cpx;
  int s;
  int e;
}sparKernel;

typedef struct
{
  float *bank;
  int s;
  int e;
}filterBank;

typedef struct
{
  double r;
  double i;
}complex_t;

class ConstQ
{
 private:
  const int minNote;
  const int maxNote;
  const int bpo;
  const int fs;
  const double thresh;

  fftw_complex *cin, *cout;
  double *fin;
  fftw_plan c2c, r2c;

  int constQBins;
  int Q;
  double minFreq;
  double maxFreq;
  double progress;

  filterBank *filterbanks;
  sparKernel *sparkernel;
  unsigned int fftLength;
  pthread_t thread;


  static void *sparseKernelThread(void *arg);
  void genSparseKernel();

  inline double noteNum2Freq(double noteNum) {
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
    return (double) 0.54 - 0.46 * cos(2*M_PI*n/(N-1));
  }

  void genFilterBanks();  
 public:
  ConstQ(int _minNote=21, int _maxNote=88, int _bpo=36, 
	 int _fs=44100, double _thresh=0.0054);

  virtual ~ConstQ();
  
  // getters
  int getConstQBins(){ return constQBins; }
  int getBPO() { return bpo; }
  int getMaxNote() { return maxNote; }
  int getMinNote() { return minNote; }
  double getKernelProgress() { return progress;};
  pthread_t getThread() { return thread; }
  unsigned int getFFTLength(){ return fftLength; }

  
  void genConstantQspec(double *samples,
			int numSample, double *constQSpec);


};
