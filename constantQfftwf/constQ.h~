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
  fftwf_complex *cpx;
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
  float r;
  float i;
}complex_t;

class ConstQ
{
 private:
  const int minNote;
  const int maxNote;
  const int bpo;
  const int fs;
  const float thresh;

  fftwf_complex *cin, *cout;
  float *fin;
  fftwf_plan c2c, r2c;

  int constQBins;
  int Q;
  double minFreq;
  double maxFreq;
  float progress;

  filterBank *filterbanks;
  sparKernel *sparkernel;
  bool isKernelUsed;
  unsigned int fftLength;
  pthread_t thread;
  
  static void *sparseKernelThread(void *arg);
  void genSparseKernel();

  void genFilterBanks();  

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

  inline float hamming(int n, int N)
  {
    return (float) 0.54 - 0.46 * cos(2*M_PI*n/(N-1));
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

  
  void genConstantQspec(float *samples,
			int numSample, float *constQSpec);


};
