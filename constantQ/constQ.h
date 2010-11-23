
#include <math.h>
#include "kiss_fft.h"

typedef struct
{
  kiss_fft_cpx *cpx;
  int numComp;
}sparKernel;

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
  int constQBins;
  sparKernel *sparkernel;

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

  int sparseKernel(int minNote, int maxNote, int bpo,
	       int fs, float thresh);

 public:
  ConstQ(int _minNote=21, int _maxNote=88, int _bpo=36, 
	 int _fs=44100, float _thresh=0.0054);

  virtual ~ConstQ();

};
