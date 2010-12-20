#include "constQ.h"
#include <math.h>

class Chroma: public ConstQ
{
private:
  
  inline float gausswin(float n, float N, float alpha = 2.5) { 
    return exp(-0.5*pow(alpha*(2 * n - N + 1) /(N-1),2));
  }

public:
  Chroma(int _minNote=21, int _maxNote=88, int _bpo=36, 
	 int _fs=44100, float _thresh=0.0054)
    :ConstQ(_minNote, _maxNote, _bpo, _fs, _thresh) {};
  virtual ~Chroma(){};

  void genPitchSpec(float *constQSpec, float *pitchSpec);
  void chromaFold(float *pitchSpec, float *chroma, int numBins = 12);

  void genChroma(float *samples,
		 int numSample, float *chroma);
};
