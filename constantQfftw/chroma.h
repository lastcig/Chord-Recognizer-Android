#include "constQ.h"
#include <math.h>

class Chroma: public ConstQ
{
private:
  
  inline double gausswin(double n, double N, double alpha = 2.5) { 
    return exp(-0.5*pow(alpha*(2 * n - N + 1) /(N-1),2));
  }

public:
  Chroma(int _minNote=21, int _maxNote=88, int _bpo=36, 
	 int _fs=44100, double _thresh=0.0054)
    :ConstQ(_minNote, _maxNote, _bpo, _fs, _thresh) {};
  virtual ~Chroma(){};

  void genPitchSpec(double *constQSpec, double *pitchSpec);
  void chromaFold(double *pitchSpec, double *chroma, int numBins = 12);

  void genChroma(double *samples,
		 int numSample, double *chroma);
};
