#include "chroma.h"
#include <iostream>  
#include <cstring>
#include <time.h>

void Chroma::genPitchSpec(float *constQSpec, float *pitchSpec)
{
  int bpp = getBPO()/12; 	// bins per pitch

  float gaussFilt[bpp];
  
  for(int n=0; n<bpp; n++){
    gaussFilt[n] = gausswin(n, bpp);
  }

  int offset = bpp/2;
  for(int b=0; b<getConstQBins(); b++){
    pitchSpec[(b+offset)/bpp] += constQSpec[b] * gaussFilt[(b+offset)%bpp];
  }

}

void Chroma::chromaFold(float *pitchSpec, float *chroma, int numBins)
{
  int numPitches = getMaxNote()-getMinNote() + 1;
  for(int p=0; p<numPitches; p++){
    int chromaIndex = (p+getMinNote())%numBins;

    chroma[chromaIndex] += pitchSpec[p];
  }
}

void Chroma::genChroma(kiss_fft_scalar *samples,
		       int numSample, float *chroma)
{
  int numPitches = getMaxNote()-getMinNote() + 1;

  float constQSpec[getConstQBins()];
  memset(constQSpec, 0, sizeof(float)*getConstQBins());
  genConstantQspec(samples, numSample, constQSpec);

  float pitchSpec[numPitches];
  memset(pitchSpec, 0, sizeof(float)*numPitches);
  genPitchSpec(constQSpec, pitchSpec);

  chromaFold(pitchSpec, chroma);

}


int main(void){
  Chroma *chroma = new Chroma(21, 88, 36, 44100);
  time_t prev;
  time(&prev);
  while(1){
    time_t now;
    time(&now);
    if ( now-prev > 0){
      prev = now;
      std::cout << chroma->getKernelProgress()<<"\n";
    }
    if (chroma->getKernelProgress() == 1){
      std::cout << chroma->getKernelProgress()<<"\n";
      break;
    }
  }

  pthread_join(chroma->getThread(), NULL);
  // generate sine wave
  double pi=
    3.141592653589793238462643383279502884197169399375105820974944;

  float *signal = new float[4096];
  for (int i=0; i<4096; i++)
    signal[i] = sin(2* pi* 440 * i/float(44100));

  float c[12] ={0};
  chroma->genChroma(signal, 4096, c);
  

  for(int i = 0; i<12; i++)
    std::cout << c[i] << "\n";



  delete chroma;
  return 0;
}
