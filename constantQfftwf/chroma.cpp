#include "chroma.h"
#include <iostream>  
#include <cstring>
#include <string>
#include <sstream>
#include <time.h>

void Chroma::genPitchSpec(float *constQSpec, float *pitchSpec)
{
  //int bpp = getBPO()/12; 	// bins per pitch

  // float gaussFilt[bpp];
  
  // for(int n=0; n<bpp; n++){
  //   gaussFilt[n] = gausswin(n, bpp);
  // }

  // int offset = bpp/2;
  int idx = 0;
  for(int b=0; b<getConstQBins(); b+=3){
    pitchSpec[idx] += constQSpec[b];
    idx++;
  }

}

void Chroma::chromaFold(float *pitchSpec, float *chroma, int numBins)
{
  int numPitches = getMaxNote()-getMinNote() + 1;
  for(int p=0; p<numPitches; p++){
    int chromaIndex = (p+getMinNote())%numBins;

    chroma[chromaIndex] += pitchSpec[p];
  }
  
  float sum = 0;
  for (int i = 0; i < 12; i++)
      sum += chroma[i];
  if (sum)
    for (int i = 0; i< 12; i++)
      chroma[i]/=sum;
}

void Chroma::genChroma(float *samples,
		       int numSample, float *chroma)
{
  int numPitches = getMaxNote()-getMinNote() + 1;

  float *constQSpec = new float[getConstQBins()];
  memset(constQSpec, 0, sizeof(float)*getConstQBins());
  genConstantQspec(samples, numSample, constQSpec);
  
  float *pitchSpec;
  if(getBPO()==12){
    pitchSpec = constQSpec;
  }else{
    pitchSpec = new float[numPitches];
    memset(pitchSpec, 0, sizeof(float)*numPitches);
    genPitchSpec(constQSpec, pitchSpec);
    delete []constQSpec;
  }

  chromaFold(pitchSpec, chroma);
  delete []pitchSpec;
}

template <class T>
bool from_string(T& t, 
                 const std::string& s, 
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}


int main(int argc, char *argv[])
{
  Chroma *chroma = new Chroma(21, 88, 36, 44100);
  // time_t prev;
  // time(&prev);
  // while(1){
  //   time_t now;
  //   time(&now);
  //   if ( now-prev > 0){
  //     prev = now;
  //     std::cout << chroma->getKernelProgress()<<"\n";
  //   }
  //   if (chroma->getKernelProgress() == 1){
  //     std::cout << chroma->getKernelProgress()<<"\n";
  //     break;
  //   }
  // }

  // if(chroma->getThread() != NULL)
  //  pthread_join(chroma->getThread(), NULL);

  float freq;
  from_string<float>(freq, argv[1], std::dec);

  std::cout << "Freq : " << freq << " hz\n";
  // generate sine wave

  float *signal = new float[4096];
  for (int i=0; i<4096; i++)
    signal[i] = sin(2*M_PI* freq * i/float(44100));

  float c[12] ={0};

  chroma->genChroma(signal, 4096, c);

  for(int i = 0; i<12; i++)
    std::cout << c[i] << "\n";

  delete chroma;
  return 0;
}
