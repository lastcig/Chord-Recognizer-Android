/**
 * @file   constQ.cpp
 * @author Taemin Cho <taemin@taemin-linux>
 * @date   Fri Nov 26 06:20:38 2010
 * 
 * @brief  Generate Constant Q spectrum from Audio segment
 * 
 * 
 */

//#include <jni.h>
#include "constQ.h"
#include <iostream>  
#include <cstring>

/** 
 * Constructor of ConstQ class
 * 
 * Generate sparseKernel under given parameters
 * 
 * @param _minNote : lowest pitch to analyze
 * @param _maxNote : maximum pitch to analyze
 * @param _bpo : bins per octave
 * @param _fs : sample rate
 * @param _thresh : sparseKernel thereshold
 */

ConstQ::ConstQ(int _minNote, int _maxNote, int _bpo, 
	       int _fs, float _thresh)
  :minNote(_minNote), maxNote(_maxNote),
   bpo(_bpo), fs(_fs), thresh(_thresh)
{
  genSparseKernel();
}

/** 
 * Destructor 
 * 
 * clear all allocated memory 
 *
 */

ConstQ::~ConstQ()
{
  free(fft);
  free(fftr);

  for (int k=0; k<constQBins; k++){
    delete []sparkernel[k].cpx;
  }

  delete []sparkernel;
}

/** 
 * Generate sparseKernel
 * 
 * @param minNote 
 * @param maxNote 
 * @param bpo 
 * @param fs 
 * @param thresh 
 * 
 * @return 
 */

void ConstQ::genSparseKernel()
{
  // get frequencies of min and max Note
  minFreq = noteNum2Freq((float)minNote);
  maxFreq = noteNum2Freq((float)maxNote);

  // get Q value and the number of total constant Q bins
  Q = 1./(pow(2.,(1./bpo)) - 1);
  constQBins = ceil(bpo * log2(maxFreq/minFreq));
  
  // memory alloc for sparkernel bins
  sparkernel = new sparKernel[constQBins];
  
  fftLength =
    nextPowerOf2((unsigned int) ceil(Q * fs/minFreq));

  // fft setup
  fft = kiss_fft_alloc(fftLength, 0, NULL, NULL);
  fftr = kiss_fftr_alloc(fftLength, 0, NULL, NULL);

  for (int k=constQBins; k > 0; k--){
    double centerFreq = (minFreq * pow(2,(float(k-1)/bpo)));

    int upperBound = ceil(centerFreq * pow(2,  1./12) * fftLength/fs);
    int lowerBound = floor(centerFreq * pow(2, -1./12) * fftLength/fs);

    sparkernel[k-1].s = lowerBound;
    sparkernel[k-1].e = upperBound;
    sparkernel[k-1].cpx = new kiss_fft_cpx[upperBound - lowerBound];
  }
  pthread_create(&thread, NULL, ConstQ::sparseKernelThread, this);

}

void *ConstQ::sparseKernelThread(void *arg)
{
  ConstQ *This = (ConstQ*) arg;

  unsigned int _fftLength = This->fftLength;
  int _constQBins = This->constQBins;
  int _Q = This->Q;
  double _minFreq = This->minFreq;
  double _maxFreq = This->maxFreq;
  int _fs = This->fs;
  float _thresh = This->thresh;
  int _bpo = This->bpo;

  kiss_fft_cpx tempKernel[_fftLength];
  memset(tempKernel, 0, _fftLength * sizeof(kiss_fft_cpx));

  kiss_fft_cpx specKernel[_fftLength];

  for (int k=_constQBins; k > 0; k--){
    This->progress = 1.0/_constQBins *(_constQBins - k);

    double centerFreq = (_minFreq * pow(2,(float(k-1)/_bpo)));
    int length = ceil(_Q * _fs / centerFreq);
    for(int n=0; n<length; n++){
      kiss_fft_cpx cpx = This->e(2*pi*_Q*n/length);
      double hwin = This->hamming(n, length) / length;
      tempKernel[n].r = (kiss_fft_scalar) hwin * cpx.r;
      tempKernel[n].i = (kiss_fft_scalar) hwin * cpx.i;
    }

    kiss_fft(This->fft, tempKernel, specKernel);

    int upperBound = This->sparkernel[k-1].e;
    int lowerBound = This->sparkernel[k-1].s;

    // cut under the threshold, Conjugate and scaling
    for(int i=lowerBound; i<upperBound; i++){
      if (This->mag(specKernel[i]) > _thresh){
	This->sparkernel[k-1].cpx[i-lowerBound].r = specKernel[i].r / _fftLength;
	This->sparkernel[k-1].cpx[i-lowerBound].i = - specKernel[i].i / _fftLength; // Conjugate
      } else {
	This->sparkernel[k-1].cpx[i-lowerBound].r = 0;
	This->sparkernel[k-1].cpx[i-lowerBound].i = 0;
      }
    }

  }

  This->progress = 1.0;
}


/** 
 * generate Constant Q spectrum from audio segment
 * 
 * @param samples : pointer of segmented audio signal
 * @param numSample : number of audio samples in the segment
 * @param constQSpec : 
 */


void ConstQ::genConstantQspec(kiss_fft_scalar *samples, 
			      int numSample, float *constQSpec)
{
  // alloc input memory for zero padding 
  kiss_fft_scalar timeData[fftLength];
  memset(timeData, 0, fftLength * sizeof(kiss_fft_scalar));
  std::memcpy(timeData, samples, numSample * sizeof(kiss_fft_scalar));  
  
  // alloc output memory
  kiss_fft_cpx freqData[fftLength/2 + 1];
  
  // calculate fft
  kiss_fftr(fftr, timeData, freqData);
    
  for(int n=0; n<constQBins; n++){
    float real = 0;
    float imag = 0;
    for(int i = sparkernel[n].s; i< sparkernel[n].e; i++){
      real += (sparkernel[n].cpx[i-sparkernel[n].s].r * freqData[i].r) 
	- (sparkernel[n].cpx[i-sparkernel[n].s].i * freqData[i].i);
      imag += (sparkernel[n].cpx[i-sparkernel[n].s].r * freqData[i].i) 
	+ (sparkernel[n].cpx[i-sparkernel[n].s].i * freqData[i].r);
    }
    constQSpec[n] = sqrt(real*real + imag*imag);
  }
}
