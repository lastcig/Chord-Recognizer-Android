/**
 * @file   constQ.cpp
 * @author Taemin Cho <taemin@taemin-linux>
 * @date   Fri Nov 26 06:20:38 2010
 * 
 * @brief  Generate Constant Q spectrum from Audio segment
 * 
 * 
 */

#include "constQ_neon.h"
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
  std::cout<<"constructed\n";
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
  fftwf_free(cout);
  fftwf_free(fin);
  fftwf_destroy_plan(r2c);

  for (int k=0; k<constQBins; k++){
    fftwf_free(sparkernel[k].cpx);
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
  anglestep = 2*M_PI*Q;

  // memory alloc for sparkernel bins
  sparkernel = new sparKernel[constQBins];
  
  fftLength =
    nextPowerOf2((unsigned int) ceil(Q * fs/minFreq));

  // fftw setup
  cin = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * fftLength);
  memset(cin, 0, sizeof(fftwf_complex) * fftLength);
  cout = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * fftLength);
  c2c = fftwf_plan_dft_1d(fftLength, cin, cout, FFTW_FORWARD, FFTW_ESTIMATE);

  fin = (float *) fftwf_malloc(sizeof(float) * fftLength);
  r2c = fftwf_plan_dft_r2c_1d(fftLength, fin, cout, FFTW_ESTIMATE);

  centerFreqz = new float[constQBins];

  for (int k=constQBins; k > 0; k--){

    centerFreqz[k-1] = (float) (minFreq * pow(2,(float(k-1)/bpo)));

    sparkernel[k-1].s = floor(centerFreqz[k] * pow(2, -1./12) * fftLength/fs);
    sparkernel[k-1].e = ceil(centerFreqz[k] * pow(2,  1./12) * fftLength/fs);

    std::cout<<"k: "<< k << "\n";

    sparkernel[k-1].cpx = (fftwf_complex *) 
      fftwf_malloc(sizeof(fftwf_complex) * (sparkernel[k-1].e - sparkernel[k-1].s));

    //new complex_t[sparkernel[k-1].e - sparkernel[k-1].s];

  }

  pthread_create(&thread, NULL, ConstQ::sparseKernelThread, this);

  // std::complex <float> *tempKernel;
  // tempKernel = reinterpret_cast<std::complex <float> *> (cin); 

  // std::complex <float> *specKernel;
  // specKernel = reinterpret_cast<std::complex <float> *> (cout); 

  // std::complex <float> *sparKernel;

  // for (int k=constQBins; k > 0; k--){
  //   progress = 1.0/constQBins * (constQBins - k);

  //   float centerFreq = (minFreq * pow(2,(float(k-1)/bpo)));
  //   int length = ceil(Q * fs / centerFreq);

  //   for(int n=0; n<length; n++){
  //     std::complex <float> cpx(0, (2*M_PI*Q*n/length));
  //     cpx = exp(cpx);
  //     float hwin = HAMMING(n, length) / length;
  //     tempKernel[n] = hwin * cpx;
  //   }

  //   fftwf_execute(c2c);

  //   int upperBound = ceil(centerFreq * pow(2,  1./12) * fftLength/fs);
  //   int lowerBound = floor(centerFreq * pow(2, -1./12) * fftLength/fs);

  //   sparkernel[k-1].s = lowerBound;
  //   sparkernel[k-1].e = upperBound;
  //   sparkernel[k-1].cpx = (fftwf_complex *) 
  //     fftwf_malloc(sizeof(fftwf_complex) * (upperBound - lowerBound));

  //   sparKernel = 
  //     reinterpret_cast<std::complex <float> *> (sparkernel[k-1].cpx);

  //   // cut under the threshold, Conjugate and scaling
  //   for(int i=lowerBound; i<upperBound; i++){
  //     if (abs(specKernel[i]) > thresh)
  // 	sparKernel[i-lowerBound] = conj(specKernel[i] / (float)fftLength);
  //     else
  // 	sparKernel[i-lowerBound] = 0;
  //   }

  // }
  
  // fftwf_destroy_plan(c2c); // no more needed
  // progress = 1.0;

}

void *ConstQ::sparseKernelThread(void *arg)
{
  ConstQ *This = (ConstQ*) arg;

  float _fftLength = (float) This->fftLength;
  int _constQBins = This->constQBins;
  int _Q = This->Q;
  float _anglestep = This->anglestep;
  float *_centerFreqz = This->centerFreqz;
  int _fs = This->fs;
  float _thresh = This->thresh;

  // std::complex <float> *tempKernel
  //   = reinterpret_cast<std::complex <float> *> (This->cin); 

  complex_t *tempKernel = (complex_t *) This->cin;

  //float *tempKernel = (float *) This->cin;

  // complex_t *specKernel = (complex_t*) This->cout; 

  std::cout<<"thread\n";
  std::complex <float> *specKernel;
  specKernel = reinterpret_cast<std::complex <float> *> (This->cout); 

  for (int k=_constQBins; k > 0; k--){
    This->progress = 1.0/_constQBins * (_constQBins - k);

    int length = ceil(_Q * _fs / _centerFreqz[k-1]);

    for(int n=0; n<length; n++){
      float hwin = HAMMING(n, length) / length;
      tempKernel[n].r = hwin * cos(_anglestep*n/length);
      tempKernel[n].i = hwin * sin(_anglestep*n/length);
    }

    fftwf_execute(This->c2c);

    int upperBound = This->sparkernel[k-1].e;
    int lowerBound = This->sparkernel[k-1].s;

    std::complex <float> *sparKernel;
    sparKernel = 
      reinterpret_cast<std::complex <float> *> (This->sparkernel[k-1].cpx);

    // cut under the threshold, Conjugate and scaling
    for(int i=lowerBound; i<upperBound; i++){
      if (abs(specKernel[i]) > _thresh)
    	sparKernel[i-lowerBound] = conj(specKernel[i] / _fftLength);
      else
    	sparKernel[i-lowerBound] = 0;
    }

    // complex_t *kernel = (complex_t *) This->sparkernel[k-1].cpx;

    // // cut under the threshold, Conjugate and scaling
    // for(int i=lowerBound; i<upperBound; i++){
    //   int idx = i-lowerBound;
    //   if (ABS(specKernel[i]) > _thresh){
    // 	kernel[idx].r = specKernel[i].r / _fftLength;
    // 	kernel[idx].i = - specKernel[i].i / _fftLength; // Conjugate
    //   }
    //   else{
    // 	kernel[idx].r = 0;
    // 	kernel[idx].i = 0;
    //   }
    // }

  }
  
  fftwf_destroy_plan(This->c2c); // no more needed
  fftwf_free(This->cin);
  delete []_centerFreqz;
  This->progress = 1.0;

  return NULL;
}


/** 
 * generate Constant Q spectrum from audio segment
 * 
 * @param samples : pointer of segmented audio signal
 * @param numSample : number of audio samples in the segment
 * @param constQSpec : 
 */


void ConstQ::genConstantQspec(float *samples, 
			      int numSample, float *constQSpec)
{
  // clean input memory for zero padding 
  memset(fin, 0, fftLength * sizeof(float));
  std::memcpy(fin, samples, numSample * sizeof(float));  
  
  // fft calculate
  fftwf_execute(r2c);

  std::complex <float> *freqData;
  freqData = 
    reinterpret_cast<std::complex <float> *> (cout);
    
  for(int n=0; n<constQBins; n++){
    std::complex <float> *kernel;
    kernel = 
      reinterpret_cast<std::complex <float> *> (sparkernel[n].cpx);

    std::complex <float> temp(0, 0);

    for(int i = sparkernel[n].s; i< sparkernel[n].e; i++){
      temp += kernel[i-sparkernel[n].s] * freqData[i];
    }

    constQSpec[n] = abs(temp);
  }

}
