#include "constQ.h"

#include <iostream>  

ConstQ::ConstQ(int _minNote, int _maxNote, int _bpo, 
	       int _fs, float _thresh)
  :minNote(_minNote), maxNote(_maxNote),
   bpo(_bpo), fs(_fs), thresh(_thresh)
{
  sparseKernel(minNote, maxNote, bpo, fs, thresh);
}

ConstQ::~ConstQ()
{
  for (int k=0; k<constQBins; k++){
    delete []sparkernel[k].cpx;
    //    std::cout << k <<"clear mem\n";
  }

  delete []sparkernel;
}

int ConstQ::sparseKernel(int minNote, int maxNote, int bpo,
		  int fs, float thresh)
{
  double minFreq(noteNum2Freq((float)minNote));
  double maxFreq(noteNum2Freq((float)maxNote));

  float Q = 1./(pow(2.,(1./bpo)) - 1);
  constQBins = ceil(bpo * log2(maxFreq/minFreq));
  
  // memory alloc for sparkernel bins
  sparkernel = new sparKernel[constQBins];
  
  unsigned int fftLength =
    nextPowerOf2((unsigned int) ceil(Q * fs/minFreq));

  // fft setup
  fft = kiss_fft_alloc(fftLength, 0, NULL, NULL);

  kiss_fft_cpx *tempKernel = new kiss_fft_cpx[fftLength];
  // error check
  if (tempKernel==NULL){
    std::cout << "memory alloc failure\n";
    return -1;
  }

  kiss_fft_cpx *specKernel = new kiss_fft_cpx[fftLength];
  // error check
  if (specKernel==NULL){
    std::cout << "memory alloc failure\n";
    return -1;
  }
  
  for (int k=constQBins; k > 0; k--){
    double centerFreq = (minFreq * pow(2,(float(k-1)/bpo)));
    int length = ceil(Q * fs / centerFreq);
    for(int n=0; n<length; n++){
      kiss_fft_cpx cpx = e(2*pi*Q*n/length);
      double hwin = hamming(n, length) / length;
      tempKernel[n].r = (kiss_fft_scalar) hwin * cpx.r;
      tempKernel[n].i = (kiss_fft_scalar) hwin * cpx.i;
    }

    kiss_fft(fft, tempKernel, specKernel);

    int upperBound = ceil(centerFreq * pow(2,  1./12) * fftLength/fs);
    int lowerBound = floor(centerFreq * pow(2, -1./12) * fftLength/fs);

    sparkernel[k-1].numComp = upperBound;
    sparkernel[k-1].cpx = new kiss_fft_cpx[upperBound];

    // cut under the threshold, Conjugate and scaling
    for(int i=lowerBound; i<upperBound; i++){
      if (mag(specKernel[i]) > thresh){
	sparkernel[k-1].cpx[i].r = specKernel[i].r / fftLength;
	sparkernel[k-1].cpx[i].i = - specKernel[i].i / fftLength; // Conjugate
      }
      if (k == 1)
	std::cout << sparkernel[k-1].cpx[i].r << " " << sparkernel[k-1].cpx[i].i << "\n";
    }


  }

  delete []tempKernel;
  delete []specKernel;

}

// create  (JNIEnv *, jobject, jint numSamples)
// {
// 	KissFFT* fft = new KissFFT();
// 	fft->config = kiss_fftr_alloc(numSamples,0,NULL,NULL);
// 	fft->spectrum = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * numSamples);
// 	fft->numSamples = numSamples;
// 	return (jlong)fft;
// 	return 0;
// }

// JNIEXPORT void JNICALL Java_com_badlogic_gdx_audio_analysis_KissFFT_destroy(JNIEnv *, jobject, jlong handle)
// {
// 	KissFFT* fft = (KissFFT*)handle;
// 	free(fft->config);
// 	free(fft->spectrum);
// 	free(fft);
// }


int main(void){
  ConstQ *constQ = new ConstQ();
  //  std::cout << constQ->bpo;
  
  delete constQ;
  return 0;
}



