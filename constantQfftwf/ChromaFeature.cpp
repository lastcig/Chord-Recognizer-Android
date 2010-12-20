#include "chroma.h"

JNIEXPORT jlong JNICALL Java_com_taemincho_feature_Chroma_create
(JNIEnv *, jobject, jint minNote, jint maxNote, jint bpo, jint fs, jfloat thresh){
  Chroma *chroma = new Chroma(minNote, maxNote, bpo, fs, thresh);
  return (jlong)chroma;
  return 0;
}

JNIEXPORT void JNICALL Java_com_taemincho_feature_Chroma_destroy
(JNIEnv *, jobject, jlong handle){
  Chroma *This = (Chroma *)handle;
  delete This;
}

JNIEXPORT void JNICALL Java_com_taemincho_feature_Chroma_genChroma
(JNIEnv *, jobject, jlong handle, jfloatArray samples, jint numSample, jfloatArray chroma){
  kiss_fft_scalar *source = 
    (kiss_fft_scalar*)env->GetFloatArrayElements(samples, NULL);
  float *dest = (float*)env->GetFloatArrayElements(chroma, NULL);
  Chroma *This = (Chroma *)handle;
  This->genChroma(source, (int)numSample, dest);
  
}
