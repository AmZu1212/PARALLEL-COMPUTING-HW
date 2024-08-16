#include <immintrin.h> // for AVX
#include <iostream>

using namespace std;


int main(){

     float array[8] = {0, 1, 2, 3, 4, 5, 6, 7};
     float onesArray[8] = {1, 1, 1, 1, 1, 1, 1, 1};
     float result;
     __m256 oneVector = _mm256_loadu_ps(onesArray);
     oneVector = _mm256_hadd_ps(oneVector, oneVector);
     oneVector = _mm256_hadd_ps(oneVector, oneVector);
     oneVector = _mm256_hadd_ps(oneVector, oneVector);
     result = _mm256_cvtss_f32(oneVector);
     cout << result << endl;

     for (size_t i = 0; i < 8; i++)
     {
          cout << array[i]  << " | ";
     }
     cout << endl;

     __m256 myVector = _mm256_loadu_ps(array);
     myVector = _mm256_mul_ps(myVector, myVector);
     float *newArray = (float*) &myVector; 

     __m256 myVector2 = _mm256_loadu_ps(array);
     myVector2 = _mm256_add_ps(myVector2, myVector2);
     float *newArray2 = (float*) &myVector2; 



     __m256 myVector3 = _mm256_loadu_ps(array);
     myVector3 = _mm256_sub_ps(myVector3, myVector3);
     float *newArray3 = (float*) &myVector3;



     __m256 myVector4 = _mm256_loadu_ps(array);
     myVector4 = _mm256_sqrt_ps(myVector4);
     float *newArray4 = (float*) &myVector4;



     for (size_t i = 0; i < 8; i++)
     {
          cout << newArray[i]  << " | ";
     }
     cout << endl;

     for (size_t i = 0; i < 8; i++)
     {
          cout << newArray2[i] << " | ";
     }
     cout << endl;

     for (size_t i = 0; i < 8; i++)
     {
          cout << newArray3[i]  << " | ";
     }
     cout << endl;

     for (size_t i = 0; i < 8; i++)
     {
          cout << newArray4[i]  << " | ";
     }
     cout << endl;

    return 0;
}