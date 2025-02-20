#ifndef FFT_H
#define FFT_H

#include <complex>
#include <vector>
#include <random>
#include <ctime>

typedef std::complex<double> Complex;

class FFT {
public:
    static void transform(std::vector<Complex>& data, int isInverse);

private:
    static void fftRadix2(std::vector<Complex>& data, int isInverse);
    static void fftMixedRadix(std::vector<Complex>& data, int isInverse);
};

// Функция для генерации случайных чисел (метод суммы 12)
class RandomGenerator {
public:
    static double Sum12();
};

#endif // FFT_H
