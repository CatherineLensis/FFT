// FFT.h
#ifndef FFT_H
#define FFT_H
#include <complex>
#include <vector>
#include <random>
#include <ctime>
#include <map>

typedef std::complex<double> Complex;

class FFT {
public:
    // Основная функция БПФ для смешанных радиксов
    static void fft(std::vector<std::complex<double>>& data);

private:
    // Факторизация числа на радиксы 2, 3, 5
    static std::map<int, int> factorize(int N);

    // БПФ для радикса-2
    static void fft_radix2(std::vector<std::complex<double>>& data);

    // БПФ для радикса-3
    static void fft_radix3(std::vector<std::complex<double>>& data);

    // БПФ для радикса-5
    static void fft_radix5(std::vector<std::complex<double>>& data);

    // Рекурсивная функция для смешанных радиксов
    static void fft_mixed_radix(std::vector<std::complex<double>>& data, const std::map<int, int>& factors);

    // Перестановка данных для смешанных радиксов
    static void mixed_radix_permute(std::vector<std::complex<double>>& data, int N, const std::map<int, int>& factors);
};

#endif // FFT_H
