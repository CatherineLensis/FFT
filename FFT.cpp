#include "FFT.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

const double PI = 3.14159265358979323846;

std::map<int, int> FFT::factorize(int N) {
    std::map<int, int> factors;
    while (N % 2 == 0) {
        factors[2]++;
        N /= 2;
    }
    while (N % 3 == 0) {
        factors[3]++;
        N /= 3;
    }
    while (N % 5 == 0) {
        factors[5]++;
        N /= 5;
    }
    if (N != 1) {
        throw std::invalid_argument("N must be a product of 2, 3, and 5 only");
    }
    return factors;
}

void FFT::fft_radix2(std::vector<std::complex<double>>& data) {
    int N = data.size();
    if (N <= 1) return;

    std::vector<std::complex<double>> even(N / 2);
    std::vector<std::complex<double>> odd(N / 2);
    for (int i = 0; i < N / 2; i++) {
        even[i] = data[2 * i];
        odd[i] = data[2 * i + 1];
    }

    fft_radix2(even);
    fft_radix2(odd);

    for (int k = 0; k < N / 2; k++) {
        std::complex<double> t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        data[k] = even[k] + t;
        data[k + N / 2] = even[k] - t;
    }
}

void FFT::fft_radix3(std::vector<std::complex<double>>& data) {
    int N = data.size();
    if (N <= 1) return;

    int N1 = N / 3;
    std::vector<std::complex<double>> f0(N1), f1(N1), f2(N1);
    for (int i = 0; i < N1; i++) {
        f0[i] = data[i * 3];
        f1[i] = data[i * 3 + 1];
        f2[i] = data[i * 3 + 2];
    }

    fft_radix3(f0);
    fft_radix3(f1);
    fft_radix3(f2);

    std::complex<double> w = std::polar(1.0, -2 * PI / 3);
    for (int k = 0; k < N1; k++) {
        std::complex<double> wk = std::polar(1.0, -2 * PI * k / N);
        std::complex<double> t1 = wk * f1[k];
        std::complex<double> t2 = wk * wk * f2[k];
        data[k] = f0[k] + t1 + t2;
        data[k + N1] = f0[k] + w * t1 + w * w * t2;
        data[k + 2 * N1] = f0[k] + w * w * t1 + w * t2;
    }
}

void FFT::fft_radix5(std::vector<std::complex<double>>& data) {
    int N = data.size();
    if (N <= 1) return;

    int N1 = N / 5;
    std::vector<std::complex<double>> f[5];
    for (int i = 0; i < 5; i++) {
        f[i].resize(N1);
        for (int j = 0; j < N1; j++) {
            f[i][j] = data[j * 5 + i];
        }
        fft_radix5(f[i]);
    }

    for (int k = 0; k < N1; k++) {
        std::complex<double> wk = std::polar(1.0, -2 * PI * k / N);
        std::complex<double> w[5] = {1.0, wk, wk * wk, wk * wk * wk, wk * wk * wk * wk};
        std::complex<double> t[5];
        for (int i = 0; i < 5; i++) {
            t[i] = w[i] * f[i][k];
        }
        data[k] = f[0][k] + t[1] + t[2] + t[3] + t[4];
        data[k + N1] = f[0][k] + w[1] * t[1] + w[2] * t[2] + w[3] * t[3] + w[4] * t[4];
        data[k + 2 * N1] = f[0][k] + w[2] * t[1] + w[4] * t[2] + w[1] * t[3] + w[3] * t[4];
        data[k + 3 * N1] = f[0][k] + w[3] * t[1] + w[1] * t[2] + w[4] * t[3] + w[2] * t[4];
        data[k + 4 * N1] = f[0][k] + w[4] * t[1] + w[3] * t[2] + w[2] * t[3] + w[1] * t[4];
    }
}

void FFT::mixed_radix_permute(std::vector<std::complex<double>>& data, int N, const std::map<int, int>& factors) {
    // Простая перестановка для демонстрации; для полной корректности нужна обобщенная версия
    std::vector<std::complex<double>> temp = data;
    int stride = 1;
    for (auto [radix, power] : factors) {
        int size = std::pow(radix, power);
        for (int i = 0; i < N; i++) {
            int idx = 0, t = i;
            for (int j = 0; j < power; j++) {
                idx = idx * radix + (t % radix);
                t /= radix;
            }
            idx = idx * (N / size) + (i / size);
            data[i] = temp[idx];
        }
        stride *= size;
    }
}

void FFT::fft_mixed_radix(std::vector<std::complex<double>>& data, const std::map<int, int>& factors) {
    int N = data.size();
    if (N == 1) return;

    if (factors.at(2) > 0) {
        int size = 1 << factors.at(2);
        std::vector<std::vector<std::complex<double>>> sub_data(N / size, std::vector<std::complex<double>>(size));
        for (int i = 0; i < N / size; i++) {
            for (int j = 0; j < size; j++) {
                sub_data[i][j] = data[i + j * (N / size)];
            }
            fft_radix2(sub_data[i]);
        }
        for (int i = 0; i < N / size; i++) {
            for (int j = 0; j < size; j++) {
                data[i + j * (N / size)] = sub_data[i][j];
            }
        }
    } else if (factors.at(3) > 0) {
        int size = std::pow(3, factors.at(3));
        std::vector<std::vector<std::complex<double>>> sub_data(N / size, std::vector<std::complex<double>>(size));
        for (int i = 0; i < N / size; i++) {
            for (int j = 0; j < size; j++) {
                sub_data[i][j] = data[i + j * (N / size)];
            }
            fft_radix3(sub_data[i]);
        }
        for (int i = 0; i < N / size; i++) {
            for (int j = 0; j < size; j++) {
                data[i + j * (N / size)] = sub_data[i][j];
            }
        }
    } else if (factors.at(5) > 0) {
        int size = std::pow(5, factors.at(5));
        std::vector<std::vector<std::complex<double>>> sub_data(N / size, std::vector<std::complex<double>>(size));
        for (int i = 0; i < N / size; i++) {
            for (int j = 0; j < size; j++) {
                sub_data[i][j] = data[i + j * (N / size)];
            }
            fft_radix5(sub_data[i]);
        }
        for (int i = 0; i < N / size; i++) {
            for (int j = 0; j < size; j++) {
                data[i + j * (N / size)] = sub_data[i][j];
            }
        }
    }
}

void FFT::fft(std::vector<std::complex<double>>& data) {
    int N = data.size();
    auto factors = factorize(N);
    mixed_radix_permute(data, N, factors);
    fft_mixed_radix(data, factors);
}
