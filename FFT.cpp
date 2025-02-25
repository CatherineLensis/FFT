#include "FFT.h"
#include <cmath>
#include <iostream>

#define M_PI 3.14159265358979323846

double RandomGenerator::Sum12() {
    static std::random_device rd;
    static std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    double sum = 0.0;
    for (int i = 0; i < 12; i++) {
        sum += dist(rng);
    }
    return sum / 12.0;
}

void FFT::transform(std::vector<Complex>& data, int isInverse) {
    size_t n = data.size();

    // Проверка на кратность 2, 3, 5
    size_t temp = n;
    while (temp > 1) {
        if (temp % 2 == 0) temp /= 2;
        else if (temp % 3 == 0) temp /= 3;
        else if (temp % 5 == 0) temp /= 5;
        else {
            std::cerr << "FFT не поддерживает длину " << n << " (должна быть кратна 2, 3 или 5)\n";
            return;
        }
    }

    if ((n & (n - 1)) == 0) {
        fftRadix2(data, isInverse);
    }
    else {
        fftMixedRadix(data, isInverse);
    }

    // Масштабирование для обратного преобразования
    if (isInverse) {
        double scale = 1.0 / static_cast<double>(n);
        for (auto& x : data) x *= scale;
    }
}

void FFT::fftRadix2(std::vector<Complex>& data, int isInverse) {
    size_t n = data.size();
    if (n <= 1) return;

    // Перестановка бит
    size_t logN = std::log2(n);
    for (size_t i = 0, j = 0; i < n; ++i) {
        if (i < j) std::swap(data[i], data[j]);
        size_t mask = n >> 1;
        while (j & mask) {
            j ^= mask;
            mask >>= 1;
        }
        j ^= mask;
    }

    // Алгоритм Cooley-Tukey
    for (size_t len = 2; len <= n; len *= 2) {
        double angle = (isInverse ? 2.0 : -2.0) * M_PI / len; // Поменяли знак
        Complex wLen(cos(angle), sin(angle));
        for (size_t i = 0; i < n; i += len) {
            Complex w(1);
            for (size_t j = 0; j < len / 2; ++j) {
                Complex u = data[i + j];
                Complex v = data[i + j + len / 2] * w;
                data[i + j] = u + v;
                data[i + j + len / 2] = u - v;
                w *= wLen;
            }
        }
    }
}

void FFT::fftMixedRadix(std::vector<Complex>& data, int isInverse) {
    size_t n = data.size();
    if (n <= 1) return;

    if (n % 5 == 0) {
        size_t m = n / 5;
        std::vector<Complex> x0(m), x1(m), x2(m), x3(m), x4(m);

        // Разделение на подмассивы
        for (size_t i = 0; i < m; ++i) {
            x0[i] = data[i * 5];
            x1[i] = data[i * 5 + 1];
            x2[i] = data[i * 5 + 2];
            x3[i] = data[i * 5 + 3];
            x4[i] = data[i * 5 + 4];
        }

        // Рекурсивное применение БПФ
        fftMixedRadix(x0, isInverse);
        fftMixedRadix(x1, isInverse);
        fftMixedRadix(x2, isInverse);
        fftMixedRadix(x3, isInverse);
        fftMixedRadix(x4, isInverse);

        // Комбинирование результатов (бабочка радикса-5)
        for (size_t k = 0; k < m; ++k) {
            double angle = (isInverse ? 2.0 : -2.0) * M_PI * k / n;
            Complex w(cos(angle), sin(angle));
            Complex w2 = w * w;
            Complex w3 = w2 * w;
            Complex w4 = w2 * w2;

            Complex t0 = x0[k];
            Complex t1 = x1[k] * w;
            Complex t2 = x2[k] * w2;
            Complex t3 = x3[k] * w3;
            Complex t4 = x4[k] * w4;

            // Формулы для радикса-5
            Complex a0 = t0 + t1 + t2 + t3 + t4;
            Complex a1 = t0 + t1 * Complex(0.309016994374947, -0.9510565162951535) +
                t2 * Complex(-0.809016994374947, -0.587785252292473) +
                t3 * Complex(-0.809016994374947, 0.587785252292473) +
                t4 * Complex(0.309016994374947, 0.9510565162951535);
            Complex a2 = t0 + t1 * Complex(-0.809016994374947, -0.587785252292473) +
                t2 * Complex(0.309016994374947, 0.9510565162951535) +
                t3 * Complex(0.309016994374947, -0.9510565162951535) +
                t4 * Complex(-0.809016994374947, 0.587785252292473);
            Complex a3 = t0 + t1 * Complex(-0.809016994374947, 0.587785252292473) +
                t2 * Complex(0.309016994374947, -0.9510565162951535) +
                t3 * Complex(0.309016994374947, 0.9510565162951535) +
                t4 * Complex(-0.809016994374947, -0.587785252292473);
            Complex a4 = t0 + t1 * Complex(0.309016994374947, 0.9510565162951535) +
                t2 * Complex(-0.809016994374947, 0.587785252292473) +
                t3 * Complex(-0.809016994374947, -0.587785252292473) +
                t4 * Complex(0.309016994374947, -0.9510565162951535);

            if (isInverse) {
                a1 = std::conj(a1);
                a2 = std::conj(a2);
                a3 = std::conj(a3);
                a4 = std::conj(a4);
            }

            data[k] = a0;
            data[k + m] = a1;
            data[k + 2 * m] = a2;
            data[k + 3 * m] = a3;
            data[k + 4 * m] = a4;
        }
    }
    else if (n % 3 == 0) {
        size_t m = n / 3;
        std::vector<Complex> x0(m), x1(m), x2(m);

        // Разделение на подмассивы
        for (size_t i = 0; i < m; ++i) {
            x0[i] = data[i * 3];
            x1[i] = data[i * 3 + 1];
            x2[i] = data[i * 3 + 2];
        }

        // Рекурсивное применение БПФ
        fftMixedRadix(x0, isInverse);
        fftMixedRadix(x1, isInverse);
        fftMixedRadix(x2, isInverse);

        // Комбинирование результатов
        for (size_t k = 0; k < m; ++k) {
            double angle = (isInverse ? 2.0 : -2.0) * M_PI * k / n;
            Complex w(cos(angle), sin(angle));
            Complex w2 = w * w;

            Complex t0 = x0[k];
            Complex t1 = x1[k] * w;
            Complex t2 = x2[k] * w2;

            Complex omega = Complex(-0.5, isInverse ? 0.8660254037844386 : -0.8660254037844386);
            Complex omega2 = omega * omega;

            data[k] = t0 + t1 + t2;
            data[k + m] = t0 + t1 * omega + t2 * omega2;
            data[k + 2 * m] = t0 + t1 * omega2 + t2 * omega;
        }
    }
    else if (n % 2 == 0) {
        size_t m = n / 2;
        std::vector<Complex> x0(m), x1(m);

        // Разделение на подмассивы
        for (size_t i = 0; i < m; ++i) {
            x0[i] = data[i * 2];
            x1[i] = data[i * 2 + 1];
        }

        // Рекурсивное применение БПФ
        fftMixedRadix(x0, isInverse);
        fftMixedRadix(x1, isInverse);

        // Комбинирование результатов
        for (size_t k = 0; k < m; ++k) {
            double angle = (isInverse ? 2.0 : -2.0) * M_PI * k / n;
            Complex w(cos(angle), sin(angle));

            Complex t0 = x0[k];
            Complex t1 = x1[k] * w;

            data[k] = t0 + t1;
            data[k + m] = t0 - t1;
        }
    }
}
