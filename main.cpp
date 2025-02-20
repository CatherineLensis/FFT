#include "FFT.h"
#include <iostream>

double computeRMSE(const std::vector<Complex>& original, const std::vector<Complex>& restored) {
    double error = 0.0;
    size_t N = original.size();

    for (size_t i = 0; i < N; ++i) {
        double diff_real = original[i].real() - restored[i].real();
        double diff_imag = original[i].imag() - restored[i].imag();
        error += diff_real * diff_real + diff_imag * diff_imag;
    }

    return std::sqrt(error / N);
}

int main() {
    setlocale(LC_ALL, "Russian");
    const size_t N = 9; // Пример длины, кратной 3 (можно менять на 8, 25 и т.д.)

    std::vector<Complex> data(N);

    for (size_t i = 0; i < N; ++i) {
        data[i] = Complex(RandomGenerator::Sum12(), RandomGenerator::Sum12());
    }
    std::vector<Complex> original_data = data;

    std::cout << "Исходные данные:\n";
    for (const auto& num : data) {
        std::cout << num << "\n";
    }

    FFT::transform(data, 0); // Прямое преобразование
    std::cout << "\nПосле прямого FFT:\n";
    for (const auto& num : data) {
        std::cout << num << "\n";
    }

    FFT::transform(data, 1); // Обратное преобразование
    std::cout << "\nПосле обратного FFT:\n";
    for (const auto& num : data) {
        std::cout << num << "\n";
    }

    double error = computeRMSE(original_data, data);
    std::cout << "\nСреднеквадратичная ошибка (RMSE): " << error << "\n";

    return 0;
}