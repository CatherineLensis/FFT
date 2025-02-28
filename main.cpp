#include "FFT.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

// Класс RandomGenerator (пример реализации Sum12)
class RandomGenerator {
public:
    static double Sum12() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dis(-1.0, 1.0);
        double sum = 0.0;
        for (int i = 0; i < 12; i++) {
            sum += dis(gen);
        }
        return sum / 12.0; // Пример нормализации
    }
};

// Функция для вычисления RMSE
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
    const size_t N = 2700; // 2^2 * 3^3 * 5^2

    std::vector<Complex> data(N);

    // Заполнение случайными данными
    for (size_t i = 0; i < N; ++i) {
        data[i] = Complex(RandomGenerator::Sum12(), RandomGenerator::Sum12());
    }
    std::vector<Complex> original_data = data;

    // Вывод исходных данных (ограничиваем до первых 5 для читаемости)
    std::cout << "Исходные данные (первые 5):\n";
    for (size_t i = 0; i < 5 && i < N; ++i) {
        std::cout << data[i] << "\n";
    }

    // Прямое преобразование
    FFT::transform(data, false);
    std::cout << "\nПосле прямого FFT (первые 5):\n";
    for (size_t i = 0; i < 5 && i < N; ++i) {
        std::cout << data[i] << "\n";
    }

    // Обратное преобразование
    FFT::transform(data, true);
    std::cout << "\nПосле обратного FFT (первые 5):\n";
    for (size_t i = 0; i < 5 && i < N; ++i) {
        std::cout << data[i] << "\n";
    }

    // Вычисление и вывод ошибки
    double error = computeRMSE(original_data, data);
    std::cout << "\nСреднеквадратичная ошибка (RMSE): " << error << "\n";

    return 0;
}
