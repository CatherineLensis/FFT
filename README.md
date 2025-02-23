# Быстрое преобразование Фурье (FFT) на C++

Этот проект реализует класс быстрого прямого и обратного преобразования Фурье (Fast Fourier Transform, FFT) для комплексных чисел с длиной входных данных, кратной 2, 3 или 5. Реализация включает поддержку радикс-2 и смешанного радикса (radix-3), а также демонстрирует точность алгоритма путем сравнения исходных данных с восстановленными после прямого и обратного преобразования.

## Описание

Проект состоит из трех основных файлов:
- **`FFT.h`**: Заголовочный файл с определением класса `FFT` и вспомогательного класса `RandomGenerator`.
- **`FFT.cpp`**: Реализация методов класса `FFT`, включая алгоритмы `fftRadix2` (для длин, кратных 2) и `fftMixedRadix` (для длин, кратных 3 или 5).
- **`main.cpp`**: Тестовая программа, которая генерирует случайные комплексные числа, выполняет прямое и обратное преобразование и вычисляет среднеквадратичную ошибку (RMSE).

### Основные возможности
- Поддержка длин входных данных, кратных 2, 3 или 5 (например, 8, 9, 25).
- Реализация алгоритма Cooley-Tukey для радикс-2 и смешанного радикса для радикс-3.
- Генерация случайных комплексных чисел с помощью метода суммы 12 равномерных распределений.
- Вычисление RMSE для оценки точности восстановления данных.

### Требования
- Компилятор C++ с поддержкой C++11 или выше (например, g++, clang++).
- Стандартная библиотека C++ (`<complex>`, `<vector>`, `<random>` и т.д.).

## Установка и запуск

1. Склонируйте репозиторий:
   ```bash
   git clone <URL_репозитория>
   cd <имя_папки>
   ```
2. Скомпилируйте проект:
   ```bash
   g++ -o fft main.cpp FFT.cpp -std=c++11
   ```
3. Запустите программу
   ```bash
   ./fft
   ```
Пример вывода
   ```text
   Исходные данные:
   (0.0452289,-0.129608)
   (0.0421805,0.213138)
   (-0.0391098,-0.108629)
   (-0.0650625,-0.0380168)
   (-0.0110111,0.0431392)
   (0.248634,0.445509)
   (-0.113791,-0.0403525)
   (-0.0912338,-0.0731161)
   (0.0967341,0.0743496)
   
   После прямого FFT:
   (0.11257,0.386413)
   (-0.086214,-0.330596)
   (0.84444,0.41581)
   (-0.454235,-0.187928)
   (0.414108,-0.660683)
   (-0.523043,-0.169207)
   (-0.0592082,-0.822418)
   (0.0821415,0.593408)
   (0.0765014,-0.391272)
   
   После обратного FFT:
   (0.0452289,-0.129608)
   (0.0421805,0.213138)
   (-0.0391098,-0.108629)
   (-0.0650625,-0.0380168)
   (-0.0110111,0.0431392)
   (0.248634,0.445509)
   (-0.113791,-0.0403525)
   (-0.0912338,-0.0731161)
   (0.0967341,0.0743496)
   
   Среднеквадратичная ошибка (RMSE): 8.66976e-17
   ```
Ошибка RMSE на уровне (10^{-17}) подтверждает высокую точность алгоритма, близкую к машинной точности double.

## Структура кода
-**`FFT::transform`**: Основной метод, который выбирает подходящий алгоритм (радикс-2 или смешанный радикс) и выполняет масштабирование для обратного преобразования.
-**`FFT::fftRadix2`**: Реализация БПФ для длин, являющихся степенью 2, с использованием алгоритма Cooley-Tukey.
-**`FFT::fftMixedRadix`**: Реализация для длин, кратных 3 или 5, с рекурсивным разбиением и комбинированием результатов.
-**`RandomGenerator::Sum12`**: Генерация псевдослучайных чисел методом суммы 12 равномерных распределений.
-**`computeRMSE`**: Функция для вычисления среднеквадратичной ошибки между двумя векторами комплексных чисел.
-**`main`**: Тестовая функция, выполняющая генерацию данных, прямое и обратное преобразование, а также вывод результатов.
## Ограничения
Поддержка радикс-5 в `fftMixedRadix` пока не полностью реализована (требуется доработка комбинации результатов).
Длина входных данных должна быть строго кратна 2, 3 или 5, иначе программа выведет ошибку.
Возможные улучшения
Завершить реализацию радикс-5 в `fftMixedRadix`.
Добавить поддержку произвольных комбинаций радиксов (например, (2^a \cdot 3^b \cdot 5^c)) через динамическое разбиение.
Оптимизировать производительность для больших массивов (например, кэш-дружественные алгоритмы).

Автор: Ekaterina Lazko

Дата: Февраль 2025

---


