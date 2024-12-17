#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdlib>

// Функция для обращения матрицы методом Гаусса–Жордана.
// A - входная матрица размера N×N
// Ainv - выходная матрица (буфер памяти такой же размерности), куда будет записана обратная
// N - размерность
bool invert_matrix(double* A, double* Ainv, int N) {
    // Создаем увеличенную матрицу размером N x 2N
    // В левой части будет A, в правой единичная матрица, после преобразований правая часть станет Ainv
    std::vector<double> aug(N * 2 * N, 0.0);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            aug[i * (2 * N) + j] = A[i * N + j];
        }
        aug[i * (2 * N) + (N + i)] = 1.0;
    }

    // Прямой ход метода Гаусса–Жордана
    for (int i = 0; i < N; i++) {
        // Находим строку с максимальным по модулю ведущим элементом для численной устойчивости
        double maxEl = std::fabs(aug[i * (2 * N) + i]);
        int maxRow = i;
        for (int k = i + 1; k < N; k++) {
            double val = std::fabs(aug[k * (2 * N) + i]);
            if (val > maxEl) {
                maxEl = val;
                maxRow = k;
            }
        }

        // Если ведущий элемент ~0, матрица необратима
        if (std::fabs(aug[maxRow * (2 * N) + i]) < 1e-14) {
            return false; // необратима
        }

        // Меняем местами i-ую и maxRow-строки
        if (maxRow != i) {
            for (int c = 0; c < 2 * N; c++) {
                std::swap(aug[i * (2 * N) + c], aug[maxRow * (2 * N) + c]);
            }
        }

        // Нормируем ведущий элемент до 1, делим всю строку
        double diagVal = aug[i * (2 * N) + i];
        for (int c = 0; c < 2 * N; c++) {
            aug[i * (2 * N) + c] /= diagVal;
        }

        // Обнуляем столбец i во всех строках, кроме i-ой
        for (int r = 0; r < N; r++) {
            if (r != i) {
                double factor = aug[r * (2 * N) + i];
                for (int c = 0; c < 2 * N; c++) {
                    aug[r * (2 * N) + c] -= factor * aug[i * (2 * N) + c];
                }
            }
        }
    }

    // Теперь правая половина матрицы aug должна быть Ainv
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Ainv[i * N + j] = aug[i * (2 * N) + (N + j)];
        }
    }

    return true;
}

int main() {
    int N = 1000; // Размер матрицы
    std::vector<double> A(N*N), Ainv(N*N);

    // Инициализируем матрицу случайными значениями
    srand(123);
    for (int i = 0; i < N * N; i++) {
        A[i] = (double)rand() / RAND_MAX;
    }

    auto start = std::chrono::high_resolution_clock::now();
    bool ok = invert_matrix(A.data(), Ainv.data(), N);
    auto end = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(end - start).count();

    if (!ok) {
        std::cout << "Матрица необратима.\n";
    } else {
        std::cout << "Обращение матрицы завершено.\n";
    }

    std::cout << "Размер: " << N << "x" << N << ", время: " << elapsed << " с\n";

    return 0;
}
