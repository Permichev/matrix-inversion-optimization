#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

bool invert_matrix(double* A, double* Ainv, int N) {
    std::vector<double> aug(N * 2 * N, 0.0);

    // Инициализация матрицы aug
    #pragma clang loop vectorize(enable)
    for (int i = 0; i < N; i++) {
        #pragma clang loop vectorize(enable)
        for (int j = 0; j < N; j++) {
            aug[i*(2*N) + j] = A[i*N + j];
        }
        aug[i*(2*N) + (N+i)] = 1.0;
    }

    // Метод Гаусса–Жордана
    for (int i = 0; i < N; i++) {
        double maxEl = std::fabs(aug[i*(2*N) + i]);
        int maxRow = i;
        for (int k = i+1; k < N; k++) {
            double val = std::fabs(aug[k*(2*N) + i]);
            if (val > maxEl) {
                maxEl = val;
                maxRow = k;
            }
        }

        if (std::fabs(aug[maxRow*(2*N) + i]) < 1e-14) {
            return false; // матрица необратима
        }

        if (maxRow != i) {
            #pragma clang loop vectorize(enable)
            for (int c = 0; c < 2*N; c++) {
                std::swap(aug[i*(2*N) + c], aug[maxRow*(2*N) + c]);
            }
        }

        double diagVal = aug[i*(2*N) + i];
        #pragma clang loop vectorize(enable)
        for (int c = 0; c < 2*N; c++) {
            aug[i*(2*N) + c] /= diagVal;
        }

        // Параллелизация по строкам при обнулении столбца i
        #pragma omp parallel for
        for (int r = 0; r < N; r++) {
            if (r != i) {
                double factor = aug[r*(2*N) + i];
                #pragma clang loop vectorize(enable)
                for (int c = 0; c < 2*N; c++) {
                    aug[r*(2*N) + c] -= factor * aug[i*(2*N) + c];
                }
            }
        }
    }

    // Копируем результат в Ainv
    #pragma clang loop vectorize(enable)
    for (int i = 0; i < N; i++) {
        #pragma clang loop vectorize(enable)
        for (int j = 0; j < N; j++) {
            Ainv[i*N + j] = aug[i*(2*N) + (N+j)];
        }
    }

    return true;
}

int main() {
    int N = 2000; // Большой размер для измерения производительности
    std::vector<double> A(N*N), Ainv(N*N);

    srand(123);
    for (int i = 0; i < N*N; i++) {
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
