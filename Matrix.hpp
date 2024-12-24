#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
template <typename T>
class Matrix {
protected:
    size_t rows_, cols_;
    std::vector<T> data_;

    // SVD Decomposition (Jacobi Method)
    std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> svd() const {
        size_t m = rows_;
        size_t n = cols_;
        size_t min_dim = std::min(m, n);
        
        // Initialize matrices with correct dimensions
        Matrix<T> U = *this;                    // m x n
        Matrix<T> V(n, n);                      // n x n
        Matrix<T> S(m, n);                // m x n
        
        // Initialize V as identity matrix
        for (size_t i = 0; i < n; ++i) {
            V.set(i, i, T(1));
        }
        
        // Jacobi iteration
        const T epsilon = 1e-10;
        bool converge = false;
        int max_iterations = 1000;
        
        while (!converge && max_iterations-- > 0) {
            converge = true;
            // Only iterate over the minimum dimension
            for (size_t p = 0; p < min_dim; ++p) {
                for (size_t q = p + 1; q < n; ++q) {
                    T alpha = 0.0, beta = 0.0, gamma = 0.0;
                    
                    // Calculate only over actual rows
                    for (size_t i = 0; i < m; ++i) {
                        alpha += U.get(i, p) * U.get(i, p);
                        beta += U.get(i, q) * U.get(i, q);
                        gamma += U.get(i, p) * U.get(i, q);
                    }
                    
                    if (std::abs(gamma) > epsilon) {
                        converge = false;
                        T zeta = (beta - alpha) / (2.0 * gamma);
                        T t = std::copysign(1.0 / (std::abs(zeta) + std::sqrt(1 + zeta * zeta)), zeta);
                        T c = 1.0 / std::sqrt(1 + t * t);
                        T s = c * t;
                        
                        // Update U
                        for (size_t i = 0; i < m; ++i) {
                            T up = U.get(i, p), uq = U.get(i, q);
                            U.set(i, p, c * up - s * uq);
                            U.set(i, q, s * up + c * uq);
                        }
                        
                        // Update V
                        for (size_t i = 0; i < n; ++i) {
                            T vp = V.get(i, p), vq = V.get(i, q);
                            V.set(i, p, c * vp - s * vq);
                            V.set(i, q, s * vp + c * vq);
                        }
                    }
                }
            }
        }
        
        // Compute and normalize singular values
        for (size_t j = 0; j < min_dim; ++j) {
            T sigma = 0.0;
            for (size_t i = 0; i < m; ++i) {
                sigma += U.get(i, j) * U.get(i, j);
            }
            sigma = std::sqrt(std::max(sigma, epsilon));
            S.set(j, j, sigma);
            
            // Normalize corresponding column in U
            if (sigma > epsilon) {
                for (size_t i = 0; i < m; ++i) {
                    U.set(i, j, U.get(i, j) / sigma);
                }
            }
        }
        
        return {U, S, V};
    }
public:
    // Constructors
    Matrix(size_t rows, size_t cols)
        : rows_(rows), cols_(cols), data_(rows * cols, T(0)) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Matrix dimensions must be positive.");
        }
    }

    Matrix(size_t rows, size_t cols, const std::vector<T>& values)
        : rows_(rows), cols_(cols), data_(values) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Matrix dimensions must be positive.");
        }
        if (values.size() != rows * cols) {
            throw std::invalid_argument("Initial values do not match matrix dimensions.");
        }
    }

    // Accessors
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }

    virtual T get(size_t row, size_t col) const {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range("Index out of range.");
        }
        return data_[row * cols_ + col];
    }

    virtual void set(size_t row, size_t col, T value) {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range("Index out of range.");
        }
        data_[row * cols_ + col] = value;
    }

    // Basic operations
    Matrix<T> transpose() const {
        Matrix<T> result(cols_, rows_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.set(j, i, get(i, j));
            }
        }
        return result;
    }

    Matrix<T> gram_matrix() const {
        return this->transpose() * (*this);
    }

    size_t rank() const {
        if (rows_ == 0 || cols_ == 0) {
            return 0;
        }

        Matrix<T> A(*this);
        size_t row = 0, col = 0;

        while (row < A.rows_ && col < A.cols_) {
            size_t max_row = row;
            for (size_t i = row + 1; i < A.rows_; ++i) {
                if (std::abs(A.get(i, col)) > std::abs(A.get(max_row, col))) {
                    max_row = i;
                }
            }

            if (A.get(max_row, col) == T(0)) {
                ++col;
            } else {
                if (max_row != row) {
                    for (size_t j = 0; j < A.cols_; ++j) {
                        std::swap(A.data_[row * A.cols_ + j], A.data_[max_row * A.cols_ + j]);
                    }
                }

                for (size_t i = row + 1; i < A.rows_; ++i) {
                    T factor = A.get(i, col) / A.get(row, col);
                    for (size_t j = col; j < A.cols_; ++j) {
                        A.data_[i * A.cols_ + j] -= factor * A.get(row, j);
                    }
                }

                ++row;
                ++col;
            }
        }

        size_t rank = 0;
        for (size_t i = 0; i < A.rows_; ++i) {
            bool non_zero_row = false;
            for (size_t j = 0; j < A.cols_; ++j) {
                if (std::abs(A.get(i, j)) > T(0)) {
                    non_zero_row = true;
                    break;
                }
            }
            if (non_zero_row) {
                ++rank;
            }
        }

        return rank;
    }

    Matrix<T> pseudo_inverse() const {
        T eps = 1e-10;
        auto [U, S, V] = svd();
        // Создаем диагональную матрицу с обратными сингулярными значениями
        size_t m = rows_;
        size_t n = cols_;
        Matrix<T> S_inv(n, m);
        
        for (size_t i = 0; i < std::min(n,m); ++i) {
            T s = S.get(i, i);
            if (s > eps) {
                S_inv.set(i, i, T(1) / s);
            }
        }
        
        // Для псевдообратной матрицы: A⁺ = V * Σ⁺ * Uᵀ
        // где Σ⁺ - диагональная матрица с обратными ненулевыми сингулярными значениями
        Matrix<T> V_T = V.transpose();
        Matrix<T> U_T = U.transpose();
        //V * S_inv * Ut;
        return V*S_inv*U_T;
    }

    // Overloaded operators
    Matrix<T> operator+(const Matrix<T>& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Matrix dimensions must match for addition.");
        }

        Matrix<T> result(rows_, cols_);
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] + other.data_[i];
        }
        return result;
    }

    Matrix<T> operator-(const Matrix<T>& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Matrix dimensions must match for subtraction.");
        }

        Matrix<T> result(rows_, cols_);
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] - other.data_[i];
        }
        return result;
    }

    Matrix<T> operator*(const Matrix<T>& other) const {
        if (cols_ != other.rows_) {
            throw std::invalid_argument("Matrix dimensions must match for multiplication.");
        }

        Matrix<T> result(rows_, other.cols_);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < other.cols_; ++j) {
                T sum = 0;
                for (size_t k = 0; k < cols_; ++k) {
                    sum += get(i, k) * other.get(k, j);
                }
                result.set(i, j, sum);
            }
        }
        return result;
    }

    Matrix<T> operator*(const T& scalar) const {
        Matrix<T> result(rows_, cols_);
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] * scalar;
        }
        return result;
    }
    Matrix<T> operator/(const T& scalar) const {
        Matrix<T> result(rows_, cols_);
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] / scalar;
        }
        return result;
    }
};


// SquareMatrix class
template <typename T>
class SquareMatrix : public Matrix<T> {
private:
    bool has_converged(const Matrix<T>& old_matrix, const Matrix<T>& new_matrix, double tolerance) const {
        size_t n = this->rows();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                T diff = std::abs(old_matrix.get(i, j) - new_matrix.get(i, j));
                if (diff > tolerance) {
                    return false;
                }
            }
        }
        return true;
    }
protected:
    void qr_decomposition(SquareMatrix<T>& Q, SquareMatrix<T>& R) const{
        size_t n = this->rows();
        size_t m = this->cols();
        if (n != m) {
            throw std::invalid_argument("QR decomposition is only valid for square matrices.");
        }

        // Инициализация матриц Q и R
        Q = Matrix<T>(n, n);
        R = Matrix<T>(n, n);

        // Ортогонализация Грама-Шмидта для столбцов
        for (size_t k = 0; k < n; ++k) {
            // Копируем k-й столбец в Q
            for (size_t i = 0; i < n; ++i) {
                Q.set(i, k, this->get(i, k));
            }

            // Вычитание проекций
            for (size_t j = 0; j < k; ++j) {
                T dot_product = T(0);
                for (size_t i = 0; i < n; ++i) {
                    dot_product += Q.get(i, j) * Q.get(i, k);
                }
                for (size_t i = 0; i < n; ++i) {
                    Q.set(i, k, Q.get(i, k) - dot_product * Q.get(i, j));
                }
            }

            // Нормализация
            T norm = T(0);
            for (size_t i = 0; i < n; ++i) {
                norm += Q.get(i, k) * Q.get(i, k);
            }
            norm = std::sqrt(norm);
            for (size_t i = 0; i < n; ++i) {
                Q.set(i, k, Q.get(i, k) / norm);
            }

            // Заполнение матрицы R
            for (size_t j = k; j < n; ++j) {
                T dot_product = T(0);
                for (size_t i = 0; i < n; ++i) {
                    dot_product += Q.get(i, k) * this->get(i, j);
                }
                R.set(k, j, dot_product);
            }
        }
    } // Perform QR decomposition

public:
    SquareMatrix(size_t size) : Matrix<T>(size, size) {}
    SquareMatrix(size_t size, const std::vector<T>& values) : Matrix<T>(size, size, values) {}
    SquareMatrix<T> minor(size_t row, size_t col) const {
        if (this->rows_ <= 1 || this->cols_ <= 1) {
            throw std::invalid_argument("Matrix must be at least 2x2 to have a minor.");
        }

        // Создаем новый минор размером (rows_-1)x(cols_-1)
        SquareMatrix<T> result(this->rows_ - 1);

        size_t mi = 0;
        for (size_t i = 0; i < this->rows_; ++i) {
            if (i == row) continue; // Пропускаем строку `row`
            size_t mj = 0;
            for (size_t j = 0; j < this->cols_; ++j) {
                if (j == col) continue; // Пропускаем столбец `col`
                result.set(mi, mj, this->get(i, j));
                ++mj;
            }
            ++mi;
        }

        return result;
    }
    T determinant() const {
        if (this->rows_ != this->cols_) {
            throw std::invalid_argument("Matrix must be square.");
        }

        // Копируем данные матрицы, чтобы избежать изменений оригинала
        SquareMatrix<T> A = *this;
        T det = 1;

        for (size_t i = 0; i < A.rows_; ++i) {
            // Находим максимальный элемент для устойчивости (метод с выбором главного элемента)
            size_t max_row = i;
            for (size_t j = i + 1; j < A.rows_; ++j) {
                if (std::abs(A.get(j, i)) > std::abs(A.get(max_row, i))) {
                    max_row = j;
                }
            }

            if (A.get(max_row, i) == T(0)) {
                return T(0); // Если элемент на диагонали ноль, определитель равен нулю
            }

            // Меняем строки, если нужно, и меняем знак определителя
            if (max_row != i) {
                for (size_t j = 0; j < A.cols_; ++j) {
                    std::swap(A.data_[i * A.cols_ + j], A.data_[max_row * A.cols_ + j]);
                }
                det *= -1; // Меняем знак детерминанта при обмене строк
            }

            // Приводим нижнюю часть матрицы к нулям
            for (size_t j = i + 1; j < A.rows_; ++j) {
                T factor = A.get(j, i) / A.get(i, i);
                for (size_t k = i; k < A.cols_; ++k) {
                    A.data_[j * A.cols_ + k] -= A.get(i, k) * factor;
                }
            }

            det *= A.get(i, i); // Умножаем на диагональный элемент
        }

        return det;
    }
    SquareMatrix<T> inverse() const {
        if (this->rows_ != this->cols_) {
            throw std::invalid_argument("Matrix must be square.");
        }

        T det = this->determinant();
        if (det == T(0)) {
            throw std::invalid_argument("Matrix is singular, cannot compute the inverse.");
        }

        SquareMatrix<T> adj_A = this->adj();
        SquareMatrix<T> result(this->rows_);

        // Делим каждое значение матрицы алгебраических дополнений на детерминант
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                result.set(i, j, adj_A.get(i, j) / det);
            }
        }

        return result;
    }
    SquareMatrix<T> power(size_t exponent) const {
        if (this->rows_ != this->cols_) {
            throw std::invalid_argument("Matrix must be square.");
        }

        SquareMatrix<T> result(this->rows_);
        for (size_t i = 0; i < this->rows_; ++i) {
            result.set(i, i, T(1)); // Инициализируем как единичную матрицу
        }

        if (exponent == 0) {
            return result; // Если степень 0, возвращаем единичную матрицу
        }

        SquareMatrix<T> base = *this;

        // Если степень отрицательная, сначала находим обратную матрицу
        if (exponent < 0) {
            base = this->inverse();
            exponent = -exponent;
        }

        while (exponent > 0) {
            if (exponent % 2 == 1) {
                result = result * base;
            }
            base = base * base; // Возводим в квадрат
            exponent /= 2;
        }

        return result;
    } // Matrix exponentiation
    bool is_symmetric() const {
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (this->get(i, j) != this->get(j, i)) {
                    return false;
                }
            }
        }
        return true;
    }
    SquareMatrix<T> adj() const {
        if (this->rows_ != this->cols_) {
            throw std::invalid_argument("Matrix must be square.");
        }

        SquareMatrix<T> result(this->rows_);

        // Вычисляем алгебраические дополнения для каждой позиции
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                // Получаем минор для позиции (i, j)
                SquareMatrix<T> minor_matrix = this->minor(i, j);
                // Алгебраическое дополнение — это минор с учетом знака
                result.set(i, j, std::pow(-1, i + j) * minor_matrix.determinant());
            }
        }
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                result.set(j, i, this->get(i, j));
            }
        }
        // Возвращаем транспонированную матрицу
        return result;
    }

    std::vector<T> solve(const Matrix<T>& B) const {
        if (this->rows_ != this->cols_) {
            throw std::invalid_argument("Matrix must be square.");
        }
        if (B.rows() != this->rows_ || B.cols() != 1) {
            throw std::invalid_argument("Vector size must match the number of rows in the matrix.");
        }

        size_t n = this->rows_;
        Matrix<T> augmented_matrix(n, n + 1);

        // Заполняем расширенную матрицу (A|B)
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                augmented_matrix.set(i, j, this->get(i, j));
            }
            augmented_matrix.set(i, n, B.get(i, 0)); // добавляем столбец B
        }

        // Прямой ход метода Гаусса
        for (size_t i = 0; i < n; ++i) {
            // Находим строку с максимальным элементом в столбце i
            size_t max_row = i;
            for (size_t j = i + 1; j < n; ++j) {
                if (std::abs(augmented_matrix.get(j, i)) > std::abs(augmented_matrix.get(max_row, i))) {
                    max_row = j;
                }
            }

            // Если элемент на главной диагонали ноль, матрица вырождена
            if (augmented_matrix.get(max_row, i) == T(0)) {
                throw std::invalid_argument("Matrix is singular, cannot compute the solution.");
            }

            // Меняем текущую строку с строкой с максимальным элементом
            for (size_t j = 0; j <= n; ++j) {
                T temp = augmented_matrix.get(i, j);
                augmented_matrix.set(i, j, augmented_matrix.get(max_row, j));
                augmented_matrix.set(max_row, j, temp);
            }

            // Приводим строки ниже текущей к нулю
            for (size_t j = i + 1; j < n; ++j) {
                T factor = augmented_matrix.get(j, i) / augmented_matrix.get(i, i);
                for (size_t k = i; k <= n; ++k) {
                    augmented_matrix.set(j, k, augmented_matrix.get(j, k) - factor * augmented_matrix.get(i, k));
                }
            }
        }

        // Обратный ход
        std::vector<T> solution(n);
        for (int i = n - 1; i >= 0; --i) {
            solution[i] = augmented_matrix.get(i, n);
            for (size_t j = i + 1; j < n; ++j) {
                solution[i] -= augmented_matrix.get(i, j) * solution[j];
            }
            solution[i] /= augmented_matrix.get(i, i);
        }

        return solution;
    }
    Matrix<T> operator/(const Matrix<T>& other) const{
        if (this->rows_ != this->cols_ || other.rows_ != other.cols_) {
            throw std::invalid_argument("Both matrices must be square.");
        }
        
        // Проверяем, что у другой матрицы есть обратная
        T det = other.determinant();
        if (det == T(0)) {
            throw std::invalid_argument("Matrix is singular, cannot divide.");
        }

        // Находим обратную матрицу для other
        Matrix<T> inv_other = other.inverse();

        // Умножаем текущую матрицу на обратную матрицу
        Matrix<T> result = *this * inv_other;

        return result;
    } // div matr/matr
    // Eigenvalues
    std::vector<T> eigenvalues() const {
        size_t n = this->rows();
        if (n != this->cols()) {
            throw std::invalid_argument("Matrix must be square to compute eigenvalues.");
        }

        // Инициализация текущей матрицы
        Matrix<T> current_matrix = *this;

        // QR-итерации
        const double tolerance = 1e-9;
        const size_t max_iterations = 1000;
        for (size_t iter = 0; iter < max_iterations; ++iter) {
            SquareMatrix<T> Q(n, n), R(n, n);
            current_matrix.qr_decomposition(Q, R);

            // Обновление матрицы: A_new = R * Q
            Matrix<T> next_matrix(n, n);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    T sum = T(0);
                    for (size_t k = 0; k < n; ++k) {
                        sum += R.get(i, k) * Q.get(k, j);
                    }
                    next_matrix.set(i, j, sum);
                }
            }

            // Проверка сходимости
            if (has_converged(current_matrix, next_matrix, tolerance)) {
                current_matrix = next_matrix;
                break;
            }

            current_matrix = next_matrix;
        }

        // Извлечение собственных значений из диагонали
        std::vector<T> eigenvalues;
        for (size_t i = 0; i < n; ++i) {
            eigenvalues.push_back(current_matrix.get(i, i));
        }

        return eigenvalues;
    }


    SquareMatrix<T> sin() const {
        SquareMatrix<T> result = *this;  // Начинаем с самой матрицы (x)
        SquareMatrix<T> term = *this;    // Первая степень: x
        T factorial = 1;
        size_t max_terms = 10;     // Число слагаемых в ряду

        for (size_t n = 1; n <= max_terms; ++n) {
            term = term * (*this) * (*this);  // Умножаем на x^2
            factorial *= (2 * n) * (2 * n + 1); // (2n)(2n+1)
            if (n % 2 == 0) {
                result = result - term / factorial; // Четные члены с минусом
            } else {
                result = result + term / factorial; // Нечетные члены с плюсом
            }
        }
        return result;
    }

    SquareMatrix<T> cos() const {
        SquareMatrix<T> result(this->rows()); // Начинаем с единичной матрицы
        for (int i = 0; i < this->rows(); i++) result.set(i,i, T(1));
        SquareMatrix<T> term(this->rows());   // Первая степень: 1
        for (int i = 0; i < this->rows(); i++) term.set(i,i, T(1));
        T factorial = 1;
        size_t max_terms = 10;     // Число слагаемых в ряду

        for (size_t n = 1; n <= max_terms; ++n) {
            term = term * (*this) * (*this);  // Умножаем на x^2
            factorial *= (2 * n - 1) * (2 * n); // (2n-1)(2n)
            if (n % 2 == 0) {
                result = result + term / factorial; // Четные члены с плюсом
            } else {
                result = result - term / factorial; // Нечетные члены с минусом
            }
        }
        return result;
    }

    
    SquareMatrix<T> tan() const {
        return this->sin() / this->cos();
    }

    SquareMatrix<T> exp() const {
        SquareMatrix<T> result(this->rows()); // Начинаем с единичной матрицы
        for (int i = 0; i < this->rows(); i++) result.set(i,i, T(1));
        SquareMatrix<T> term(this->rows());   // Первая степень: 1
        for (int i = 0; i < this->rows(); i++) term.set(i,i, T(1));
        T factorial = 1;
        size_t max_terms = 10;  // Число слагаемых в ряду

        for (size_t n = 1; n <= max_terms; ++n) {
            term = term * (*this);  // Умножаем на x^n
            factorial *= n;         // factorial(n)
            result = result + term / factorial;
        }
        return result;
    }

    SquareMatrix<T> log() const {
        //if (*this == identity(rows_)) {
        //    return identity(rows_);
        //}
         SquareMatrix<T> id(this->rows()); // Начинаем с единичной матрицы
        for (int i = 0; i < this->rows(); i++) id.set(i,i, T(1));
        SquareMatrix<T> result = *this - id;  // A - I
        SquareMatrix<T> term = result;
        T factorial = 1;
        size_t max_terms = 10;  // Число слагаемых в ряду

        for (size_t n = 1; n <= max_terms; ++n) {
            term = term * result;  // Умножаем на (A - I)^n
            factorial *= n;        // factorial(n)
            result = result + term / T(n);  // Делим на n
        }
        return result;
    }
};


// IdentityMatrix class
template <typename T>
class IdentityMatrix : public SquareMatrix<T> {
public:
    IdentityMatrix(size_t size) : SquareMatrix<T>(size) {
        for (size_t i = 0; i < size; ++i) {
            this->set(i, i, T(1));
        }
    }

    T get(size_t row, size_t col) const override {
        return (row == col) ? T(1) : T(0);
    }

    void set(size_t row, size_t col, T value) override {
        if (value != T(1) || row != col) {
            throw std::logic_error("Cannot modify IdentityMatrix.");
        }
    }
};
