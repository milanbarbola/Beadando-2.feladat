#include <iostream>
#include <vector>

template<typename T>
class Vector {
private:
    std::vector<T> data;

public:
    Vector() = default;
    explicit Vector(int size) : data(size) {}
    Vector(std::initializer_list<T> values) : data(values) {}

    // Accessor for the vector elements
    T& operator[](int index) {
        return data[index];
    }

    // Const accessor for the vector elements
    const T& operator[](int index) const {
        return data[index];
    }

    // Vector addition
    Vector<T> operator+(const Vector<T>& other) const {
        Vector<T> result(data.size());
        for (int i = 0; i < data.size(); ++i) {
            result[i] = data[i] + other[i];
        }
        return result;
    }

    // Vector subtraction
    Vector<T> operator-(const Vector<T>& other) const {
        Vector<T> result(data.size());
        for (int i = 0; i < data.size(); ++i) {
            result[i] = data[i] - other[i];
        }
        return result;
    }

    // Scalar multiplication
    Vector<T> operator*(T scalar) const {
        Vector<T> result(data.size());
        for (int i = 0; i < data.size(); ++i) {
            result[i] = data[i] * scalar;
        }
        return result;
    }

    // Scalar division
    Vector<T> operator/(T scalar) const {
        Vector<T> result(data.size());
        for (int i = 0; i < data.size(); ++i) {
            result[i] = data[i] / scalar;
        }
        return result;
    }

    // Print the vector elements
    void print() const {
        for (const auto& value : data) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
};

template<typename T>
class Matrix {
private:
    std::vector<T> data;
    int size;

public:
    Matrix(int size) : size(size), data(size * size) {}

    // Accessor for the matrix elements
    T& operator()(int row, int col) {
        return data[row * size + col];
    }

    // Const accessor for the matrix elements
    const T& operator()(int row, int col) const {
        return data[row * size + col];
    }

    // Matrix addition
    Matrix<T> operator+(const Matrix<T>& other) const {
        Matrix<T> result(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }

    // Matrix subtraction
    Matrix<T> operator-(const Matrix<T>& other) const {
        Matrix<T> result(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return result;
    }

    // Scalar multiplication
    Matrix<T> operator*(T scalar) const {
        Matrix<T> result(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result(i, j) = (*this)(i, j) * scalar;
            }
        }
        return result;
    }

    // Scalar division
    Matrix<T> operator/(T scalar) const {
        Matrix<T> result(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result(i, j) = (*this)(i, j) / scalar;
            }
        }
        return result;
    }

    // Matrix multiplication
    Matrix<T> operator*(const Matrix<T>& other) const {
        Matrix<T> result(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                T sum = 0;
                for (int k = 0; k < size; ++k) {
                    sum += (*this)(i, k) * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    // Compute the determinant of the matrix
    T determinant() const {
        if (size == 1) {
            return (*this)(0, 0);
        } else if (size == 2) {
            return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
        } else {
            T det = 0;
            for (int i = 0; i < size; ++i) {
                Matrix<T> submatrix(size - 1);
                for (int j = 1; j < size; ++j) {
                    for (int k = 0; k < size; ++k) {
                        if (k < i) {
                            submatrix(j - 1, k) = (*this)(j, k);
                        } else if (k > i) {
                            submatrix(j - 1, k - 1) = (*this)(j, k);
                        }
                    }
                }
                det += (*this)(0, i) * submatrix.determinant() * (i % 2 == 0 ? 1 : -1);
            }
            return det;
        }
    }

    // Compute the inverse of the matrix
    Matrix<T> inverse() const {
        Matrix<T> result(size);
        T det = determinant();
        if (det == 0) {
            std::cout << "Matrix is not invertible." << std::endl;
            return result;
        }

        Matrix<T> adjoint(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                Matrix<T> submatrix(size - 1);
                for (int m = 0; m < size; ++m) {
                    for (int n = 0; n < size; ++n) {
                        if (m < i && n < j) {
                            submatrix(m, n) = (*this)(m, n);
                        } else if (m < i && n > j) {
                            submatrix(m, n - 1) = (*this)(m, n);
                        } else if (m > i && n < j) {
                            submatrix(m - 1, n) = (*this)(m, n);
                        } else if (m > i && n > j) {
                            submatrix(m - 1, n - 1) = (*this)(m, n);
                        }
                    }
                }
                adjoint(j, i) = submatrix.determinant() * ((i + j) % 2 == 0 ? 1 : -1);
            }
        }

        result = adjoint / det;
        return result;
    }

    // Print the matrix elements
    void print() const {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    // Matrix-vector multiplication
    Vector<T> operator*(const Vector<T>& vec) const {
        Vector<T> result(size);
        for (int i = 0; i < size; ++i) {
            T sum = 0;
            for (int j = 0; j < size; ++j) {
                sum += (*this)(i, j) * vec[j];
            }
            result[i] = sum;
        }
        return result;
    }
};


int main() {
    Vector<double> v1({1.0, 2.0, 3.0});
    Vector<double> v2({4.0, 5.0, 6.0});

    Vector<double> v3 = v1 + v2;
    Vector<double> v4 = v1 - v2;
    Vector<double> v5 = v1 * 2.0;
    Vector<double> v6 = v2 / 2.0;

    Matrix<double> m1(4);
    m1(0, 0) = 1; m1(0, 1) = 2; m1(0, 2) = 3; m1(0, 3) = 4;
    m1(1, 0) = 4; m1(1, 1) = 5; m1(1, 2) = 6; m1(1, 3) = 6;
    m1(2, 0) = 7; m1(2, 1) = 8; m1(2, 2) = 10; m1(2, 3) = 6;
    m1(3, 0) = 7; m1(3, 1) = 8; m1(3, 2) = 9; m1(3, 3) = 6;

    Matrix<double> m2(4);
    m2(0, 0) = 9; m2(0, 1) = 8; m2(0, 2) = 7; m2(0, 3) = 7;
    m2(1, 0) = 6; m2(1, 1) = 5; m2(1, 2) = 4; m2(1, 3) = 8;
    m2(2, 0) = 3; m2(2, 1) = 2; m2(2, 2) = 1; m2(2, 3) = 7;
    m2(3, 0) = 5; m2(3, 1) = 2; m2(3, 2) = 10; m2(3, 3) = 7;

    Matrix<double> m3(3);
    m3(0, 0) = 1; m3(0, 1) = 2; m3(0, 2) = 3;
    m3(1, 0) = 4; m3(1, 1) = 5; m3(1, 2) = 6;
    m3(2, 0) = 7; m3(2, 1) = 8; m3(2, 2) = 9;

    

    Matrix<double> m4 = m1 + m2;
    Matrix<double> m5 = m1 - m2;
    Matrix<double> m6 = m1 * 2;
    Matrix<double> m7 = m2 / 2;
    Matrix<double> m8 = m1 * m2;
    Matrix<double> m9 = m1.inverse();
    Vector<double> m10 = m3 * v1;

    std::cout << "Determinant of m1: " << m1.determinant() << std::endl;
    std::cout <<"" << std::endl;

    std::cout << "Determinant of m2: " << m2.determinant() << std::endl;
    std::cout <<"" << std::endl;

    std::cout << "Inverse of m1" << std::endl;
    m9.print();
    std::cout<< "" << std::endl;

    std::cout << "Matrix vector multiplication" << std:: endl;
    m10.print();
    std::cout<< "" << std::endl;
    
    std::cout<< "Testing other functinalities of the classes" << std::endl;
    m4.print(); // Here you can test the other functionalities like the sum (m4), subtraction (m5), multiplication with scalar (m6), etc.

    return 0;
}
