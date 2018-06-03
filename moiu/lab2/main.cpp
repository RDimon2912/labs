#include <bits/stdc++.h>

const double eps = (double) 1e-7;

template <typename T>
struct Matrix {
    size_t n_, m_;
    std::vector<std::vector<T>> matrix_;

    Matrix() : n_(0), m_(0) {}
    Matrix(size_t n) : n_(n), m_(n), matrix_(n, std::vector<T>(n, 0)) {
        for (size_t i = 0; i < n; i ++) {
            matrix_[i][i] = 1.;
        }
    }
    Matrix(size_t n, size_t m) : n_(n), m_(m), matrix_(n, std::vector<T>(m, 0)) {};

    Matrix(const Matrix& t) : n_(t.n_), m_(t.m_), matrix_(t.matrix_) {}
    Matrix(Matrix&& t) : n_(t.n_), m_(t.m_), matrix_(std::move(t.matrix_)) {}
    Matrix(const std::vector<std::vector<T>>& matrix)
            : n_(matrix.size())
            , m_(matrix.front().size())
            , matrix_(std::move(matrix))
    {
    }

    Matrix& operator =(const Matrix& t) {
        n_ = t.n_;
        m_ = t.m_;
        matrix_ = t.matrix_;
        return *this;
    }

    Matrix&operator= (Matrix&& t) {
        n_ = t.n_;
        m_ = t.m_;
        matrix_ = std::move(t.matrix_);
        return *this;
    }


    Matrix& operator *= (const Matrix& t) {
        assert(m_ == t.n_);
        std::vector<std::vector<T>> result(n_, std::vector<T>(t.m_, 0));
        for (size_t i = 0; i < n_; i ++) {
            for (size_t j = 0; j < t.m_; j ++) {
                for (size_t k = 0; k < m_; k ++) {
                    result[i][j] += matrix_[i][k] * t.matrix_[k][j];
                }
            }
        }
        m_ = t.m_;
        matrix_ = result;
        return *this;
    }

    template <typename U>
    Matrix&operator *= (const U& t) {
        for (size_t i = 0; i < n_; i ++) {
            for (size_t j = 0; j < m_; j ++) {
                matrix_[i][j] *= t;
            }
        }
        return *this;
    }

    template<typename U>
    friend Matrix operator *(const Matrix& a, const U& b) {
        Matrix res(a);
        res *= b;
        return res;
    }

    Matrix& operator +=(const Matrix& t) {
        assert(n_ == t.n_ && m_ == t.m_);
        for (size_t i = 0; i < n_; i ++) {
            for (size_t j = 0; j < m_; j ++) {
                matrix_[i][j] += t(i, j);
            }
        }
        return *this;
    }

    friend Matrix operator +(const Matrix& a, const Matrix& b) {
        Matrix res(a); res += b; return res;
    }

    Matrix& operator -=(const Matrix& t) {
        assert(n_ == t.n_ && m_ == t.m_);
        for (size_t i = 0; i < n_; i ++) {
            for (size_t j = 0; j < m_; j ++) {
                matrix_[i][j] -= t(i, j);
            }
        }
        return *this;
    }

    friend Matrix operator -(const Matrix& a, const Matrix& b) {
        Matrix res(a); res -= b; return res;
    }

    T& operator()(size_t i, size_t j) {
        return matrix_[i][j];
    }

    const T& operator()(size_t i, size_t j) const  {
        return matrix_[i][j];
    }

    T& operator()(size_t i) {
        assert(n_ == 1 || m_ == 1);
        if (n_ == 1) {
            return matrix_[0][i];
        }
        return matrix_[i][0];
    }

    const T& operator()(size_t i) const  {
        assert(n_ == 1 || m_ == 1);
        if (n_ == 1) {
            return matrix_[0][i];
        }
        return matrix_[i][0];
    }

    friend std::istream& operator >>(std::istream& is, Matrix & a) {
        is >> a.n_ >> a.m_;
        a.matrix_.assign(a.n_, std::vector<T>(a.m_));
        for (size_t i = 0; i < a.n_; i ++) {
            for (size_t j = 0; j < a.m_; j ++) {
                is >> a.matrix_[i][j];
            }
        }
        return is;
    }

    friend std::ostream& operator << (std::ostream& os, const Matrix & a) {
        os << std::fixed << std::setprecision(3);
        for (size_t i = 0; i < a.n_; i ++) {
            for (size_t j = 0; j < a.m_; j ++) {
                os << (abs(a.matrix_[i][j]) < 1e-7 ? 0 : a.matrix_[i][j]) << " ";
            }
            os << "\n";
        }
        return os;
    }

    T Det() const {
        assert(n_ == m_);
        if (n_ == 1) {
            return matrix_[0][0];
        }
        T ans = 0.;
        for (size_t j = 0; j < m_; j ++) {
            ans += Adj(0, j) * matrix_[0][j];
        }
        return ans;
    }

    T Minor(size_t i, size_t j) const  {
        assert(i < n_ && j < m_);
        Matrix b(n_ - 1, m_ - 1);
        for (size_t i_ = 0; i_ < b.n_; i_ ++) {
            for (size_t j_ = 0; j_ < b.m_; j_ ++) {
                size_t x = i_ >= i ? i_ + 1 : i_;
                size_t y = j_ >= j ? j_ + 1 : j_;
                b(i_, j_) = matrix_[x][y];
            }
        }
        return b.Det();
    }

    T Adj(size_t i, size_t j) const {
        return Minor(i, j) * ( (i + j) & 1 ? -1. : 1. );
    }

    Matrix GetAdjMatrix() const {
        Matrix a(n_, m_);
        for (size_t i = 0; i < n_; i ++) {
            for (size_t j = 0; j < m_; j ++) {
                a(i, j) = Adj(i, j);
            }
        }
        return a;
    }

    Matrix GetT() const {
        Matrix t(m_, n_);
        for (size_t i = 0; i < n_; i ++) {
            for (size_t j = 0; j < m_; j ++) {
                t(j, i) = matrix_[i][j];
            }
        }
        return t;
    }

    Matrix GetInv() const {
        T det = Det();
        assert(abs(det) > 1e-7);
        Matrix ans = GetAdjMatrix().GetT();
        ans *= 1 / det;
        return ans;
    }

    bool IsIdentity() const {
        bool ok = true;
        for (size_t i = 0; i < n_; i ++) {
            for (size_t j = 0; j < m_; j ++) {
                if (i == j) ok &= abs(matrix_[i][j] - 1.) < 1e-7;
                else ok &= abs(matrix_[i][j]) < 1e-7;
            }
        }
        return ok;
    }

    Matrix Row(size_t i) {
        Matrix res(1, m_);
        for (size_t j = 0; j < m_; j++) {
            res(0, j) = matrix_[i][j];
        }
        return res;
    }

    Matrix Col(size_t j) {
        Matrix res(n_, 1);
        for (size_t i = 0; i < n_; i++) {
            res(i, 0) = matrix_[i][j];
        }
        return res;
    }

    void SwapRows(size_t i, size_t j) {
        std::swap(matrix_[i], matrix_[j]);
    }

    size_t Rows() const {
        return n_;
    }

    size_t Cols() const {
        return m_;
    }

    void PopRow() {
        matrix_.pop_back();
        n_ = matrix_.size();
    }

    void PopCol() {
        for (size_t i = 0; i < n_; i ++) {
            matrix_[i].pop_back();
            m_ = matrix_[i].size();
        }
    }
};

template<typename T>
bool LUP(const Matrix<T> &A, Matrix<T>&L, Matrix<T>& U,  Matrix<T> &P) {
    if (fabs(A.Det()) < 1e-7) {
        return false;
    }
    //n - размерность исходной матрицы
    const size_t n = A.Rows();

    Matrix<T> C = A, E(n);

    //загружаем в матрицу P единичную матрицу
    P = Matrix<T>(n);
    L = Matrix<T>(n);
    U = Matrix<T>(n);

    for( int i = 0; i < n; i++ ) {
        //поиск опорного элемента
        T pivotValue = 0;
        int pivot = -1;
        for( int row = i; row < n; row++ ) {
            if( fabs(C( row , i )) > pivotValue ) {
                pivotValue = fabs(C(row , i ));
                pivot = row;
            }
        }
        if( fabs(pivotValue) > 1e-7 ) {
            //меняем местами i-ю строку и строку с опорным элементом
            P.SwapRows(pivot, i);
            C.SwapRows(pivot, i);
            for( int j = i+1; j < n; j++ ) {
                C( j , i ) /= C( i , i );
                for( int k = i+1; k < n; k++ )
                    C( j , k ) -= C(j , i ) * C( i , k );
            }
        }
    }

    for (size_t i = 0; i < n; i ++) {
        for (size_t j = 0; j < i; j ++) {
            L(i, j) = C(i, j);
        }
    }
    U = C - L + E;
    return true;
}

template <typename T>
bool GetInv(const Matrix<T>& A, Matrix<T>& inv) {
    Matrix<T> L, U, P;
    if (!LUP(A, L, U, P)) {
        return false;
    }

    Matrix<T> a = A;
    inv = Matrix<T>(A.Rows());
    int n = (int) a.Rows();
    for (int i = n - 1; i >= 0; i --) {

        for (int j = i; j >= 0; j --) {
            T res = L(j, i);
            for (int k = j + 1; k < n; k ++) {
                res -= U(j, k) * inv(k, i);
            }
            inv(j, i) = res / U(j, j);
        }

        for (int j = i - 1; j >= 0; j --) {
            T res = U(i, j);
            for (int k = j + 1; k < n; k ++) {
                res -= L(k, j) * inv(i, k);
            }
            inv(i, j) = res / L(j, j);
        }
    }

    inv *= P;

    return true;
}



class SimplexSolver {
public:
    SimplexSolver() = default;

    void MainStage(
            const Matrix<double>& A_matrix,
            const Matrix<double>& c_row_vector,
            const Matrix<double>& x_vector,
            const std::vector<size_t>& J_basis_vector)
    {
        n = A_matrix.Rows();
        m = A_matrix.Cols();
        A = A_matrix;
        c_row = c_row_vector;
        x = x_vector;
        J_basis = J_basis_vector;

        Stage1();
        Stage2();

        int iteration = 10;
        while (iteration--) {
            std::cout << std::endl << "Start iteration!!!" << std::endl;

            Stage3();
            Stage4();

            if (CheckPlan()) {
                std::cout << "Answer found!" << std::endl;
                std::cout << "x: " << x << std::endl;
                std::cout << "J_basis: ";
                std::copy(J_basis.begin(), J_basis.end(), std::ostream_iterator<size_t>(std::cout, " "));
                std::cout << std::endl;
                return;
            }

            Stage5();
            Stage6();
            Stage7();

            if (!IsLimited()) {
                std::cout << "Not limited!" << std::endl;
                return;
            }

            Stage8();
        }
    }

private:
    // input
    size_t n, m;
    Matrix<double> A;
    Matrix<double> b;
    Matrix<double> c_row;
    Matrix<double> x;
    std::vector<size_t> J_basis;

    // stages
    Matrix<double> A_basis;
    Matrix<double> c_basis;


    void Stage1() {
        A_basis = Matrix<double>(n, J_basis.size());

        for (size_t i = 0; i < J_basis.size(); ++i) {
            for (size_t j = 0; j < n; j ++) {
                A_basis(j, i) = A(j, J_basis.at(i));
            }
        }

        std::cout << "A_basis:" << std::endl
                  << A_basis << std::endl;

        c_basis = Matrix<double>(1, J_basis.size());
        for (size_t i = 0; i < J_basis.size(); ++i) {
            c_basis(0, i) = c_row(0, J_basis.at(i));
        }
    }

    Matrix<double> A_basis_inversed;
    void Stage2() {
        const auto det = A_basis.Det();
        if (std::abs(det) < eps) {
            throw std::runtime_error("Determinant is zero!");
        }
        GetInv(A_basis, A_basis_inversed);
    }

    Matrix<double> u_row;
    void Stage3() {
        u_row = c_basis * A_basis_inversed;

        std::cout << "u' = " << u_row << std::endl;
    }

    Matrix<double> delta_row;
    void Stage4() {
        delta_row = u_row * A - c_row;
        std::cout << "delta' = " << delta_row << std::endl;
    }

    bool CheckPlan() {
        for (size_t i = 0; i < delta_row.Cols(); ++i) {
            if (IsBasisIndex(i)) {
                continue;
            }

            if (delta_row(0, i) < 0) {
                return false;
            }
        }

        return true;
    }

    size_t j0;
    void Stage5() {
        for (size_t i = 0; i < m; ++i) {
            if (!IsBasisIndex(i) && delta_row(0, i) < 0) {
                j0 = i;

                std::cout << "j0 = " << j0 << std::endl; // zero based indices
                return;
            }
        }

        throw std::runtime_error("Could not found non-basis index");
    }

    Matrix<double> z;
    void Stage6() {
        z = A_basis_inversed * A.Col(j0);

        std::cout << "z:" << std::endl
                  << z << std::endl;
    }

    Matrix<double> theta;
    void Stage7() {
        theta = Matrix<double>(1, J_basis.size());

        for (size_t i = 0; i < J_basis.size(); ++i) {
            if (z(i, 0) <= eps) {
                theta(0, i) = INFINITY;
            } else {
                theta(0, i) = x(0, J_basis.at(i)) / z(i, 0);
            }
        }

        std::cout << "theta: " << theta << std::endl;
    }

    size_t theta_min_index;
    bool IsLimited() {
        theta_min_index = 0;
        for (size_t i = 0; i < theta.Cols(); ++i) {
            if (theta(0, theta_min_index) > theta(0, i)) {
                theta_min_index = i;
            }
        }

        std::cout << "theta_min_index = " << theta_min_index << std::endl;

        return theta(0, theta_min_index) != INFINITY;
    }

    void Stage8() {
        size_t before_replace_index = J_basis.at(theta_min_index);

        Matrix<double> x_delta = z * theta(0, theta_min_index);
        for (size_t i = 0; i < J_basis.size(); ++i) {
            x(J_basis.at(i)) -= x_delta(i);
        }

        x(j0) = theta(theta_min_index);
        std::cout << "x = " << x << std::endl;

        J_basis.at(theta_min_index) = j0;
        std::cout << "J_basis = ";
        for (size_t i = 0; i < J_basis.size(); ++i) {
            std::cout << J_basis.at(i) << " ";
        }
        std::cout << std::endl;

        for (size_t i = 0; i < c_basis.Cols(); ++i) {
            c_basis(i) = c_row(J_basis.at(i));
        }

        for (size_t i = 0; i < A_basis.Rows(); i ++) {
            A_basis(i, theta_min_index) = A(i, J_basis.at(theta_min_index));
        }
        GetInv(A_basis, A_basis_inversed);
    }


    // helpers
    bool IsBasisIndex(size_t i) {
        for (size_t j = 0; j < J_basis.size(); ++j) {
            if (J_basis.at(j) == i) {
                return true;
            }
        }
        return false;
    }
};

int test_1() {
    SimplexSolver simplexSolver;
    Matrix<double> A({
            {0, 1,  4,  1, 0,  -3, 5,  0},
            {1, -1, 0,  1, 0,  0,  1,  0},
            {0, 7,  -1, 0, -1, 3,  8,  0},
            {1, 1,  1,  1, 0,  3,  -3, 1}
    });

    Matrix<double> b ({{
            6, 10, -2, 15
    }});

    Matrix<double> c ({{
            -5, -2, 3, -4, -6, 0, -1, -5
    }});

    Matrix<double> x ({{
            4, 0, 0, 6, 2, 0, 0, 5
    }});

    std::vector<size_t> J_basis {
            0, 3, 4, 7
    };

    simplexSolver.MainStage(A, c, x, J_basis);

    return 0;
}

int test_2() {
    SimplexSolver simplexSolver;

    Matrix<double> A({
             {0, 1,  4,  1, 0,  -3, 1,  0},
             {1, -1, 0,  1, 0,  0,  0,  0},
             {0, 7,  -1, 0, -1, 3,  -1,  0},
             {1, 1,  1,  1, 0,  3,  -1, 1}
     });

    Matrix<double> b ({{
           6, 10, -2, 15
   }});

    Matrix<double> c ({{
           -5, -2, 3, -4, -6, 0, 1, -5
   }});

    Matrix<double> x ({{
           10, 0, 1.5, 0, 0.5, 0, 0, 3.5
   }});

    std::vector<size_t> J_basis {
            0, 2, 4, 7
    };

    simplexSolver.MainStage(A, c, x, J_basis);

    return 0;

}

int test_3() {
    SimplexSolver simplexSolver;
    Matrix<double> A({
                             {0, 1,  4,  1, 0,  -8, 1,  5},
                             {0, -1, 0,  -1, 0,  0,  0,  0},
                             {0, 2,  -1, 0, -1, 3,  -1,  0},
                             {1, 1,  1,  1, 0,  3,  1, 1}
                     });

    Matrix<double> b ({{
                               36, -11, 10, 20
                       }});

    Matrix<double> c ({{
                               -5, -2, 3, -4, -6, 0, 1, -5
                       }});

    Matrix<double> x ({{
                               4, 5, 0, 6, 0, 0, 0, 5
                       }});

    std::vector<size_t> J_basis {
            0, 1, 3, 7
    };

    simplexSolver.MainStage(A, c, x, J_basis);

    return 0;
}

int test_4() {
    SimplexSolver simplexSolver;
    Matrix<double> A({
                             {0, 1,  1, 1,    0,  -8, 1,    5},
                             {0, -1, 0, -7.5, 0,  0,  0,    2},
                             {0, 2,  1, 0,    -1, 3,  -1.4, 0},
                             {1, 1,  1, 1,    0,  3,  1,    1}
                     });

    Matrix<double> b ({{
                               15, -45, 1.8, 19
                       }});

    Matrix<double> c ({{
                               -6, -9, -5, 2, -6, 0, 1, 3
                       }});

    Matrix<double> x ({{
                               4, 0, 6, 6, 0, 0, 3, 0
                       }});

    std::vector<size_t> J_basis {
            0, 2, 3, 6
    };

    simplexSolver.MainStage(A, c, x, J_basis);

    return 0;
}

int test_5() {
    SimplexSolver simplexSolver;
    Matrix<double> A({
                             {0, 1,  1,  -7.5, 0,  0, 0,  2},
                             {0, 2, 1, 0, -1, 3, -1.5, 0},
                             {1, -1, 1, -1, 0, 3, 1, 1},
                     });

    Matrix<double> b ({{
                               6, 1.5, 10
                       }});

    Matrix<double> c ({{
                               -6, -9, -5, 2, -6, 0, 1, 3
                       }});

    Matrix<double> x ({{
                               4, 0, 6, 0, 4.5, 0, 0, 0
                       }});

    std::vector<size_t> J_basis {
            0, 2, 4
    };

    simplexSolver.MainStage(A, c, x, J_basis);

    return 0;
}


int main() {
//    test_1();
//    test_2();
//    test_3();
//    test_4();
    test_5();
}