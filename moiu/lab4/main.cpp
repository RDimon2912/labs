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


std::tuple<Matrix<double>, std::vector<size_t>> DualSimplexMethod(
        const Matrix<double>& A,
        const Matrix<double>& b,
        const Matrix<double>& c,
        std::vector<size_t> j_basis)
{
    while (true) {
        Matrix<double> A_basis(A.Rows(), j_basis.size());
        for (size_t i = 0; i < j_basis.size(); ++i) {
            for (size_t j = 0; j < A.Rows(); j ++){
                A_basis(j, i) = A(j, j_basis[i]);
            }
        }

        std::cout << "A_basis:" << std::endl << A_basis << std::endl;

        Matrix<double> A_basis_inversed;
        GetInv(A_basis, A_basis_inversed);
        std::cout << "A_basis_inversed:" << std::endl << A_basis_inversed << std::endl;

        Matrix<double> y(1, j_basis.size());
        {
            Matrix<double> c_basis(1, j_basis.size());
            for (size_t i = 0; i < j_basis.size(); ++i) {
                c_basis(i) = c(j_basis[i]);
            }

            y = c_basis * A_basis_inversed;
        }

        Matrix<double> theta = y * A - c.GetT();

        std::cout << "y: " << y << std::endl;
        std::cout << "theta: " << theta << std::endl;
        std::cout << "J_basis: ";
        for (const auto& e : j_basis) {
            std::cout << e << " ";
        }
        std::cout << std::endl;

        Matrix<double> cappa_basis = A_basis_inversed * b;
        std::cout << "cappa_basis: " << std::endl << cappa_basis << std::endl;

        size_t k = 0;
        for (size_t i = 0; i < cappa_basis.Rows(); ++i) {
            if (cappa_basis(k) > cappa_basis(i)) {
                k = i;
            }
        }

        if (cappa_basis(k) >= 0) {
            return std::make_tuple(cappa_basis, j_basis);
        }

        Matrix<double> mu = (A_basis_inversed.Row(k) * A).GetT();
        std::cout << "mu: " << mu << std::endl;

        {
            size_t min_index = 0;
            for (size_t i = 0; i < mu.Rows(); ++i) {
                if (mu(min_index) > mu(i)) {
                    min_index = i;
                }
            }

            if (mu(min_index) >= 0) {
                throw std::runtime_error("Task isn't joint!");
            }
        }

        size_t sigma_min_index = 0;
        double sigma_min = INFINITY;

        {
            std::vector<size_t> J_not_basis(A.Cols());
            iota(J_not_basis.begin(), J_not_basis.end(), 0);
            J_not_basis.erase(
                    remove_if(
                            J_not_basis.begin(), J_not_basis.end(),
                            [&j_basis](const auto index) {
                                return find(j_basis.begin(), j_basis.end(), index) !=
                                       j_basis.end();
                            }
                    ),
                    J_not_basis.end()
            );

            for (size_t i = 0; i < J_not_basis.size(); ++i) {
                size_t index = J_not_basis[i];
                double sigma = (mu(index) >= 0 ? INFINITY : -theta(index) /
                                                            mu(index));

                if (sigma < sigma_min) {
                    sigma_min = sigma;
                    sigma_min_index = index;
                }
            }
        }

        std::cout << "sigma_min: " << sigma_min << std::endl;
        std::cout << "sigma_min_index: " << sigma_min_index << std::endl;

        y +=  A_basis_inversed.Row(k) * sigma_min;

        theta += (mu * sigma_min).GetT();

        j_basis[k] = sigma_min_index;
    }
}

template <typename T>
std::string PrintContainer(T&& container) {
    std::string result = "{";
    for (const auto& e : container) {
        result += std::to_string(e);
        result += " ";
    }

    result.pop_back();
    result += "}";
    return result;
}

int test_1() {
    Matrix<double> A ({
                              {-2, -1, 1, -7, 0, 0, 0,  2},
                              {4,  2,  1, 0,  1, 5, -1, -5},
                              {1,  1,  0, -1, 0, 3, -1, 1}
                      });

    Matrix<double> b ({{-2, 4, 3}}); b = b.GetT();
    Matrix<double> c ({{2, 2, 1, -10, 1, 4, -2, -3}}); c = c.GetT();

    Matrix<double> y_beg ({{1, 1, 1}}); y_beg = y_beg.GetT();

    std::vector<size_t> j_basis {1, 4, 6};

    Matrix<double> result_cappa;
    std::vector<size_t> result_j_basis;
    std::tie(result_cappa, result_j_basis) = DualSimplexMethod(A, b, c, j_basis);

    std::cout << "Result cappa: " << result_cappa << std::endl;
    std::cout << "Result j_basis: " << PrintContainer(result_j_basis) << std::endl;
    return 0;
}

int test_2() {
    Matrix<double> A ({
                              {-2, -1, 1, -7, 0, 0, 0,  2},
                              {4,  2,  1, 0,  1, 5, -1, -5},
                              {1,  1,  0, -1, 0, 3, -1, 1}
                      });

    Matrix<double> b ({{-2, -4, -2}}); b = b.GetT();
    Matrix<double> c ({{5, 2, 3, -16, 1, 3, -3, -12}}); c = c.GetT();

    Matrix<double> y_beg ({{1, 2, -1}}); y_beg = y_beg.GetT();

    std::vector<size_t> j_basis {0, 1, 2};

    Matrix<double> result_cappa;
    std::vector<size_t> result_j_basis;
    std::tie(result_cappa, result_j_basis) = DualSimplexMethod(A, b, c, j_basis);

    std::cout << "Result cappa: " << result_cappa << std::endl;
    std::cout << "Result j_basis: " << PrintContainer(result_j_basis) << std::endl;
    return 0;
}

int main() {
//    test_1();
    test_2();

}
