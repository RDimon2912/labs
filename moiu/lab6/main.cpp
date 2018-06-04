#include <bits/stdc++.h>

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

    T Sum() const {
        T result = 0;
        for (size_t i = 0; i < n_; i ++) {
            for (size_t j = 0; j < m_; j ++) {
                result += matrix_[i][j];
            }
        }
        return result;
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

    Matrix Row(size_t i) const {
        Matrix res(1, m_);
        for (size_t j = 0; j < m_; j++) {
            res(0, j) = matrix_[i][j];
        }
        return res;
    }

    Matrix Col(size_t j) const {
        Matrix res(n_, 1);
        for (size_t i = 0; i < n_; i++) {
            res(i, 0) = matrix_[i][j];
        }
        return res;
    }

    Matrix TopRows(size_t count) const {
        assert(n_ >= count);
        Matrix res(count, m_);
        for (size_t i = 0; i < count; i ++) {
            res.matrix_[i] = matrix_[i];
        }
        return res;
    }

    Matrix BottomRows(size_t count) const  {
        assert(n_ >= count);
        Matrix res(count, m_);
        for (size_t i = 0; i < count; i ++) {
            res.matrix_[i] = matrix_[n_ - count + i];
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
//    if (fabs(A.Det()) < 1e-7) {
//        return false;
//    }
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


template <typename C, typename T>
size_t FindIndex(const C& c, const T& t) {
    return std::find(std::begin(c), std::end(c), t) - std::begin(c);
};

template <typename C, typename T>
const T* FindPtr(const C& c, const T& t) {
    auto result = std::find(std::begin(c), std::end(c), t);
    return result != std::end(c) ? &*result : nullptr;
}

template <typename C, typename T>
void Erase(C& c, const T& t) {
    c.erase(std::remove(std::begin(c), std::end(c), t), std::end(c));
};

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

Matrix<double> TakeCols(const Matrix<double>& from, const std::vector<size_t>& cols) {
    Matrix<double> result(from.Rows(), cols.size());
    for (size_t i = 0; i < cols.size(); ++i) {
        for (size_t j = 0; j < from.Rows(); j ++) {
            result(j, i) = from(j, cols[i]);
        }
    }

    return result;
}

Matrix<double> SolveTask(
        const Matrix<double>& A,
        const Matrix<double>& b,
        const Matrix<double>& c,
        const Matrix<double>& D,
        Matrix<double> x,
        std::vector<size_t> J_r,
        std::vector<size_t> J_r_extend
) {
    bool skip = false;

    while (true) {
        // step 1
        std::cout << std::endl << std::endl << "Start iteration" << std::endl;
        std::cout << D.Cols() << " " << x.Rows() << std::endl;
        Matrix<double> c_beg = D * x + c;
        std::cout << "c_beg: " << c_beg << std::endl;

        Matrix<double> A_r(A.Rows(), J_r.size());
        Matrix<double> c_op(1, J_r.size());
        for (size_t i = 0; i < J_r.size(); ++i) {
            for (size_t j = 0; j < A.Rows(); j ++) {
                A_r(j, i) = A(j, J_r[i]);
            }
            c_op(i) = c_beg(J_r[i]);
        }

        Matrix<double> A_r_extend = TakeCols(A, J_r_extend);

        std::cout << "A_r: " << std::endl << A_r << std::endl;
        std::cout << "A_r_extend: " << std::endl << A_r_extend << std::endl;
        std::cout << "c_op: " << c_op << std::endl;

        Matrix<double> A_r_inversed;
        GetInv(A_r, A_r_inversed);
        Matrix<double> u = (c_op * -1) * A_r_inversed;
        std::cout << "u': " << u << std::endl;

        Matrix<double> delta = u * A + c_beg.GetT();
        std::cout << "delta: " << delta << std::endl;

        // step 2
        std::vector<size_t> J_not_star;
        for (size_t i = 0; i < A.Cols(); ++i) {
            bool ok = true;
            for (const auto value : J_r_extend) {
                if (i == value) {
                    ok = false;
                }
            }

            if (ok) {
                J_not_star.push_back(i);
            }
        }

        size_t j0 = 0;
        for (size_t i : J_not_star) {
            if (delta(j0) > delta(i)) {
                j0 = i;
            }
        }
        std::cout << "j0: " << j0 << std::endl;

        if (delta(j0) >= 0) {
            return x;
        }

        // step 3
        Matrix<double> H(J_r_extend.size() + J_r.size(), J_r_extend.size() + J_r.size());
        for (size_t i = 0; i < J_r_extend.size(); ++i) {
            for (size_t j = 0; j < J_r_extend.size(); ++j) {
                H(i, j) = D(J_r_extend[i], J_r_extend[j]);
            }
        }
        Matrix<double> A_r_extend_t = A_r_extend.GetT();

        for (size_t i = 0; i < J_r_extend.size(); i ++) {
            for (size_t j = 0; j < J_r.size(); j ++) {
                H(i, J_r_extend.size() + j) = A_r_extend_t(i, j);
            }
        }
        for (size_t i = 0; i < J_r.size(); i ++) {
            for (size_t j = 0; j < J_r_extend.size(); j ++) {
                H(i + J_r_extend.size(), j) = A_r_extend(i, j);
            }
        }

        std::cout << "H: " << std::endl << H << std::endl;

        Matrix<double> bb(J_r_extend.size() + J_r.size(), 1);
        for (size_t i = 0; i < J_r_extend.size(); ++i) {
            bb(i) = D(J_r_extend[i], j0);
        }
        for (size_t i = 0; i < J_r.size(); ++i) {
            bb(i + J_r_extend.size()) = A(i, j0);
        }
        std::cout << "bb: " << bb << std::endl;

        Matrix<double> H_inversed;
        GetInv(H, H_inversed);
        Matrix<double> vector = H_inversed * -1 * bb;
        std::cout << "vector: " << vector << std::endl;

        Matrix<double> l(A.Cols(), 1);
        l(j0) = 1;
        for (size_t i = 0; i < J_r_extend.size(); ++i) {
            l(J_r_extend[i]) = vector(i);
        }
        std::cout << "l: " << l << std::endl;

        Matrix<double> theta(J_r_extend.size() + 1, 1);
        double theta_min = INFINITY;
        size_t j_star_index = -1;
        for (size_t i = 0; i < J_r_extend.size(); ++i) {
            size_t j = J_r_extend[i];
            if (l(j) < 0) {
                theta(i) = -x(j) / l(j);
                if (theta(i) < theta_min) {
                    theta_min = theta(i);
                    j_star_index = j;
                }
            }
        }
        double delta_theta = (l.GetT() * D * l)(0, 0);
        if (delta_theta > 0) {
            theta(J_r_extend.size()) = abs(delta(j0)) / delta_theta;
            if (theta(J_r_extend.size()) < theta_min) {
                theta_min = theta(J_r_extend.size());
                j_star_index = j0;
            }
        }

        std::cout << "theta: " << theta << std::endl;
        std::cout << "theta_min: " << theta_min << std::endl;
        std::cout << "j_star_index: " << j_star_index << std::endl;
        if (theta_min == INFINITY){
            throw std::runtime_error("Not limited!");
        }

        x += l * theta_min;
        std::cout << "x: " << x << std::endl;

        if (j_star_index == j0) { // case a)
            J_r_extend.push_back(j0);
            continue;
        }

        // case b)
        if (FindPtr(J_r_extend, j_star_index) != nullptr && FindPtr(J_r, j_star_index) == nullptr) {
            Erase(J_r_extend, j_star_index);
            delta(j0) += theta_min * delta_theta;
            continue;
        }

        auto s = FindIndex(J_r, j_star_index);
        Matrix<double> e_s(1, A.Rows());
        e_s(s) = 1;

        size_t j_r = -1;
        bool ok = false;
        for (const auto j : J_r_extend) {
            if (FindPtr(J_r, j)) {
                continue;
            }
            GetInv(A_r, A_r_inversed);

            Matrix<double> result = e_s * A_r_inversed * A.Col(j);

            if (std::abs(result(0, 0)) > 1e-6) {
                j_r = j;
                ok = true;
            }
        }

        if (ok) { // case c)
            Erase(J_r, j_star_index);
            Erase(J_r_extend, j_star_index);
            J_r.push_back(j_r);
            delta(j0) += theta_min * delta_theta;
        } else { // case d)
            Erase(J_r, j_star_index);
            Erase(J_r_extend, j_star_index);
            J_r.push_back(j0);
            J_r_extend.push_back(j0);
        }
    }
}

int test_1() {
    Matrix<double> A ({
            {1, 2, 0, 1, 0, 4,  -1, -3},
            {1, 3, 0, 0, 1, -1, -1, 2},
            {1, 4, 1, 0, 0, 2,  -2, 0}
    });

    Matrix<double> b ({{4, 5, 6}}); b = b.GetT();

    Matrix<double> B ({
            {1,  1, -1, 0, 3,  4,  -2, 1},
            {2,  6, 0,  0, 1,  -5, 0,  -1},
            {-1, 2, 0,  0, -1, 1,  1,  1}
    });

    Matrix<double>  d ({{7, 3, 3}}); d = d.GetT();

    Matrix<double> D = B.GetT() * B;
    std::cout << "D: " << std::endl << D << std::endl;

    Matrix<double>  c = (d.GetT() * -1 * B).GetT();
    std::cout << "c: " << c << std::endl;

    Matrix<double>  x ({{0, 0, 6, 4, 5, 0, 0, 0}}); x = x.GetT();
    std::cout << "x: " << x << std::endl;
    std::vector<size_t> J_r {3 - 1, 4 - 1, 5 - 1};
    std::vector<size_t> J_r_extend {3 - 1, 4 - 1, 5 - 1};

    auto answer = SolveTask(A, b, c, D, x, J_r, J_r_extend);
    std::cout << "Answer: " << answer << std::endl;

    auto result = c.GetT() * answer + answer.GetT() * 0.5 * D * answer;
    std::cout << "result: " << result << std::endl;

    return 0;
}

int test_2() {
    Matrix<double> A ({
              {11, 0, 0, 1, 0, -4,  -1, 1},
              {1, 1, 0, 0, 1, -1, -1, 1},
              {1, 1, 1, 0, 1, 2,  -2, 1}
      });

    Matrix<double> b ({{8, 2, 5}}); b = b.GetT();

    Matrix<double> B ({
              {1,  -1, 0, 3,  -1, 5,  -2, 1},
              {2,  5, 0,  0, -1,  4, 0,  0},
              {-1, 3, 0,  5, 4, -1,  -2,  1}
      });

    Matrix<double>  d ({{6, 10, 9}}); d = d.GetT();

    Matrix<double> D = B.GetT() * B;
    std::cout << "D: " << std::endl << D << std::endl;

    Matrix<double>  c = (d.GetT() * -1 * B).GetT();
    std::cout << "c: " << c << std::endl;

    Matrix<double>  x ({{0.7273, 1.2727, 3.0, 0, 0, 0, 0, 0}}); x = x.GetT();
    std::cout << "x: " << x << std::endl;
    std::vector<size_t> J_r {1 - 1, 2 - 1, 3 - 1};
    std::vector<size_t> J_r_extend {1 - 1, 2 - 1, 3 - 1};

    auto answer = SolveTask(A, b, c, D, x, J_r, J_r_extend);
    std::cout << "Answer: " << answer << std::endl;

    auto result = c.GetT() * answer + answer.GetT() * 0.5 * D * answer;
    std::cout << "result: " << result << std::endl;

    return 0;
}

int test_3() {
    Matrix<double> A ({
                              {2, -3, 1, 1, 3, 0, 1, 2},
                              {-1, 3, 1, 0, 1, 4, 5, -6},
                              {1, 1, -1, 0, 1, -2, 4, 8}
                      });

    Matrix<double> b ({{8, 4, 14}}); b = b.GetT();

    Matrix<double> B ({
                              {1, 0, 0, 3, -1, 5, 0, 1},
                              {2, 5, 0, 0, 0, 4, 0, 0},
                              {-1, 9, 0, 5, 2, -1, -1, 5}
                      });


    Matrix<double> D = B.GetT() * B;
    std::cout << "D: " << std::endl << D << std::endl;

    Matrix<double>  c({{-13, -217, 0, -117, -27, -71, 18, -99}}); c = c.GetT();
    std::cout << "c: " << c << std::endl;

    Matrix<double>  x ({{0, 2, 0, 0, 4, 0, 0, 1}}); x = x.GetT();
    std::cout << "x: " << x << std::endl;
    std::vector<size_t> J_r {2 - 1, 5 - 1, 8 - 1};
    std::vector<size_t> J_r_extend {2 - 1, 5 - 1, 8 - 1};

    auto answer = SolveTask(A, b, c, D, x, J_r, J_r_extend);
    std::cout << "Answer: " << answer << std::endl;

    auto result = c.GetT() * answer + answer.GetT() * 0.5 * D * answer;
    std::cout << "result: " << result << std::endl;

    return 0;
}

int test_4() {
    Matrix<double> A ({
                              {0, 2, 1, 4, 3, 0, -5, -10},
                              {-1, 3, 1, 0, 1, 3, -5, -6},
                              {1, 1, 1, 0, 1, -2, -5, 8}
                      });

    Matrix<double> b ({{6, 4, 14}}); b = b.GetT();

    Matrix<double> D(8); D(2, 2) = 0.0;

    std::cout << "D: " << std::endl << D << std::endl;

    Matrix<double>  c({{1, 3, -1, 3, 5, 2, -2, 0}}); c = c.GetT();
    std::cout << "c: " << c << std::endl;

    Matrix<double>  x ({{0, 2, 0, 0, 4, 0, 0, 1}}); x = x.GetT();
    std::cout << "x: " << x << std::endl;
    std::vector<size_t> J_r {2 - 1, 5 - 1, 8 - 1};
    std::vector<size_t> J_r_extend {2 - 1, 5 - 1, 8 - 1};

    auto answer = SolveTask(A, b, c, D, x, J_r, J_r_extend);
    std::cout << "Answer: " << answer << std::endl;

    auto result = c.GetT() * answer + answer.GetT() * 0.5 * D * answer;
    std::cout << "result: " << result << std::endl;

    return 0;
}

int test_5() {
    Matrix<double> A ({
              {0, 0, 1, 5, 2, 0, -5, -4},
              {1, 1, -1, 0, 1, -1, -1, -1},
              {1, 1, 1, 0, 1, 2, 5, 8}
      });

    Matrix<double> b ({{15, -1, 9}}); b = b.GetT();

    Matrix<double> D(8, 8);

    std::cout << "D: " << std::endl << D << std::endl;

    Matrix<double>  c({{1, -3, 4, 3, 5, 6, -2, 0}}); c = c.GetT();
    std::cout << "c: " << c << std::endl;

    Matrix<double>  x ({{4, 0, 5, 2, 0, 0, 0, 0}}); x = x.GetT();
    std::cout << "x: " << x << std::endl;
    std::vector<size_t> J_r {1 - 1, 3 - 1, 4 - 1};
    std::vector<size_t> J_r_extend {1 - 1, 3 - 1, 4 - 1};

    auto answer = SolveTask(A, b, c, D, x, J_r, J_r_extend);
    std::cout << "Answer: " << answer << std::endl;

    auto result = c.GetT() * answer + answer.GetT() * 0.5 * D * answer;
    std::cout << "result: " << result << std::endl;

    return 0;
}

int main() {
//    return test_1();
//    return test_2();
//    return test_3();
//    return test_4();
    return test_5();
}
