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

    Matrix TopRows(size_t count) {
        assert(n_ >= count);
        Matrix res(count, m_);
        for (size_t i = 0; i < count; i ++) {
            res.matrix_[i] = matrix_[i];
        }
        return res;
    }

    Matrix BottomRows(size_t count) {
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


using PlanIndex = std::pair<size_t, size_t>;

void DoSolveNorthWestTask(
        const Matrix<double>& workCost,
        const std::vector<PlanIndex>& basis_indices,
        Matrix<double>& dstWorkStores, Matrix<double>& dstWorkConsumers)
{
    Matrix<double> newMatrix(basis_indices.size() + 1, basis_indices.size() + 1);
    Matrix<double> b(newMatrix.Rows(), 1);

    for (size_t i = 0; i < basis_indices.size(); ++i) {
        const auto& index = basis_indices[i];
        newMatrix(i, index.first) = 1;
        newMatrix(i, workCost.Rows() + index.second) = 1;
        b(i) = workCost(index.first, index.second);
    }
    newMatrix(newMatrix.Rows() - 1, 0) = 1;
    std::cout << "new_matrix: " << std::endl << newMatrix << std::endl;
    std::cout << "b: " << std::endl << b << std::endl;

    Matrix<double> newMatrix_inverted;
    GetInv(newMatrix, newMatrix_inverted);

    Matrix<double> answer = newMatrix_inverted * b;
    dstWorkStores = answer.TopRows(workCost.Rows());
    dstWorkConsumers = answer.BottomRows(workCost.Cols());

    for (size_t i = 0; i < basis_indices.size(); ++i) {
        const auto& index = basis_indices[i];

        assert(
                abs(workCost(index.first, index.second) -
                    dstWorkStores(index.first) -
                    dstWorkConsumers(index.second)) < 1e-9

        );
    }
}

std::vector<PlanIndex> FindSmartCycle(
        const PlanIndex& start,
        const std::vector<PlanIndex>& basis_indices)
{
    auto allGraph = basis_indices;
    allGraph.push_back(start);

    std::vector<PlanIndex> path;

    std::function<bool(const PlanIndex& index, const PlanIndex& prev)> dfs;
    dfs = [&](const PlanIndex& current, const PlanIndex& prev) {
        if (!path.empty() && path.front() == current) {
            return true;
        }
        path.push_back(current);

        std::vector<PlanIndex> up_down;
        std::vector<PlanIndex> left_right;

        copy_if(allGraph.begin(), allGraph.end(), back_inserter(up_down), [&](const PlanIndex& c) {
            return c.second == current.second;
        });
        sort(up_down.begin(), up_down.end());

        copy_if(allGraph.begin(), allGraph.end(), back_inserter(left_right), [&](const PlanIndex& c) {
            return c.first == current.first;
        });
        sort(left_right.begin(), left_right.end());

        auto up_down_iter = lower_bound(up_down.begin(), up_down.end(), current);
        assert(up_down_iter != up_down.end());
        auto left_right_iter = lower_bound(left_right.begin(), left_right.end(), current);
        assert(left_right_iter != left_right.end());

        std::vector<PlanIndex> nextNodes;
        if (up_down_iter != up_down.begin()) {
            nextNodes.push_back(*std::prev(up_down_iter));
        }
        if (std::next(up_down_iter) != up_down.end()) {
            nextNodes.push_back(*std::next(up_down_iter));
        }
        if (left_right_iter != left_right.begin()) {
            nextNodes.push_back(*std::prev(left_right_iter));
        }
        if (std::next(left_right_iter) != left_right.end()) {
            nextNodes.push_back(*std::next(left_right_iter));
        }

        for (const auto& next : nextNodes) {
            if (next == prev) {
                continue;
            }

            if (dfs(next, current)) {
                return true;
            }
        }

        path.pop_back();
        return false;
    };

    assert(dfs(start, {-1, -1}));

    // remove unnecessary nodes
    std::vector<PlanIndex> result;
    result.push_back(path.front());

    auto check = [](const auto& a, const auto& b, const auto& c) {
        if (a.first == b.first && b.first == c.first) {
            return false;
        }
        if (a.second == b.second && b.second == c.second) {
            return false;
        }

        return true;
    };
    for (size_t i = 1; i + 1 < path.size(); ++i) {
        if (!check(path[i - 1], path[i], path[i + 1])) {
            continue;
        }

        result.push_back(path[i]);
    }
    if (check(result.front(), path.back(), result.back())) {
        result.push_back(path.back());
    }

    return result;
}

Matrix<double> DoTransportTask(Matrix<double> stores, Matrix<double> consumers, Matrix<double> cost) {
    auto storesSum = stores.Sum();
    auto consumersSum = consumers.Sum();

    if (storesSum < consumersSum) {
        Matrix<double> newStores(stores.Rows() + 1, 1);
        for (size_t i =0; i < stores.Rows(); i ++) {
            newStores(i) = stores(i);
        }
        newStores(newStores.Rows() - 1) = consumersSum - storesSum;

        Matrix<double> newCost(cost.Rows() + 1, cost.Cols());
        for (size_t i = 0; i < cost.Rows(); i ++) {
            for (size_t j = 0; j < cost.Cols(); j ++) {
                newCost(i, j) = cost(i, j);
            }
        }
        stores = newStores;
        cost = newCost;
    } else if (consumersSum < storesSum) {
        Matrix<double> newConsumers(consumers.Rows() + 1, 1);
        for (size_t i = 0; i < consumers.Rows(); i ++) {
            newConsumers(i) = consumers(i);
        }
        newConsumers(newConsumers.Rows() - 1) = storesSum - consumersSum;

        Matrix<double> newCost(cost.Rows(), cost.Cols() + 1);
        for (size_t i = 0; i < cost.Rows(); i ++) {
            for (size_t j = 0; j < cost.Cols(); j ++) {
                newCost(i, j) = cost(i, j);
            }
        }

        consumers = newConsumers;
        cost = newCost;
    }

    std::cout << "stores: " << stores << std::endl;
    std::cout << "consumers: " << consumers << std::endl;
    std::cout << "cost: " << cost << std::endl;

    Matrix<double> result(cost.Rows(), cost.Cols());
    Matrix<double> stores_remain = stores;
    Matrix<double> consumers_remain = consumers;

    std::vector<PlanIndex> basis_indices;

    // Fill result
    {
        size_t row = 0;
        size_t column = 0;

        do {
            const auto add = std::min(consumers_remain(column), stores_remain(row));
            result(row, column) = add;

            consumers_remain(column) -= add;
            stores_remain(row) -= add;

            basis_indices.emplace_back(row, column);
            if (stores_remain(row) < 1e-9) {
                ++row;
            } else if (consumers_remain(column) < 1e-9) {
                ++column;
            }
        } while (row < result.Rows() && column < result.Cols());
    }

    // check degenerate
    assert(basis_indices.size() == stores.Rows() + consumers.Rows() - 1);
    std::cout << "result: " << std::endl << result << std::endl;

    while (true) {
        Matrix<double> workCost(cost.Rows(), cost.Cols());
        for (const auto& basis_index : basis_indices) {
            workCost(basis_index.first, basis_index.second) =
                    cost(basis_index.first, basis_index.second);
        }

        std::cout << "workCost: " << std::endl << workCost << std::endl;

        std::vector<PlanIndex> non_basis_indices;
        for (size_t i = 0; i < workCost.Rows(); ++i) {
            for (size_t j = 0; j < workCost.Cols(); ++j) {
                PlanIndex index {i, j};
                if (find(basis_indices.begin(), basis_indices.end(), index) == basis_indices.end()) {
                    non_basis_indices.push_back(index);
                }
            }
        }

        Matrix<double> workStores;
        Matrix<double> workConsumers;
        DoSolveNorthWestTask(workCost, basis_indices, workStores, workConsumers);

        std::cout << "wordStores: " << std::endl << workStores << std::endl;
        std::cout << "wordConsumers: " << std::endl << workConsumers << std::endl;

        Matrix<double> delta(workCost.Rows(), workCost.Cols());
        auto min_delta = non_basis_indices.front();
        for (const auto& index : non_basis_indices) {
            const auto i = index.first;
            const auto j = index.second;
            delta(i, j) = cost(i, j) - workStores(i) - workConsumers(j);

            if (delta(i, j) < delta(min_delta.first, min_delta.second)) {
                min_delta = index;
            }
        }

        std::cout << "delta: " << std::endl << delta << std::endl;

        if (delta(min_delta.first, min_delta.second) >= 0) {
            return result;
        }

        std::cout << "min_delta: {" << min_delta.first << ", " << min_delta.second << "}, "
             << "value: " << delta(min_delta.first, min_delta.second) << std::endl;

        auto cycle = FindSmartCycle(min_delta, basis_indices);
        std::cout << "Cycle: " << std::endl;
        for (const auto& index : cycle) {
            std::cout << "{" << index.first << ", " << index.second << "}, ";
        }
        std::cout << std::endl;
        assert(cycle.size() % 2 == 0);

        std::vector<PlanIndex> plusNode, minusNode;
        for (size_t i = 0; i < cycle.size(); ++i) {
            if (i % 2 == 0) {
                plusNode.push_back(cycle[i]);
            } else {
                minusNode.push_back(cycle[i]);
            }
        }

        size_t minus_node_index = 0;
        for (size_t i = 1; i < minusNode.size(); ++i) {
            if (result(minusNode[i].first, minusNode[i].second) <
                result(minusNode[minus_node_index].first, minusNode[minus_node_index].second))
            {
                minus_node_index = i;
            }
        }
        auto theta = result(minusNode[minus_node_index].first, minusNode[minus_node_index].second);
        std::cout << "minus_node_index: " << minus_node_index << std::endl;
        std::cout << minusNode[minus_node_index].first << " " << minusNode[minus_node_index].second << std::endl;
        std::cout << "theta: " << theta << std::endl;

        for (const auto& index : plusNode) {
            result(index.first, index.second) += theta;
        }
        for (const auto& index : minusNode) {
            result(index.first, index.second) -= theta;
        }

        std::cout << "new_result: " << std::endl << result << std::endl;

        auto remove_iter = find(basis_indices.begin(), basis_indices.end(), minusNode[minus_node_index]);
        assert(remove_iter != basis_indices.end());
        basis_indices.erase(remove_iter);
        basis_indices.push_back(min_delta);

        std::cout << "new_basis: " << std::endl;
        for (const auto& index : basis_indices) {
            std::cout << "{" << index.first << ", " << index.second << "}, ";
        }
        std::cout << std::endl;
    }
    return cost;
}


int test_1() {
    Matrix<double> stores ({{20, 30, 25}}); stores = stores.GetT();
    Matrix<double> consumers ({{10, 10, 10, 10, 10}}); consumers = consumers.GetT();

    Matrix<double> cost ({
            {2,  8, -5, 7,  10},
            {11, 5, 8,  -8, -4},
            {1,  3, 7,  4,  2}
    });

    auto answer = DoTransportTask(stores, consumers, cost);
    std::cout << "answer: " << std::endl;
    std::cout << answer << std::endl;

    double all_cost = 0;
    for (size_t i = 0; i < cost.Rows(); ++i) {
        for (size_t j = 0; j < cost.Cols(); ++j) {
            all_cost += cost(i, j) * answer(i, j);
        }
    }

    std::cout << "all_cost: " << all_cost << std::endl;

    return 0;
}

int test_2() {
    Matrix<double> stores ({{20, 11, 18, 27}}); stores = stores.GetT();
    Matrix<double> consumers ({{11, 4, 10, 12, 8, 9, 10, 4}}); consumers = consumers.GetT();

    Matrix<double> cost({
            {-3, 6, 7,  12, 6, -3, 2,  16},
            {4,  3, 7,  10, 0, 1,  -3, 7},
            {19, 3, 2,  7,  3, 7,  8,  15},
            {1,  4, -7, -3, 9, 13, 17, 22}
    });

    auto answer = DoTransportTask(stores, consumers, cost);
    std::cout << "answer: " << std::endl;
    std::cout << answer << std::endl;

    double all_cost = 0;
    for (size_t i = 0; i < cost.Rows(); ++i) {
        for (size_t j = 0; j < cost.Cols(); ++j) {
            all_cost += cost(i, j) * answer(i, j);
        }
    }

    std::cout << "all_cost: " << all_cost << std::endl;

    return 0;
}


int test_3() {
    Matrix<double> stores ({{15, 12, 18, 20}}); stores = stores.GetT();
    Matrix<double> consumers ({{5, 5, 10, 4, 6, 20, 10, 5}}); consumers = consumers.GetT();

    Matrix<double> cost({
                                {-3, 10, 70, -3, 7, 4, 2, -20},
                                {3, 5, 8, 8, 0, 1, 7, -10},
                                {-15, 1, 0, 0, 13, 5, 4, 5},
                                {1, -5, 9, -3, -4, 7, 16, 25}
                        });

    auto answer = DoTransportTask(stores, consumers, cost);
    std::cout << "answer: " << std::endl;
    std::cout << answer << std::endl;

    double all_cost = 0;
    for (size_t i = 0; i < cost.Rows(); ++i) {
        for (size_t j = 0; j < cost.Cols(); ++j) {
            all_cost += cost(i, j) * answer(i, j);
        }
    }

    std::cout << "all_cost: " << all_cost << std::endl;

    return 0;
}

int test_4() {
    Matrix<double> stores ({{53, 20, 45, 38}}); stores = stores.GetT();
    Matrix<double> consumers ({{15, 31, 10, 3, 18}}); consumers = consumers.GetT();

    Matrix<double> cost({
            {3, 0, 3, 1, 6},
            {2, 4, 10, 5, 7},
            {-2, 5, 3, 2, 9},
            {1, 3, 5, 1, 9}
    });

    auto answer = DoTransportTask(stores, consumers, cost);
    std::cout << "answer: " << std::endl;
    std::cout << answer << std::endl;

    double all_cost = 0;
    for (size_t i = 0; i < cost.Rows(); ++i) {
        for (size_t j = 0; j < cost.Cols(); ++j) {
            all_cost += cost(i, j) * answer(i, j);
        }
    }

    std::cout << "all_cost: " << all_cost << std::endl;

    return 0;
}


int test_5() {
    Matrix<double> stores ({{13, 5, 7, 9, 10}}); stores = stores.GetT();
    Matrix<double> consumers ({{20, 5, 6, 11}}); consumers = consumers.GetT();

    Matrix<double> cost({
            {2, 6, 8, -3},
            {3, 2, 12, 4},
            {7, 2, 5, 7},
            {9, 2, 14, 9},
            {8, 7, 8, 8}
    });

    auto answer = DoTransportTask(stores, consumers, cost);
    std::cout << "answer: " << std::endl;
    std::cout << answer << std::endl;

    double all_cost = 0;
    for (size_t i = 0; i < cost.Rows(); ++i) {
        for (size_t j = 0; j < cost.Cols(); ++j) {
            all_cost += cost(i, j) * answer(i, j);
        }
    }

    std::cout << "all_cost: " << all_cost << std::endl;

    return 0;
}

int test_6() {
    Matrix<double> stores ({{7, 3, 7, 3, 7}}); stores = stores.GetT();
    Matrix<double> consumers ({{10, 10, 4, 3}}); consumers = consumers.GetT();

    Matrix<double> cost({
            {1, 1, -1, -1},
            {0, 0, 2, 6},
            {5, 4, 7, 6},
            {7, 8, 5, 7},
            {2, 5, 10, 2}
    });

    auto answer = DoTransportTask(stores, consumers, cost);
    std::cout << "answer: " << std::endl;
    std::cout << answer << std::endl;

    double all_cost = 0;
    for (size_t i = 0; i < cost.Rows(); ++i) {
        for (size_t j = 0; j < cost.Cols(); ++j) {
            all_cost += cost(i, j) * answer(i, j);
        }
    }

    std::cout << "all_cost: " << all_cost << std::endl;

    return 0;
}

int main() {
//    return test_1();
//    return test_2();
//    return test_3();
//    return test_4();
//    return test_5();
    return test_6();
}
