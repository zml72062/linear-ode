#include "base.h"


static void make_identity(GiNaC::matrix& _matrix) {
    int rows = _matrix.rows();
    for (int i = 0; i < rows; i++)
        _matrix(i, i) = 1;
}

static GiNaC::matrix subs(const GiNaC::matrix& _matrix, const GiNaC::ex& _rule) {
    return GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)_matrix).subs(_rule, GiNaC::subs_options::algebraic));
}

static GiNaC::matrix append_column(const GiNaC::matrix& _matrix, const GiNaC::matrix& _col) {
    int r = _matrix.rows(), c = _matrix.cols();
    GiNaC::matrix newmat(r, c + 1);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            newmat(i, j) = _matrix(i, j);
        newmat(i, c) = _col(i, 0);
    }
    return newmat;
}

static GiNaC::matrix append_row(const GiNaC::matrix& _matrix, const GiNaC::matrix& _row) {
    int r = _matrix.rows(), c = _matrix.cols();
    GiNaC::matrix newmat(r + 1, c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            newmat(i, j) = _matrix(i, j);
    for (int j = 0; j < c; j++)
        newmat(r, j) = _row(0, j);
    return newmat;
}


#define SUBMATRIX(m, r, nr, c, nc) GiNaC::ex_to<GiNaC::matrix>(GiNaC::sub_matrix((m), (r), (nr), (c), (nc)))
#define NORMAL(m) GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(m)).normal())


void base_diffeq::initialize(const GiNaC::matrix& _matrix) {
    if (_matrix.cols() != _matrix.rows())
        throw std::runtime_error("Input not a square matrix!");
    raw_coeff = _matrix;
    coeff = raw_coeff;
    transform = GiNaC::matrix(_matrix.rows(), _matrix.rows());
    make_identity(transform);
}


moser_t base_diffeq::get_moser_struct(const GiNaC::ex& x0) {
    moser_t info(N());
    GiNaC::symbol t("t");
    GiNaC::matrix coeff_ = subs(coeff, x == x0 + t);
    for (int i = 0; i < N(); i++)
        for (int j = 0; j < N(); j++)
            coeff_(i, j) = coeff_(i, j).series(t, 0);
    
    info.negative_power = 0;
    for (int i = 0; i < N(); i++) {
        for (int j = 0; j < N(); j++) {
            int ndegree = coeff_(i, j).ldegree(t);
            if (-ndegree > info.negative_power)
                info.negative_power = -ndegree;
        }
    }

    if (info.negative_power <= 1)
        return info;
    
    for (int i = 0; i < N(); i++) {
        for (int j = 0; j < N(); j++) {
            info.A0(i, j) = coeff_(i, j).coeff(t, -info.negative_power).normal();
            info.A1(i, j) = coeff_(i, j).coeff(t, -info.negative_power + 1).normal();
        }
    }
    return info;
}


void base_diffeq::moser_reduction_one_step(const GiNaC::ex& x0) {
    // step 1: rank reduction of A0
    auto moser_info = get_moser_struct(x0);
    if (moser_info.moser_rank() == 1) // already regular
        return;

    std::vector<int> is_indep_col(N(), 0);
    GiNaC::matrix indep_cols(N(), 0);
    int r = moser_info.A0.rank(), n_indep_cols = 0, i = 0;
    while (n_indep_cols != r) {
        auto next_col = SUBMATRIX(moser_info.A0, 0, N(), i++, 1);
        auto new_indep_cols = append_column(indep_cols, next_col);
        if ((int)new_indep_cols.rank() > n_indep_cols) {
            indep_cols = new_indep_cols;
            n_indep_cols++;
            is_indep_col[i - 1] = 1;
        }
    }

    GiNaC::matrix permutation(N(), N());
    int m1 = 0, m2 = r;
    for (int i = 0; i < N(); i++) {
        if (is_indep_col[i])
            permutation(i, m1++) = 1;
        else
            permutation(i, m2++) = 1;
    }
    
    auto newA0 = NORMAL(moser_info.A0.mul(permutation));
    auto newA0_coeff = SUBMATRIX(newA0, 0, N(), 0, r);
    auto newA0_bias  = SUBMATRIX(newA0, 0, N(), r, N() - r);
    GiNaC::matrix vars(r, N() - r);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < N() - r; j++)
            vars(i, j) = GiNaC::symbol("temp_" + std::to_string(i) + "_" + std::to_string(j));
    // there must exist a unique solution, since 
    // rank(newA0_coeff) = rank(newA0_coeff, newA0_bias) = r
    auto combini = NORMAL(newA0_coeff.solve(vars, newA0_bias));

    GiNaC::matrix lincomb(N(), N());
    for (int i = 0; i < N(); i++)
        lincomb(i, i) = 1;
    for (int i = 0; i < r; i++)
        for (int j = r; j < N(); j++)
            lincomb(i, j) = -combini(i, j - r);

    auto new_transform = NORMAL(permutation.mul(lincomb));
    update(new_transform);
    coeff = NORMAL(coeff);
    transform = NORMAL(transform);

    // step 2: decide Moser reducibility
    moser_info = get_moser_struct(x0);
    r = moser_info.A0.rank();

    GiNaC::matrix moser_discriminant(N(), N());
    for (int i = 0; i < N(); i++) {
        for (int j = 0; j < r; j++)
            moser_discriminant(i, j) = moser_info.A0(i, j);
        for (int j = r; j < N(); j++)
            moser_discriminant(i, j) = moser_info.A1(i, j);
    }

    GiNaC::symbol lambda("lambda");
    for (int i = r; i < N(); i++)
        moser_discriminant(i, i) += lambda; 

    /**
     * In principle we should test whether 
     * 
     *      moser_discriminant.determinant() == 0
     * 
     * for an arbitrary lambda, but this would be slow symbolically.
     * 
     * Instead, we check this for (N + 1) different lambda values,
     * which suffices for an at-most-N-degree polynomial.
     */
    for (int i = 0; i < N() + 1; i++) {
        if (!(bool)(subs(moser_discriminant, lambda == i).determinant().expand() == 0))
            // not Moser reducible
            return;
    }
    
    // step 3: deal with the first case, where the first r rows of Moser 
    //         discriminant are not linearly independent
    int subrank = SUBMATRIX(moser_discriminant, 0, r, 0, N()).rank();
    if (subrank < r) {
        GiNaC::matrix T(N(), N());
        make_identity(T);
        for (int i = 0; i < r; i++)
            T(i, i) = (x - x0);
        update(T);
        coeff = NORMAL(coeff);
        transform = NORMAL(transform);
        return;
    }

    // step 4: deal with the second case, where the first r rows of Moser
    //         discriminant are linearly independent
    for (int i = r; i < N(); i++)
        moser_discriminant(i, i) -= lambda; 

    std::vector<int> is_indep_row(N(), 0);
    for (int i = 0; i < r; i++)
        is_indep_row[i] = 1;
    GiNaC::matrix indep_rows = SUBMATRIX(moser_discriminant, 0, r, 0, N());
    int h = N() - moser_discriminant.rank(), n_indep_rows = r, j = r;
    while (n_indep_rows != N() - h) {
        auto next_row = SUBMATRIX(moser_discriminant, j++, 1, 0, N());
        auto new_indep_rows = append_row(indep_rows, next_row);
        if ((int)new_indep_rows.rank() > n_indep_rows) {
            indep_rows = new_indep_rows;
            n_indep_rows++;
            is_indep_row[j - 1] = 1;
        }
    }
    GiNaC::matrix inv_permutation(N(), N());
    int n1 = 0, n2 = N() - h;
    for (int i = 0; i < N(); i++) {
        if (is_indep_row[i])
            inv_permutation(n1++, i) = 1;
        else
            inv_permutation(n2++, i) = 1;
    }

    auto newD = NORMAL(inv_permutation.mul(moser_discriminant));
    auto newD_coeff = SUBMATRIX(newD, 0, N() - h, 0, N()).transpose();
    auto newD_bias  = SUBMATRIX(newD, N() - h, h, 0, N()).transpose();
    GiNaC::matrix varsD(N() - h, h);
    for (int i = 0; i < N() - h; i++)
        for (int j = 0; j < h; j++)
            varsD(i, j) = GiNaC::symbol("tempD_" + std::to_string(i) + "_" + std::to_string(j)); 
    auto combiniD = NORMAL(newD_coeff.solve(varsD, newD_bias));
    GiNaC::matrix lincombD(N(), N());
    for (int i = 0; i < N(); i++)
        lincombD(i, i) = 1;
    for (int i = N() - h; i < N(); i++)
        for (int j = 0; j < N() - h; j++)
            lincombD(i, j) = -combiniD(j, i - (N() - h));
    update(lincombD.mul(inv_permutation).inverse());

    GiNaC::matrix T(N(), N());
    make_identity(T);
    for (int i = 0; i < r; i++)
        T(i, i) = (x - x0);
    for (int i = N() - h; i < N(); i++)
        T(i, i) = (x - x0);
    update(T);
    coeff = NORMAL(coeff);
    transform = NORMAL(transform);
}


void base_diffeq::moser_reduction(const GiNaC::ex& x0) {
    std::cerr << "DiffEqSolver: start Moser reduction\n";
    GiNaC::ex moser_rank = get_moser_struct(x0).moser_rank();
    int reduction_step = 1;
    while (true) {
        moser_reduction_one_step(x0);
        std::cerr << "DiffEqSolver: perform Moser reduction for " << reduction_step++ << " step\r";
        GiNaC::ex new_moser_rank = get_moser_struct(x0).moser_rank();
        if (new_moser_rank < moser_rank)
            moser_rank = new_moser_rank;
        else if (new_moser_rank == moser_rank)
            break;
        else
            throw std::logic_error("Moser reduction runs into an impossible branch!");
    }
    std::cerr << "\nDiffEqSolver: finish Moser reduction\n";
}


bool base_diffeq::is_analytic(const GiNaC::ex& x0) {
    return get_moser_struct(x0).negative_power == 0;
}


bool base_diffeq::is_first_order(const GiNaC::ex& x0) {
    return get_moser_struct(x0).negative_power == 1;
}


bool base_diffeq::is_regular(const GiNaC::ex& x0) {
    moser_reduction(x0);
    return get_moser_struct(x0).moser_rank() == 1;
}


void base_diffeq::update(const GiNaC::matrix& new_transform) {
    transform = transform.mul(new_transform);
    coeff = new_transform.inverse().mul(coeff.mul(new_transform).sub(GiNaC::ex_to<GiNaC::matrix>(new_transform.diff(x))));
}


void base_diffeq::clear_all_reduction() {
    coeff = raw_coeff;
    transform = GiNaC::ex_to<GiNaC::matrix>(GiNaC::unit_matrix(N()));
}

