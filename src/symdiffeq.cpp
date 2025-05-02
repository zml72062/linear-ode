#include "symdiffeq.h"
#include "symsolver.h"
#include "jordan.h"


static GiNaC::matrix subs(const GiNaC::matrix& _matrix, const GiNaC::ex& _rule) {
    return GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)_matrix).subs(_rule, GiNaC::subs_options::algebraic));
}

static void make_identity(GiNaC::matrix& _matrix) {
    int rows = _matrix.rows();
    for (int i = 0; i < rows; i++)
        _matrix(i, i) = 1;
}


#define NORMAL(mat) GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(mat)).normal())


symdiffeq::symdiffeq(const GiNaC::matrix& _coeff, const GiNaC::symbol& _var, const GiNaC::symbol& _eps)
    : eps("eps") {
    initialize(subs(_coeff, GiNaC::lst{_var == x, _eps == eps}));
}


symdiffeq symdiffeq::from_string(const std::string& _str) {
    GiNaC::symbol x("x"), eps("eps");
    GiNaC::symtab tab;
    tab["x"] = x;
    tab["eps"] = eps;
    GiNaC::parser parser(tab);
    auto result = parser(_str);
    int rows = result.nops();
    GiNaC::matrix coeff(rows, rows);
    int i = 0, j;
    for (auto r = result.begin(); r != result.end(); ++r, i++) {
        if ((int)r->nops() != rows)
            throw std::runtime_error("Input not a square matrix!");
        j = 0;
        for (auto c = r->begin(); c != r->end(); ++c, j++) {
            coeff(i, j) = *c;
        }
    }
    return symdiffeq(coeff, x, eps);
}


void symbolic_regular_struct_t::initialize(const GiNaC::matrix& _A0) {
    // step 1: compute eigenvalues of A0
    A0 = _A0;
    int N = _A0.rows();
    mult = std::vector<int>(N);
    block_lambda = std::vector<int>(N);
    block_size = std::vector<int>(N);
    
    jordan j(A0);
    lambdas = j.eigenvalues(mult.data());

    // step 2: initial Jordan decomposition of A0
    J = j.jordan_blocks(lambdas, mult.data(), &num_blocks, block_lambda.data(), block_size.data());
    P = j.jordan_transform(lambdas, num_blocks, block_lambda.data(), block_size.data());

    // step 3: initial grouping of eigenvalues
    group_eigenvalues();

    // step 4: save initial Jordan block structure
    block_structure = std::vector<int>(N, 0);
    block_structure_idx = std::vector<int>(N, 0);
    for (int i = 1; i < N; i++)
        if (J(i - 1, i) == 1) {
            block_structure[i] = block_structure[i - 1];
            block_structure_idx[i] = block_structure_idx[i - 1];
        } else {
            block_structure[i] = i;
            block_structure_idx[i] = block_structure_idx[i - 1] + 1;
        }
}


void symbolic_regular_struct_t::group_eigenvalues() {
    eiggroups = std::vector<std::vector<std::pair<int, int>>>();
    int N = A0.rows();
    std::vector<int> searched(N, 0);
    for (int i = 0; i < N; i++) {
        if (!searched[i]) {
            searched[i] = 1;
            eiggroups.push_back({{i, 0}});
            for (int j = i + 1; j < N; j++) {
                if (!searched[j]) {
                    GiNaC::ex diff = (J(j, j) - J(i, i)).expand().normal();
                    if (GiNaC::ex_to<GiNaC::numeric>(diff).is_integer()) {
                        int intdiff = GiNaC::ex_to<GiNaC::numeric>(diff).to_int();
                        eiggroups.back().push_back({j, intdiff});
                        searched[j] = 1;
                    }
                }
            }
        }
    }

    for (auto& group: eiggroups) {
        // sort each group of eigenvalues, such that the largest
        // eigenvalue is at the front
        std::sort(group.begin(), group.end(), [](const auto& p1, const auto& p2) {
            return p1.second > p2.second;
        });
        int base = group[0].second;
        for (auto& pair: group)
            pair.second -= base;
    }
}


void symbolic_regular_struct_t::update() {
    // step 1: Jordan decomposition of J
    jordan j(J);
    lambdas = j.eigenvalues(mult.data());
    auto J1 = j.jordan_blocks(lambdas, mult.data(), &num_blocks, block_lambda.data(), block_size.data());
    auto P1 = j.jordan_transform(lambdas, num_blocks, block_lambda.data(), block_size.data());
    J = J1;
    P = P1;

    // step 2: initial grouping of eigenvalues
    group_eigenvalues();

    // step 3: save Jordan block structure
    int N = A0.rows();
    block_structure[0] = 0;
    block_structure_idx[0] = 0;
    for (int i = 1; i < N; i++)
        if (J(i - 1, i) == 1) {
            block_structure[i] = block_structure[i - 1];
            block_structure_idx[i] = block_structure_idx[i - 1];
        } else {
            block_structure[i] = i;
            block_structure_idx[i] = block_structure_idx[i - 1] + 1;
        }
}


int symdiffeq::regular_reduction_one_step(const GiNaC::ex& x0) {
    // step 1: apply existing P to the full A(x)
    update(reg_struct.P);

    // step 2: check if there are eigenvalues differing by integers from each other
    std::set<int> front_indices;
    bool need = false;
    for (auto& group: reg_struct.eiggroups) {
        front_indices = std::set<int>();
        for (auto& pair: group) {
            if (pair.second == 0)
                front_indices.insert(pair.first);
            else {
                need = true;
                break;
            }
        }
        if (need)
            break;
    }

    if (!need)
        return 0;

    std::vector<int> front_indices_seq;
    int N = reg_struct.block_structure.size();
    while (front_indices.size() > 0) {
        int bl = reg_struct.block_structure[*front_indices.begin()];
        int k;
        for (k = bl; k < N; k++) {
            if (reg_struct.block_structure[k] == bl) {
                front_indices_seq.push_back(k);
                front_indices.erase(k);
            } else break;
        }
    }

    // step 3: permutation transformation
    GiNaC::matrix Q(N, N);
    int sz = front_indices_seq.size();
    for (int i = 0; i < sz; i++) {
        Q(front_indices_seq[i], i) = 1;
        front_indices.insert(front_indices_seq[i]);
    }
    int i = sz;
    for (int j = 0; j < N; j++)
        if (front_indices.find(j) == front_indices.end())
            Q(j, i++) = 1;
    update(Q);

    // step 4: symbolic reduction
    GiNaC::matrix Px(N, N);
    make_identity(Px);
    for (int i = 0; i < sz; i++)
        Px(i, i) = (x - x0);
    update(Px);

    // step 5: update J
    reg_struct.J = NORMAL(NORMAL(Q.inverse()).mul(NORMAL(reg_struct.J.mul(Q))));
    GiNaC::matrix A0_(N, N);
    GiNaC::symbol t("t");
    for (int i = 0; i < sz; i++)
        for (int j = sz; j < N; j++)
            A0_(i, j) = subs(coeff(i, j), x == x0 + t).series(t, 1).coeff(t, -1);

    for (int i = 0; i < sz; i++)
        reg_struct.J(i, i) = (reg_struct.J(i, i) - 1).normal();
    
    for (int i = 0; i < sz; i++)
        for (int j = sz; j < N; j++)
            reg_struct.J(i, j) = A0_(i, j);
    
    // step 6: update reg_struct
    reg_struct.update();

    return 1;
}


void symdiffeq::regular_reduction(const GiNaC::ex& x0) {
    while (regular_reduction_one_step(x0));
}


std::pair<GiNaC::matrix, GiNaC::lst> symdiffeq::solve(const GiNaC::ex& x0, unsigned order) {
    clear_all_reduction();
    moser_reduction(x0);
    
    GiNaC::symbol t("t");
    GiNaC::matrix coeff_ = subs(coeff, x == x0 + t);
    if (is_analytic(x0)) {
        symsolver_analytic solver(coeff_, t, order);
        auto out = subs(NORMAL(subs(transform, x == x0 + t).mul(solver.solution())), t == x - x0);
        auto radius_eqs = GiNaC::ex_to<GiNaC::lst>(((GiNaC::ex)(solver.radius_eqs())).subs(t == x));
        return {out, radius_eqs};
    }

    if (!(bool)(get_moser_struct(x0).moser_rank() == 1))
        throw std::runtime_error("not implemented: x0 is an irregular singularity");

    // solve for regular singularity
    GiNaC::matrix A0(N(), N()), A_ana(N(), N());
    for (int i = 0; i < N(); i++)
        for (int j = 0; j < N(); j++)
            A0(i, j) = coeff_(i, j).series(t, 1).coeff(t, -1).normal();

    reg_struct.initialize(A0);
    regular_reduction(x0);
    // now we can remove spurious singular part of coefficients
    for (int i = 0; i < N(); i++) {
        for (int j = 0; j < N(); j++) {
            coeff_(i, j) = coeff(i, j).subs(x == x0 + t, GiNaC::subs_options::algebraic).normal();
            auto pre_expansion = GiNaC::series_to_poly(coeff_(i, j).series(t, 0));
            pre_expansion = pre_expansion - pre_expansion.coeff(t, 0);
            A_ana(i, j) = (coeff_(i, j) - pre_expansion).normal();
        }
    }
    
    symsolver_regular solver(A_ana, t, order, reg_struct);
    auto out = subs(NORMAL(subs(transform, x == x0 + t).mul(solver.solution())), t == x - x0);
    auto radius_eqs = GiNaC::ex_to<GiNaC::lst>(((GiNaC::ex)(solver.radius_eqs())).subs(t == x));
    return {out, radius_eqs};
}


std::pair<GiNaC::matrix, GiNaC::lst> symdiffeq::solve(const GiNaC::ex& x0, const GiNaC::matrix& Y0, unsigned order) {
    clear_all_reduction();
    moser_reduction(x0);
    
    if (is_analytic(x0)) {
        GiNaC::symbol t("t");
        symsolver_analytic solver(subs(coeff, x == x0 + t), t, Y0, order);
        auto out = subs(NORMAL(subs(transform, x == x0 + t).mul(solver.solution())), t == x - x0);
        auto radius_eqs = GiNaC::ex_to<GiNaC::lst>(((GiNaC::ex)(solver.radius_eqs())).subs(t == x));
        return {out, radius_eqs};
    }

    throw std::runtime_error("not implemented: x0 is a singularity");
}


