#include "diffeq.h"
#include "solver.h"


static void make_identity(GiNaC::matrix& _matrix) {
    int rows = _matrix.rows();
    for (int i = 0; i < rows; i++)
        _matrix(i, i) = 1;
}

static GiNaC::matrix subs(const GiNaC::matrix& _matrix, const GiNaC::ex& _rule) {
    return GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)_matrix).subs(_rule, GiNaC::subs_options::algebraic));
}

static int ca_mat_jordan_blocks_with_known_eigenvalues(ca_vec_t lambda, slong * mult, slong * num_blocks, slong * block_lambda, slong * block_size, ca_mat_t A, ca_ctx_t ctx) {
    slong i, j, k, n;
    slong * ranks, * diagram;
    slong ranks_len, rank;
    fmpz * exp;
    slong exp_i;
    int status = 1;

    n = ca_mat_nrows(A);

    if (n != ca_mat_ncols(A))
    {
        /* matrix must be square */
        return 0;
    }

    ranks = (slong *)flint_malloc(sizeof(slong) * (n + 1));
    diagram = (slong *)flint_malloc(sizeof(slong) * n);

    if (status)
    {
        exp = mult;

        *num_blocks = 0;

        for (i = 0; status && i < ca_vec_length(lambda, ctx); i++)
        {
            exp_i = exp[i];

            if (exp_i == 1)
            {
                block_lambda[*num_blocks] = i;
                block_size[*num_blocks] = 1;
                *num_blocks += 1;
            }
            else
            {
                ca_mat_t B, C;

                ca_mat_init(B, n, n, ctx);
                ca_mat_init(C, n, n, ctx);

                for (j = 0; j < n; j++)
                    for (k = 0; k < n; k++)
                        if (j == k)
                            ca_sub(ca_mat_entry_ptr(B, j, j), ca_mat_entry_ptr(A, j, j), ca_vec_entry_ptr(lambda, i), ctx);
                        else
                            ca_set(ca_mat_entry_ptr(B, j, k), ca_mat_entry_ptr(A, j, k), ctx);

                ca_mat_set(C, B, ctx);

                status &= ca_mat_rank(&rank, C, ctx);

                ranks_len = 2;
                ranks[0] = n;
                ranks[1] = rank;

                j = 0;
                while (status && (ranks[j] > ranks[j + 1] && ranks[j + 1] + exp_i > n))
                {
                    ca_mat_mul(C, B, C, ctx);
                    status &= ca_mat_rank(&rank, C, ctx);
                    ranks[ranks_len] = rank;
                    j++;
                    ranks_len++;
                }

                if (status)
                {
                    /* Ferrer's diagram of an integer partition */
                    for (j = 0; j < ranks_len - 1; j++)
                        diagram[j] = ranks[j] - ranks[j + 1];

                    /* Transpose Ferrer's diagram */
                    for (j = 1; j <= diagram[0]; j++)
                    {
                        slong c = 0;

                        for (k = 0; k < ranks_len - 1; k++)
                            c += (diagram[k] >= j);

                        block_lambda[*num_blocks] = i;
                        block_size[*num_blocks] = c;
                        *num_blocks += 1;
                    }
                }

                ca_mat_clear(B, ctx);
                ca_mat_clear(C, ctx);
            }
        }
    }

    flint_free(ranks);
    flint_free(diagram);

    return status;
}


#define NORMAL(m) GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(m)).normal())


diffeq::diffeq(const GiNaC::matrix& _coeff, const GiNaC::symbol& _var) {
    initialize(subs(_coeff, _var == x));
}


diffeq diffeq::from_string(const std::string& _str) {
    GiNaC::symbol x("x");
    GiNaC::symtab tab;
    tab["x"] = x;
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
    return diffeq(coeff, x);
}


void regular_struct_t::initialize(const GiNaC::matrix& _A0) {
    // step 1: compute eigenvalues of A0
    int N = _A0.rows();
    A0.initialize(ctx);
    J.initialize(ctx, N, N);
    P.initialize(ctx, N, N);
    A0.import_from_ginac(_A0);
    lambdas.initialize(ctx, N);
    mult = std::vector<ulong>(N);
    block_lambda = std::vector<slong>(N);
    block_size = std::vector<slong>(N);
    if (!ca_mat_eigenvalues(lambdas.ptr(), mult.data(), A0.ptr(), ctx.ctx))
        throw std::runtime_error("failed to compute eigenvalues of A0!");

    // step 2: initial Jordan decomposition of A0
    if (!ca_mat_jordan_blocks_with_known_eigenvalues(lambdas.ptr(), (slong*)mult.data(), &num_blocks, block_lambda.data(), block_size.data(), A0.ptr(), ctx.ctx))
        throw std::runtime_error("failed to Jordan decompose A0!");
    ca_mat_set_jordan_blocks(J.ptr(), lambdas.ptr(), num_blocks, block_lambda.data(), block_size.data(), ctx.ctx);
    if (!ca_mat_jordan_transformation(P.ptr(), lambdas.ptr(), num_blocks, block_lambda.data(), block_size.data(), A0.ptr(), ctx.ctx)) 
        throw std::runtime_error("failed to Jordan decompose A0!");

    // step 3: initial grouping of eigenvalues
    group_eigenvalues();

    // step 4: save initial Jordan block structure
    block_structure = std::vector<int>(N, 0);
    block_structure_idx = std::vector<int>(N, 0);
    for (int i = 1; i < N; i++)
        if (ca_check_is_one(J(i - 1, i), ctx.ctx) == T_TRUE) {
            block_structure[i] = block_structure[i - 1];
            block_structure_idx[i] = block_structure_idx[i - 1];
        } else {
            block_structure[i] = i;
            block_structure_idx[i] = block_structure_idx[i - 1] + 1;
        }
}


void regular_struct_t::update() {
    // step 1: Jordan decomposition of J
    if (!ca_mat_jordan_blocks_with_known_eigenvalues(lambdas.ptr(), (slong*)mult.data(), &num_blocks, block_lambda.data(), block_size.data(), J.ptr(), ctx.ctx))
        throw std::runtime_error("failed to Jordan decompose J!");
    if (!ca_mat_jordan_transformation(P.ptr(), lambdas.ptr(), num_blocks, block_lambda.data(), block_size.data(), J.ptr(), ctx.ctx)) 
        throw std::runtime_error("failed to Jordan decompose J!");
    ca_mat_set_jordan_blocks(J.ptr(), lambdas.ptr(), num_blocks, block_lambda.data(), block_size.data(), ctx.ctx);
    
    // step 2: initial grouping of eigenvalues
    group_eigenvalues();

    // step 3: save Jordan block structure
    int N = ca_mat_nrows(A0.ptr());
    block_structure[0] = 0;
    block_structure_idx[0] = 0;
    for (int i = 1; i < N; i++)
        if (ca_check_is_one(J(i - 1, i), ctx.ctx) == T_TRUE) {
            block_structure[i] = block_structure[i - 1];
            block_structure_idx[i] = block_structure_idx[i - 1];
        } else {
            block_structure[i] = i;
            block_structure_idx[i] = block_structure_idx[i - 1] + 1;
        }
}


void regular_struct_t::group_eigenvalues() {
    eiggroups = std::vector<std::vector<std::pair<int, int>>>();
    int N = ca_mat_nrows(A0.ptr());
    flint::ca_number diff(ctx);
    std::vector<int> searched(N, 0);
    for (int i = 0; i < N; i++) {
        if (!searched[i]) {
            searched[i] = 1;
            eiggroups.push_back({{i, 0}});
            for (int j = i + 1; j < N; j++) {
                if (!searched[j]) {
                    diff.call(ca_sub, J(j, j), J(i, i));
                    if (ca_check_is_integer(diff.ptr(), ctx.ctx) == T_TRUE) {
                        fmpz_t intdiff;
                        fmpz_init(intdiff);
                        if (!ca_get_fmpz(intdiff, diff.ptr(), ctx.ctx)) {
                            fmpz_clear(intdiff);
                            throw std::runtime_error("failed to convert eigenvalue difference to integer!");
                        }
                        int intdiff_c = fmpz_get_si(intdiff);
                        eiggroups.back().push_back({j, intdiff_c});
                        searched[j] = 1;
                        fmpz_clear(intdiff);
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


int diffeq::regular_reduction_one_step(const GiNaC::ex& x0, unsigned digits) {
    // step 1: apply existing P to the full A(x)
    update(reg_struct.P.export_to_ginac(digits));

    // step 2: check if there are eigenvalues differing by integers from each other
    std::set<int> front_indices;
    bool merge = false, need = false;
    int merge_eigval = -1, modify_eigval = -1;
    for (auto& group: reg_struct.eiggroups) {
        front_indices = std::set<int>();
        for (auto& pair: group) {
            if (pair.second == 0)
                front_indices.insert(pair.first);
            else {
                need = true;
                modify_eigval = reg_struct.block_lambda[reg_struct.block_structure_idx[*front_indices.begin()]];
                if (pair.second == -1) {
                    merge = true;
                    merge_eigval = reg_struct.block_lambda[reg_struct.block_structure_idx[pair.first]];
                }
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
    flint::ca_matrix Q(reg_struct.ctx, N, N);
    int sz = front_indices_seq.size();
    for (int i = 0; i < sz; i++) {
        ca_set_si(Q(front_indices_seq[i], i), 1, reg_struct.ctx.ctx);
        front_indices.insert(front_indices_seq[i]);
    }
    int i = sz;
    for (int j = 0; j < N; j++) {
        if (front_indices.find(j) == front_indices.end()) {
            ca_set_si(Q(j, i++), 1, reg_struct.ctx.ctx);
        }
    }
    update(Q.export_to_ginac(digits));

    // step 4: symbolic reduction
    GiNaC::matrix Px(N, N);
    make_identity(Px);
    for (int i = 0; i < sz; i++)
        Px(i, i) = (x - x0);
    update(Px);
    coeff = NORMAL(coeff);
    transform = NORMAL(transform);

    // step 5: update J
    flint::ca_matrix Q_inv(reg_struct.ctx, N, N), tempJQ(reg_struct.ctx, N, N);
    ca_mat_inv(Q_inv.ptr(), Q.ptr(), reg_struct.ctx.ctx);
    ca_mat_mul(tempJQ.ptr(), reg_struct.J.ptr(), Q.ptr(), reg_struct.ctx.ctx);
    ca_mat_mul(reg_struct.J.ptr(), Q_inv.ptr(), tempJQ.ptr(), reg_struct.ctx.ctx);

    GiNaC::matrix A0_(N, N);
    GiNaC::symbol t("t");
    for (int i = 0; i < sz; i++)
        for (int j = sz; j < N; j++)
            A0_(i, j) = subs(coeff(i, j), x == x0 + t).series(t, 1).coeff(t, -1);

    if (merge)
        for (int i = 0; i < sz; i++)
            ca_set(reg_struct.J(i, i), reg_struct.lambdas[merge_eigval], reg_struct.ctx.ctx);
    else
        for (int i = 0; i < sz; i++)
            ca_sub_si(reg_struct.J(i, i), reg_struct.J(i, i), 1, reg_struct.ctx.ctx);

    flint::ca_number number(reg_struct.ctx);
    for (int i = 0; i < sz; i++) {
        for (int j = sz; j < N; j++) {
            gq2f(A0_(i, j), number.ptr(), reg_struct.ctx.ctx);
            ca_set(reg_struct.J(i, j), number.ptr(), reg_struct.ctx.ctx);
        }
    }
    
    // step 6: merge or modify eigenvalues
    if (merge) {
        slong new_length = ca_vec_length(reg_struct.lambdas.ptr(), reg_struct.ctx.ctx) - 1;
        ca_set(reg_struct.lambdas[modify_eigval], reg_struct.lambdas[new_length], reg_struct.ctx.ctx);
        ca_vec_set_length(reg_struct.lambdas.ptr(), new_length, reg_struct.ctx.ctx);
        reg_struct.mult[merge_eigval] += sz;
        reg_struct.mult[modify_eigval] = reg_struct.mult[new_length];
    } else {
        ca_sub_si(reg_struct.lambdas[modify_eigval], reg_struct.lambdas[modify_eigval], 1, reg_struct.ctx.ctx);
    }

    // step 7: update reg_struct
    reg_struct.update();

    return 1;
}


void diffeq::regular_reduction(const GiNaC::ex& x0, unsigned digits) {
    std::cerr << "DiffEqSolver: start regular reduction\n";
    int reduction_step = 1;
    do {
        std::cerr << "DiffEqSolver: perform regular reduction for " << reduction_step++ << " step\r";
    } while (regular_reduction_one_step(x0, digits));
    std::cerr << "\nDiffEqSolver: finish regular reduction\n";
}


std::pair<GiNaC::matrix, GiNaC::ex> diffeq::solve(const GiNaC::ex& x0, unsigned order, unsigned digits) {
    clear_all_reduction();
    moser_reduction(x0);
    
    GiNaC::symbol t("t");
    GiNaC::matrix coeff_ = subs(coeff, x == x0 + t);
    if (is_analytic(x0)) {
        solver_analytic solver(coeff_, t, order);
        auto out = GiNaC::ex_to<GiNaC::matrix>(
            subs(transform, x == x0 + t).mul(solver.solution(digits)).expand().subs(t == x - x0)
        );
        auto radius = solver.radius();
        return {out, radius};
    }

    if (!(bool)(get_moser_struct(x0).moser_rank() == 1))
        throw std::runtime_error("not implemented: x0 is an irregular singularity");

    // solve for regular singularity
    GiNaC::matrix A0(N(), N()), A_ana(N(), N());
    for (int i = 0; i < N(); i++)
        for (int j = 0; j < N(); j++)
            A0(i, j) = coeff_(i, j).series(t, 1).coeff(t, -1);

    reg_struct.initialize(A0);
    regular_reduction(x0, digits);
    // now we can remove spurious singular part of coefficients
    for (int i = 0; i < N(); i++) {
        for (int j = 0; j < N(); j++) {
            coeff_(i, j) = coeff(i, j).subs(x == x0 + t, GiNaC::subs_options::algebraic).normal();
            auto numer_denom = coeff_(i, j).numer_denom();
            auto numer = numer_denom[0], denom = numer_denom[1];  
            auto ldeg = denom.ldegree(t);
            auto pre_expansion = GiNaC::series_to_poly(coeff_(i, j).series(t, 0));
            pre_expansion = pre_expansion - pre_expansion.coeff(t, 0);
            A_ana(i, j) = (GiNaC::quo((numer - pre_expansion * denom).expand(), GiNaC::pow(t, ldeg), t, false) / GiNaC::quo(denom, GiNaC::pow(t, ldeg), t, false)).normal();
        }
    }
    
    solver_regular solver(A_ana, t, order, reg_struct);
    auto out =  GiNaC::ex_to<GiNaC::matrix>(
        subs(transform, x == x0 + t).mul(solver.solution(digits)).expand().subs(t == x - x0)
    );
    auto radius = solver.radius();

    return {out, radius};
}


std::pair<GiNaC::matrix, GiNaC::ex> diffeq::solve(const GiNaC::ex& x0, const GiNaC::matrix& Y0, unsigned order, unsigned digits) {
    clear_all_reduction();
    moser_reduction(x0);
    
    if (is_analytic(x0)) {
        GiNaC::symbol t("t");
        solver_analytic solver(subs(coeff, x == x0 + t), t, Y0, order);
        auto out = GiNaC::ex_to<GiNaC::matrix>(
            subs(transform, x == x0 + t).mul(solver.solution(digits)).expand().subs(t == x - x0)
        );
        auto radius = solver.radius();
        return {out, radius};
    }

    throw std::runtime_error("not implemented: x0 is a singularity");
}


