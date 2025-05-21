#include "solver.h"


solver_analytic::solver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, int order)
    : var(x), coeff(order + 1), sol(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");
    for (int i = 0; i <= order; i++) {
        coeff[i].initialize(ctx, N, N);
        sol[i].initialize(ctx, N, N);
    }
    rho = GiNaC::pow(10, int(GiNaC::Digits));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            flint::ca_rational r(ctx);
            flint::ca_poly     r_expand(ctx);
            r.import_from_ginac(A(i, j).normal(), x);
            GiNaC::ex rho_ij = r.radius_convergence(unsigned(GiNaC::Digits));
            if (rho_ij < rho)
                rho = rho_ij;
            if (rho_ij < GiNaC::pow(10, -6)) 
                throw std::runtime_error("A is not analytic at x0!");
            r.series(r_expand, order + 1);
            int len = r_expand.ptr()->length;
            for (int k = 0; k <= order; k++) {
                if (k < len)
                    ca_set(coeff[k](i, j), r_expand.coeff(k), ctx.ctx);
                else
                    ca_set_si(coeff[k](i, j), 0, ctx.ctx);
            }
        }
    }

    for (int i = 0; i < N; i++)
        ca_set_si(sol[0](i, i), 1, ctx.ctx);
}


solver_analytic::solver_analytic(const GiNaC::matrix& A, const GiNaC::ex& x, const GiNaC::matrix& Y0, int order)
    : var(x), coeff(order + 1), sol(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");
    if ((int)Y0.rows() != N)
        throw std::runtime_error("Y0 and A have inconsistent shapes!");

    for (int i = 0; i <= order; i++) {
        coeff[i].initialize(ctx, N, N);
        sol[i].initialize(ctx, N, Y0.cols());
    }
    rho = GiNaC::pow(10, int(GiNaC::Digits));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            flint::ca_rational r(ctx);
            flint::ca_poly     r_expand(ctx);
            r.import_from_ginac(A(i, j).normal(), x);
            GiNaC::ex rho_ij = r.radius_convergence(unsigned(GiNaC::Digits));
            if (rho_ij < rho)
                rho = rho_ij;
            if (rho_ij < GiNaC::pow(10, -6)) 
                throw std::runtime_error("A is not analytic at x0!");
            r.series(r_expand, order + 1);
            int len = r_expand.ptr()->length;
            for (int k = 0; k <= order; k++) {
                if (k < len)
                    ca_set(coeff[k](i, j), r_expand.coeff(k), ctx.ctx);
                else
                    ca_set_si(coeff[k](i, j), 0, ctx.ctx);
            }
        }
    }

    int C = Y0.cols();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < C; j++)
            gq2f(Y0(i, j), sol[0](i, j), ctx.ctx);
}


GiNaC::ex solver_analytic::radius() {
    return rho;
}


GiNaC::matrix solver_analytic::solution(unsigned dig) {
    int N = sol[0].ptr()->r, C = sol[0].ptr()->c;
    int order = coeff.size() - 1;
    for (int k = 0; k < order; k++) {
        std::cerr << "DiffEqSolver: solving order " << k << "\r";
        sol[k + 1].call(ca_mat_zero);
        flint::ca_matrix temp(ctx, N, C);
        for (int j = 0; j <= k; j++) {
            temp.call(ca_mat_mul, coeff[k - j].ptr(), sol[j].ptr());
            sol[k + 1].call(ca_mat_add, sol[k + 1].ptr(), temp.ptr());
        }
        sol[k + 1].call(ca_mat_div_si, sol[k + 1].ptr(), k + 1);
    }
    std::cerr << "\nDiffEqSolver: solving done\n";

    GiNaC::matrix M(N, C);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < C; j++) {
            M(i, j) = 0;
            for (int k = 0; k <= order; k++) {
                M(i, j) += (fc2g(sol[k](i, j), dig, ctx.ctx) * GiNaC::pow(var, k));
            }
        }
    }

    return M;
}


solver_regular::solver_regular(const GiNaC::matrix& A, const GiNaC::ex& x, int order, regular_t& rstruct)
    : var(x), prstruct(&rstruct), coeff(order + 1), solutions(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");
    for (int i = 0; i <= order; i++) {
        coeff[i].initialize(prstruct->ctx, N, N);
        solutions[i].initialize(prstruct->ctx, N, N);
    }
    rho = GiNaC::pow(10, int(GiNaC::Digits));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            flint::ca_rational r(prstruct->ctx);
            flint::ca_poly     r_expand(prstruct->ctx);
            r.import_from_ginac(A(i, j).normal(), x);
            GiNaC::ex rho_ij = r.radius_convergence(unsigned(GiNaC::Digits));
            if (rho_ij < rho)
                rho = rho_ij;
            if (rho_ij < GiNaC::pow(10, -6)) 
                throw std::runtime_error("A is not analytic at x0!");
            r.series(r_expand, order + 1);
            int len = r_expand.ptr()->length;
            for (int k = 0; k < order; k++) {
                if (k < len)
                    ca_set(coeff[k + 1](i, j), r_expand.coeff(k), prstruct->ctx.ctx);
                else
                    ca_set_si(coeff[k + 1](i, j), 0, prstruct->ctx.ctx);
            }
            ca_set(coeff[0](i, j), prstruct->J(i, j), prstruct->ctx.ctx);
        }
    }
}


void solver_regular::fill_in_solution(int psol, int puppersol) {
    int order = coeff.size() - 1, N = ca_mat_nrows(solutions[0].ptr());
    for (int i = 0; i < N; i++)
        ca_set_si(solutions[0](i, psol), i == psol ? 1 : 0, prstruct->ctx.ctx);
    
    flint::ca_number diag(prstruct->ctx);
    diag.call(ca_set, prstruct->J(psol, psol));

    for (int i = 1; i <= order; i++) {
        std::cerr << "DiffEqSolver: solving solution " << psol << " for order " << i << "\r";
        diag.call(ca_add_si, diag.ptr(), 1);

        flint::ca_matrix temp(prstruct->ctx, N, 1), tempY(prstruct->ctx, N, 1), prod(prstruct->ctx, N, 1);
        if (puppersol != -1)
            for (int j = 0; j < N; j++)
                ca_neg(temp(j, 0), solutions[i](j, puppersol), prstruct->ctx.ctx);

        for (int l = 1; l <= i; l++) {
            for (int j = 0; j < N; j++)
                ca_set(tempY(j, 0), solutions[i - l](j, psol), prstruct->ctx.ctx);
            prod.call(ca_mat_mul, coeff[l].ptr(), tempY.ptr());
            temp.call(ca_mat_add, temp.ptr(), prod.ptr());
        }

        flint::ca_matrix inv(prstruct->ctx, N, N), invmat(prstruct->ctx, N, N);
        inv.call(ca_mat_set_ca, diag.ptr());
        inv.call(ca_mat_sub, inv.ptr(), coeff[0].ptr());
        if (ca_mat_inv(invmat.ptr(), inv.ptr(), prstruct->ctx.ctx) != T_TRUE)
            throw GiNaC::pole_error("singular coefficient encountered!", 1);
        temp.call(ca_mat_mul, invmat.ptr(), temp.ptr());
        for (int j = 0; j < N; j++)
            ca_set(solutions[i](j, psol), temp(j, 0), prstruct->ctx.ctx);
    }
    std::cerr << "\nDiffEqSolver: solving solution " << psol << " done\n";
}


GiNaC::ex solver_regular::radius() {
    return rho;
}


GiNaC::matrix solver_regular::solution(unsigned dig) {
    int order = coeff.size() - 1, N = ca_mat_nrows(solutions[0].ptr());
    fill_in_solution(0, -1);
    for (int i = 1; i < N; i++) {
        if (ca_check_is_one(prstruct->J(i - 1, i), prstruct->ctx.ctx) == T_TRUE)
            fill_in_solution(i, i - 1);
        else
            fill_in_solution(i, -1);
    }

    GiNaC::matrix sol(N, N), finalsol(N, N);
    int blockstart = 0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            sol(i, j) = 0;
            for (int k = 0; k <= order; k++) {
                sol(i, j) += (fc2g(solutions[k](i, j), dig, prstruct->ctx.ctx) * GiNaC::pow(var, k));
            }
        }

        if (j == 0 || ca_check_is_one(prstruct->J(j - 1, j), prstruct->ctx.ctx) != T_TRUE) {
            for (int i = 0; i < N; i++)
                finalsol(i, j) = sol(i, j) * GiNaC::pow(var, fc2g(prstruct->J(j, j), dig, prstruct->ctx.ctx));
            blockstart = j;
        } else {
            for (int i = 0; i < N; i++) {
                finalsol(i, j) = sol(i, j);
                int k = j;
                while (k != blockstart) {
                    k--;
                    finalsol(i, j) += sol(i, k) * GiNaC::pow(GiNaC::log(var), j - k) / GiNaC::tgamma(j - k + 1);
                }
                finalsol(i, j) *= GiNaC::pow(var, fc2g(prstruct->J(j, j), dig, prstruct->ctx.ctx));
            }
        }
    }
    return finalsol;
}

