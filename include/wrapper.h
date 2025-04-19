#ifndef WRAPPER_H
#define WRAPPER_H


/**
 * wrapper.h - C++ wrappers for FLINT types.
 */


#include "interface.h"
#include <flint/ca_mat.h>
#include <flint/ca_vec.h>


#define NOCOPYMOVE(name)          \
public:                           \
    name(const name&) = delete;   \
    name(name&&)      = delete;   


namespace flint {

/**
 * RAII wrapper for `ca_ctx_t`.
 */
struct ca_ctx { NOCOPYMOVE(ca_ctx)
    ca_ctx()  { ca_ctx_init(ctx);  }
    ~ca_ctx() { ca_ctx_clear(ctx); }
    ca_ctx_t ctx;
};


/**
 * RAII wrapper for `ca_t`.
 */
class ca_number { NOCOPYMOVE(ca_number)
public:
    ca_number(): pctx(nullptr) { }
    ca_number(ca_ctx& ctx): pctx(&ctx) {
        ca_init(number, ctx.ctx);
    }
    ~ca_number() {
        if (pctx)
            ca_clear(number, pctx->ctx);
    }

    /**
     * Explicitly initialize the internal data.
     */
    void initialize(ca_ctx& ctx) {
        if (pctx)
            ca_clear(number, pctx->ctx);
        pctx = &ctx;
        ca_init(number, ctx.ctx);
    }
    
    /**
     * Import data from a GiNaC rational number.
     * 
     * @param gq a GiNaC rational number (can be complex)
     * 
     * Notice that this method should only be called after the `ca_number` 
     * object is fully initialized, i.e., either if the object was 
     * instantiated with the constructor taking a `ca_ctx&` argument, or 
     * if `initialize()` method has been explicitly called.
     */
    void import_from_ginac(const GiNaC::ex& gq) {
        gq2f(gq, number, pctx->ctx);
    }

    /**
     * Export data numerically to GiNaC.
     * 
     * @param dig number of digits to keep
     * @param rational whether to export a rational approximation (default 
     * `false`)
     * 
     * Notice that this method triggers numerical evaluation of the internal
     * FLINT number, which may be costly especially when `dig` is large, and 
     * potentially introduces side effects. Thus it is suggested that this 
     * method be only called at the very end of computation.
     */
    GiNaC::ex export_to_ginac(unsigned dig, bool rational = false) {
        if (rational)
            return fc2gq(number, dig, pctx->ctx);
        else
            return fc2g(number, dig, pctx->ctx);
    }

    /**
     * Inherit function calls on `ca_t` via template argument matching.
     */
    template <typename Func, typename... Args>
    void call(Func func, Args... args) {
        func(number, args..., pctx->ctx);
    }

    ca_struct* ptr() {
        return number;
    }
private:
    ca_t    number;
    ca_ctx* pctx;
};


/**
 * RAII wrapper for `ca_poly_t`.
 */
class ca_poly { NOCOPYMOVE(ca_poly)
public:
    ca_poly(): pctx(nullptr) { }
    ca_poly(ca_ctx& ctx): pctx(&ctx) {
        ca_poly_init(poly, ctx.ctx);
    }
    ~ca_poly() {
        if (pctx)
            ca_poly_clear(poly, pctx->ctx);
    }

    /**
     * Explicitly initialize the internal data.
     */
    void initialize(ca_ctx& ctx) {
        if (pctx)
            ca_poly_clear(poly, pctx->ctx);
        pctx = &ctx;
        ca_poly_init(poly, ctx.ctx);
    }
    
    /**
     * Import data from a GiNaC polynomial with rational coefficients.
     * 
     * @param gpoly a GiNaC polynomial with rational coefficients (can be complex)
     * @param var variable in the polynomial
     * 
     * Notice that this method should only be called after the `ca_poly` object
     * is fully initialized, i.e., either if the object was instantiated with 
     * the constructor taking a `ca_ctx&` argument, or if `initialize()` method
     * has been explicitly called.
     */
    void import_from_ginac(const GiNaC::ex& gpoly, const GiNaC::ex& var) {
        gqp2f(gpoly, var, poly, pctx->ctx);
    }

    ca_ptr coeff(int i) {
        return ca_poly_coeff_ptr(poly, i);
    }

    /**
     * Inherit function calls on `ca_poly_t` via template argument matching.
     */
    template <typename Func, typename... Args>
    void call(Func func, Args... args) {
        func(poly, args..., pctx->ctx);
    }

    ca_poly_struct* ptr() {
        return poly;
    }

    friend class ca_rational;
private:
    ca_poly_t   poly;
    ca_ctx*     pctx;
};


/**
 * Rational functions are treated as the division of two polynomials.
 */
class ca_rational { NOCOPYMOVE(ca_rational)
public:
    ca_rational()  { }
    ca_rational(ca_ctx& ctx): numer(ctx), denom(ctx) { }
    ~ca_rational() { }

    /**
     * Explicitly initialize the internal data.
     */
    void initialize(ca_ctx& ctx) {
        numer.initialize(ctx);
        denom.initialize(ctx);
    }

    /**
     * Import data from a GiNaC rational function.
     * 
     * @param grat a GiNaC rational function (can be complex)
     * @param var variable in the rational function
     * 
     * Notice that this method should only be called after the `ca_rational` 
     * object is fully initialized, i.e., either if the object was instantiated 
     * with the constructor taking a `ca_ctx&` argument, or if `initialize()` 
     * method has been explicitly called.
     */
    void import_from_ginac(const GiNaC::ex& grat, const GiNaC::ex& var) {
        auto gnumer_denom = grat.normal().numer_denom(),
             gnumer = gnumer_denom[0],
             gdenom = gnumer_denom[1];
        
        numer.import_from_ginac(gnumer, var);
        denom.import_from_ginac(gdenom, var);
    }

    /**
     * Give series expansion around origin up to a designated order.
     * 
     * @param result mutable reference to a fully initialized `ca_poly`
     * object, where the result is stored
     * @param order the desired expansion order
     * 
     * Notice that this method assumes that the origin is not a pole of the
     * rational function. Otherwise, the result may contain infinite or 
     * undefined coefficients.
     */
    void series(ca_poly& result, unsigned order) {
        if (result.pctx != numer.pctx || result.pctx != denom.pctx) {
            std::cerr << "Warning: result polynomial is in a different context "
                         "from that of the rational function. Result would not "
                         "be written." << std::endl;
            return;
        }
        ca_poly_div_series(result.poly, numer.poly, denom.poly, order + 1, result.pctx->ctx);
    }

    /**
     * Find the pole closest to origin and return its distance to origin. 
     * This value marks the radius of convergence of the Taylor expansion
     * around origin.
     * 
     * @param dig number of digits to keep
     */
    GiNaC::ex radius_convergence(unsigned dig) {
        int length = denom.ptr()->length;
        if (length == 0) {
            std::cerr << "Warning: empty denominator" << std::endl;
            return 0;
        }
        if (length == 1) {
            if (ca_check_is_zero(denom.coeff(0), denom.pctx->ctx) == T_TRUE) {
                std::cerr << "Warning: denominator has a pole at origin" << std::endl;
                return 0;
            } else {
                return GiNaC::pow(10, int(GiNaC::Digits));
            }
        }
        int degree = length - 1;
        ca_vec_t roots;
        std::vector<ulong> mults(degree);
        ca_vec_init(roots, degree, denom.pctx->ctx);
        if (ca_poly_roots(roots, mults.data(), denom.ptr(), denom.pctx->ctx) == 0) {
            std::cerr << "Warning: fail to find poles" << std::endl;
            ca_vec_clear(roots, denom.pctx->ctx);
            return 0;
        }
        GiNaC::ex closest = GiNaC::pow(10, int(GiNaC::Digits)), dist;
        int nroots = ca_vec_length(roots, denom.pctx->ctx);
        for (int i = 0; i < nroots; i++) {
            if ((dist = GiNaC::abs(fc2g(ca_vec_entry(roots, i), dig, denom.pctx->ctx))) < closest)
                closest = dist;
        }
        ca_vec_clear(roots, denom.pctx->ctx);

        if (closest < GiNaC::pow(10, -6))
            std::cerr << "Warning: denominator has a pole at origin" << std::endl;
        return closest;
    }
private:
    ca_poly numer;
    ca_poly denom;
};


/**
 * RAII wrapper for `ca_vec_t`.
 */
class ca_vec { NOCOPYMOVE(ca_vec)
public:
    ca_vec(): pctx(nullptr) { }
    ca_vec(ca_ctx& ctx, int n): pctx(&ctx) {
        ca_vec_init(vector, n, ctx.ctx);
    }
    ~ca_vec() {
        if (pctx)
            ca_vec_clear(vector, pctx->ctx);
    }

    /**
     * Explicitly initialize the internal data.
     */
    void initialize(ca_ctx& ctx, int n) {
        if (pctx)
            ca_vec_clear(vector, pctx->ctx);
        pctx = &ctx;
        ca_vec_init(vector, n, ctx.ctx);
    }

    ca_ptr operator[](int i) {
        return ca_vec_entry_ptr(vector, i);
    }
    
    /**
     * Inherit function calls on `ca_vec_t` via template argument matching.
     */
    template <typename Func, typename... Args>
    void call(Func func, Args... args) {
        func(vector, args..., pctx->ctx);
    }

    ca_vec_struct* ptr() {
        return vector;
    }
private:
    ca_vec_t vector;
    ca_ctx*  pctx;
};


/**
 * RAII wrapper for `ca_mat_t`.
 */
class ca_matrix { NOCOPYMOVE(ca_matrix)
public:
    ca_matrix(): pctx(nullptr) { }
    ca_matrix(ca_ctx& ctx, int r = 1, int c = 1): pctx(&ctx) {
        ca_mat_init(matrix, r, c, ctx.ctx);
    }
    ~ca_matrix() {
        if (pctx)
            ca_mat_clear(matrix, pctx->ctx);
    }

    /**
     * Explicitly initialize the internal data.
     */
    void initialize(ca_ctx& ctx, int r = 1, int c = 1) {
        if (pctx)
            ca_mat_clear(matrix, pctx->ctx);
        pctx = &ctx;
        ca_mat_init(matrix, r, c, ctx.ctx);
    }
    
    /**
     * Import data from a GiNaC matrix of rational numbers.
     * 
     * @param gmat a GiNaC matrix whose elements are all rational numbers
     * (can be complex)
     * 
     * Notice that this method should only be called after the `ca_matrix` 
     * object is fully initialized, i.e., either if the object was 
     * instantiated with the constructor taking a `ca_ctx&` argument, or 
     * if `initialize()` method has been explicitly called.
     */
    void import_from_ginac(const GiNaC::matrix& gmat) {
        int r = gmat.rows(), c = gmat.cols();
        initialize(*pctx, r, c);

        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                gq2f(gmat(i, j), ca_mat_entry(matrix, i, j), pctx->ctx);
    }

    /**
     * Export data numerically to GiNaC.
     * 
     * @param dig number of digits to keep
     * @param rational whether to export a rational approximation (default 
     * `false`)
     * 
     * Notice that this method triggers numerical evaluation of the internal
     * FLINT matrix, which may be costly especially when `dig` is large, and 
     * potentially introduces side effects. Thus it is suggested that this 
     * method be only called at the very end of computation.
     */
    GiNaC::matrix export_to_ginac(unsigned dig, bool rational = false) {
        int r = ca_mat_nrows(matrix), c = ca_mat_ncols(matrix);
        GiNaC::matrix result(r, c);

        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                if (rational)
                    result(i, j) = fc2gq(ca_mat_entry(matrix, i, j), dig, pctx->ctx);
                else
                    result(i, j) = fc2g(ca_mat_entry(matrix, i, j), dig, pctx->ctx);
        
        return result;
    }

    ca_ptr operator()(int i, int j) {
        return ca_mat_entry_ptr(matrix, i, j);
    }

    /**
     * Inherit function calls on `ca_mat_t` via template argument matching.
     */
    template <typename Func, typename... Args>
    void call(Func func, Args... args) {
        func(matrix, args..., pctx->ctx);
    }

    ca_mat_struct* ptr() {
        return matrix;
    }
private:
    ca_mat_t matrix;
    ca_ctx*  pctx;
};


} // namespace flint


#endif // WRAPPER_H

