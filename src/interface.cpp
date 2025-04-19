#include "interface.h"
#include <flint/fmpq.h>
#include <flint/qqbar.h>
#include <sstream>


void gq2f(const GiNaC::ex& gq, ca_t f, ca_ctx_t ctx) {
    if (GiNaC::ex_to<GiNaC::numeric>(gq).is_real()) {
        std::ostringstream ostr;
        ostr << gq.normal();
        std::size_t pkt, E;
        int power = 0;
        std::string ostrstr;
        ostrstr = ostr.str();
        fmpq_t fq;
        fmpq_init(fq);

        DETERMINE:
        if (((pkt = ostrstr.find('.')) == std::string::npos) &
            ((E = ostrstr.find('E')) == std::string::npos))
            fmpq_set_str(fq, ostrstr.c_str(), 10);
        else if (E == std::string::npos)
            fmpq_set_str(fq, (ostrstr.substr(0, pkt) + ostrstr.substr(pkt + 1) 
                 + "/1" + std::string(ostrstr.size() - pkt - 1, '0')).c_str(), 10);
        else {
            power = std::stoi(ostrstr.substr(E + 1));
            ostrstr = ostrstr.substr(0, E);
            goto DETERMINE;
        }

        ca_set_fmpq(f, fq, ctx);
        if (power != 0) {
            ca_t ten_power;
            ca_init(ten_power, ctx);
            gq2f(GiNaC::pow(10, power), ten_power, ctx);
            ca_mul(f, f, ten_power, ctx);
            ca_clear(ten_power, ctx);
        }
        fmpq_clear(fq);
    } else {
        auto real = gq.real_part(), imag = gq.imag_part();
        ca_t realq, imagq;
        qqbar_t qqr, qqi, qq;
        ca_init(realq, ctx);
        ca_init(imagq, ctx);
        gq2f(real, realq, ctx);
        gq2f(imag, imagq, ctx);
        qqbar_init(qqr);
        qqbar_init(qqi);
        qqbar_init(qq);
        ca_get_qqbar(qqr, realq, ctx);
        ca_get_qqbar(qqi, imagq, ctx);
        qqbar_set_re_im(qq, qqr, qqi);
        ca_set_qqbar(f, qq, ctx);
        qqbar_clear(qqr);
        qqbar_clear(qqi);
        qqbar_clear(qq);
        ca_clear(realq, ctx);
        ca_clear(imagq, ctx);
    }
}


void gqp2f(const GiNaC::ex& gqp, const GiNaC::ex& var, ca_poly_t f, ca_ctx_t ctx) {
    auto gqpex = gqp.expand();
    int d = gqpex.degree(var);
    for (int i = 0; i <= d; i++) {
        ca_t coeff;
        ca_init(coeff, ctx);
        gq2f(gqpex.coeff(var, i), coeff, ctx);
        ca_poly_set_coeff_ca(f, i, coeff, ctx);
        ca_clear(coeff, ctx);
    }
}


static GiNaC::ex fq2g(const ca_t fc, ca_ctx_t ctx) {
    ctx->options[CA_OPT_PRINT_FLAGS] = CA_PRINT_REPR;
    char* str = ca_get_str(fc, ctx);
    std::string strcpp(str);
    flint_free(str);
    ctx->options[CA_OPT_PRINT_FLAGS] = CA_PRINT_DEFAULT;

    try {
        return GiNaC::parser()(strcpp);
    } catch (GiNaC::parse_error& err) {
        fmpq_t fcq;
        fmpq_init(fcq);
        if (ca_get_fmpq(fcq, fc, ctx)) {
            char* str = fmpq_get_str(NULL, 10, fcq);
            std::string strcpp(str);
            flint_free(str);
            fmpq_clear(fcq);
            return GiNaC::parser()(strcpp);
        }
        fmpq_clear(fcq);
        throw GiNaC::parse_error("Cannot find a suitable representation");
    }
}


GiNaC::ex fc2g(const ca_t fc, unsigned dig, ca_ctx_t ctx) {
    try {
        if (ca_check_is_rational(fc, ctx) == T_TRUE)
            return fq2g(fc, ctx);
    } catch (GiNaC::parse_error& err) { }

    ca_t real, imag;
    ca_init(real, ctx);
    ca_init(imag, ctx);
    ca_re(real, fc, ctx);
    ca_im(imag, fc, ctx);
    
    try {
        if (ca_check_is_rational(real, ctx) == T_TRUE && ca_check_is_rational(imag, ctx) == T_TRUE) {
            GiNaC::ex result = fq2g(real, ctx) + GiNaC::I * fq2g(imag, ctx);
            ca_clear(real, ctx);
            ca_clear(imag, ctx);
            return result;
        }
    } catch (GiNaC::parse_error& err) { }

    ca_clear(real, ctx);
    ca_clear(imag, ctx);

    ctx->options[CA_OPT_PRINT_FLAGS] = (CA_PRINT_N | (CA_PRINT_DIGITS * dig));
    char* str = ca_get_str(fc, ctx);
    std::string strcpp(str);
    flint_free(str);
    ctx->options[CA_OPT_PRINT_FLAGS] = CA_PRINT_DEFAULT;
    return GiNaC::parser()(strcpp);
}


static GiNaC::ex to_fraction_pos(ca_ptr fvalue, unsigned digits, ca_ctx_t ctx) {
    if (ca_check_is_rational(fvalue, ctx) == T_TRUE)
        return fc2g(fvalue, digits, ctx);
    
    GiNaC::ex denom_bound = GiNaC::pow(10, (int)digits), 
              numer_prev = 1, denom_prev = 0, numer, denom = 1,
              n_value = fc2g(fvalue, 2 * digits, ctx);
    ca_t intp, residue;
    ca_init(intp, ctx);
    ca_init(residue, ctx);
    ca_floor(intp, fvalue, ctx);
    ca_sub(residue, fvalue, intp, ctx);
    numer = fc2g(intp, digits, ctx);

    while (GiNaC::pow(denom, 2) <= denom_bound && GiNaC::abs(n_value - numer / denom) >= 1 / denom_bound) {
        ca_inv(residue, residue, ctx);
        ca_floor(intp, residue, ctx);
        ca_sub(residue, residue, intp, ctx);

        GiNaC::ex new_term = fc2g(intp, digits, ctx),
                  numer_new = new_term * numer + numer_prev,
                  denom_new = new_term * denom + denom_prev;
        numer_prev = numer;
        denom_prev = denom;
        numer = numer_new;
        denom = denom_new;
    }

    ca_clear(intp, ctx);
    ca_clear(residue, ctx);
    
    return (numer / denom).normal();
}


GiNaC::ex fc2gq(ca_ptr fc, unsigned dig, ca_ctx_t ctx) {
    GiNaC::ex real_part, imag_part;
    ca_t real, imag, real_sgn, imag_sgn;
    ca_init(real, ctx);
    ca_init(imag, ctx);
    ca_init(real_sgn, ctx);
    ca_init(imag_sgn, ctx);

    ca_re(real, fc, ctx);
    ca_im(imag, fc, ctx);
    ca_sgn(real_sgn, real, ctx);
    ca_sgn(imag_sgn, imag, ctx);

    if (ca_check_is_neg_one(real_sgn, ctx) == T_TRUE) {
        ca_t neg_real;
        ca_init(neg_real, ctx);
        ca_neg(neg_real, real, ctx);
        real_part = -to_fraction_pos(neg_real, dig, ctx);
        ca_clear(neg_real, ctx);
    } else {
        real_part = to_fraction_pos(real, dig, ctx);
    }

    if (ca_check_is_neg_one(imag_sgn, ctx) == T_TRUE) {
        ca_t neg_imag;
        ca_init(neg_imag, ctx);
        ca_neg(neg_imag, imag, ctx);
        imag_part = -to_fraction_pos(neg_imag, dig, ctx);
        ca_clear(neg_imag, ctx);
    } else {
        imag_part = to_fraction_pos(imag, dig, ctx);
    }

    ca_clear(real, ctx);
    ca_clear(imag, ctx);
    ca_clear(real_sgn, ctx);
    ca_clear(imag_sgn, ctx);

    return real_part + GiNaC::I * imag_part;
}

