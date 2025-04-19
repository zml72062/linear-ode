#include <flint/flint.h>
#include <flint/ca.h>
#include <ginac/ginac.h>


extern "C" {

truth_t ca_check_is_zero(const ca_t x, ca_ctx_t ctx) {
    ctx->options[CA_OPT_PRINT_FLAGS] = (CA_PRINT_N | (CA_PRINT_DIGITS * GiNaC::Digits));
    char* str = ca_get_str(x, ctx);
    std::string strcpp(str);
    flint_free(str);
    ctx->options[CA_OPT_PRINT_FLAGS] = CA_PRINT_DEFAULT;
    bool result = (GiNaC::abs(GiNaC::parser()(strcpp)) < GiNaC::pow(10, -((int)GiNaC::Digits / 2)));
    if (result)
        return T_TRUE;
    else
        return T_FALSE;
}

}

