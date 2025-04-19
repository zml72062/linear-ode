#ifndef INTERFACE_H
#define INTERFACE_H


/**
 * interface.h - Interfaces between GiNaC and FLINT.
 */

#include <ginac/ginac.h>
#include <flint/ca.h>
#include <flint/ca_poly.h>


/**
 * gq2f: GiNaC rational number to FLINT
 * 
 * @param gq a GiNaC rational number (can be complex)
 * @param f pointer to a FLINT `ca_struct`
 * @param ctx pointer to a FLINT `ca_ctx_struct`
 */
void gq2f(const GiNaC::ex& gq, ca_t f, ca_ctx_t ctx);


/**
 * gqp2f: GiNaC polynomial with rational coefficients to FLINT
 * 
 * @param gqp a GiNaC polynomial with rational coefficients (can be complex)
 * @param var variable in the polynomial
 * @param f pointer to a FLINT `ca_poly_struct`
 * @param ctx pointer to a FLINT `ca_ctx_struct`
 */
void gqp2f(const GiNaC::ex& gqp, const GiNaC::ex& var, ca_poly_t f, ca_ctx_t ctx);


/**
 * fc2g: FLINT complex number to GiNaC
 * 
 * @param fc a FLINT complex number
 * @param dig number of digits to keep
 * @param ctx pointer to a FLINT `ca_ctx_struct`
 * 
 * Notice that this method triggers numerical evaluation of the FLINT number,
 * which may be costly especially when `dig` is large, and potentially 
 * introduces side effects. Thus it is suggested that this method be only
 * called at the very end of computation.
 */
GiNaC::ex fc2g(const ca_t fc, unsigned dig, ca_ctx_t ctx);


/**
 * fc2gq: GiNaC rational approximation of FLINT complex number
 * 
 * @param fc a FLINT complex number
 * @param dig number of digits to keep
 * @param ctx pointer to a FLINT `ca_ctx_struct`
 * 
 * Notice that this method triggers numerical evaluation of the FLINT number,
 * which may be costly especially when `dig` is large, and potentially 
 * introduces side effects. Thus it is suggested that this method be only
 * called at the very end of computation.
 */
GiNaC::ex fc2gq(ca_ptr fc, unsigned dig, ca_ctx_t ctx);



#endif // INTERFACE_H

