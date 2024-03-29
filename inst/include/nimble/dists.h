/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */

// Based on Writing R Extensions 6.6.1 and mgcv's mgcv.h
#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
/* If we are compiling with a version of R before FCONE and the explicit supplying of extra arguments 
   was introduced, then FCONE has to be defined */ 
#ifndef FCONE
#define FCONE
#endif


#ifndef __DISTS
#define __DISTS

#include "Utils.h"

using std::max;

bool R_IsNA_ANY(double*, int);
bool R_IsNaN_ANY(double*, int);
bool ISNAN_ANY(double*, int);
bool R_FINITE_ANY(double*, int);

extern "C" {
// BLAS/LAPACK routines
  // prototypes for dgemm_, dtrmv_, dtrsv_, dtrsm_, dtrmm_, dpotr_ not needed when #include "R_ext/lapack.h"
  // NIMBLE C wrappers called from R
  SEXP C_dmnorm_chol(SEXP, SEXP, SEXP, SEXP, SEXP); 
  SEXP C_rmnorm_chol(SEXP, SEXP, SEXP);
  SEXP C_dmvt_chol(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
  SEXP C_rmvt_chol(SEXP, SEXP, SEXP, SEXP); 
  SEXP C_dlkj_corr_cholesky(SEXP, SEXP, SEXP, SEXP); 
  SEXP C_rlkj_corr_cholesky(SEXP, SEXP); 
  SEXP C_dwish_chol(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_rwish_chol(SEXP, SEXP, SEXP);
  SEXP C_dinvwish_chol(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_rinvwish_chol(SEXP, SEXP, SEXP);
  SEXP C_dcat(SEXP, SEXP, SEXP);
  SEXP C_rcat(SEXP, SEXP);
  SEXP C_dmulti(SEXP, SEXP, SEXP, SEXP);
  SEXP C_rmulti(SEXP, SEXP);
  SEXP C_ddirch(SEXP, SEXP, SEXP);
  SEXP C_rdirch(SEXP);
  SEXP C_dinterval(SEXP, SEXP, SEXP, SEXP);
  SEXP C_rinterval(SEXP, SEXP, SEXP);
  SEXP C_dexp_nimble(SEXP, SEXP, SEXP);
  SEXP C_rexp_nimble(SEXP, SEXP);
  SEXP C_pexp_nimble(SEXP, SEXP, SEXP, SEXP);
  SEXP C_qexp_nimble(SEXP, SEXP, SEXP, SEXP);
  SEXP C_ddexp(SEXP, SEXP, SEXP, SEXP);
  SEXP C_rdexp(SEXP, SEXP, SEXP);
  SEXP C_pdexp(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_qdexp(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_dinvgamma(SEXP, SEXP, SEXP, SEXP);
  SEXP C_rinvgamma(SEXP, SEXP, SEXP);
  SEXP C_pinvgamma(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_qinvgamma(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_dsqrtinvgamma(SEXP, SEXP, SEXP, SEXP);
  SEXP C_rsqrtinvgamma(SEXP, SEXP, SEXP);
  SEXP C_dcar_normal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_dcar_proper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_rcar_proper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
}

// NOTE: R CMD SHLIB seems to handle C++ code without using wrapping the functions in 'extern "C"'; note that some of these functions have a bit of C++ syntax

// core scalar d/r functions provided by NIMBLE to extend R
double dcat(double, double*, int, int);
double rcat(double*, int);
double dmulti(double*, double, double*, int, int);
void rmulti(int*, double, double*, int);
double ddirch(double*, double*, int, int);
void rdirch(double*, double*, int);

double dmnorm_chol(double*, double*, double*, int, double, int, int);
void rmnorm_chol(double *, double*, double*, int, double);
double dmvt_chol(double*, double*, double*, double, int, double, int, int);
void rmvt_chol(double *, double*, double*, double, int, double);
double dwish_chol(double*, double*, double, int, double, int, int);
void rwish_chol(double*, double*, double, int, double, int);
double dinvwish_chol(double*, double*, double, int, double, int, int);
void rinvwish_chol(double*, double*, double, int, double, int);
double dlkj_corr_cholesky(double*, double, int, int);
void rlkj_corr_cholesky(double*, double, int);

double dinterval(double, double, double*, int, int);
double rinterval(double, double*, int);

double dcar_normal(double*, double*, double*, double*, double, int, int, int, int, int);
void rcar_normal(int, double*, double*, double*, double, int, int);

double dcar_proper(double*, double*, double*, double*, double*, double*, double, double, double*, int, int, int);
void rcar_proper(double*, double*, double*, double*, double*, double*, double, double, double*, int, int);


// SHOULD BE IN nimDists.h 
// Chris comment on above line: I don't think this is the case...

extern "C" {
  SEXP C_rt_nonstandard(SEXP, SEXP, SEXP, SEXP);
  SEXP C_dt_nonstandard(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_pt_nonstandard(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_qt_nonstandard(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
}

double dt_nonstandard(double, double, double, double, int);
double rt_nonstandard(double, double, double);
double pt_nonstandard(double, double, double, double, int, int);
double qt_nonstandard(double, double, double, double, int, int);

double dconstraint(double, double, int);
double rconstraint(double);

double rexp_nimble(double);
double dexp_nimble(double, double, int);
double pexp_nimble(double, double, int, int);
double qexp_nimble(double, double, int, int);

double rdexp(double, double);
double ddexp(double, double, double, int);
double pdexp(double, double, double, int, int);
double qdexp(double, double, double, int, int);

double rinvgamma(double, double);
double dinvgamma(double, double, double, int);
double pinvgamma(double, double, double, int, int);
double qinvgamma(double, double, double, int, int);

double rsqrtinvgamma(double, double);
double dsqrtinvgamma(double, double, double, int);

double rflat();
double dflat(double, int);
double rhalfflat();
double dhalfflat(double, int);
#endif
