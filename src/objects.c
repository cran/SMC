/*
 *  objects.c
 *  SMC
 *
 *  Created by Gopi Goswami on Tue May 22 2006.
 *  Copyright (C) 2006 Gopika R. Goswami
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  For a copy of the GNU General Public License please write to the
 *  Free Software Foundation, Inc.
 *  59 Temple Place, Suite 330.
 *  Boston, MA  02111-1307 USA.
 *
 *  For bugs in the code please contact:
 *  goswami@stat.harvard.edu
 *
 *
 *  SYNOPSIS
 *
 *
 *
 *  DESCRIPTION
 *
 * 
 *
 */

#include <assert.h>
#include <R.h>           
#include <Rinternals.h>
#include <Rmath.h>
#include "utils.h"
#include "objects.h"

#ifndef DEBUG_OBJECTS
#define DEBUG_OBJECTS 0
#endif

#if (!DEBUG_OBJECTS)
#define DEBUG(_x) ((void) 0)
#endif

/* The static functions of this file */
static SEXP
getListElement (SEXP list, char *str);

static SEXP
printList (SEXP list);

static int
gather_time_details (Sampler *ss);

static int
args_list1_init (ArgsList1 *al);

static int
args_list2_init (ArgsList2 *al);

static int
args_list3_init (ArgsList3 *al);

static int
sampler_register_streams_N_log_weights (Sampler *ss, SEXP currentStreams,
                                        SEXP currentLogWeights);

static int
sampler_register_resamp (Sampler *ss, SEXP aList);

static int
sampler_register_acceptance_rates (Sampler *ss, SEXP aList);

static double
sampler_adjust_log_weights (int nn, double *lw, double *aw);

static int
sampler_register_summary (Sampler *ss, SEXP summ);

static SEXP
propagate_func_user_Rfunc (Sampler *ss, int currentPeriod, int nStreamsToGenerate,
                           SEXP lag1Streams, SEXP lag1LogWeights);

static Rboolean
resampCriterion_func_user_Rfunc (Sampler *ss, int currentPeriod, SEXP currentStreams,
                                 SEXP currentLogWeights);

static Rboolean
resampCriterion_func_builtin (Sampler *ss, int currentPeriod, SEXP currentStreams,
                              SEXP currentLogWeights);

static SEXP
resamp_func_user_Rfunc (Sampler *ss, int currentPeriod, SEXP currentStreams,
                        SEXP currentLogWeights);

static SEXP
resamp_func_builtin_PPW (Sampler *ss, int currentPeriod, SEXP currentStreams,
                         SEXP currentLogWeights);

static SEXP
summary_func_user_Rfunc (Sampler *ss, int currentPeriod, SEXP currentStreams,
                         SEXP currentLogWeights);

static SEXP
summary_func_builtin (Sampler *ss, int currentPeriod, SEXP currentStreams,
                      SEXP currentLogWeights);

static SEXP
MHUpdate_func_user_Rfunc (Sampler *ss, int currentPeriod, int nMHSteps,
                          SEXP currentStreams, SEXP lag1Streams,
                          SEXP lag1LogWeights);

static SEXP
one_iter_without_MH (Sampler *ss, SEXP notRequired);

static SEXP
one_iter_with_MH (Sampler *ss, SEXP notRequired);

static SEXP
register_summary (Sampler *ss, SEXP draws);

static SEXP
register_propUniqueStreamIds (Sampler *ss, SEXP draws);

static SEXP
register_streams (Sampler *ss, SEXP draws);

static SEXP
register_log_weights (Sampler *ss, SEXP draws);

static SEXP
register_acceptance_rates (Sampler *ss, SEXP draws);

static SEXP
register_all (Sampler *ss, SEXP draws);

static SEXP
make_draws (Sampler *ss);


/*
 * The following has been borrowed from src/main/deriv.c from the R
 * distribution.
 */
/* static SEXP */
/* lang5(SEXP s, SEXP t, SEXP u, SEXP v, SEXP w) */
/* { */
/*     PROTECT(s); */
/*     s = LCONS(s, list4(t, u, v, w)); */
/*     UNPROTECT(1); */
/*     return s; */
/* } */


/*
 * The following has been taken from the "Manual" section of R
 * website, it gets the list element named str, or returns NULL, if
 * not found. It has been modified a little bit.
 */
static SEXP
getListElement (SEXP list, char *str)
{
        SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
        int ii, found = 0, nn = length(list);
     
        for (ii = 0; ii < nn; ii++) {
                if (strcmp(CHAR(STRING_ELT(names, ii)), str) == 0) {
                        elmt = VECTOR_ELT(list, ii); ++found; break;
                }
        }
        if (found == 0) {
                char *errMsg1, errMsg2[MAX_LINE_LENGTH];

                errMsg1 = (char *) R_alloc(MAX_LINE_LENGTH, sizeof(errMsg1));
                for (ii = 0; ii < nn; ii++) {
                        errMsg1 = strcat(errMsg1, CHAR(STRING_ELT(names, ii)));
                        errMsg1 = strcat(errMsg1,
                                         ((ii == (nn - 1)) ? "" : ", "));
                }
                        
                sprintf(errMsg2,
                        "No element called \"%s\" found in the SEXP list, " \
                        "check your C code; given elements are:\n%s\n",
                        str, errMsg1);
                error(errMsg2);
        }
        return elmt;
}


static SEXP
printList (SEXP list)
{
        SEXP names = getAttrib(list, R_NamesSymbol);
        int ii, nn = length(list);

        Rprintf("The given list has the following components:\n");
        for (ii = 0; ii < nn; ++ii)
                Rprintf("%s\n", CHAR(STRING_ELT(names, ii)));
        return 0;
}


static int
gather_time_details (Sampler *ss)
{
        SEXP SEXPTmp;
        double *doublesTmp;
        
        PROTECT(SEXPTmp = eval(ss->procTimeFuncCall, ss->procTimeFuncEnv));
        doublesTmp           = REAL(SEXPTmp);
        ss->timeDetails->usr = doublesTmp[0];
        ss->timeDetails->sys = doublesTmp[1];
        UNPROTECT(1);
        return 0;
}


/*******************************************************************************/


/*
 * The following args_list*_init functions use PROTECT statements
 * during initialization. They return the nProtected to be kept track
 * of by the its caller.
 */
static int
args_list1_init (ArgsList1 *al)
{
        int nComps = 0, comp = 0, nProtected = 0, nProtectedCaller = 0;
        SEXP names;

        al->posCurrentPeriod      = nComps++;
        al->posNStreamsToGenerate = nComps++;
        al->posLag1Streams        = nComps++;
        al->posLag1LogWeights     = nComps++;

        PROTECT(al->argsList = allocVector(VECSXP, nComps)); ++nProtectedCaller;
        PROTECT(names        = allocVector(STRSXP, nComps)); ++nProtected;

        for (comp = 0; comp < nComps; ++comp)
                SET_VECTOR_ELT(al->argsList, comp, R_NilValue);
        comp = 0;
        SET_STRING_ELT(names, comp, mkChar("currentPeriod")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("nStreamsToGenerate")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("lag1Streams")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("lag1LogWeights")); ++comp;
        
        setAttrib(al->argsList, R_NamesSymbol, names);
        UNPROTECT(nProtected);
        return nProtectedCaller;
}


static int
args_list2_init (ArgsList2 *al)
{
        int nComps = 0, comp = 0, nProtected = 0, nProtectedCaller = 0;
        SEXP names;

        al->posCurrentPeriod     = nComps++;
        al->posCurrentStreams    = nComps++;
        al->posCurrentLogWeights = nComps++;

        PROTECT(al->argsList = allocVector(VECSXP, nComps)); ++nProtectedCaller;
        PROTECT(names        = allocVector(STRSXP, nComps)); ++nProtected;

        for (comp = 0; comp < nComps; ++comp)
                SET_VECTOR_ELT(al->argsList, comp, R_NilValue);
        comp = 0;
        SET_STRING_ELT(names, comp, mkChar("currentPeriod")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("currentStreams")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("currentLogWeights")); ++comp;
        
        setAttrib(al->argsList, R_NamesSymbol, names);
        UNPROTECT(nProtected);
        return nProtectedCaller;
}


static int
args_list3_init (ArgsList3 *al)
{
        int nComps = 0, comp = 0, nProtected = 0, nProtectedCaller = 0;
        SEXP names;

        al->posCurrentPeriod  = nComps++;
        al->posNMHSteps       = nComps++;
        al->posCurrentStreams = nComps++;
        al->posLag1Streams    = nComps++;
        al->posLag1LogWeights = nComps++;

        PROTECT(al->argsList = allocVector(VECSXP, nComps)); ++nProtectedCaller;
        PROTECT(names        = allocVector(STRSXP, nComps)); ++nProtected;

        for (comp = 0; comp < nComps; ++comp)
                SET_VECTOR_ELT(al->argsList, comp, R_NilValue);
        comp = 0;
        SET_STRING_ELT(names, comp, mkChar("currentPeriod")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("nMHSteps")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("currentStreams")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("lag1Streams")); ++comp;
        SET_STRING_ELT(names, comp, mkChar("lag1LogWeights")); ++comp;
        
        setAttrib(al->argsList, R_NamesSymbol, names);
        UNPROTECT(nProtected);
        return nProtectedCaller;
}


Sampler *
sampler_new (SEXP opts)
{
        Sampler *ss;
        SEXP SEXPTmp;
        
        ss                    = (Sampler *) R_alloc(1, sizeof(struct Sampler));
        ss->nStreams          = INTEGER(getListElement(opts, "nStreams"))[0];
        ss->nPeriods          = INTEGER(getListElement(opts, "nPeriods"))[0];
        ss->nStreamsPreResamp = INTEGER(getListElement(opts, "nStreamsPreResamp"))[0];
        ss->dimPerPeriod      = INTEGER(getListElement(opts, "dimPerPeriod"))[0];
        ss->dimSummPerPeriod  = INTEGER(getListElement(opts, "dimSummPerPeriod"))[0];
        ss->returnStreams     = LOGICAL(getListElement(opts, "returnStreams"))[0];
        ss->returnLogWeights  = LOGICAL(getListElement(opts, "returnLogWeights"))[0];
        ss->nMHSteps          = INTEGER(getListElement(opts, "nMHSteps"))[0];
        ss->verboseLevel      = INTEGER(getListElement(opts, "verboseLevel"))[0];

        ss->printEstTimeAt = 10; ss->printEstTimeNTimes = 10;
        /* FIXME: The setting for ss->printDotAt, is it all right? */
        ss->printInitialDotsWhen = ss->printEstTimeAt / 10;
        ss->printDotAt = 0; ss->nDotsPerLine = 20;
        ss->eachDotWorth = (int) ceil((ss->nPeriods - ss->printEstTimeAt + 1.0) / \
                                      (ss->printEstTimeNTimes * ss->nDotsPerLine));
        ss->nProtected = 0;

        /* The user provided functions */
        ss->propagateFunc     = getListElement(opts, "propagateFunc");
        ss->propagateArgsList = (ArgsList1 *) R_alloc(1, sizeof(struct ArgsList1));
        ss->nProtected       += args_list1_init(ss->propagateArgsList);        

        ss->resampCriterionFunc     = getListElement(opts, "resampCriterionFunc");
        ss->resampCriterionArgsList = (ArgsList2 *) R_alloc(1, sizeof(struct ArgsList2));
        ss->nProtected             += args_list2_init(ss->resampCriterionArgsList);

        ss->resampFunc     = getListElement(opts, "resampFunc");
        ss->resampArgsList = (ArgsList2 *) R_alloc(1, sizeof(struct ArgsList2));
        ss->nProtected    += args_list2_init(ss->resampArgsList);

        ss->summaryFunc     = getListElement(opts, "summaryFunc");
        ss->summaryArgsList = (ArgsList2 *) R_alloc(1, sizeof(struct ArgsList2));
        ss->nProtected     += args_list2_init(ss->summaryArgsList);

        ss->MHUpdateFunc     = getListElement(opts, "MHUpdateFunc");
        ss->MHUpdateArgsList = (ArgsList3 *) R_alloc(1, sizeof(struct ArgsList3));
        ss->nProtected      += args_list3_init(ss->MHUpdateArgsList);
        
        SEXPTmp = getListElement(opts, "doCallFunc");
        PROTECT(ss->doCallFuncCall = lang4(SEXPTmp, R_NilValue,
                                           R_NilValue, R_NilValue));
        ++(ss->nProtected);
        ss->doCallFuncEnv = getListElement(opts, "doCallFuncEnv");
        
        SEXPTmp = getListElement(opts, "procTimeFunc");
        PROTECT(ss->procTimeFuncCall = lang1(SEXPTmp));
        ++(ss->nProtected);
        ss->procTimeFuncEnv = getListElement(opts, "procTimeFuncEnv");

        ss->timeDetails = (TimeDetails *) R_alloc(1, sizeof(struct TimeDetails));

        PROTECT(ss->SEXPCurrentPeriod = allocVector(INTSXP, 1)); ++(ss->nProtected);
        PROTECT(ss->SEXPNStreamsToGenerate = allocVector(INTSXP, 1)); ++(ss->nProtected);
        PROTECT(ss->SEXPNMHSteps = allocVector(INTSXP, 1)); ++(ss->nProtected);

        ss->dotsList = getListElement(opts, "dotsList");

        ss->SEXPCurrentStreams      = R_NilValue;
        PROTECT(ss->SEXPLag1Streams = allocMatrix(REALSXP, ss->nStreams, ss->dimPerPeriod));
        ++(ss->nProtected);

        ss->SEXPCurrentLogWeights           = R_NilValue;
        PROTECT(ss->SEXPCurrentAdjWeights   = allocVector(REALSXP, ss->nStreamsPreResamp));
        ++(ss->nProtected);
        PROTECT(ss->SEXPLag1LogWeights      = allocVector(REALSXP, ss->nStreams));
        ++(ss->nProtected);
        PROTECT(ss->SEXPLag1AdjWeights      = allocVector(REALSXP, ss->nStreams));
        ++(ss->nProtected);
        PROTECT(ss->SEXPAcceptanceRates     = allocVector(REALSXP, ss->nStreams));
        ++(ss->nProtected);
        PROTECT(ss->SEXPSummary             = allocVector(REALSXP, ss->dimSummPerPeriod));
        ++(ss->nProtected);
        PROTECT(ss->SEXPPropUniqueStreamIds = allocVector(REALSXP, 1));
        ++(ss->nProtected);
        
        ss->scratch_RC                  = (ResampleContext *) R_alloc(1, sizeof(struct ResampleContext));
        ss->scratch_RC->streamIds       = (int *) R_alloc(ss->nStreams, sizeof(int));
        ss->scratch_RC->uniqueStreamIds = (int *) R_alloc(ss->nStreams, sizeof(int));
        ss->scratch_RC->partialSum      = (double *) R_alloc(ss->nStreamsPreResamp, sizeof(double));
        return ss;
}


int
sampler_print (Sampler *ss)
{
        PRINT_STUB_INT(ss->nStreams);
        PRINT_STUB_INT(ss->nPeriods);
        PRINT_STUB_INT(ss->dimPerPeriod);
        PRINT_STUB_INT(ss->dimSummPerPeriod);
        PRINT_STUB_INT(ss->verboseLevel);
        PRINT_STUB_INT(ss->printDotAt);
        PRINT_STUB_INT(ss->nDotsPerLine);
        PRINT_STUB_INT(ss->eachDotWorth);
        PRINT_STUB_INT(ss->returnStreams);
        PRINT_STUB_INT(ss->returnLogWeights);
        return 0;
}


int
sampler_register_streams_N_log_weights (Sampler *ss, SEXP currentStreams,
                                        SEXP currentLogWeights)
{
        int ii, jj, count = 0;
        int ns = ss->nStreams, dPP = ss->dimPerPeriod;
        double *destS  = REAL(ss->SEXPLag1Streams);
        double *srcS   = REAL(currentStreams);
        double *destLW = REAL(ss->SEXPLag1LogWeights);
        double *srcLW  = REAL(currentLogWeights);

        /*
         * Note: The matrices, srcS and destS are stored in
         * (column-wise) vec form internally, and we are copying them
         * by iterating row-wise. But, that's fine, since count \in \{
         * 0, 1, ..., ns \times dpp \}, which achives the goal of
         * copying all the elements.
         */
        for (ii = 0; ii < ns; ++ii) {
                for (jj = 0; jj < dPP; ++jj) {
                        destS[count] = srcS[count]; ++count; 
                }
                destLW[ii] = srcLW[ii];
        }
        return 0;
}


static int
sampler_register_resamp (Sampler *ss, SEXP aList)
{
        if (ss->resampDone == TRUE)
                REAL(ss->SEXPPropUniqueStreamIds)[0] = REAL(getListElement(aList, "propUniqueStreamIds"))[0];
        else
                REAL(ss->SEXPPropUniqueStreamIds)[0] = NA_REAL;
        return 0;
}


static int
sampler_register_acceptance_rates (Sampler *ss, SEXP aList)
{
        int ii;
        int ns = ss->nStreams;
        double *dest = REAL(ss->SEXPAcceptanceRates);

        if (ss->resampDone == TRUE) {
                double *src  = REAL(getListElement(aList, "acceptanceRates"));
                
                for (ii = 0; ii < ns; ++ii) dest[ii] = src[ii];
        }
        else {
                for (ii = 0; ii < ns; ++ii) dest[ii] = NA_REAL;
        }
        return 0;
}


static double
sampler_adjust_log_weights (int nn, double *lw, double *aw)
{
        int jj;
        double mlw = R_NegInf, sumaw = 0.0;

        for (jj = 0; jj < nn; ++jj) mlw = MAX(mlw, lw[jj]);
        for (jj = 0; jj < nn; ++jj) {
                aw[jj] = exp(lw[jj] - mlw);
                sumaw += aw[jj];
        }
        return sumaw;
}


static int
sampler_register_summary (Sampler *ss, SEXP summ)
{
        int ii, dspp = ss->dimSummPerPeriod;
        double *dest = REAL(ss->SEXPSummary), *src = REAL(summ);

        for (ii = 0; ii < dspp; ++ii) dest[ii] = src[ii];
        return 0;
}


static SEXP
propagate_func_user_Rfunc (Sampler *ss, int currentPeriod, int nStreamsToGenerate,
                           SEXP lag1Streams, SEXP lag1LogWeights)
{
        ArgsList1 *al = ss->propagateArgsList;
        
        INTEGER(ss->SEXPCurrentPeriod)[0]      = currentPeriod + 1;
        INTEGER(ss->SEXPNStreamsToGenerate)[0] = nStreamsToGenerate;

        SET_VECTOR_ELT(al->argsList, al->posCurrentPeriod, ss->SEXPCurrentPeriod);
        SET_VECTOR_ELT(al->argsList, al->posNStreamsToGenerate, ss->SEXPNStreamsToGenerate);
        SET_VECTOR_ELT(al->argsList, al->posLag1Streams, lag1Streams);
        SET_VECTOR_ELT(al->argsList, al->posLag1LogWeights, lag1LogWeights);

        SETCADR(ss->doCallFuncCall, ss->propagateFunc);
        SETCADDR(ss->doCallFuncCall, al->argsList);
        SETCADDDR(ss->doCallFuncCall, ss->dotsList);
        return eval(ss->doCallFuncCall, ss->doCallFuncEnv);
}


static Rboolean
resampCriterion_func_user_Rfunc (Sampler *ss, int currentPeriod, SEXP currentStreams,
                                 SEXP currentLogWeights)
{        
        ArgsList2 *al = ss->resampCriterionArgsList;
        SEXP SEXPTmp;
        Rboolean res;
        
        INTEGER(ss->SEXPCurrentPeriod)[0] = currentPeriod + 1;

        SET_VECTOR_ELT(al->argsList, al->posCurrentPeriod, ss->SEXPCurrentPeriod);
        SET_VECTOR_ELT(al->argsList, al->posCurrentStreams, currentStreams);
        SET_VECTOR_ELT(al->argsList, al->posCurrentLogWeights, currentLogWeights);

        SETCADR(ss->doCallFuncCall, ss->resampCriterionFunc);
        SETCADDR(ss->doCallFuncCall, al->argsList);
        SETCADDDR(ss->doCallFuncCall, ss->dotsList);         
        PROTECT(SEXPTmp = eval(ss->doCallFuncCall, ss->doCallFuncEnv));
        res = LOGICAL(SEXPTmp)[0];
        UNPROTECT(1);
        return res;
}


static Rboolean
resampCriterion_func_builtin (Sampler *ss, int currentPeriod, SEXP currentStreams,
                              SEXP currentLogWeights)
{
        int ii, nspr = ss->nStreamsPreResamp;
        double *sclw = REAL(currentLogWeights);
        double *scaw = REAL(ss->SEXPCurrentAdjWeights);
        double weight, sumWeights = 0.0, sumWeightsSq = 0.0;
        double meanWeights, meanWeightsSq, CVWeightsSq;

        /*
         * Compute the coefficient of variation of the weights and
         * hence the effective sample size.
         */
        sampler_adjust_log_weights(nspr, sclw, scaw);
        for (ii = 0; ii < nspr; ++ii) {
                weight        = scaw[ii];
                sumWeights   += weight;
                sumWeightsSq += SQR(weight);
        }
        meanWeights   = sumWeights / nspr;
        meanWeightsSq = SQR(meanWeights);
        CVWeightsSq   = (sumWeightsSq / nspr - meanWeightsSq) / meanWeightsSq;
        PHONY(PRINT_STUB_DOUBLE(CVWeightsSq));
        /*
         * The builtin resampCriterion: resample as and when the
         * CVWeightsSq >= 1.0.
         */
        return ((CVWeightsSq >= 1.0) ? (TRUE) : (FALSE));
}


static SEXP
resamp_func_user_Rfunc (Sampler *ss, int currentPeriod, SEXP currentStreams,
                        SEXP currentLogWeights)
{        
        ArgsList2 *al = ss->resampArgsList;
        
        INTEGER(ss->SEXPCurrentPeriod)[0] = currentPeriod + 1;

        SET_VECTOR_ELT(al->argsList, al->posCurrentPeriod, ss->SEXPCurrentPeriod);
        SET_VECTOR_ELT(al->argsList, al->posCurrentStreams, currentStreams);
        SET_VECTOR_ELT(al->argsList, al->posCurrentLogWeights, currentLogWeights);

        SETCADR(ss->doCallFuncCall, ss->resampFunc);
        SETCADDR(ss->doCallFuncCall, al->argsList);
        SETCADDDR(ss->doCallFuncCall, ss->dotsList);         
        return eval(ss->doCallFuncCall, ss->doCallFuncEnv);
}


/*
 * The following returns a R list with the following components:
 * currentStreams
 * currentLogWeights
 * propUniqueStreamIds
 */
static SEXP
resamp_func_builtin_PPW (Sampler *ss, int currentPeriod, SEXP currentStreams,
                         SEXP currentLogWeights)
{
        ResampleContext *rc = ss->scratch_RC;
        int nspr = ss->nStreamsPreResamp, dpp = ss->dimPerPeriod;
        int ns = ss->nStreams, *sids = rc->streamIds, ii, jj, kk;
        int nusids, *usids = rc->uniqueStreamIds;
        int nComps = 0, nProtected = 0;
        double *ps = rc->partialSum;
        double sum, uu;
        SEXP resampCurrentStreams, resampCurrentLogWeights, resampPropUniqueStreamIds;
        SEXP retList, names;
        double *rcs, *rclw;
        double *scs  = REAL(currentStreams);
        double *sclw = REAL(currentLogWeights);
        double *scaw = REAL(ss->SEXPCurrentAdjWeights);
        void *vmax = vmaxget( );

        PROTECT(resampCurrentStreams    = allocMatrix(REALSXP, ns, dpp));
        ++nComps; ++nProtected;
        PROTECT(resampCurrentLogWeights = allocVector(REALSXP, ns));
        ++nComps; ++nProtected;
        rcs  = REAL(resampCurrentStreams);
        rclw = REAL(resampCurrentLogWeights);

        sampler_adjust_log_weights(nspr, sclw, scaw);
        ps[0] = scaw[0];
        for (jj = 1; jj < nspr; ++jj) {
                ps[jj] = ps[jj - 1] + scaw[jj];
        }
        sum = ps[nspr - 1]; nusids = 0;
        /* resample the streams with probability proportional to their
         * weights */
        for (jj = 0; jj < ns; ++jj) {
                uu = runif(0, sum);
                for (ii = 0; ii < nspr; ++ii) {
                        if (uu <= ps[ii]) { sids[jj] = ii; break; }
                }
                /* copying the resampled stream */
                for (kk = 0; kk < dpp; ++kk)
                        rcs[kk * ns + jj] = scs[kk * nspr + sids[jj]];
                /* making the resampled logWeights = 0 */
                rclw[jj] = 0;
                /* find the unique stream and register it */
                if (utils_is_int_in_iarray(sids[jj], nusids, usids) == FALSE) {
                        usids[nusids] = sids[jj]; ++nusids;
                }
        }
        rc->nUniqueStreamIds    = nusids;
        rc->propUniqueStreamIds = nusids / ((double) nspr);
        PROTECT(resampPropUniqueStreamIds = allocVector(REALSXP, 1));
        ++nComps; ++nProtected;
        REAL(resampPropUniqueStreamIds)[0] = rc->propUniqueStreamIds;

        PROTECT(retList = allocVector(VECSXP, nComps)); ++nProtected;
        PROTECT(names   = allocVector(STRSXP, nComps)); ++nProtected;
        nComps = 0;
        SET_VECTOR_ELT(retList, nComps, resampCurrentStreams);
        SET_STRING_ELT(names,   nComps, mkChar("currentStreams"));
        ++nComps;
        SET_VECTOR_ELT(retList, nComps, resampCurrentLogWeights);
        SET_STRING_ELT(names,   nComps, mkChar("currentLogWeights"));
        ++nComps;
        SET_VECTOR_ELT(retList, nComps, resampPropUniqueStreamIds);
        SET_STRING_ELT(names,   nComps, mkChar("propUniqueStreamIds"));
        setAttrib(retList, R_NamesSymbol, names);
        UNPROTECT(nProtected);
        vmaxset(vmax);
        return retList;        
}


static SEXP
summary_func_user_Rfunc (Sampler *ss, int currentPeriod, SEXP currentStreams,
                         SEXP currentLogWeights)
{        
        ArgsList2 *al = ss->summaryArgsList;
        
        INTEGER(ss->SEXPCurrentPeriod)[0] = currentPeriod + 1;

        SET_VECTOR_ELT(al->argsList, al->posCurrentPeriod, ss->SEXPCurrentPeriod);
        SET_VECTOR_ELT(al->argsList, al->posCurrentStreams, currentStreams);
        SET_VECTOR_ELT(al->argsList, al->posCurrentLogWeights, currentLogWeights);

        SETCADR(ss->doCallFuncCall, ss->summaryFunc);
        SETCADDR(ss->doCallFuncCall, al->argsList);
        SETCADDDR(ss->doCallFuncCall, ss->dotsList);         
        return eval(ss->doCallFuncCall, ss->doCallFuncEnv);
}


static SEXP
summary_func_builtin (Sampler *ss, int currentPeriod, SEXP currentStreams,
                      SEXP currentLogWeights)
{
        int dspp = ss->dimSummPerPeriod, ns = ss->nStreams;
        int ii, jj, start, nProtected = 0;
        double *scs  = REAL(currentStreams);
        double *sclw = REAL(currentLogWeights);
        double *scaw = REAL(ss->SEXPCurrentAdjWeights);
        double sumcaw, *summ;
        SEXP SEXPSumm;
        void *vmax = vmaxget( );

        PROTECT(SEXPSumm = allocVector(REALSXP, dspp)); ++nProtected;
        summ = REAL(SEXPSumm);        
        /*
         * Note: here dimSummPerPeriod == dimPerPeriod and we only
         * provide the weighted average of each of the dimensions
         */
        sumcaw = sampler_adjust_log_weights(ns, sclw, scaw);
        for (ii = 0; ii < dspp; ++ii) {
                summ[ii] = 0.0; start = ii * ns;
                for (jj = 0; jj < ns; ++jj)
                        summ[ii] += scaw[jj] * scs[start + jj];
                summ[ii] /= sumcaw;
        }
        UNPROTECT(nProtected);
        vmaxset(vmax);
        return SEXPSumm;
}


static SEXP
MHUpdate_func_user_Rfunc (Sampler *ss, int currentPeriod, int nMHSteps,
                          SEXP currentStreams, SEXP lag1Streams,
                          SEXP lag1LogWeights)
{
        ArgsList3 *al = ss->MHUpdateArgsList;
        
        INTEGER(ss->SEXPCurrentPeriod)[0] = currentPeriod + 1;
        INTEGER(ss->SEXPNMHSteps)[0]      = nMHSteps;

        SET_VECTOR_ELT(al->argsList, al->posCurrentPeriod, ss->SEXPCurrentPeriod);
        SET_VECTOR_ELT(al->argsList, al->posNMHSteps, ss->SEXPNMHSteps);
        SET_VECTOR_ELT(al->argsList, al->posCurrentStreams, currentStreams);
        SET_VECTOR_ELT(al->argsList, al->posLag1Streams, lag1Streams);
        SET_VECTOR_ELT(al->argsList, al->posLag1LogWeights, lag1LogWeights);
        
        SETCADR(ss->doCallFuncCall, ss->MHUpdateFunc);
        SETCADDR(ss->doCallFuncCall, al->argsList);
        SETCADDDR(ss->doCallFuncCall, ss->dotsList);         
        return eval(ss->doCallFuncCall, ss->doCallFuncEnv);
}


static SEXP
one_iter_without_MH (Sampler *ss, SEXP notRequired)
{
        int nProtected = 0, cp = ss->currentPeriod, nstg = ss->nStreamsPreResamp;
        SEXP l1s, l1lw, cs, clw;
        SEXP retList1, retList2, summ;

        /* propagate */
        l1s  = ss->SEXPLag1Streams;
        l1lw = ss->SEXPLag1LogWeights;
        PROTECT(retList1 = (*(ss->propagate))(ss, cp, nstg, l1s, l1lw));
        ++nProtected;
        PHONY(printList(retList1););
        cs             = getListElement(retList1, "currentStreams");
        clw            = getListElement(retList1, "currentLogWeights");
        ss->resampDone = (*(ss->resampCriterion))(ss, cp, cs, clw);
        /* resamp? */
        if (ss->resampDone == TRUE) {
                /* resamp */
                PROTECT(retList2 = (*(ss->resamp))(ss, cp, cs, clw)); ++nProtected;
                cs  = getListElement(retList2, "currentStreams");
                clw = getListElement(retList2, "currentLogWeights");
        }
        else {
                retList2 = R_NilValue;
        }
        /* summary computation */
        PROTECT(summ = (*(ss->summary))(ss, cp, cs, clw)); ++nProtected;
        sampler_register_resamp(ss, retList2);
        sampler_register_streams_N_log_weights(ss, cs, clw);
        sampler_register_summary(ss, summ);
        UNPROTECT(nProtected);
        return 0;
}


static SEXP
one_iter_with_MH (Sampler *ss, SEXP notRequired)
{
        int nProtected = 0, cp = ss->currentPeriod;
        int nstg = ss->nStreamsPreResamp, nmhs = ss->nMHSteps;
        SEXP l1s, l1lw, cs, clw;
        SEXP retList1, retList2, retList3, summ;
        
        /* propagate */
        l1s  = ss->SEXPLag1Streams;
        l1lw = ss->SEXPLag1LogWeights;
        PROTECT(retList1 = (*(ss->propagate))(ss, cp, nstg, l1s, l1lw));
        ++nProtected;
        PHONY(printList(retList1););
        cs             = getListElement(retList1, "currentStreams");
        clw            = getListElement(retList1, "currentLogWeights");
        ss->resampDone = (*(ss->resampCriterion))(ss, cp, cs, clw);
        
        /* resamp? */
        if (ss->resampDone == TRUE) {
                SEXP l1s  = ss->SEXPLag1Streams;
                SEXP l1lw = ss->SEXPLag1LogWeights;
                
                /* resamp */
                PROTECT(retList2 = (*(ss->resamp))(ss, cp, cs, clw)); ++nProtected;
                cs  = getListElement(retList2, "currentStreams");
                clw = getListElement(retList2, "currentLogWeights");
                
                /* MHUpdate only after resamp */
                PROTECT(retList3 = (*(ss->MHUpdate))(ss, cp, nmhs, cs, l1s, l1lw));
                ++nProtected;
                cs = getListElement(retList3, "currentStreams");
        }
        else {
                retList2 = R_NilValue;
                retList3 = R_NilValue;
        }
        /* summary computation */
        PROTECT(summ = (*(ss->summary))(ss, cp, cs, clw)); ++nProtected;
        sampler_register_resamp(ss, retList2);
        sampler_register_acceptance_rates(ss, retList3);
        sampler_register_streams_N_log_weights(ss, cs, clw);
        sampler_register_summary(ss, summ);        
        UNPROTECT(nProtected);
        return 0;
}


static SEXP
register_summary (Sampler *ss, SEXP draws)
{
        int ii, dspp = ss->dimSummPerPeriod;
        int start = dspp * ss->currentPeriod;
        SEXP summary = VECTOR_ELT(draws, ss->summaryPos);
        double *dest = REAL(summary), *src = REAL(ss->SEXPSummary);

        /* summary is a ss->dimSummPeriod \times ss->nPeriods matrix */
        for (ii = 0; ii < dspp; ++ii) dest[start + ii] = src[ii];
        return R_NilValue;
}


static SEXP
register_propUniqueStreamIds (Sampler *ss, SEXP draws)
{
        double *dest = REAL(VECTOR_ELT(draws, ss->propUniqueStreamIdsPos));
        
        dest[ss->currentPeriod] = REAL(ss->SEXPPropUniqueStreamIds)[0];
        
        return R_NilValue;
}


static SEXP
register_streams (Sampler *ss, SEXP draws)
{
        int ii, nn = ss->nStreams * ss->dimPerPeriod;
        int start = nn * ss->currentPeriod;
        double *dest = REAL(VECTOR_ELT(draws, ss->streamsPos));
        double *src  = REAL(ss->SEXPLag1Streams);

        for (ii = 0; ii < nn; ++ii) dest[start + ii] = src[ii];
        return R_NilValue;
}


static SEXP
register_log_weights (Sampler *ss, SEXP draws)
{
        int ii, ns = ss->nStreams;
        int start = ns * ss->currentPeriod;
        double *dest = REAL(VECTOR_ELT(draws, ss->logWeightsPos));
        double *src  = REAL(ss->SEXPLag1LogWeights);

        for (ii = 0; ii < ns; ++ii) dest[start + ii] = src[ii];
        return R_NilValue;
}


static SEXP
register_acceptance_rates (Sampler *ss, SEXP draws)
{
        int ii, ns = ss->nStreams;
        int start = ns * ss->currentPeriod;
        double *dest = REAL(VECTOR_ELT(draws, ss->acceptanceRatesPos));
        double *src  = REAL(ss->SEXPAcceptanceRates);

        for (ii = 0; ii < ns; ++ii) dest[start + ii] = src[ii];
        return R_NilValue;
}


static SEXP
register_all (Sampler *ss, SEXP draws)
{
        register_summary(ss, draws);
        register_propUniqueStreamIds(ss, draws);
        if (ss->returnStreams == TRUE)    register_streams(ss, draws);
        if (ss->returnLogWeights == TRUE) register_log_weights(ss, draws);
        if (ss->nMHSteps > 0)             register_acceptance_rates(ss, draws);
        return R_NilValue;
}


int
sampler_init (Sampler *ss)
{
        ss->propagate = propagate_func_user_Rfunc;
        if (ss->resampCriterionFunc == R_NilValue)
                ss->resampCriterion = resampCriterion_func_builtin;
        else
                ss->resampCriterion = resampCriterion_func_user_Rfunc;
        if (ss->resampFunc == R_NilValue)
                ss->resamp = resamp_func_builtin_PPW;
        else
                ss->resamp = resamp_func_user_Rfunc;
        if (ss->summaryFunc == R_NilValue)
                ss->summary = summary_func_builtin;
        else
                ss->summary = summary_func_user_Rfunc;
        if (ss->MHUpdateFunc == R_NilValue) {
                ss->MHUpdate       = NULL;
                ss->samplerOneIter = one_iter_without_MH;
        }
        else {
                ss->MHUpdate       = MHUpdate_func_user_Rfunc;
                ss->samplerOneIter = one_iter_with_MH;
        }
        /*
         * The most common case is to just return summary, and hence
         * avoid the 'if' in the most common case and call
         * register_summary.
         */
        if ((ss->returnStreams == TRUE) ||
            (ss->returnLogWeights == TRUE))
                ss->registerThisDraw = register_all;
        else
                ss->registerThisDraw = register_summary;
        return 0;
}


static SEXP
make_draws (Sampler *ss)
{
        int ii, period, ns = ss->nStreams, np = ss->nPeriods;
        int nComps = 0, nProtected = 0;
        Rboolean initialEstDone = FALSE;
        double timeStartUsr, timeEndUsr, timeToFinish;
        double timeStartSys, timeEndSys;
        /*
         * draws is a list containing streams (if returnStreams ==
         * TRUE), log_weights (if returnLogWeights == TRUE), summary
         * and propUniqueStreamIds
         */
        SEXP draws, summary, propUniqueStreamIds, names;

        gather_time_details(ss);
        timeStartUsr = ss->timeDetails->usr;
        timeStartSys = ss->timeDetails->sys;
        ++nComps; /* for the summary component is always there */
        ++nComps; /* for the propUniqueStreamIds component is always there */
        if (ss->returnStreams == TRUE) ++nComps;
        if (ss->returnLogWeights == TRUE) ++nComps;
        if (ss->nMHSteps > 0) ++nComps;
        PROTECT(draws = allocVector(VECSXP, nComps)); ++(ss->nProtected);
        PROTECT(names = allocVector(STRSXP, nComps)); ++nProtected;
        nComps = 0;
        ss->summaryPos = nComps;
        PROTECT(summary = allocMatrix(REALSXP, ss->dimSummPerPeriod, np)); 
        SET_VECTOR_ELT(draws, nComps, summary);
        UNPROTECT(1);
        SET_STRING_ELT(names, nComps, mkChar("summary")); ++nComps;
        ss->propUniqueStreamIdsPos = nComps;
        PROTECT(propUniqueStreamIds = allocVector(REALSXP, np));
        SET_VECTOR_ELT(draws, nComps, propUniqueStreamIds);
        UNPROTECT(1);
        SET_STRING_ELT(names, nComps, mkChar("propUniqueStreamIds")); ++nComps;
        if (ss->returnStreams == TRUE) {
                SEXP streams, streamsDim;
                int ii, nn, *intsTmp, dims[3];

                ss->streamsPos = nComps;
                dims[0] = ns;
                dims[1] = ss->dimPerPeriod;
                dims[2] = np;
                nn = dims[0] * dims[1] * dims[2];
                PROTECT(streams    = allocVector(REALSXP, nn)); 
                PROTECT(streamsDim = allocVector(INTSXP, 3)); 
                intsTmp = INTEGER(streamsDim);
                for (ii = 0; ii < 3; ++ii) intsTmp[ii] = dims[ii];
                setAttrib(streams, R_DimSymbol, streamsDim);
                SET_VECTOR_ELT(draws, nComps, streams);
                UNPROTECT(2); 
                SET_STRING_ELT(names, nComps, mkChar("streams")); ++nComps;
        }
        if (ss->returnLogWeights == TRUE) {
                SEXP logWeights;
                
                ss->logWeightsPos  = nComps;
                PROTECT(logWeights = allocMatrix(REALSXP, ns, np)); 
                SET_VECTOR_ELT(draws, nComps, logWeights);
                UNPROTECT(1);
                SET_STRING_ELT(names, nComps, mkChar("logWeights")); ++nComps;
        }
        if (ss->nMHSteps > 0) {
                SEXP acceptanceRates;

                ss->acceptanceRatesPos  = nComps;
                PROTECT(acceptanceRates = allocMatrix(REALSXP, ns, np));
                SET_VECTOR_ELT(draws, nComps, acceptanceRates);
                UNPROTECT(1);
                SET_STRING_ELT(names, nComps, mkChar("acceptanceRates"));
        }
        setAttrib(draws, R_NamesSymbol, names);
        UNPROTECT(nProtected);
        
        GetRNGstate( );
        period = 0;             /* fake initialization: for gcc -Wall */
        for (ii = 0; ii < np; ++ii) {
                ss->currentPeriod = ii;
                (*(ss->samplerOneIter))(ss, R_NilValue);
                (*(ss->registerThisDraw))(ss, draws);
                if (ss->verboseLevel >= 1) {
                        period = ii + 1;
                        if (initialEstDone == FALSE) {
                                if ((period % ss->printInitialDotsWhen) == 0) Rprintf(".");
                        }
                        else {
                                if (period == ss->printDotAt) {
                                        Rprintf(".");
                                        ss->printDotAt += ss->eachDotWorth;
                                }
                        }
                        if (period == ss->printEstTimeAt) {
                                gather_time_details(ss);
                                timeEndUsr = ss->timeDetails->usr;
                                timeToFinish = ((timeEndUsr - timeStartUsr) / period * \
                                                (np - period));
                                Rprintf("\n[Time to finish (est): %8.0f secs, " \
                                        "this period: %8d]", ceil(timeToFinish), period);
                                if (initialEstDone == FALSE) initialEstDone = TRUE;
                                ss->printEstTimeAt += ss->eachDotWorth * ss->nDotsPerLine;
                                ss->printDotAt = period + 1;
                        }
                }
        }
        PutRNGstate( );
        if (ss->verboseLevel >= 1) {
                gather_time_details(ss);
                timeEndUsr = ss->timeDetails->usr;
                timeEndSys = ss->timeDetails->sys;
                Rprintf("\n[Total time: %10.0f secs (usr), %5.0f secs (sys), " \
                        "this period: %10d]\n", (timeEndUsr - timeStartUsr),
                        (timeEndSys - timeStartSys), period);
        }
        return draws;        
}


SEXP
sampler_run (Sampler *ss)
{
        int nComps = 0;
        SEXP draws, samplerObj, names;

        PHONY(sampler_print(ss););
        draws              = make_draws(ss); ++nComps;
        PROTECT(samplerObj = allocVector(VECSXP, nComps));
        PROTECT(names      = allocVector(STRSXP, nComps)); 
        ss->nProtected    += 2;
        
        /* fill in the samplerObj */
        SET_VECTOR_ELT(samplerObj, 0, draws);
        SET_STRING_ELT(names, 0, mkChar("draws"));
        setAttrib(samplerObj, R_NamesSymbol, names);
        PHONY(PRINT_STUB_INT(ss->nProtected););
        UNPROTECT(ss->nProtected);
        return samplerObj;
}


SEXP
SMCMainC (SEXP argsList)
{
        Sampler *ss;
        SEXP samplerObj;
        
        ss = sampler_new(argsList);
        sampler_init(ss);
        PROTECT(samplerObj = sampler_run(ss));        
        UNPROTECT(1);
        return samplerObj;
}


