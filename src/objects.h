/*
 *  objects.h
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
 */

#ifndef OBJECTS_H
#define OBJECTS_H

#include <R.h>
#include <Rinternals.h>
#include "R_ext/Boolean.h"

/*
 * The forward declaration of all the typedefs
 */
typedef struct ResampleContext ResampleContext;
typedef struct TimeDetails TimeDetails;
typedef struct ArgsList1 ArgsList1;
typedef struct ArgsList2 ArgsList2;
typedef struct ArgsList3 ArgsList3;
typedef struct Sampler Sampler;

typedef SEXP     (*SamplerUtil) (Sampler *, SEXP);
typedef SEXP     (*FuncPtr1) (Sampler *, int, int, SEXP, SEXP);
typedef Rboolean (*FuncPtr2_1) (Sampler *, int, SEXP, SEXP);
typedef SEXP     (*FuncPtr2_2) (Sampler *, int, SEXP, SEXP);
typedef SEXP     (*FuncPtr3) (Sampler *, int, int, SEXP, SEXP, SEXP);


struct ResampleContext
{
        int *streamIds, nUniqueStreamIds, *uniqueStreamIds;
        double propUniqueStreamIds, *partialSum;
};

struct TimeDetails
{
        double usr, sys;
};

struct ArgsList1 
{
        int posCurrentPeriod, posNStreamsToGenerate;
        int posLag1Streams, posLag1LogWeights;
        SEXP argsList;
};
        
struct ArgsList2 
{
        int posCurrentPeriod;
        int posCurrentStreams, posCurrentLogWeights;
        SEXP argsList;
};
        
struct ArgsList3 
{
        int posCurrentPeriod, posNMHSteps;                       
        int posCurrentStreams, posLag1Streams, posLag1LogWeights;
        SEXP argsList;
};
        
struct Sampler {
        int nStreams, nStreamsPreResamp, thisStream;
        int nPeriods, currentPeriod, dimPerPeriod, dimSummPerPeriod;
        int nMHSteps;

        Rboolean returnStreams, returnLogWeights, resampDone;
        
        int verboseLevel, printEstTimeAt, printEstTimeNTimes;
        int printInitialDotsWhen, printDotAt, eachDotWorth, nDotsPerLine;

        int nProtected;
        
        SEXP propagateFunc;
        ArgsList1 *propagateArgsList;        
        FuncPtr1 propagate;
        
        SEXP resampCriterionFunc;
        ArgsList2 *resampCriterionArgsList;
        FuncPtr2_1 resampCriterion;
        
        SEXP resampFunc;
        ArgsList2 *resampArgsList;
        FuncPtr2_2 resamp;
        
        SEXP summaryFunc;
        ArgsList2 *summaryArgsList;
        FuncPtr2_2 summary;
        
        SEXP MHUpdateFunc;
        ArgsList3 *MHUpdateArgsList;
        FuncPtr3 MHUpdate;

        SEXP doCallFuncCall, doCallFuncEnv, dotsList;        
        SEXP procTimeFuncCall, procTimeFuncEnv;
        TimeDetails *timeDetails;
        
        SamplerUtil samplerOneIter, registerThisDraw;

        SEXP SEXPCurrentPeriod, SEXPNStreamsToGenerate, SEXPNMHSteps;        
        /*
         * SEXPCurrentStreams & SEXPLag1Streams are of dimension:
         * nStreams \times dimPerPeriod
         */
        SEXP SEXPCurrentStreams, SEXPLag1Streams;

        /*
         * l(og) weights and norm(alized) weights are of dimension:
         * nStreams
         */
        SEXP SEXPCurrentLogWeights, SEXPCurrentAdjWeights;
        SEXP SEXPLag1LogWeights, SEXPLag1AdjWeights;
        SEXP SEXPAcceptanceRates;

        /*
         * SEXPSummary is of dimension:
         * dimSummPerPeriod
         */
        SEXP SEXPSummary, SEXPPropUniqueStreamIds;
        
        ResampleContext *scratch_RC;
        int summaryPos, propUniqueStreamIdsPos;
        int streamsPos, logWeightsPos, acceptanceRatesPos;
};

extern Sampler *
sampler_new (SEXP opts);

extern int
sampler_print (Sampler *ss);

extern int
sampler_init (Sampler *ss);

extern SEXP
sampler_run (Sampler *ss);

extern SEXP
SMCMainC (SEXP argsList);

#endif /* OBJECTS_H */

