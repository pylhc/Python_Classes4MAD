#include "integration.h"

/* Performs calculation of auxilliary variables after integration, event finding
   are complete. This is the only place where memory is allocated outside of
   an "Init" function. */
void AuxVarCalc( IData *GS ) {
  int i;
  
  assert(GS);

  if( (GS->nAuxVars > 0) && (GS->pointsIdx > 0 ) ) {
    
    /* Allocate space for aux variable calculations if this is the first
       time we have done aux calculations */
    if( GS->gAuxPoints == NULL ) {
      GS->gAuxPoints = (double **)PyMem_Malloc(GS->pointsIdx*sizeof(double *));
      assert(GS->gAuxPoints);
      
      for( i = 0; i < GS->pointsIdx; i++ ) {
	GS->gAuxPoints[i] = (double *)PyMem_Malloc(GS->nAuxVars*sizeof(double));
	assert(GS->gAuxPoints[i]);
      }
    }
    /* Otherwise, realloc space for the additional pointers, PyMem_Malloc
       space for the new aux points */
    else {
      double **temp = NULL;
      temp = (double **)PyMem_Realloc(GS->gAuxPoints,GS->pointsIdx*sizeof(double *));
      assert(temp);
      GS->gAuxPoints = temp;
      for( i = GS->auxIdx; i < GS->pointsIdx; i++ ) {
	GS->gAuxPoints[i] = (double *)PyMem_Malloc(GS->nAuxVars*sizeof(double));
	assert(GS->gAuxPoints[i]);
      }
    }
    
   /* Perform calculation on each integrated point */
    for( i = GS->auxIdx; i < GS->pointsIdx; i++ ) {
      if( GS->nExtInputs > 0 ) {
	FillCurrentExtInputValues(GS, GS->gTimeV[i]);
      }
      auxvars((unsigned) GS->phaseDim, (unsigned) GS->paramDim, GS->gTimeV[i], 
	      GS->gPoints[i], GS->gParams, GS->gAuxPoints[i], 
	      (unsigned) GS->extraSpaceSize, GS->gExtraSpace, 
	      (unsigned) GS->nExtInputs, GS->gCurrentExtInputVals);
    }
    
    GS->auxIdx = GS->pointsIdx;    
  }
}


/* Takes IData struct, (global), time, point, and phase space dimension
   as input.
   Saves the point and time in the global points/times arrays.
   Increments the global point counter */
void OutputPoint( IData *GS, double t, double *y ) {
  int i, j;
  int skipPoint = 0;
  int insertSpecPoints = 0;
  int addSpecIdx = 0;

  /* If there is space left in the points/times vectors */
  if( GS->timeIdx < GS->maxPts && GS->calcSpecTimes != 0 && GS->specTimesLen > 0) {

    /* Calculate how many specific-time points should be inserted
       between the previous point and this point */
    if( GS->direction > 0) {
      for( i = GS->specTimesIdx; (i < GS->specTimesLen &&  GS->gSpecTimes[i] <= t); i++ ) {
	if( GS->gSpecTimes[i] < t ) {
	  insertSpecPoints++;
	}
	/* We will not insert but will skip over multiple spectimes identical to t */
	addSpecIdx++;
      }
    }
    else {
      for( i = GS->specTimesIdx; (i < GS->specTimesLen &&  GS->gSpecTimes[i] >= t); i++ ) {
	if( GS->gSpecTimes[i] > t ) {
	  insertSpecPoints++;
	}
	/* We will not insert but will skip over multiple spectimes identical to t */
	addSpecIdx++;
      }
    }

    /* If t is a duplicate time, then insertSpecPoints can be > 0 only once, not 
       on subsequent t's that are identical. */

    /* Insert these points */
    for( j = 0; (GS->timeIdx < GS->maxPts && j < insertSpecPoints); j++ ) {
      for( i = 0; i < GS->phaseDim; i++ ) {
	GS->gPoints[GS->pointsIdx][i] = GS->cContSolFun(i, GS->gSpecTimes[GS->specTimesIdx+j]);
      }
      GS->pointsIdx++;
      GS->gTimeV[GS->timeIdx] = GS->gSpecTimes[GS->specTimesIdx+j];
      GS->timeIdx++;
    }

    /* Skip ahead the right amount in the spec times array */
    GS->specTimesIdx += addSpecIdx;
  }
 
/*   if( GS->timeIdx < GS->maxPts ) { */
/*      If no events were found, we may assume that the point */
/*        to be output is different from the previous point  */
/*     if( GS->eventFound == 0 ) */
      
/*       skipPoint = 0; */
/*      If this is the first point, output it  */
/*     else if( GS->timeIdx == 0 ) */
/*       skipPoint = 0; */
/*      Otherwise, check that this point doesn't duplicate  */
/*        the last point (we assume that integration points */
/*        differ)  */
/*     else if ( t - GS->gTimeV[GS->timeIdx-1] == 0 ) */
/*       skipPoint = 1; */
/*     else */
/*       skipPoint = 0; */

/*      If we don't need to skip the point because it's a  */
/*        duplicate, output it */ 
/*     if( skipPoint == 0 ) { */
/*        Save the point */ 
/*       for(i = 0; i < GS->phaseDim; i++) { */
/* 	GS->gPoints[GS->pointsIdx][i] = y[i]; */
/*       } */
/*       GS->pointsIdx++; */

/*        Save the time */
/*       GS->gTimeV[GS->timeIdx] = t;  */
/*       GS->timeIdx++; */
/*     } */
/*   } */

    if( GS->timeIdx < GS->maxPts ) {
      if( (GS->timeIdx == 0) || (t - GS->gTimeV[GS->timeIdx-1] != 0 ) ) {
	for( i = 0; i < GS->phaseDim; i++ ) {
	  GS->gPoints[GS->pointsIdx][i] = y[i];
	}
	GS->pointsIdx++;
	GS->gTimeV[GS->timeIdx] = t;
	GS->timeIdx++;
      }
    }
}

/* Takes IData struct (global), time, point, phaseDim, phase space dimension
   and current index into refinement point buffer as input.
   Saves the point and time in the global refinement points/times arrays
   for later incorporation into the trajectory. */
void BufferRefinePoint( IData *GS, double t, double *y, int idx ) {
  int i;
  
  if( (idx >= 0) && (idx < GS->refine) ) {
    for( i = 0; i < GS->phaseDim; i++ ) {
      GS->gRefinePoints[idx][i] = y[i];
    }
    GS->gRefineTimes[idx] = t;
  }
}

/* Takes IData struct (global), time, point, phaseDim, phase space dimension
   and current index into event point buffer as input.
   Saves the point and time in the global refinement points/times arrays
   for later incorporation into the trajectory. */
void BufferEventPoint( IData *GS, double t, double *y, int idx ) {
  int i;
  
  if( (idx >= 0) && (idx < GS->haveActive) ) {
    for( i = 0; i < GS->phaseDim; i++ ) {
      GS->gEventPointBuf[idx][i] = y[i];
    }
    GS->gEventTimeBuf[idx] = t;
  }

}


void SavePoints(IData *GS, int TotEvents, int foundTerm) {
  int i, j, term = 0;
  
  /* Negative TotEvents means a terminal event was found */
  if( TotEvents < 0 ) {
    TotEvents = -TotEvents;
    term = 1;
  }
  term = foundTerm;

  /* !!!!!!!!! Mergesort on output the contents of the refine, event buffers */
  /* Remember to account for terminal events being last! */
  
  /* If no events were found, just output refinement points, if necessary,
     and the last integration point */
  if ( TotEvents == 0 ) {
    if ( GS->refine > 0 ) {
      for( i = 0; i < GS->refine; i++ ) {
	OutputPoint(GS, GS->gRefineTimes[i], GS->gRefinePoints[i]);
      }
    }
    // Output the last point (not an event point if no terminal events) 
    OutputPoint(GS, GS->lastTime, GS->lastPoint);
  } 
  /* If events were found, output is more complicated */
  else if ( TotEvents > 0 ) { 
    i = 0; j = 0; 
    /* If there are no refinement points, just output event points */
    if ( GS->refine <= 0 ) { 
      for( i = 0; i < TotEvents; i++ ) { 
 	OutputPoint(GS, GS->gEventTimeBuf[i], GS->gEventPointBuf[i]); 
      } 
    } 
    else {
      int k = 0; 
      /* Merge the event points and refinement points */
      while ( i + j < TotEvents + GS->refine ) { 
	/* If there are no more event points to include, or if the current refinement
	   point is earlier than the current event point, save the current 
	   refinement point */
 	if ( (i >= TotEvents) ||  
 	     (( j < GS->refine) && (GS->gRefineTimes[j] < GS->gEventTimeBuf[i])) ) { 
 	  OutputPoint(GS, GS->gRefineTimes[j], GS->gRefinePoints[j]); 
 	  j++; 
 	} 
	/* If there are no more refinement points to include, or if the current event
	   point is earlier than the current refinement point, save the current 
	   event point */
	else if ( (j >= GS->refine) ||  
 		    (( i < TotEvents) && (GS->gEventTimeBuf[i] < GS->gRefineTimes[j])) ) { 
 	  OutputPoint(GS, GS->gEventTimeBuf[i], GS->gEventPointBuf[i]); 
 	  i++; 
 	} 
 	else { /* The current event time and the current refine time are 
 		  identical. Only save one. */ 
 	  OutputPoint(GS, GS->gRefineTimes[j], GS->gRefinePoints[j]);
 	  i++; j++; 
 	}
      }
    } 

    if( term == 0 ) {
      /* If there were no terminal events, we need to output the last integrated
	 point, too */
      OutputPoint(GS, GS->lastTime, GS->lastPoint);
    }
  }
}

void FillCurrentExtInputValues( IData *GS, double t ) {
  
  int i;

  for( i = 0; i < GS->nExtInputs; i++ ) {
    if( GS->gExtInputLens[i] > 0 ) {
      GS->gCurrentExtInputVals[i] = GetCurrentExtInputValue( GS, t, i );
    }
    else {
      GS->gCurrentExtInputVals[i] = 0;
    }
  }
}


double GetCurrentExtInputValue( IData *GS, double t, int idx ) {
  int i, curidx = 0;


  /* If there is only one time point for the external input, then return its 
     associated value iff t is that exact time. */
  if( GS->gExtInputLens[idx] == 1 && GS->gExtInputTimes[idx][0] == t ) {
    return GS->gExtInputVals[idx][0];
  }
  
  /* Forward integration case */
  if( GS->direction > 0 ) {
    /* If the given time is outside the valid time domain for the external input,
       return 0.0 */
    if( t < GS->gExtInputTimes[idx][0] || t > GS->gExtInputTimes[idx][GS->gExtInputLens[idx]-1] ) {
      return 0.0;
    }
    
    curidx = FindCurrentExtInputIndex(GS, t, idx); 
/*     //    for( i = 0; i < GS->gExtInputLens[idx]; i++ ) { */
/*     //  if( GS->gExtInputTimes[idx][i]  >= t ) { */
/*     //	curidx = i - 1; */
/*     //	break; */
/*     // } */
/*     //} */
  }
  /* Backwards case */
  else if( GS->direction < 0 ) {
    /* If the given time is outside the valid time domain for the external input,
       return 0.0 */
    if( t > GS->gExtInputTimes[idx][0] || t < GS->gExtInputTimes[idx][GS->gExtInputLens[idx]-1] ) {
      return 0.0;
    }

    curidx = FindCurrentExtInputIndex(GS, t, idx); 
/*     //    for( i = 0; i < GS->gExtInputLens[idx]; i++ ) { */
/*     //  if( GS->gExtInputTimes[idx][i]  <= t ) { */
/*     //	curidx = i - 1; */
/*     //	break; */
/*     // } */
/*     //} */
  }
  /* If direction not set, return 0.0 */
  else
    return 0.0;

  /* Record the current index */
  GS->gCurrentExtInputIndex[idx] = curidx;

  /* If we are right on t for which ext input value exists */
  if( GS->gExtInputTimes[idx][curidx+1] == t ) {
    return GS->gExtInputVals[idx][curidx+1];
  }
  /* Since we assume the time, value arrays are oriented correctly for the given direction,
     the interpolation proceeds identically in both cases. */
  else {
    double temp = t - GS->gExtInputTimes[idx][curidx];
    /* Linear interpolation. We assume each time in gExtInputTimes is distinct */
    return (GS->gExtInputVals[idx][curidx+1] - GS->gExtInputVals[idx][curidx]) / 
      (GS->gExtInputTimes[idx][curidx+1] - GS->gExtInputTimes[idx][curidx]) * temp + GS->gExtInputVals[idx][curidx];
  }
}

int FindCurrentExtInputIndex( IData *GS, double t, int idx ) {
  int i, curidx = 0, newidx = 0;

  curidx = GS->gCurrentExtInputIndex[idx];

  /* We assume direction != 0 */

  /* Forward case */
  if( GS->direction > 0 ) {
    /* If we are not ahead of ourselves */
    if( GS->gExtInputTimes[idx][curidx] < t ) {
      for( i = curidx; i < GS->gExtInputLens[idx]; i++ ) {
	if( GS->gExtInputTimes[idx][i] >= t ) {
	  newidx = i - 1;
	  break;
	}
      }
    }
    else { 
      for( i = curidx; i >= 0; i-- ) {
	if( GS->gExtInputTimes[idx][i] < t ) {
	  newidx = i;
	  break;
	}
      }
    }
  }
  else {
    /* If we are not behind ourselves */
    if( GS->gExtInputTimes[idx][curidx] >= t ) {
      for( i = curidx; i < GS->gExtInputLens[idx]; i++ ) {
	if( GS->gExtInputTimes[idx][i] < t ) {
	  newidx = i - 1;
	  break;
	}
      }
    }
    else {
      for( i = curidx; i >= 0; i-- ) {
	if( GS->gExtInputTimes[idx][i] >= t ) {
	  newidx = i;
	  break;
	}
      }
    }
  }
   
  return newidx;
}
