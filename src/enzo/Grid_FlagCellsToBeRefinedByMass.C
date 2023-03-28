/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY MASS)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::FlagCellsToBeRefinedByMass(int level, int method, int RestrictFlag)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* error check */

  int ThisFlaggingMethod = CellFlaggingMethod[method];
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  if (MassFlaggingField == NULL && ThisFlaggingMethod == 2) {
    fprintf(stderr, "Mass Flagging Field is undefined.\n");
    return -1;
  }
 
  if (ParticleMassFlaggingField == NULL && ThisFlaggingMethod == 4) {
    fprintf(stderr, "Mass Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Compute the ModifiedMinimumMass */
 
  float ModifiedMinimumMassForRefinement =
    MinimumMassForRefinement[method]*
//      (MinimumMassForRefinementLevelExponent[method]*(log10(level+1e1)-1) +1) * POW(8, -level); //1
//      POW(log(level + exp(1.0)), MinimumMassForRefinementLevelExponent[method]) * POW(8, -level); //2
      POW(RefineBy, level*MinimumMassForRefinementLevelExponent[method]); //3
  if (ProblemType == 28)
    ModifiedMinimumMassForRefinement = 0;
 
  /* Flag points */

  float *ffield;
  if (ThisFlaggingMethod == 2)
    ffield = MassFlaggingField;
  else if (ThisFlaggingMethod == 4)
    ffield = ParticleMassFlaggingField;
  else {
    ENZO_VFAIL("Unrecognized mass refinement flagging method (%"ISYM")\n", 
	    method)
  }

  /* If RestrictFlag is true, then cells will only be flagged if they
     are flagged by mass.  This happens when using must-refine
     particles and must happen after all of the other flagging
     methods. */

  if (RestrictFlag) {
    for (i = 0; i < size; i++)
      FlaggingField[i] = ((ffield[i] > ModifiedMinimumMassForRefinement) &&
			  (FlaggingField[i] > 0)) ? 1 : 0;
  } else {
    for (i = 0; i < size; i++)
      FlaggingField[i] += (ffield[i] > ModifiedMinimumMassForRefinement) ? 1 : 0;
  }

  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i] > 0)
      NumberOfFlaggedCells++;
 
  /* remove MassFlaggingField. */
 
  if (ThisFlaggingMethod == 2) {
    delete [] MassFlaggingField;
    MassFlaggingField = NULL;
  } else if (ThisFlaggingMethod == 4) {

    delete [] ParticleMassFlaggingField;
    ParticleMassFlaggingField = NULL;
  }
 
  return NumberOfFlaggedCells;
 
}
