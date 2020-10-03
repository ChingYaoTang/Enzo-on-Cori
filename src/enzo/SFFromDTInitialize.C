/***********************************************************************
/
/  INITIALIZE SIMULATION WITH THE STRUCTURE OF DT
/
/  written by: ChingYao Tang
/  date:       July, 2020
/
/  PURPOSE: Initializes simulation with the structure of the flow driven by stochastic forcing
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/


#include "preincludes.h"
#include <string.h>
#include <stdio.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "phys_constants.h"

void MHDCTSetupFieldLabels();

int SFFromDTInitialize(FILE *fptr, FILE *Outfptr, 
             HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";

  /* declarations */

  char line[MAX_LINE_LENGTH];
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  char *SFFromDTFile;
  int ret;

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("pointer dummy = %p\n", dummy);
    printf("pointer SFFromDTFile = %p\n", SFFromDTFile);
  }
  /* set default parameters specifying the random force field */


  /* set other default parameters */


  /* read input from file */
  if (debug) printf("SFFromDTInitialize: reading problem-specific parameters.\n");

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    if (sscanf(line, "SFFromDTFile = %s", dummy) == 1) {
      SFFromDTFile = dummy;
      ret++;
    }
    /* If the dummy char space was used, then make another. */
    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "SFFromDT") && line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR )
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("pointer dummy = %p\n", dummy);
    printf("pointer SFFromDTFile = %p\n", SFFromDTFile);
  }


  if (MetaData.TopGridRank != 3) {
      fprintf(stderr,"SFFromDTInitialize: Only 3D is availabe at this point.\n");
      return FALSE;
  }

  if (DualEnergyFormalism != 1) {
      fprintf(stderr,"SFFromDTInitialize: Only DualEnergyFormalism=1 is allowed.\n");
      return FALSE;
  }

  /* Begin grid initialization */
/*
  HierarchyEntry *CurrentGrid; // all level 0 grids on this processor first
  CurrentGrid = &TopGrid;
  while (CurrentGrid != NULL) {
      if (CurrentGrid->GridData->SFFromDTInitializeGrid(SFFromDTFile) == FAIL) {
          fprintf(stderr, "Error in SFFromDTInitializeGrid.\n");
          return FAIL;
      }
      CurrentGrid = CurrentGrid->NextGridThisLevel;
  }
*/
  if (TopGrid.GridData->SFFromDTInitializeGrid(SFFromDTFile) == FAIL) {
    ENZO_FAIL("Error in SFFromDTInitializeGrid.");
  }


  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;

  if ( UseMHD ) {
    DataLabel[count++] = BxName;
    DataLabel[count++] = ByName;
    DataLabel[count++] = BzName;
  }
  if ( HydroMethod == MHD_RK ){
    DataLabel[count++] = PhiName;
  }
  MHDCTSetupFieldLabels();

  if (MultiSpecies) {
    DataLabel[count++] =  ElectronName;
    DataLabel[count++] =  HIName;
    DataLabel[count++] =  HIIName;
    DataLabel[count++] =  HeIName;
    DataLabel[count++] =  HeIIName;
    DataLabel[count++] =  HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[count++] =  HMName;
      DataLabel[count++] =  H2IName;
      DataLabel[count++] =  H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[count++] =  DIName;
      DataLabel[count++] =  DIIName;
      DataLabel[count++] =  HDIName;
    }
  }  // if Multispecies

  for (int i = 0; i < count; i++)
    DataUnits[i] = NULL;

/* Write parameters to parameter output file */

  if (debug) printf("SFFromDTInitialize: writing parameters to output file.\n");

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "(output)SFFromDTFile    = %s\n\n", SFFromDTFile);

  }
  delete [] dummy;
  return SUCCESS;
}
