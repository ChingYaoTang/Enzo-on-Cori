/***********************************************************************
/
/  GRID CLASS: SFFromDTInitializeGrid
/
/  written by: ChingYao Tang
/  date:       July, 2020
/
/  PURPOSE: Initializes grid for a star formation simulation with DT structure
/
/  RETURNS: FAIL or SUCCESS
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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
              float *TemperatureUnits, float *TimeUnits,
              float *VelocityUnits, FLOAT Time);

int grid::SFFromDTInitializeGrid(char *SFFromDTFile)
{
  /* declarations */

  int DensNum, GENum, TENum;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, B1Num, B2Num, B3Num, PhiNum;
  
  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[DensNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++] = TotalEnergy;
  FieldType[GENum = NumberOfBaryonFields++] = InternalEnergy;

  int ivel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;

  if (UseMHD) {
    FieldType[B1Num = NumberOfBaryonFields++] = Bfield1;
    FieldType[B2Num = NumberOfBaryonFields++] = Bfield2;
    FieldType[B3Num = NumberOfBaryonFields++] = Bfield3;
  
    if (HydroMethod == MHD_RK) {
      FieldType[PhiNum = NumberOfBaryonFields++] = PhiField;
    }
    if (UsePoissonDivergenceCleaning) {
      FieldType[NumberOfBaryonFields++] = Phi_pField;
    }
  }

  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
        VelocityUnits;

  int size = 1;

  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  printf("size = %d\n", size);

  /* Allocate fields. */
  this->AllocateGrids();


  /* set baryon fields by reading the  */
  hid_t    file, space, dset;
  herr_t   status;
  hsize_t  dims[3];
  int      ndims;

  file = H5Fopen(SFFromDTFile, H5F_ACC_RDONLY, H5P_DEFAULT);

  int index, n;
  dset   = H5Dopen(file, "/Density");
  space  = H5Dget_space(dset);
  ndims  = H5Sget_simple_extent_dims(space, dims, NULL);
  double *RBuff = new double[dims[0]*dims[1]*dims[2]]; 
  printf("size of RBuff = %d\n", dims[0]*dims[1]*dims[2]);
  printf("point RBuff = %p\n", RBuff);

  /* density */
  status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
  for(int k=0; k<GridDimension[2]; k++)
    for(int j=0; j<GridDimension[1]; j++)
      for(int i=0; i<GridDimension[0]; i++){
        index                       = i + GridDimension[0]*(j + GridDimension[1]*k);
        n                           = (i*dims[1]+j)*dims[2]+k;
        BaryonField[DensNum][index] = RBuff[n];
        
      }
  status = H5Dclose (dset);
  status = H5Sclose (space);
  printf("Allocate D\n");


  /* internal energy*/
  dset   = H5Dopen(file, "/GasEnergy");
  status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
#ifdef BYHAND
// sum over GE of all cells
  double *avg_GE = new double(0.0);
  for(int kk=0; kk<=dims[2]; kk++)
    for(int jj=0; jj<=dims[1]; jj++)
      for(int ii=0; ii<=dims[0]; ii++){
        n        = (ii*dims[1]+jj)*dims[2]+kk;
// divide by number of cell => average GE  
        *avg_GE += (RBuff[n]/(dims[0]*dims[1]*dims[2]));
      }
  printf("Avg GE = %e\n", *avg_GE);
#endif
// allocate avg GE to all cells
  for(int k=0; k<GridDimension[2]; k++)
    for(int j=0; j<GridDimension[1]; j++)
      for(int i=0; i<GridDimension[0]; i++){
        index                     = i + GridDimension[0]*(j + GridDimension[1]*k);
#ifdef BYHAND
        BaryonField[GENum][index] = *avg_GE;
#else
        n                         = (i*dims[1]+j)*dims[2]+k;
        BaryonField[GENum][index] = RBuff[n];
#endif
      }
  status = H5Dclose (dset);
  printf("Allocate GE\n");


  /* velocity */
#ifdef BYHAND
// reckon the avg velocity
  double *Vel_sq = new double[dims[0]*dims[1]*dims[2]]; 
  double *avg_KE = new double(0.0);
#endif
// Vel_sq = Vx^2
  dset   = H5Dopen(file, "/x-velocity");
  status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
  printf("Read in Vx\n");
  for(int k=0; k<GridDimension[2]; k++)
    for(int j=0; j<GridDimension[1]; j++)
      for(int i=0; i<GridDimension[0]; i++){
        index     = i + GridDimension[0]*(j + GridDimension[1]*k);
        n         = (i*dims[1]+j)*dims[2]+k;
#ifdef BYHAND
        Vel_sq[n] = (RBuff[n]*RBuff[n]);
//        printf("Vx^2 at %d = %g\n", n, Vel_sq[n]);
#endif
// alllocate original data
        BaryonField[ivel][index] = RBuff[n];
      }
  status = H5Dclose (dset);

// + Vy^2
  dset   = H5Dopen(file, "/y-velocity");
  status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
  printf("Read in Vy\n");
  for(int k=0; k<GridDimension[2]; k++)
    for(int j=0; j<GridDimension[1]; j++)
      for(int i=0; i<GridDimension[0]; i++){
        index      = i + GridDimension[0]*(j + GridDimension[1]*k);
        n          = (i*dims[1]+j)*dims[2]+k;
#ifdef BYHAND
        Vel_sq[n] += (RBuff[n]*RBuff[n]);
//        printf("Vx^2+Vy^2 at %d = %g\n", n, Vel_sq[n]);
#endif
// alllocate original data
        BaryonField[ivel+1][index] = RBuff[n];
      }
  status = H5Dclose (dset);

// + Vz^2
  dset   = H5Dopen(file, "/z-velocity");
  status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
  printf("Read in Vz\n");
  for(int k=0; k<GridDimension[2]; k++)
    for(int j=0; j<GridDimension[1]; j++)
      for(int i=0; i<GridDimension[0]; i++){
        index      = i + GridDimension[0]*(j + GridDimension[1]*k);
        n          = (i*dims[1]+j)*dims[2]+k;
#ifdef BYHAND
        Vel_sq[n] += (RBuff[n]*RBuff[n]);
        if(Vel_sq[n]>1e2)
            printf("V^2 at (%d,%d,%d) = %e\n", i, j, k, Vel_sq[n]);
#endif
// alllocate original data
        BaryonField[ivel+2][index] = RBuff[n];
      }
  status = H5Dclose (dset);

#ifdef BYHAND
  printf("Avg KE = %e\n", *avg_KE);
// avg_KE = 0.5*Sum(Vel_sq)/(dim^3)
  for(int kk=0; kk<=dims[2]; kk++)
    for(int jj=0; jj<=dims[1]; jj++)
      for(int ii=0; ii<=dims[0]; ii++){
        n        = (ii*dims[1]+jj)*dims[2]+kk;
        *avg_KE = *avg_KE + (Vel_sq[n]/(dims[0]*dims[1]*dims[2]*2));
        if(*avg_KE>1e4)
          printf("ke at (%d,%d,%d) = %e\n", ii, jj, kk, *avg_KE);
      }
  printf("Avg KE = %g\n", *avg_KE);

// rescale velocity fields on all cells
  for(int k=0; k<GridDimension[2]; k++)
    for(int j=0; j<GridDimension[1]; j++)
      for(int i=0; i<GridDimension[0]; i++){
        index = i + GridDimension[0]*(j + GridDimension[1]*k);
        BaryonField[ivel][index]   *= (*avg_GE)/(*avg_KE); 
        BaryonField[ivel+1][index] *= (*avg_GE)/(*avg_KE);
        BaryonField[ivel+2][index] *= (*avg_GE)/(*avg_KE);
      }
  printf("Coefficient = %g\n", (*avg_GE)/(*avg_KE));
#endif
  printf("Allocate Vel\n");

  /* total energy*/
#ifdef BYHAND
 printf("TE=GE+0.5*V^2");
#else
  dset   = H5Dopen(file, "/TotalEnergy");
  status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
#endif
  for(int k=0; k<GridDimension[2]; k++)
    for(int j=0; j<GridDimension[1]; j++)
      for(int i=0; i<GridDimension[0]; i++){
        index      = i + GridDimension[0]*(j + GridDimension[1]*k);
#ifdef BYHAND        
        BaryonField[TENum][index] = BaryonField[GENum][index] + 0.5*( pow(BaryonField[ivel][index],2) 
                                                                    + pow(BaryonField[ivel+1][index],2) 
                                                                    + pow(BaryonField[ivel+2][index],2) );
#else
        n          = (i*dims[1]+j)*dims[2]+k;
        BaryonField[TENum][index] = RBuff[n];
#endif
      }
  printf("Allocate TE\n");
#ifdef BYHAND  
  delete avg_GE; 
  delete Vel_sq; 
  delete avg_KE; 
#else
  status = H5Dclose (dset);
#endif

//  if (HydroMethod == MHD_RK) {
//      for (int i = 0; i < size; i++) {
//          BaryonField[iBx  ][i]  = ;
//      }
//  }
  
  /* If doing multi-species (HI, etc.), set these. */
  if (MultiSpecies > 0) {
    dset   = H5Dopen(file, "/Electron_Density");
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
    for(int k=0; k<GridDimension[2]; k++)
      for(int j=0; j<GridDimension[1]; j++)
        for(int i=0; i<GridDimension[0]; i++){
          index = i + GridDimension[0]*(j + GridDimension[1]*k);
          n     = (i*dims[1]+j)*dims[2]+k;
          BaryonField[DeNum][index] = RBuff[n];
        
        }
    status = H5Dclose (dset);

    dset   = H5Dopen(file, "/HI_Density");
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
    for(int k=0; k<GridDimension[2]; k++)
      for(int j=0; j<GridDimension[1]; j++)
        for(int i=0; i<GridDimension[0]; i++){
          index = i + GridDimension[0]*(j + GridDimension[1]*k);
          n     = (i*dims[1]+j)*dims[2]+k;
          BaryonField[HINum][index] = RBuff[n];
        
        }
    status = H5Dclose (dset);

    dset   = H5Dopen(file, "/HII_Density");
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
    for(int k=0; k<GridDimension[2]; k++)
      for(int j=0; j<GridDimension[1]; j++)
        for(int i=0; i<GridDimension[0]; i++){
          index = i + GridDimension[0]*(j + GridDimension[1]*k);
          n     = (i*dims[1]+j)*dims[2]+k;
          BaryonField[HIINum][index] = RBuff[n];
        
        }
    status = H5Dclose (dset);

    dset   = H5Dopen(file, "/HeI_Density");
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
    for(int k=0; k<GridDimension[2]; k++)
      for(int j=0; j<GridDimension[1]; j++)
        for(int i=0; i<GridDimension[0]; i++){
          index = i + GridDimension[0]*(j + GridDimension[1]*k);
          n     = (i*dims[1]+j)*dims[2]+k;
          BaryonField[HeINum][index] = RBuff[n];
        
        }
    status = H5Dclose (dset);

    dset   = H5Dopen(file, "/HeII_Density");
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
    for(int k=0; k<GridDimension[2]; k++)
      for(int j=0; j<GridDimension[1]; j++)
        for(int i=0; i<GridDimension[0]; i++){
          index = i + GridDimension[0]*(j + GridDimension[1]*k);
          n     = (i*dims[1]+j)*dims[2]+k;
          BaryonField[HeIINum][index] = RBuff[n];
        
        }
    status = H5Dclose (dset);

    dset   = H5Dopen(file, "/HeIII_Density");
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
    for(int k=0; k<GridDimension[2]; k++)
      for(int j=0; j<GridDimension[1]; j++)
        for(int i=0; i<GridDimension[0]; i++){
          index = i + GridDimension[0]*(j + GridDimension[1]*k);
          n     = (i*dims[1]+j)*dims[2]+k;
          BaryonField[HeIIINum][index] = RBuff[n];
        
        }
    status = H5Dclose (dset);

  
    if (MultiSpecies > 1) {
      dset   = H5Dopen(file, "/HM_Density");
      status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
      for(int k=0; k<GridDimension[2]; k++)
        for(int j=0; j<GridDimension[1]; j++)
          for(int i=0; i<GridDimension[0]; i++){
            index = i + GridDimension[0]*(j + GridDimension[1]*k);
            n     = (i*dims[1]+j)*dims[2]+k;
            BaryonField[HMNum][index] = RBuff[n];
        
          }
      status = H5Dclose (dset);
          
      dset   = H5Dopen(file, "/H2I_Density");
      status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
      for(int k=0; k<GridDimension[2]; k++)
        for(int j=0; j<GridDimension[1]; j++)
          for(int i=0; i<GridDimension[0]; i++){
            index = i + GridDimension[0]*(j + GridDimension[1]*k);
            n     = (i*dims[1]+j)*dims[2]+k;
            BaryonField[H2INum][index] = RBuff[n];
        
          }
      status = H5Dclose (dset);

      dset   = H5Dopen(file, "/H2II_Density");
      status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
      for(int k=0; k<GridDimension[2]; k++)
        for(int j=0; j<GridDimension[1]; j++)
          for(int i=0; i<GridDimension[0]; i++){
            index = i + GridDimension[0]*(j + GridDimension[1]*k);
            n     = (i*dims[1]+j)*dims[2]+k;
            BaryonField[H2IINum][index] = RBuff[n];
        
          }
      status = H5Dclose (dset);
    }
    if (MultiSpecies > 2) {
      dset   = H5Dopen(file, "/DI_Density");
      status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
      for(int k=0; k<GridDimension[2]; k++)
        for(int j=0; j<GridDimension[1]; j++)
          for(int i=0; i<GridDimension[0]; i++){
            index = i + GridDimension[0]*(j + GridDimension[1]*k);
            n     = (i*dims[1]+j)*dims[2]+k;
            BaryonField[DINum][index] = RBuff[n];
        
          }
      status = H5Dclose (dset);

      dset   = H5Dopen(file, "/DII_Density");
      status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
      for(int k=0; k<GridDimension[2]; k++)
        for(int j=0; j<GridDimension[1]; j++)
          for(int i=0; i<GridDimension[0]; i++){
            index = i + GridDimension[0]*(j + GridDimension[1]*k);
            n     = (i*dims[1]+j)*dims[2]+k;
            BaryonField[DIINum][index] = RBuff[n];
        
          }
      status = H5Dclose (dset);

      dset   = H5Dopen(file, "/HDI_Density");
      status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, RBuff);
      for(int k=0; k<GridDimension[2]; k++)
        for(int j=0; j<GridDimension[1]; j++)
          for(int i=0; i<GridDimension[0]; i++){
            index = i + GridDimension[0]*(j + GridDimension[1]*k);
            n     = (i*dims[1]+j)*dims[2]+k;
            BaryonField[HDINum][index] = RBuff[n];
        
          }
      status = H5Dclose (dset);

    }
  }
  printf("Allocate Multi-species\n");
  
  
  status = H5Fclose (file);
  delete RBuff; 
  
  return SUCCESS;
}
