/***********************************************************************
/
/  CREATES STAR PARTICLES FROM EXISTING PARTICLES
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  NOTES:  negative types mark particles that have just been before 
/          and not been converted into a star particle.
/
************************************************************************/
#include <stdlib.h>
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

void InsertStarAfter(Star * &Node, Star * &NewNode);

int grid::FindNewStarParticles(int level)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0)
    return SUCCESS;

  int i;
  Star *NewStar, *cstar;
  bool exists;

  for (i = 0; i < NumberOfParticles; i++)
    if (ParticleType[i] == -PARTICLE_TYPE_SINGLE_STAR ||
	ParticleType[i] == -PARTICLE_TYPE_BLACK_HOLE ||
	ParticleType[i] == -PARTICLE_TYPE_CLUSTER ||
	ParticleType[i] == -PARTICLE_TYPE_COLOR_STAR ||
	ParticleType[i] == -PARTICLE_TYPE_SIMPLE_SOURCE ||
	ABS(ParticleType[i]) == PARTICLE_TYPE_MBH ||
	(StarParticleRadiativeFeedback == TRUE &&
	 ParticleType[i] == PARTICLE_TYPE_STAR)) {

      // Check if it already exists (wasn't activated on the last
      // timestep, usually because of insufficient mass)
      exists = false;
      for (cstar = Stars; cstar; cstar = cstar->NextStar)
	if (cstar->Identifier == ParticleNumber[i]) {
	  cstar->SetLevel(level);
	  exists = true;
	  break;
	}

      if (!exists) {
	NewStar = new Star(this, i, level);

	/* If using an IMF for Pop III stars, assign the mass after
	   merging new (massless) star particles to avoid unnecessary
	   calls to the IMF. */

	if (ParticleType[i] == -PARTICLE_TYPE_SINGLE_STAR)
	  if (PopIIIInitialMassFunction == FALSE){
	    NewStar->AssignFinalMass(PopIIIStarMass);	
        if (debug2) 
            printf("Grid_FindNewStarParticles: Assign FinalMass %e Msun to particle[star_id=%d, star_type=%d, star_mass=%e Msun] & insert it to the star list of this grid[id=%d]\n", NewStar->FinalMass,NewStar->Identifier, NewStar->type, NewStar->Mass, ID);
      }
    if (ParticleType[i] == -PARTICLE_TYPE_SIMPLE_SOURCE) 
	  NewStar->AssignFinalMass(PopIIIStarMass);
	InsertStarAfter(Stars, NewStar);
	NumberOfStars++;
      }

    }

  return SUCCESS;

}
