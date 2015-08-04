/*****************************************************************************
 * calc_dose.cpp:                                                            *
 *    functions calc_dose:  dose calculation algorithm                       *
 *              proc_pdose: parallel photon dose calculation process         *
 *              proc_edose: parallel electron dose calculation process       *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/06/28        *
 *                                                                           *
 *****************************************************************************/
 
// ****************************************
// includes
// ****************************************

#include <math.h>
#include <pthread.h>
#include "definitions.h"
#include "global.h"
#include "treatment_plan.h"
#include "e_beam_model.h"
#include "p_beam_model.h"
#include "multi_electron.h"
#include "ranmar.h"
#include "sobseq.h"
#include "portal_dose.h"

// ****************************************
// declare functions and global variables
// ****************************************

real etime(void);
void start_message(int, int);
real batch_message(int, int, real, ranmar_state);

// treatment plan parameters
extern treatment_plan plan;

// the accelerator head models
extern p_beam_model *plinac;    // photon beam
extern e_beam_model *elinac;    // electron beam

// simulate photon histories by KERMA approximation
void kerma_photon(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *, portal_dose *);

// simulate photon histories by multiple photon transport
void multi_photon(particle_parameters &, int, ranmar &,
#ifdef USE_SOBOL
                  real *,
#endif // USE_SOBOL
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *, portal_dose *);

// simulate electron history
void one_electron(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *, portal_dose *);

// lock shared memory for update (total_weight, beam_dose, beam_error)
pthread_mutex_t lock_dose;

// particle weight counter
volatile double total_weight;

// number of batches, repeats
int n_batch, n_repeat, n_rotate, beam_id;

// number of histories per batch
long n_history_per_batch;

// array of CPU times for each process
real *proc_time;

// portal dose matrix
portal_dose *portal;

// dose factor
real dose_factor;

// low energy correction for collision stopping power
real s_correction;

#ifdef PHSP_READ
real_3 iso_center;
long positronCount=0, negatronCount=0, photonCount=0;
long primaryCount=0, scatterCount=0;
real pETotal=0.0, eETotal=0.0;
real beamArea;
real col_cdistance = 44.999;
#endif

// ****************************************
// function proc_pdose
// ****************************************

void *
proc_pdose(void *arg)
{
   // initial process time in seconds
   real cpu_start = etime();

#ifdef USE_SOBOL
   real zeta[SOBOL_DIM]; // array to store Sobol random numbers
#endif

   // present process (thread) index
   register int i_process = *((int *) arg);

   // local particle weight counter
   double batch_weight = ZERO;

   // temporary dose matrix for present batch and process
   array_3d<double> *batch_dose;
   if ( (batch_dose = new array_3d<double>(dim, ZERO)) == NULL ) {
      xvmc_error("proc_pdose","cannot construct batch_dose",8); }

#ifdef CHECK_ENERGY
   sum_energy_type sum_energy; // test energy conservation
#endif // CHECK_ENERGY

   // variable to store the particle parameters
   particle_parameters particle;

   // initialize random sequence for the present process and beam
   ranmar rndm(plan.ini_rndm1 + i_process, plan.ini_rndm2 + beam_id,
               plan.ini_rndm3,             plan.ini_rndm4);

#ifdef USE_SOBOL
   // initalize Sobol' sequence for the present process
   sobseq sobol(SOBOL_DIM);
   int    sobol_count = 0;
   // each process takes a definite sample of the Sobol' sequence
   long ini_sobol = i_process*n_batch*n_history_per_batch/n_process;
   for (long i_p = 0; i_p < ini_sobol; ++i_p) sobol.next();
#endif // USE_SOBOL

#ifdef PHSP_WRITE
	// open a phase space file (phspWrite.bin) in /tmp
    FILE *ostrm = fopen("/home/jokim/XVMC/PHSP/phspWrite.bin","wb");
    if(ostrm==NULL) {
     	 printf("\n ERROR: opening file /home/jokim/XVMC/PHSP/phspWrite.bin");
		 exit(-1);
    }		
    /* Write Head */
    char ModeRW[20]="MODE2";
    fwrite(ModeRW, sizeof(char), 5, ostrm);
    if (strncmp(ModeRW, "MODE2", 5) != 0) {
       printf("ERROR: %s is unrecognized MODE\n", ModeRW);
      exit(-1);
    }	

    int nTotal = 0;
    int nPhoton = 0;
    float mxPhotEng = 0.0;  // Maximum Photon Energy
    float mnElectEng = 100.0; // Minimum Electron Energy
    float nSrcPart = 0.0;   // Later should be changed with total weight

    fwrite(&nTotal, sizeof(int), 1, ostrm);
    fwrite(&nPhoton, sizeof(int), 1, ostrm);

    float floatArray[20];
    floatArray[0] = mxPhotEng;
    floatArray[1] = mnElectEng;
    floatArray[2] = nSrcPart;

    fwrite(floatArray, sizeof(float), 3, ostrm);

    char charArray[200];
    charArray[7] ='\0';
    fwrite(charArray, sizeof(char), 7, ostrm); /* padding */ 
    /* Enf writing head */
#endif  // PHSP_WRITE

#ifdef PHSP_READ
    // For Multi-core processors
    char sPhspFile[256];
    sprintf(sPhspFile, "phspRead%d.bin", i_process);

    col_cdistance = plinac->get_col_C_distance();
    printf ("PHSP_READ> Open Phase Space File at %f cm: %s for CPU %d\n", col_cdistance, sPhspFile, i_process); 
    
    FILE *istrm = fopen(sPhspFile,"rb"); /* input stream */
    if(istrm==NULL) {
       printf("\n ERROR: opening input %s \n", sPhspFile);
		 exit(-1);
    }
    rewind (istrm);
    char ModeRW[20];
    if (fread(ModeRW, sizeof(char), 5, istrm) == 5) {
       ModeRW[5] = '\0';
    }
    else {
       printf ("ERROR: Can not read MODE_RW\n");
       exit(-1);
    }

    int intArray[10];
    if(fread(intArray, sizeof(int), 2, istrm)!=2){
       printf("ERROR: reading intArray from header");
       exit(-1);
    }

    int nTotal = intArray[0];
    int nPhoton = intArray[1];

    float floatArray[20];
    if (fread(floatArray, sizeof(float), 3, istrm)!=3) {
       printf("ERROR: reading floatArray from header");
       exit(-1);
    }

    float mxPhotEng = floatArray[0];
    float mnElectEng = floatArray[1];
    float nSrcPart = floatArray[2];

    // for Temporary Debug 
    printf ("PHSP_READ> nTotal = %d nPhoton= %d mxPhotEng = %E  mnElectEng = %E  nSrcPart = %E\n",nTotal, nPhoton, mxPhotEng,mnElectEng,nSrcPart);

    /* padding */
    char charArray[200];
    if (fread(charArray, sizeof(char), 7, istrm)!=7) {
       printf("ERROR: reading mode from header\n");
       exit(-1);
    }
    charArray[7] = '\0';
    long nPhsp = 0;
    int  nCycl = 0;
#endif // PHSP_READ

   // start of batch loop
   for (int i_batch = i_process; i_batch < n_batch; i_batch += n_process)
   {
      
      // initialize the dose array for this batch
      batch_weight = ZERO;
      // *batch_dose  = ZERO;
      memset((void *) batch_dose->data, 0, batch_dose->size*sizeof(double));
    
      // initialize portal dose image for this batch
      if (portal != NULL) portal->reset_batch_dose_portal(n_repeat);

#ifdef CHECK_ENERGY
      // test energy conservation
      sum_energy.init();
#endif

      // start of history per batch loop
      for (register long n = 0; n < n_history_per_batch; ++n)
      {
         // a dynamic MLC must change the leaf positions from time to time,
         // in static mode this function does nothing
         plinac->modifier->adjust(n,n_history_per_batch);

#ifdef USE_SOBOL
         // generate the next set of quasi random numbers for the
         // present process, set counter to 0
         sobol.next(); sobol_count = 0;
#endif // USE_SOBOL

#ifndef PHSP_READ
         // get particle from linac
         if ( plinac->emit(particle,
#ifdef USE_SOBOL
                           sobol, sobol_count,
#endif // USE_SOBOL
                           rndm) )
#else

        // printf("PHSP_READ> Isocenter: x = %f  y = %f  z = %f\n", iso_center.x, iso_center.y, iso_center.z); 

	int twoToThe29th = 536870912 ; //  512 * 1024 * 1024;
	int twoToThe30th = 1073741824 ; //  1024 * 1024 * 1024;

	if(nPhsp == nTotal-1) {
	    rewind (istrm);
	    if (fread(ModeRW, sizeof(char), 5, istrm) == 5) {
	       ModeRW[5] = '\0';
	    }
	    else {
	       printf ("ERROR: Can not read MODE_RW\n");
	       exit(-1);
	    }

	    if(fread(intArray, sizeof(int), 2, istrm)!=2){
	       printf("ERROR: reading intArray from header");
	       exit(-1);
	    }
	    
	    if (fread(floatArray, sizeof(float), 3, istrm)!=3) {
	       printf("ERROR: reading floatArray from header");
	       exit(-1);
	    }

	    /* padding */
	    if (fread(charArray, sizeof(char), 7, istrm)!=7) {
	       printf("ERROR: reading mode from header\n");
	       exit(-1);
	    }
	    charArray[7] = '\0';
	    
	    nCycl++;
	    // printf("--- Rewind Phsp File at nTotal = %d x %d\n", nTotal, nCycl);	  
	    nPhsp = 0;
	}
	
       float phaseArray[20];
	    int LATCH;
       if( fread(&LATCH, sizeof(int),   1, istrm) != 1){
           printf("ERROR:  Failed to read Latch\n");
           exit (-1);
       }	
	    if( fread(phaseArray, sizeof(float), 7, istrm) != 7){
           printf("ERROR: Failed to read Data\n");
           exit (-1);
       }
       nPhsp ++;
       particle.energy = fabs(phaseArray[0]);
       particle.pos.x  = phaseArray[1];
       particle.pos.y  = phaseArray[2];
       particle.dir.x  = phaseArray[3];
       particle.dir.y  = phaseArray[4];
	    particle.weight = phaseArray[5];
       particle.pos.z  = phaseArray[6];  // Last Interaction Z Position
       /*
	    printf ("%8d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
	    LATCH, phaseArray[0], phaseArray[1],phaseArray[2],
	    phaseArray[3],phaseArray[4],phaseArray[5],phaseArray[6]);
        */
   	 /* Determine the charge of the particle: */
   	 /* particleIsPositron = (egs4Latch % twoToThe29th) / twoToThe28th;
      	particleIsNegatron = (egs4Latch % twoToThe30th) / twoToThe29th; */

   	int particleIsNegatron = (((LATCH/twoToThe30th) % 2) == 1);
   	int particleIsPositron = (((LATCH/twoToThe29th) % 2) == 1);
   	if (particleIsPositron && particleIsNegatron) {
          printf("ERROR: Particle with LATCH = %d is both a "
                 "positron and a negatron.\n",  LATCH);
          printf("Energy = %f, X = %f\n", particle.energy, particle.pos.x);
   	}

	   int particleCharge = 0; // default is photon
   	if(particleIsNegatron) {
         particle.type = ELECTRON;
	      particle.energy -= 0.511034; // Only Kinetic Energy
	      eETotal += particle.energy;
	      negatronCount++;
   	}
   	else if(particleIsPositron)  {
         particle.type = POSITRON;  // position 
	      particle.energy -= 0.511034; // Only Kinetic Energy
	      eETotal += particle.energy;
	      positronCount++;
   	}
   	else {
         particle.type = PHOTON;  // photon
	      pETotal += particle.energy;
	      photonCount++;
	      if(particle.pos.z < 20.0) primaryCount++;
	      else scatterCount++;
   	}

	particle.dir.z = 1.0;
	real RT = particle.dir.x*particle.dir.x+particle.dir.y*particle.dir.y;
	if (RT < 1.0){ 
	  particle.dir.z = sqrt(1.0-RT);
	}
	if(particle.weight < 0.0) {
	   particle.dir.z *= -1;
	   particle.weight = fabs(particle.weight);
	}
/*
	printf ("  %d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %4d %4d %4d\n", 
	    particle.type, particle.energy, particle.weight,
	    particle.pos.x, particle.pos.y, particle.pos.z,
	    particle.dir.x, particle.dir.y, particle.dir.z,
	    particle.i.x, particle.i.y, particle.i.z);
*/	

#define TRACE ON  
#ifndef TRACE
   //printf("\n Read a particle");
   // particle.pos.z = 45.0;

   float phspPlane = col_cdistance; // Scoring Plane in Previous Step (it was 45.0)
		
	float R = 1.0e+10;  // Infinity 
	if (particle.dir.z != 0.0) {
      R = (100.0-iso_center.z-phspPlane)/particle.dir.z;
	   particle.pos.x = iso_center.x + particle.pos.x + particle.dir.x * R;
	   particle.pos.y = iso_center.y + particle.pos.y + particle.dir.y * R;
	   particle.pos.z = 0.0; // Set to 0.0 always for xVMC
   }

	particle.i.x = int(particle.pos.x/voxel_size.x);
	particle.i.y = int(particle.pos.y/voxel_size.y);
	particle.i.z = 0; 
#else
   particle.pos.z = col_cdistance; // Phase Space Position
   real origin_point_dist = 1.0e10; // Infinity
   if (particle.dir.z != 0.0) 
      origin_point_dist = sqrt(particle.pos.x*particle.pos.x
                              +particle.pos.y*particle.pos.y
                              +particle.pos.z*particle.pos.z);

   if ( plinac->modifier->transport(particle,rndm) )
   {
/*
	printf ("m %8d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %4d %4d %4d\n", 
	     particle.type,  particle.energy, particle.weight,
	     particle.pos.x, particle.pos.y,  particle.pos.z,
	     particle.dir.x, particle.dir.y,  particle.dir.z,
	     particle.i.x,   particle.i.y,    particle.i.z);
*/
   }
#endif
      // trace particle to the simulation grid (take beam angles into account)
      // return true if the particle hits the calculation cube

      if(plinac->trace2cube(particle.pos, particle.dir, particle.i, origin_point_dist, rndm))
#endif // PHSP_READ
         {
            // count weight,
            // just one because we want to normalize by monitor units
            batch_weight += ONE;

/*
	printf ("-> %d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %4d %4d %4d\n", 
	     particle.type, particle.energy, particle.weight,
	     particle.pos.x, particle.pos.y, particle.pos.z,
	     particle.dir.x, particle.dir.y, particle.dir.z,
	     particle.i.x,   particle.i.y,   particle.i.z);
*/
#ifdef CHECK_ENERGY
            // test energy conservation
            sum_energy.start += particle.energy*particle.weight;
#endif // CHECK_ENERGY

#ifdef USE_SOBOL
            // generate array of Sobol random numbers
            for (register int is=0; is<SOBOL_DIM; ++is)
            {
               if (sobol_count < SOBOL_DIM)
               {
                  zeta[is] = sobol.number(sobol_count); ++sobol_count;
               }
               else
               {
                  zeta[is] = rndm.number();
               }
            }
#endif // USE_SOBOL

#ifdef PHSP_WRITE  // --------------------------------------------------------------
	     float phspPlane = 45.0; // Phase Space Position

	     nTotal++;  // counts total number of particles
	     nSrcPart += particle.weight;  // counts total particle weight (initial source)
/*
	     printf ("%8d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %8d %8d %8d\n", 
		      particle.type, particle.energy, particle.weight,
				     particle.pos.x, particle.pos.y, particle.pos.z,
				     particle.dir.x, particle.dir.y, particle.dir.z,
				     particle.i.x, particle.i.y, particle.i.z);
*/				

	     if (particle.type == PHOTON) {
		     nPhoton++;
		     if (particle.energy > mxPhotEng) mxPhotEng = particle.energy;
	     }
	     else {
		     if (particle.energy < mnElectEng) mnElectEng = particle.energy;
	     }

	     floatArray[0] = particle.energy;  			// particle energy
  	     floatArray[1] = particle.pos.x;				// particle x position
  	     floatArray[2] = particle.pos.y;   			// particle y position
  	     floatArray[3] = particle.dir.x;				// particle x direction
  	     floatArray[4] = particle.dir.y;				// particle y direction
	     floatArray[5] = particle.weight;          // particle weight
	     floatArray[6] = phspPlane;  					// particle z position (New Position)
	     // floatArray[6] = particle.pos.z;  					// particle z position

	     float R = 0.0;
	     if (particle.dir.z != 0.0) {
		     R = (100.0-phspPlane)/particle.dir.z;
		     floatArray[1] = (particle.pos.x-iso_center.x) - particle.dir.x * R;
		     floatArray[2] = (particle.pos.y-iso_center.y) - particle.dir.y * R;
		     // floatArray[6] = phspPlane - particle.pos.z;
	     }

	     int LATCH = 0;  // for Photon but will be changed with LATCH later
  	     if(fwrite(&LATCH, sizeof(int),   1, ostrm) != 1){
    	     printf("ERROR: Failed to write LATCH\n");
    	     exit(-1);
  	     }

	     if(fwrite(floatArray, sizeof(float), 7, ostrm) != 7){
     	     printf("ERROR: Failed to write phsp data\n");
     	     exit (-1);
	     }

	     goto EndOfBatchLoop; // jump to Enf of Batch Loop
#endif // PHSP_WRITE  ------------------------------------------------------------


            if (particle.type == PHOTON)
            {
               // simulate photon history, with or without energy check
               if (particle.energy <= k0_cut)
               {
                  // by KERMA approximation
                  kerma_photon(particle, rndm,
#ifdef CHECK_ENERGY
                               sum_energy,
#endif // CHECK_ENERGY
                               batch_dose, portal);
               }
               else
               {
                  // or by multiple photon transport (default)
                  multi_photon(particle, n_repeat, rndm,
#ifdef USE_SOBOL
                               zeta,
#endif // USE_SOBOL
#ifdef CHECK_ENERGY
                               sum_energy,
#endif // CHECK_ENERGY
                               batch_dose, portal);
               }

            }
            else
            {
               // this is a charged particle
               one_electron(particle, rndm,
#ifdef CHECK_ENERGY
                            sum_energy,
#endif // CHECK_ENERGY
                            batch_dose, portal);
            }
#ifdef PHSP_WRITE
				EndOfBatchLoop:
				continue;
#endif // PHSP_WRITE
         }
         else
         {
            // count weight,
            // just one because we want to normalize by monitor units
            batch_weight += ONE;

#ifdef CHECK_ENERGY
            // test energy conservation
            sum_energy.start += particle.energy*particle.weight;
            if (particle.type == PHOTON)
            {
               sum_energy.ploss += particle.energy*particle.weight;
            }
            else
            {
               sum_energy.eloss += particle.energy*particle.weight;
            }
#endif // CHECK_ENERGY

         }
      }
      // end of history per batch loop
#ifdef PHSP_WRITE
				goto EndOfBatch;
#endif // PHSP_WRITE
      // lock shared memory for update
      pthread_mutex_lock(&lock_dose);

         // count total weight
         total_weight += batch_weight;

         // put "batch_dose->matrix" into "beam_dose->matrix"
         real rho,smw,dose;
         for (register int i=0; i<dim.x; ++i)
         {
            for (register int j=0; j<dim.y; ++j)
            {
               for (register int k=0; k<dim.z; ++k)
               {
                  // mass density
                  rho   = density->matrix[i][j][k];
                  switch (plan.result)
                  {
                  case DOSE:
#ifdef MONITORBSF
                  case ABS_DOSE_MBSF:	  
                     // energy dose + monitor backscatter correction
#endif	
                  case ABS_DOSE:	  
                     // energy dose
                     smw = rho;
                     break;		     
                  default:
                     // dose times stopping power ratio s_medium(rho)/s_water
                     smw   = dens_scol->matrix[i][j][k]
                           - dens_ccol->matrix[i][j][k]*s_correction;
                     break;
                  }
                  dose  = batch_dose->matrix[i][j][k]*dose_factor/batch_weight;
#ifdef PHSP_READ
                  dose *= (1.0*nTotal/nSrcPart/beamArea);
#endif	
                  // dose in vacuum is 0
                  const real rho_min = 0.0001;
                  dose = rho > rho_min ? dose/smw : ZERO;
                  beam_dose->matrix[i][j][k]  += dose;
                  beam_error->matrix[i][j][k] += dose*dose;
               }
            }
         }

         // update portal image dose
         if (portal != NULL) portal->update_batch(
             dose_factor*voxel_size.x*voxel_size.y*voxel_size.z, batch_weight);

         // check energy conservation
#ifdef CHECK_ENERGY
         tot_energy.start += sum_energy.start;
         tot_energy.edepo += sum_energy.edepo;
         tot_energy.pdepo += sum_energy.pdepo;
         tot_energy.eloss += sum_energy.eloss;
         tot_energy.ploss += sum_energy.ploss;
#endif // CHECK_ENERGY

      // unlock shared memory
      pthread_mutex_unlock(&lock_dose);

      proc_time[i_process] =
         batch_message(i_batch,i_process,cpu_start,rndm.status());

#ifdef PHSP_WRITE
		EndOfBatch:
		continue;
#endif // PHSP_WRITE

   }
   // end of batch loop

#ifdef PHSP_WRITE
	rewind(ostrm);	
	/* Write Head */
	fwrite(ModeRW, sizeof(char), 5, ostrm);
  	if (strncmp(ModeRW, "MODE2", 5) != 0) {
      printf("ERROR: %s is unrecognized MODE\n", ModeRW);
     	exit(-1);
  	}	
	fwrite(&nTotal, sizeof(int), 1, ostrm);
	fwrite(&nPhoton, sizeof(int), 1, ostrm);
	
	floatArray[0] = mxPhotEng;
	floatArray[1] = mnElectEng;
	floatArray[2] = nSrcPart;

	// for Temporary Debug 
//	printf ("nTotal = %d nPhoton= %d mxPhotEng = %E  mnElectEng = %E  nSrcPart = %E\n",nTotal, nPhoton, mxPhotEng,mnElectEng,nSrcPart);

	fwrite(floatArray, sizeof(float), 3, ostrm);
	rewind(ostrm);	

#endif // PHSWRITE

#ifdef PHSP_READ
   printf ("PHSP_READ> Close %s\n", sPhspFile);
   printf ("PHSP_READ> No of Photons = %ld (%5.2f)  No of Negatrons = %ld (%5.2f)  No of Positron = %ld (%5.2f)\n",
            photonCount, 100.0*photonCount/(photonCount+negatronCount+positronCount),
            negatronCount, 100.0*negatronCount/(photonCount+negatronCount+positronCount),
	    positronCount, 100.0*positronCount/(photonCount+negatronCount+positronCount));  
   printf ("PHSP_READ> Primary Photon = %5.2f  Scatter Photon = %5.2f\n", 
            100.0*primaryCount/(primaryCount+scatterCount), 100.0*scatterCount/(primaryCount+scatterCount));
   printf ("PHSP_READ> Avg Photon Energy = %5.2f (MeV)  Avg Electron Energy = %5.2f (MeV)\n",
            pETotal/photonCount, eETotal/(negatronCount+positronCount));
   printf ("PHSP_READ> %d Times Phase Space Cycled\n",nCycl);
   fclose(istrm);
#endif
   
   // free temporary dose matrix memory
   delete batch_dose; batch_dose = NULL;

   return(NULL);
}

// ****************************************
// function proc_edose
// ****************************************

void *
proc_edose(void *arg)
{
   // initial process time in seconds
   real cpu_start = etime();

   // present process (thread) index
   register int i_process = *((int *) arg);

   // local particle weight counter
   double batch_weight = ZERO;

   // temporary dose matrix for present batch and process
   array_3d<double> *batch_dose;
   if ( (batch_dose = new array_3d<double>(dim, ZERO)) == NULL ) {
      xvmc_error("proc_edose","cannot construct batch_dose",8); }

#ifdef CHECK_ENERGY
   sum_energy_type sum_energy; // test energy conservation
#endif // CHECK_ENERGY

   // variable to store the particle parameters
   particle_parameters particle;

   // save particle type and energy
   particle_type p_type;
   real          p_energy;

   // set this flag to indicate that the present particle (electron)
   // is a particle from the primary source (scatter foils)
   bool          primary_particle;

   // create and trace electron histories (history repetition)
   multi_electron electron;

   // initialize random sequence for the present process and beam
   ranmar rndm(plan.ini_rndm1 + i_process, plan.ini_rndm2 + beam_id,
               plan.ini_rndm3,             plan.ini_rndm4);

#ifdef USE_SOBOL
   // initalize Sobol' sequence for the present process
   sobseq sobol(SOBOL_DIM);
   int    sobol_count = 0;
   // each process takes a definite sample of the Sobol' sequence
   long ini_sobol = i_process*n_batch*n_history_per_batch/n_process;
   for (long i_p = 0; i_p < ini_sobol; ++i_p) sobol.next();
#endif // USE_SOBOL

   // start of batch loop
   for (int i_batch = i_process; i_batch < n_batch; i_batch += n_process)
   {
      // initialize the dose array for the present batch
      batch_weight = ZERO;
      // *batch_dose  = ZERO;
      memset((void *) batch_dose->data, 0, batch_dose->size*sizeof(double));

#ifdef CHECK_ENERGY
      // test energy conservation
      sum_energy.init();
#endif

      // start of history per batch loop
      for (register long n = 0; n < n_history_per_batch; ++n)
      {

         // get particle type, particle energy and set flag if the
         // particle (electron) is from the primary source
         elinac->emit(p_type,p_energy,primary_particle,rndm);

         // createelectron history for this energy
         multi_electron_create(electron, p_energy, rndm);

         // start of additional repetitions for gantry rotation
         for (register int i_rotate = 0; i_rotate < n_rotate; ++i_rotate)
         {

            // start of history repetitions
            // (different starting positions and directions)
            for (register int i_repeat = 0; i_repeat < n_repeat; ++i_repeat)
            {

#ifdef USE_SOBOL
               // generate the next set of quasi random numbers for the
               // present process, set counter to 0
               sobol.next(); sobol_count = 0;
#endif // USE_SOBOL

               // get weight, starting position, direction and voxel index
               // at the calculation cube surface
               if ( elinac->emit(particle.weight, particle.pos,
                                 particle.dir,    particle.i,
                                 primary_particle,
#ifdef USE_SOBOL
                                 sobol, sobol_count,
#endif // USE_SOBOL
                                 rndm) )
               {
                  // update particle type and energy
                  particle.type   = p_type;
                  particle.energy = p_energy;

                  // count weight
                  batch_weight += particle.weight;

#ifdef CHECK_ENERGY
                  // test energy conservation
                  sum_energy.start += particle.energy*particle.weight;
#endif // CHECK_ENERGY

                  // trace electron history through the calculation cube
                  multi_electron_trace(electron, particle,rndm,
#ifdef CHECK_ENERGY
                                 sum_energy,
#endif // CHECK_ENERGY
                                 batch_dose, NULL);
               }
               else
               {
                  // count weight
                  batch_weight += particle.weight;

#ifdef CHECK_ENERGY
                  // test energy conservation
                  sum_energy.start += p_energy*particle.weight;
                  sum_energy.eloss += p_energy*particle.weight;
#endif // CHECK_ENERGY

               }

            }
            // end of history repetitions

         }
         // end of additional repetitions for gantry rotation

      }
      // end of history per batch loop

      // lock shared memory for update
      pthread_mutex_lock(&lock_dose);

         // count total weight
         total_weight += batch_weight;

         // put "batch_dose->matrix" into "beam_dose->matrix"
         real rho,smw,dose;
         for (register int i=0; i<dim.x; ++i)
         {
            for (register int j=0; j<dim.y; ++j)
            {
               for (register int k=0; k<dim.z; ++k)
               {
                  // mass density
                  rho   = density->matrix[i][j][k];
                  switch (plan.result)
                  {
                  case DOSE:
                  case ABS_DOSE:
                     // energy dose
                     smw = rho;
                     break;
                  default:
                     // dose times stopping power ratio s_medium(rho)/s_water
                     smw   = dens_scol->matrix[i][j][k];
                     break;;
                  }
                  dose  = batch_dose->matrix[i][j][k]*dose_factor/batch_weight;
                  // dose in vacuum is 0
                  const real rho_min = 0.0001;
                  dose = rho > rho_min ? dose/smw : ZERO;
                  beam_dose->matrix[i][j][k]  += dose;
                  beam_error->matrix[i][j][k] += dose*dose;
               }
            }
         }

         // check energy conservation
#ifdef CHECK_ENERGY
         tot_energy.start += sum_energy.start;
         tot_energy.edepo += sum_energy.edepo;
         tot_energy.pdepo += sum_energy.pdepo;
         tot_energy.eloss += sum_energy.eloss;
         tot_energy.ploss += sum_energy.ploss;
#endif // CHECK_ENERGY

      // unlock shared memory
      pthread_mutex_unlock(&lock_dose);

      proc_time[i_process] =
         batch_message(i_batch,i_process,cpu_start,rndm.status());

   }
   // end of batch loop

   // free temporary dose matrix memory
   delete batch_dose; batch_dose = NULL;

   return(NULL);
}

// ****************************************
// function calc_dose
// ****************************************

void calc_dose(beam_core *this_beam)
{
   int             i_process;   // process (thread) index
   int            *p_ind;       // pointer to array of process (thread) indices
   pthread_t      *thread;      // pointer to array of processes (threads)
   void           *retval;      // process return value
   pthread_attr_t  thread_attr; // process (thread) attributes
#ifdef PHSP_READ
   iso_center.x = this_beam->iso_center.x;
   iso_center.y = this_beam->iso_center.y;
   iso_center.z = this_beam->iso_center.z;
   printf("PHSP_READ> Isocenter: x = %f  y = %f  z = %f\n", 
           iso_center.x, iso_center.y, iso_center.z); 
#endif
   xvmc_message("==============================================",1);
   if (this_beam->type == PHOTON)
   {
      xvmc_message("= Start simulation for beam:",this_beam->id,
                   "(photon beam)",0);
   }
   else
   {
      xvmc_message("= Start simulation for beam:",this_beam->id,
                   "(electron beam)",0);
   }
   xvmc_message("==============================================",0);

   // allocate memory for process array
   if ( (thread = new pthread_t[n_process]) == NULL)
   {
      xvmc_error("calc_dose","cannot allocate memory for thread array",8);
   }

   // allocate memory for process index array
   if ( (p_ind = new int[n_process]) == NULL)
   {
      xvmc_error("calc_dose",
                 "cannot allocate memory for thread index array",8);
   }

   // initialize process index array
   for (i_process = 0; i_process<n_process; ++i_process)
      p_ind[i_process] = i_process;

   // allocate memory for process CPU time array
   if ( (proc_time = new real[n_process]) == NULL)
   {
      xvmc_error("calc_dose",
                 "cannot allocate memory for process CPU time array",8);
   }

   // initialize process CPU time array
   for (i_process = 0; i_process<n_process; ++i_process)
      proc_time[i_process] = ZERO;

   // initialize thread attributes
   if (pthread_attr_init(&thread_attr))
   {
      xvmc_error("calc_dose","cannot initialize thread attributes",8);
   }

#ifdef OSF1
   // set process (thread) stack size to size of batch_dose
   // plus 2 Mbytes (DEC Unix only)
   size_t stacksize = dim.x*dim.y*dim.z*sizeof(double) + 2097152;
   if (pthread_attr_setstacksize(&thread_attr,stacksize))
   {
      xvmc_error("calc_dose","cannot set stack size",8);
   }
#endif

   // volume of one voxel
   real voxel_volume = voxel_size.x*voxel_size.y*voxel_size.z;

   // dose factor
   // beam area to calculate dose per fluence
   real beam_area = fabs(this_beam->open_x2-this_beam->open_x1)*
                    fabs(this_beam->open_y2-this_beam->open_y1);
#ifdef PHSP_READ
   beamArea = beam_area;
#endif    		    
   const real beam_area_min = 0.0001;
   if (beam_area > beam_area_min)
   {
      dose_factor = beam_area;
   }
   else
   {
      // for pencil beams we calculate dose per number of particles
      dose_factor = ONE;
   }

   // convert to 10^(-10) Gy
   dose_factor *= 1.602;

   // divide by voxel volume
   dose_factor /= voxel_volume;

   // divide by number of batches
   dose_factor /= double(this_beam->n_batch);

#ifdef MONITORBSF
   real openX1 = fabs(this_beam->open_x1);
   real openX2 = fabs(this_beam->open_x2);
   real openY1 = fabs(this_beam->open_y1);
   real openY2 = fabs(this_beam->open_y2); 
   //printf ("MONITORBSF> Field Size: X1=%4.2f  X2=%4.2f  Y1=%4.2f  Y2=%4.2f\n", 
   //         openX1, openX2, openY1, openY2);   
   // Calculate Monitor Backscatter Factor, MBSF
   real rY1 = 1.54 - 8.45E-2 * openY1 + 4.47E-5 * openY1*openY1*openY1;
   real rX1 = 0.40 - 1.87E-2 * openX1;
   real pY1 = 3.95E-2*openY1 - 3.55E-5 * openY1*openY1*openY1;
   //printf ("MONITORBSF> rY1 rX1 pY1 = %e %e %e\n", rY1, rX1, pY1);
   real rY2 = 1.54 - 8.45E-2 * openY2 + 4.47E-5 * openY2*openY2*openY2;
   real rX2 = 0.40 - 1.87E-2 * openX2;
   real pY2 = 3.95E-2*openY2 - 3.55E-5 * openY2*openY2*openY2;
   //printf ("MONITORBSF> rY2 rX2 pY2 = %e %e %e\n", rY2, rX2, pY2);
   real rY = rY1 + rY2;
   real rX = (rX1 + rX2)*(pY1 + pY2);
   real rFS = rX + rY;
   //printf ("MONITORBSF> rY rX rFS = %e %e %e\n", rY, rX, rFS);

   openX1 = 5.0;
   openX2 = 5.0;
   openY1 = 5.0;
   openY2 = 5.0; 
   // Calculate Monitor Backscatter Factor, MBSF
   rY1 = 1.54 - 8.45E-2 * openY1 + 4.47E-5 * openY1*openY1*openY1;
   rX1 = 0.40 - 1.87E-2 * openX1;
   pY1 = 3.95E-2*openY1 - 3.55E-5 * openY1*openY1*openY1;
   //printf ("MONITORBSF> rY1 rX1 pY1 = %e %e %e\n", rY1, rX1, pY1);   
   rY2 = 1.54 - 8.45E-2 * openY2 + 4.47E-5 * openY2*openY2*openY2;
   rX2 = 0.40 - 1.87E-2 * openX2;
   pY2 = 3.95E-2 *openY2- 3.55E-5 * openY2*openY2*openY2;   
   //printf ("MONITORBSF> rY2 rX2 pY2 = %e %e %e\n", rY2, rX2, pY2);
   rY = rY1 + rY2;
   rX = (rX1 + rX2)*(pY1 + pY2);
   real rREF = rX + rY;
   //printf ("MONITORBSF> rY rX rREF = %e %e %e\n", rY, rX, rREF);   
   
   real MBSF = (1.0 + rREF/100.)/(1.0 + rFS/100.0);
   // printf (DEBUG> Dose Type = %d \n", plan.result);
   if (plan.result == ABS_DOSE_MBSF && this_beam->type == PHOTON) {
      dose_factor *= (MBSF*MBSF);
      printf ("MONITORBSF> Monitor Backscatter Factor = %10.5f\n", MBSF*MBSF);
   }
#endif 

   // the following is necessary to convert dose to medium into dose to water
   if (this_beam->type == PHOTON)
   {
      // estimate secondary electron energy
      real ee = plinac->average_energy()/2.5;
      // and calculate low energy correction for collision stopping power
      if (ee > ONE) s_correction = exp(-(ee-ONE));
      else          s_correction = ONE;
   }
   else
   {
      s_correction = ZERO;
   }

   // portal dose image
   portal = this_beam->portal;

   // number of batches, repetitions, histories per batch, etc.
   beam_id             = this_beam->id;
   n_batch             = this_beam->n_batch;
   n_repeat            = this_beam->n_repeat;
   n_rotate            = this_beam->n_rotate;
   n_history_per_batch = this_beam->n_history/n_batch;

#ifdef PHSP_READ_REMOVED
    char sPhspFile[256];
    sprintf(sPhspFile, "phspRead%d.bin", i_process);
    // open a phase space file (phspRead.bin) in /tmp
    // for Temporary Debug 
    printf ("PHSP_READ> Open Phase Space File in Dose_Calc: %s for CPU %d\n", sPhspFile, i_process); 
    
    FILE *istrm = fopen(sPhspFile,"rb"); /* input stream */
    if(istrm==NULL) {
       printf("\n ERROR: opening input %s \n", sPhspFile);
		 exit(-1);
    }
    rewind (istrm);
    char ModeRW[20];
    if (fread(ModeRW, sizeof(char), 5, istrm) == 5) {
       ModeRW[5] = '\0';
    }
    else {
       printf ("ERROR: Can not read MODE_RW\n");
       exit(-1);
    }

    int intArray[10];
    if(fread(intArray, sizeof(int), 2, istrm)!=2){
       printf("ERROR: reading intArray from header");
       exit(-1);
    }

    int nTotal = intArray[0];
    int nPhoton = intArray[1];

    float floatArray[20];
    if (fread(floatArray, sizeof(float), 3, istrm)!=3) {
       printf("ERROR: reading floatArray from header");
       exit(-1);
    }

    float mxPhotEng = floatArray[0];
    float mnElectEng = floatArray[1];
    float nSrcPart = floatArray[2];

    // for Temporary Debug 
    printf ("PHSP_READ> nTotal = %d nPhoton= %d mxPhotEng = %E  mnElectEng = %E  nSrcPart = %E\n",nTotal, nPhoton, mxPhotEng,mnElectEng,nSrcPart);

    /* padding */
    char charArray[200];
    if (fread(charArray, sizeof(char), 7, istrm)!=7) {
       printf("ERROR: reading mode from header\n");
       exit(-1);
    }
    charArray[7] = '\0';
    long nPhsp = 0;
    int  nCycl = 0;
    
    fclose(istrm);

    this_beam->n_history = nTotal;
    n_history_per_batch = this_beam->n_history/n_batch;    
#endif // PHSP_READ


   // reset total particle weight
   total_weight = ZERO;

   // make the threads (processes)
   for (i_process = 0; i_process < n_process; ++i_process)
   {
      // to guarantee that the correct process index is passed to the thread
      // we use "&p_ind[i_process]" instead of "&i_process", the value of
      // "i_process" can change before the process (proc_pdose or proc_edose)
      // starts
      if (this_beam->type == PHOTON)
      {
         if (pthread_create(&thread[i_process], &thread_attr,
                             proc_pdose, &p_ind[i_process]))
         {
            xvmc_error("calc_dose","cannot create thread for photon beam",8);
         }
      }
      else
      {
         if (pthread_create(&thread[i_process], &thread_attr,
                             proc_edose, &p_ind[i_process]))
         {
            xvmc_error("calc_dose","cannot create thread for electron beam",8);
         }
      }
   }

   // join (collapse) the threads (processes) and find maximum CPU time
   real cpu_time = ZERO;
   for (i_process = 0; i_process < n_process; ++i_process)
   {
      if (pthread_join(thread[i_process], &retval))
      {
         xvmc_error("calc_dose","cannot join threads",8);
      }
      if (proc_time[i_process] > cpu_time) cpu_time = proc_time[i_process];
   }

   // free thread memory
   delete [] thread;     thread    = NULL;
   delete [] p_ind;      p_ind     = NULL;
   delete [] proc_time;  proc_time = NULL;

   // get total CPU time for this beam
   this_beam->cpu_time = cpu_time;

#ifdef CHECK_ENERGY
   double total_energy = tot_energy.edepo + tot_energy.pdepo +
                         tot_energy.eloss + tot_energy.ploss;

   energy_message("total sampled energy:            ",tot_energy.start,"MeV",1);
   energy_message("total electron energy deposit:   ",tot_energy.edepo,"MeV",0);
   energy_message("total photon energy deposit:     ",tot_energy.pdepo,"MeV",0);
   energy_message("total electron energy loss:      ",tot_energy.eloss,"MeV",0);
   energy_message("total photon energy loss:        ",tot_energy.ploss,"MeV",0);
   energy_message("total energy deposit and loss:   ",total_energy,"MeV",0);

   tot_energy.start /= total_weight;
   tot_energy.edepo /= total_weight;
   tot_energy.pdepo /= total_weight;
   tot_energy.eloss /= total_weight;
   tot_energy.ploss /= total_weight;
   total_energy     /= total_weight;

   energy_message("average sampled energy:          ",tot_energy.start,"MeV",1);
   energy_message("average electron energy deposit: ",tot_energy.edepo,"MeV",0);
   energy_message("average photon energy deposit:   ",tot_energy.pdepo,"MeV",0);
   energy_message("average electron energy loss:    ",tot_energy.eloss,"MeV",0);
   energy_message("average photon energy loss:      ",tot_energy.ploss,"MeV",0);
   energy_message("average energy deposit and loss: ",total_energy,"MeV",0);
#endif // CHECK_ENERGY

   if (this_beam->type == ELECTRON)
   {
      if (elinac->add_photons(cpu_time,n_batch))
      {
         xvmc_message("CPU time to add the photon background:",cpu_time,"s",1);
         this_beam->cpu_time += cpu_time;
      }
      else
      {
         xvmc_message("Photon background for electron beam not available!",1);
      }
   }
}
