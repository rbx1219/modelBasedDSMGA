
 #include <math.h>
 #include <iostream>
 #include <fstream>
 #include <cstdlib>
 #include "statistics.h"
 #include "dsmga.h"
 #include "global.h"

 #define MAX_GEN 200

 using namespace std;


 int main (int argc, char *argv[]) {

   if (argc != 7)
     {
       printf ("GA ell numConvergence lower upper alpha k\n");
       return -1;
     }

   int ell = atoi (argv[1]);
   int numConvergence = atoi (argv[2]);	// problem size
   int lower = atoi(argv[3]);
   int upper = atoi(argv[4]);
   BTRAP_ALPHA = atof(argv[5]);
   TRAP_K = atoi(argv[6]);

   if ( ell % TRAP_K != 0 )
    {
        printf("ell mod K != 0 \n");
        return -1;
    }
    
   
   
   //int nInitial = (int) (0.5 * ell * log((double)ell) / log(2.71828));	// initial population size

   int nInitial = 10;
   if (lower > nInitial) 
       nInitial = lower;


   int j;

   Statistics st, st_last;
   Statistics st_generation;

   int left, right, middle;

   int populationSize = nInitial/2;
   bool foundOptima;
       
   if ( lower < 0 || upper < 0) 
   {

       if (SHOW_BISECTION) printf("Bisection phase 1\n");

       do 
       {

           populationSize *= 2;

           if (SHOW_BISECTION) printf("[%d]: ", populationSize);

           foundOptima = true;
           st_generation.reset();

           for ( j = 0;  j < numConvergence;  ++j) 
           {
               DSMGA ga(ell, populationSize, 4, 1, 0, MAX_GEN, -1);
               int generation = ga.doIt(false);


               if (!ga.foundOptima()) 
               {

                   foundOptima = false;
                   st_generation.reset();
                   
                   if (SHOW_BISECTION) {
                       printf("-");
                       fflush(NULL);
                   }
                   break;
               }

               if (SHOW_BISECTION) 
               {
                   printf("+");
                   st_generation.record(generation);
                   fflush(NULL);
               }
           }	

           if (SHOW_BISECTION)
               printf("generation = %f \n", st_generation.getMean());
               
           if (foundOptima)
               st_last = st_generation;

       } while (!foundOptima);

       left = populationSize/2;
       right = populationSize;
   }
   
   else 
   {
       left = lower;
       right = upper;
   }

   
   middle = (left + right)/2; 

   if (SHOW_BISECTION) printf("Bisection phase 2\n");

   while ((right > 1.05 * left) && right > left + 2) 
   {

       middle = (left + right) / 2;

       if (SHOW_BISECTION) printf("[%d]: ", middle);

       foundOptima = true;
       st_generation.reset();

       for (j=0; j<numConvergence; j++) 
       {

           DSMGA ga(ell, middle, 2, 1, 0, MAX_GEN, -1);

           int generation = ga.doIt(false);

           if (!ga.foundOptima()) {
               foundOptima = false;
               if (SHOW_BISECTION) {
                   printf("-");
                   st_generation.reset();
                   fflush(NULL);
               }
               break;
           }

           if (SHOW_BISECTION) {
               printf("+");
               st_generation.record(generation);
               fflush(NULL);
           }
       }

       if (foundOptima) 
           right = middle;
       else 
           left = middle;

       if (SHOW_BISECTION) 
           printf("  generation = %f\n", st_generation.getMean());
           
       if (foundOptima)
           st_last = st_generation;
   }



   middle = (left + right) / 2;

   printf("---------------------\n");
   printf("population = %d\n", middle);
   printf("generation = %f \n", st_last.getMean());
   printf("nfe = %f \n", middle*st_last.getMean() );
   printf("=====================\n");
   

   return EXIT_SUCCESS;

 }


