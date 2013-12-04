
#ifndef _DSMGA_H
#define _DSMGA_H

/**
@author Tian-Li Yu
*/

#include "global.h"
#include "chromosome.h"
#include "bbi.h"
#include "statistics.h"
#include "dsmclusteringchromosome.h"

class DSMGA
{

    public:
        DSMGA (int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
            double n_pm, int n_maxGen, int n_maxFe);

        ~DSMGA ();

        void init (int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
            double n_pm, int n_maxGen, int n_maxFe);

        void initializePopulation ();
        void evaluate ();

        void selection ();

        /** tournament selection without replacement*/
        void tournamentSelection ();

        void crossOver ();
        void crossOver (Chromosome &, Chromosome &, Chromosome &, Chromosome &);
        void bbUniformXO (Chromosome &, Chromosome &, Chromosome &, Chromosome &,
            double);

        void mutation ();
        void bbMutation (Chromosome &);

        /** get BBI via DSM clustering */
        void getBBI ();
        void getBBIOneMax ();
        void getBBIMKTrap ();

        void oneRun (bool output = true);
        int doIt (bool output = true);

        bool shouldTerminate ();
        int getNextPopulation ();

        bool foundOptima ();

        int getGeneration () const
        {
            return generation;
        }

        void replacePopulation ();
    	void fullReplace ();
        void DRTR ();
        void RTR ();

        void showStatistics ();

        int findInitialPopulation();
        bool findCorrectModel(int requiredMatchedM);

        bool isSteadyState ();

		static int generation;
		static int threshold;
		static int * best;
	protected:

        BBI bbi;                                  // BB Information
        bool bbiCalculated;

        int ell;                                  // chromosome length
        int nInitial;                             // initial population size
        int nCurrent;                             // current population size
        int selectionPressure;

        double pc;                                // prob of XO
        double pm;                                // prob of Mutation
        Chromosome *population;
        Chromosome *offspring;
        int *selectionIndex;
        int maxGen;
        int maxFe;
        int repeat;
        int fe;
//        int generation;
        int bestIndex;

        Statistics stFitness;

        int getDistance (Chromosome & c1, Chromosome & c2);
        int getBBDistance (Chromosome & c1, Chromosome & c2);

        double previousFitnessMean;

};
#endif
