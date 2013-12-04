
#include "dsmga.h"
#include "chromosome.h"
#include <iostream>

int DSMGA::generation = 0;
int DSMGA::threshold = 120000000;
int *DSMGA::best = NULL;

DSMGA::DSMGA (int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
		double n_pm, int n_maxGen, int n_maxFe) {
	fe = 0;
	previousFitnessMean = -INF;
	ell = n_ell;
	nInitial = n_nInitial;
	nCurrent = nInitial;
	selectionPressure = n_selectionPressure;
	pc = n_pc;
	pm = n_pm;
	maxGen = n_maxGen;
	maxFe = n_maxFe;
	bbiCalculated = false;

	population = new Chromosome[nInitial];
	offspring = new Chromosome[nInitial];

	for (int i=0; i<nInitial; i++) {
		population[i].init(ell);
		offspring[i].init(ell);
	}
	selectionIndex = new int[nInitial];

	bbi.init (ell);

	initializePopulation ();

}


DSMGA::~DSMGA () {
	delete[]population;
	delete[]offspring;
	delete[]selectionIndex;
}


void DSMGA::getBBI () {
	if (!bbiCalculated) {
		//getBBIOneMax ();
		//getBBIMKTrap ();



		DSMClusteringChromosome dsmcc;

		dsmcc.init (ell, 0.3333, 0.3333);
		//dsmcc.createDSMmu_enfe (population, nCurrent, selectionPressure, selectionIndex);
		//dsmcc.createDSMmu (population, nCurrent, selectionPressure, selectionIndex);
		dsmcc.createDSMCC (population, nCurrent, selectionPressure, selectionIndex);


		dsmcc.speedyHillClimbingV3 ();
		bbi = dsmcc.getBBI ();


	}
	bbiCalculated = true;
}


int DSMGA::getDistance (Chromosome & c1, Chromosome & c2) {

	int i;

	int distance = 0;

	for (i = 0; i < c1.lengthLong; i++)
		distance += myBD.getHammingDistance(c1.gene[i], c2.gene[i]);

	return distance;
}


int DSMGA::getBBDistance (Chromosome & c1, Chromosome & c2) {
	getBBI ();

	int i, j;

	int distance = 0;

	for (i = 0; i < bbi.bbNum; i++) {

		bool same = true;
		for (j = 1; j <= bbi.bb[i][0]; j++)
			if (c1.getVal(bbi.bb[i][j]) != c2.getVal(bbi.bb[i][j]))
				same = false;

		if (!same)
			distance++;
	}

	return distance;
}


void DSMGA::replacePopulation () {
	int i;

	if (SHOW_POPULATION) {
		printf ("===== Current population =====\n");
		for (i = 0; i < nCurrent; i++) {
			printf ("(%d)", i);
			if (SHORT_HAND)
				population[i].shortPrintOut ();
			else
				population[i].printOut ();
			printf ("\n");
		}
	}

	if (SHOW_SELECTION_INDEX) {
		printf ("===== Selection index =====\n");
		for (i = 0; i < nCurrent; i++)
			printf ("%d ", selectionIndex[i]);
		printf ("\n");
	}
	//fullReplace();
	RTR();
	//DRTR();
}


void DSMGA::fullReplace () {
	int i;

	for (i = 0; i < nCurrent; i++)
		population[i] = offspring[i];

	bbiCalculated = false;
}


/** DRTR with window size = N */
void DSMGA::DRTR () {

	int i, j;

	int windowSize = (nCurrent < ell/20)? nCurrent : ell/20;
	//int windowSize = nCurrent;
	int *randArray = new int[windowSize];

	for (i = 0; i < nCurrent; ++i) {
		int index;
		int distance;
		int minDistance = ell + 1;                // max distance
		double minFitness = INF;
		int index_f = 1;
		double f = INF;

		myRand.uniformArray (randArray, windowSize, 0, nCurrent - 1);

		for (j = 0; j < windowSize; ++j) {


			distance = getBBDistance (offspring[i], population[randArray[j]]);

			if (distance < minDistance) {
				index = randArray[j];
				minDistance = distance;
				minFitness = population[index].getFitness ();
			} else if (distance == minDistance
					&& minFitness > population[randArray[j]].getFitness ()) {
				index = randArray[j];
				minFitness = population[index].getFitness ();
			}

			if (f > population[randArray[j]].getFitness()) {
				index_f = randArray[j];
				f = population[randArray[j]].getFitness();
			}
		}

		if (minDistance >= 2) {
			index = index_f;
		}

		if (offspring[i].getFitness () > population[index].getFitness ()) {

			population[index] = offspring[i];

			if (SHOW_REPLACEMENT) {
				printf ("Replacing (%d) with ", index);
				if (SHORT_HAND)
					offspring[i].shortPrintOut ();
				else
					offspring[i].printOut ();
				printf ("\n");
			}
		} else if (SHOW_REPLACEMENT) {
			{
				printf ("Not replacing (%d) with ", index);
				if (SHORT_HAND)
					offspring[i].shortPrintOut ();
				else
					offspring[i].printOut ();

				printf ("\n");
			}

		}
	}

	bbiCalculated = false;

	delete []randArray;

}

/** RTR with window size = N */
void DSMGA::RTR () {

	int i, j;

	int windowSize = (nCurrent < ell/20)? nCurrent : ell/20;
	//int windowSize = nCurrent;
	int *randArray = new int[windowSize];

	for (i = 0; i < nCurrent; i++) {
		int index;
		int distance;
		int minDistance = ell + 1;                // max distance
		double minFitness = INF;

		myRand.uniformArray (randArray, windowSize, 0, nCurrent - 1);

		for (j = 0; j < windowSize; j++) {
			//distance = getBBDistance (offspring[i], population[randArray[j]]);
			distance = getDistance (offspring[i], population[randArray[j]]);
			if (distance < minDistance) {
				index = randArray[j];
				minDistance = distance;
				minFitness = population[index].getFitness ();
			} else if (distance == minDistance
					&& minFitness > population[randArray[j]].getFitness ()) {
				index = randArray[j];
				minFitness = population[index].getFitness ();
			}
		}

		if (offspring[i].getFitness () > population[index].getFitness ()) {

			population[index] = offspring[i];

			if (SHOW_REPLACEMENT) {
				printf ("Replacing (%d) with ", index);
				if (SHORT_HAND)
					offspring[i].shortPrintOut ();
				else
					offspring[i].printOut ();
				printf ("\n");
			}
		} else if (SHOW_REPLACEMENT) {
			{
				printf ("Not replacing (%d) with ", index);
				if (SHORT_HAND)
					offspring[i].shortPrintOut ();
				else
					offspring[i].printOut ();

				printf ("\n");
			}

		}
	}

	bbiCalculated = false;

	delete []randArray;

}


bool DSMGA::isSteadyState () {

	if (stFitness.getNumber () <= 0)
		return false;

	if (previousFitnessMean < stFitness.getMean ()) {
		previousFitnessMean = stFitness.getMean () + 1e-6;
		return false;
	}

	return true;
}


bool DSMGA::findCorrectModel(int requiredMatchedM) {
	getBBIMKTrap();
	BBI perfectBBI(ell);
	perfectBBI = bbi;

	bbiCalculated = false;
	selection();
	getBBI();

	int numofMatchedBB = perfectBBI.getNumofMatchedBB(bbi);

	if (numofMatchedBB >= requiredMatchedM) return true;
	else return false;
}


int DSMGA::findInitialPopulation () {

	int i,j;

	getBBIMKTrap();
	BBI perfectBBI(ell);
	perfectBBI = bbi;

	bbiCalculated = false;
	selection();
	getBBI();

	int oldPopulationSize;
	int populationSize = nCurrent;
	int numofMatchedBB;

	do {

		oldPopulationSize = populationSize;
		populationSize *= 2;

		Chromosome * newPopulation = new Chromosome[populationSize];

		for (i=0; i<oldPopulationSize; i++) {
			newPopulation[i].init(ell);
			newPopulation[i] = population[i];
		}

		for (i=oldPopulationSize; i<populationSize; i++) {

			for (j=0; j<ell; j++)
				if (myRand.flip())
					newPopulation[i].setVal(j,0);
				else
					newPopulation[i].setVal(j,1);
		}

		delete []population;
		delete []selectionIndex;

		population = newPopulation;
		selectionIndex = new int[populationSize];

		nCurrent = populationSize;

		bbiCalculated = false;
		selection();
		getBBI();

		numofMatchedBB = perfectBBI.getNumofMatchedBB(bbi);

		if (populationSize > 10000) exit(-1);

	} while ( numofMatchedBB < perfectBBI.bbNum - 1 );

	int left = populationSize / 2;
	int right = populationSize;
	int middle;

	do {
		middle = (left + right) / 2;
		nCurrent = middle;

		bbiCalculated = false;
		selection();

		getBBI();

		/*
		   bbi.printOut();
		   printf("\n%d\n=========\n", middle);
		   */
		numofMatchedBB = perfectBBI.getNumofMatchedBB(bbi);
		if (numofMatchedBB < perfectBBI.bbNum - 1)
			left = middle;
		else
			right = middle;
	} while (left + 2 < right);

	printf("\n");

	return middle;

}


int DSMGA::doIt (bool output) {
	generation = 0;

	while (!shouldTerminate ()) {
		oneRun (output);
	}
	return generation;
}


void DSMGA::oneRun (bool output) {
	int i;

	selection ();
	crossOver ();
	mutation ();
	replacePopulation ();

	//Building Schemata Info

	int BestSchemata[bbi.bbNum];

	int curBB , curPiece , curLocus;
	for(curBB = 0 ; curBB < bbi.bbNum ; curBB++)
	{
		int schemataTotal = pow( 2 , bbi.bb[curBB][0] ) ;
		double *candidate = new double[schemataTotal];
		int *schemataCounts = new int[schemataTotal];
		int maxIndex = 0;


		for(int j = 0 ; j < schemataTotal ; j++)
			candidate[j] = 0.0 , schemataCounts[j] = 0;


		for( curPiece = 0 ; curPiece < nCurrent ; curPiece++)
		{
			int count = 0;
			for(curLocus = 0 ; curLocus < bbi.bb[curBB][0] ; curLocus++)
				if(population[curPiece].getVal( bbi.bb[curBB][curLocus + 1] ) == 1)
					count += pow(2 , curLocus);
			candidate[count] += population[curPiece].getFitness();
			schemataCounts[count] ++ ;
		}

		double maxFit = candidate[0] / schemataCounts[0]; 
		for(int j = 1 ; j < schemataTotal ; j++)
			if(candidate[j] / schemataCounts[j] > maxFit)
			{
				maxFit = candidate[j] / schemataCounts[j] ; 
				maxIndex = j;
			}
		BestSchemata[curBB] = maxIndex;
		delete[] candidate;
		delete[] schemataCounts;
	}

	DSMGA::best = BestSchemata;	

	//

	double max = -INF;
	stFitness.reset ();
	for (i = 0; i < nCurrent; i++) {
		double fitness = population[i].getFitness ();
		if (fitness > max) {
			max = fitness;
			bestIndex = i;
		}
		stFitness.record (fitness);
	}

	if (output)
		showStatistics ();
	generation++;
}


bool DSMGA::shouldTerminate () {
	bool
		termination = false;

	if (maxFe != -1) {
		if (fe > maxFe)
			termination = true;
	}

	if (maxGen != -1) {
		if (generation > maxGen)
			termination = true;
	}

	if (population[0].getMaxFitness() <= stFitness.getMax() )
		termination = true;

	return termination;

}


bool DSMGA::foundOptima () {
	return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA::showStatistics () {

	printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f Chromsome Length:%d\n",
			generation, stFitness.getMax (), stFitness.getMean (),
			stFitness.getMin (), population[0].getLength ());
	printf ("best chromosome:");
	population[bestIndex].printOut ();
	printf ("\n");

	if (SHOW_LINKAGE)
		bbi.printOut ();

	fflush(NULL);
}


void
DSMGA::getBBIOneMax () {
	int i;

	for (i = 0; i < ell; i++) {
		bbi.bb[i][0] = 1;
		bbi.bb[i][1] = i;
	}

	bbi.bbNum = ell;
}


void DSMGA::getBBIMKTrap() {
	int i,j;

	for (i=0; i<ell/TRAP_K; i++) {
		bbi.bb[i][0] = TRAP_K;
		for (j=1; j<=TRAP_K; j++)
			bbi.bb[i][j] = TRAP_K*i + (j-1);
	}

	bbi.bbNum = ell/TRAP_K;
}


void DSMGA::selection () {
	tournamentSelection ();
}


// tournamentSelection without replacement
void DSMGA::tournamentSelection () {
	int i, j;

	int randArray[selectionPressure * nCurrent];

	for (i = 0; i < selectionPressure; i++)
		myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0,
				nCurrent - 1);

	for (i = 0; i < nCurrent; i++) {

		int winner = 0;
		double winnerFitness = -INF;

		for (j = 0; j < selectionPressure; j++) {
			int challenger = randArray[selectionPressure * i + j];
			double challengerFitness = population[challenger].getFitness ();

			if (challengerFitness > winnerFitness) {
				winner = challenger;
				winnerFitness = challengerFitness;
			}

		}
		selectionIndex[i] = winner;
	}
}


void DSMGA::initializePopulation () {
	int i, j;

	for (i = 0; i < nInitial; i++)
		for (j = 0; j < ell; j++)
			if (myRand.flip ())
				population[i].setVal (j, 1);
			else
				population[i].setVal (j, 0);

}


void DSMGA::crossOver () {
	int i;
	getBBI ();

	if ((nCurrent & 0x1) == 0) {
		// nCurrent is even
		for (i = 0; i < nCurrent; i += 2) {
			crossOver (population[selectionIndex[i]],
					population[selectionIndex[i + 1]], offspring[i],
					offspring[i + 1]);
		}
	} else {
		for (i = 0; i < nCurrent - 1; i += 2) {
			crossOver (population[selectionIndex[i]],
					population[selectionIndex[i + 1]], offspring[i],
					offspring[i + 1]);
		}
		offspring[nCurrent - 1] = population[selectionIndex[nCurrent - 1]];
	}

}


void DSMGA::crossOver (Chromosome & p1, Chromosome & p2, Chromosome & c1, Chromosome & c2) {
	if (myRand.uniform () < pc) {
		double prob = 0.5;
		bbUniformXO (p1, p2, c1, c2, prob);
	} else {
		c1 = p1;
		c2 = p2;
	}
}


void DSMGA::mutation () {

	int i;
	getBBI ();

	for (i = 0; i < nCurrent; i++) {
		bbMutation (offspring[i]);
	}
}


void DSMGA::bbMutation (Chromosome & ch) {
	int i, j;

	for (i = 0; i < bbi.bbNum; i++) {
		if (myRand.uniform () < pm) {
			for (j = 1; j <= bbi.bb[i][0]; j++) {
				ch.setVal (bbi.bb[i][j], myRand.flip ());
			}
		}
	}
}


void DSMGA::bbUniformXO (Chromosome & p1, Chromosome & p2, Chromosome & c1, Chromosome & c2, double prob) {
	int i, j;

	for (i = 0; i < bbi.bbNum; i++) {
		if (myRand.flip (prob)) {
			for (j = 1; j <= bbi.bb[i][0]; j++) {
				c1.setVal (bbi.bb[i][j], p1.getVal(bbi.bb[i][j]));
				c2.setVal (bbi.bb[i][j], p2.getVal(bbi.bb[i][j]));
			}
		} else {
			for (j = 1; j <= bbi.bb[i][0]; j++) {
				c1.setVal (bbi.bb[i][j], p2.getVal(bbi.bb[i][j]));
				c2.setVal (bbi.bb[i][j], p1.getVal(bbi.bb[i][j]));
			}
		}
	}

}

