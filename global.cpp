/***************************************************************************
 *   Copyright (C) 2005 by Tian-Li Yu,,,                                   *
 *   tianliyu@fishlaptop.ytgroup                                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/



#include <cstdio>
#include <cstdlib>
#include "global.h"
#include "myrand.h"
#include "statistics.h"

double  BTRAP_ALPHA = 1.0;
int  TRAP_K = 5;

int  maxMemory = 0;
bool SHOW_HC = false;
bool SHOW_DSM = false;
bool SHOW_LINKAGE = true;
bool SHOW_POPULATION = false;
bool SHORT_HAND = false;
bool SHOW_SELECTION_INDEX = false;
bool SHOW_REPLACEMENT = false;
bool SHOW_MAPPING = false;
bool SHOW_COMPRESSION = false;
bool SHOW_BUSINESSMEN = false;
bool SHOW_BISECTION = true;
bool SHOW_THRESHOLD = false;


char outputFilename[100];
MyRand myRand;
BitwiseDistance myBD;
TableLookUp myTableLookUp;

void outputErrMsg (const char *errMsg)
{
  printf ("%s\n", errMsg);
  exit (1);
}

int
pow2 (int x)
{
  return (1 << x);
}

void
findMax (int number, int index[], double values[], int size)
{
  int i, j, k;
  double *max = new double[number];


  for (i = 0; i < number; i++)
    {
      index[i] = -1;
      max[i] = -1;
    }

  for (i = 0; i < size; i++)
    for (j = 0; j < number; j++)

      if (values[i] > max[j])
	{
	  for (k = number - 1; k > j; k--)
	    {
	      index[k] = index[k - 1];
	      max[k] = max[k - 1];
	    }
	  index[j] = i;
	  max[j] = values[i];

	  break;
	}

  delete[]max;
}

/** find number of maxima and return them in index, only for positive number*/
void
findMax (int number, int index[], int values[], int size)
{
  int i, j, k;
  int *max = new int[number];


  for (i = 0; i < number; i++)
    {
      index[i] = -1;
      max[i] = -1;
    }

  for (i = 0; i < size; i++)
    for (j = 0; j < number; j++)

      if (values[i] > max[j])
	{
	  for (k = number - 1; k > j; k--)
	    {
	      index[k] = index[k - 1];
	      max[k] = max[k - 1];
	    }
	  index[j] = i;
	  max[j] = values[i];

	  break;
	}

  delete[]max;
}

/** findMin--int */
void
findMin (int number, int index[], int values[], int size)
{
  int i, j, k;
  int *min = new int[number];


  for (i = 0; i < number; i++)
    {
      index[i] = -1;
      min[i] = (int)INF * 2;
    }

  for (i = 0; i < size; i++)
    for (j = 0; j < number; j++)

      if (values[i] < min[j])
	{
	  for (k = number - 1; k > j; k--)
	    {
	      index[k] = index[k - 1];
	      min[k] = min[k - 1];
	    }
	  index[j] = i;
	  min[j] = values[i];

	  break;
	}

  delete[]min;
}


/** findMin--double */
void
findMin (int number, int index[], double values[], int size)
{
  int i, j, k;
  double *min = new double[number];


  for (i = 0; i < number; i++)
    {
      index[i] = -1;
      min[i] = INF * 2;
    }

  for (i = 0; i < size; i++)
    for (j = 0; j < number; j++)

      if (values[i] < min[j])
	{
	  for (k = number - 1; k > j; k--)
	    {
	      index[k] = index[k - 1];
	      min[k] = min[k - 1];
	    }
	  index[j] = i;
	  min[j] = values[i];

	  break;
	}

  delete[]min;
}

int myDoubleCompare(const void *a, const void *b)
{
	if (*((double *)a) > *((double *)b))
		return 1;
	else if (*((double *)a) < *((double *)b))
		return -1;
	else 
		return 0;
}

// return the x which x*exp(x) = w
double lambertw(double w)
{
    double upper = 100.0;
    double lower = 0.0;

    assert(upper * exp(upper) > w);
    assert(w > 0);

    double middle;

    while ((upper - lower) / upper > 0.01) {
	middle = (upper + lower) / 2;
	
	double z = middle * exp(middle);
	if (z < w)
	    lower = middle;
	else
	    upper = middle;
    }

    return (upper + lower) / 2;
}
