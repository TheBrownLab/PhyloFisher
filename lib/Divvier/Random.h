/*
 * Random.h
 *
 *  Created on: Oct 20, 2017
 *      Author: simon
 */

#ifndef RANDOM_H_
#define RANDOM_H_

/* -----------------------------------------------------------------------
 * Name            : rngs.h  (header file for the library file rngs.c)
 * Author          : Steve Park & Dave Geyer
 * Language        : ANSI C
 * Latest Revision : 09-22-98
 * -----------------------------------------------------------------------
 */

#if !defined( _RNGS_ )
#define _RNGS_

#include <cstdlib>
#include <ctype.h>
#include <iostream>
#include <numeric>
#include <assert.h>
#include <cfloat>
#include <cmath>

double Random(void);
void   PlantSeeds(long x);
void   GetSeed(long *x);
void   PutSeed(long x);
void   SelectStream(int index);
void   TestRandom(void);

#endif

int RandInt(int begin,int end);
inline double RandDouble(double Low, double High)       { return Low + (fabs(High - Low) * Random()); }


#endif /* RANDOM_H_ */
