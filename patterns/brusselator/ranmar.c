/************************************************************************
 This random number generator originally appeared in "Toward a Universal
 Random Number Generator" by George Marsaglia and Arif Zaman.
 Florida State University Report: FSU-SCRI-87-50 (1987)

 It was later modified by F. James and published in "A Review of Pseudo-
 random Number Generators"

 Converted from FORTRAN to C by Phil Linttell, James F. Hickling
 Management Consultants Ltd, Aug. 14, 1989.

 THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
       (However, a newly discovered technique can yield
         a period of 10^600. But that is still in the development stage.)

 It passes ALL of the tests for random number generators and has a period
   of 2^144, is completely portable (gives bit identical results on all
   machines with at least 24-bit mantissas in the floating point
   representation).

 The algorithm is a combination of a Fibonacci sequence (with lags of 97
   and 33, and operation "subtraction plus one, modulo one") and an
   "arithmetic sequence" (using subtraction).

 On a Vax 11/780, this random number generator can produce a number in
    13 microseconds.
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TRUE    1
#define FALSE   0

float u[97], c, cd, cm;
int i97, j97, test;

/*int rmarin(int ij, int kl);
int ranmar(float rvec[], int len);
*/

int rmarin(int ij, int kl)
{

        float s, t;
        int i, j, k, l, m;
        int ii, jj;

        /* Change FALSE to TRUE in the next statement to test the
           random routine.*/

        test = TRUE;

        if ( ( ij < 0 || ij > 31328 ) ||
                ( kl < 0 || kl > 30081 ) )
        {
                printf ("RMARIN: The first random number seed must have a "
                        "value between 0 and 31328\n");
                printf ("        The second random number seed must have a "
                        "value between 0 and 30081\n");
                return 1;
        }

        i = (int)fmod(ij/177.0, 177.0) + 2;
        j = (int)fmod((double)ij, 177.0) + 2;
        k = (int)fmod(kl/169.0, 178.0) + 1;
        l = (int)fmod((double)kl, 169.0);

        for ( ii=0; ii<=96; ii++ )
        {
                s = (float)0.0;
                t = (float)0.5;
                for ( jj=0; jj<=23; jj++ )
                {
                        m = (int)fmod( fmod((double)i*j,179.0)*k , 179.0 );
                        i = j;
                        j = k;
                        k = m;
                        l = (int)fmod( 53.0*l+1.0 , 169.0 );
                        if ( fmod((double)l*m,64.0) >= 32)
                                s = s + t;
                        t = (float)(0.5 * t);
                }
                u[ii] = s;
        }

        c  = (float)(  362436.0 / 16777216.0);
        cd = (float)( 7654321.0 / 16777216.0);
        cm = (float)(16777213.0 / 16777216.0);

        i97 = 96;
        j97 = 32;

        test = TRUE;

        return 0;
}


int ranmar(float rvec[], int len)
{
        float uni;
        int ivec;

        if ( !test )
        {
                printf ("RANMAR: Call the initialization routine (RMARIN) "
                        "before calling RANMAR.\n");
                return 1;
        }

        for ( ivec=0; ivec < len; ivec++)
        {
                uni = u[i97] - u[j97];
                if ( uni < 0.0F )
                        uni = uni + 1.0;
                u[i97] = uni;
                i97--;
                if ( i97 < 0 )
                        i97 = 96;
                j97--;
                if ( j97 < 0 )
                        j97 = 96;
                c = c - cd;
                if ( c < 0.0F )
                        c = c + cm;
                uni = uni - c;
                if ( uni < 0.0F )
                        uni = uni + 1.0;
                rvec[ivec] = uni;
        }
        return 0;
}
