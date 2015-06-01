/*

    lse.c - logsumexp implementation in C

    author: Christoph Wehmeyer <christoph.wehmeyer@fu-berlin.de>

*/

#include "_lse.h"

// old m$ visual studio is not c99 compliant (vs2010 eg. is not)
#ifdef _MSC_VER
	#include <math.h>
	#include <float.h>
	#define INFINITY (DBL_MAX+DBL_MAX)
	#define NAN (INFINITY-INFINITY)
#endif

/* _sort()is based on examples from http://www.linux-related.de (2004) */
void _sort( double *array, int L, int R )
{
    int l, r;
    double swap;
    if( R-L > 25 ) /* use quicksort */
    {
        l = L-1;
        r = R;
        for(;;)
        {
            while( array[++l] < array[R] );
            while( array[--r] > array[R] && r > l );
            if( l >= r  ) break;
            swap = array[l];
            array[l] = array[r];
            array[r] = swap;
        }
        swap = array[l];
        array[l] = array[R];
        array[R] = swap;
        _sort( array, L, l-1 );
        _sort( array, l+1, R );
    }
    else /* use insertion sort */
    {
        for( l=L+1; l<=R; ++l )
        {
            swap = array[l];
            for( r=l-1; r>=L && swap<array[r]; --r )
                array[r+1] = array[r];
            array[r+1] = swap;
        }
    }
}

double _logsumexp( double *array, int length )
{
    int i;
    double sum=0.0;
    _sort( array, 0, length-1 );
    if( -INFINITY == array[length-1] )
        return -INFINITY;
    for( i=0; i<length-1; ++i )
        sum += exp( array[i] - array[length-1] );
    return array[length-1] + log( sum + 1.0 );
}

double _logsumexp_pair( double a, double b )
{
    if( ( -INFINITY == a ) && ( -INFINITY == b ) )
        return -INFINITY;
    if( b > a )
        return b + log( 1.0 + exp( a - b ) );
    return a + log( 1.0 + exp( b - a ) );
}
