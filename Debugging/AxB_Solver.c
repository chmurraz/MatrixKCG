#include "kcg_types.h"
#include "AxB_Solver.h"


void AxB_Solver( kcg_real A[4][4], kcg_real B[4], int rank, kcg_real X[4])
{

	kcg_int size = 4;

	kcg_real C[4][4 + 1];

	kcg_int i,j,k;

    for(  i = 0; i < rank; i++) {

		for(  j = 0; j < rank; j++) {

			C[i][j] = A[i][j];

		}

		C[i][j] = B[i];

	}


    for(  k = 0; k < rank-1; k++) {

        for(  i = k+1; i < rank; i++) {

			C[i][k] = divprotect(C[i][k], C[k][k]);

            for(  j = k+1; j < rank+1; j++) {

                C[i][j] = C[i][j] - C[i][k]*C[k][j];

			}
		}
	}



    for(  i = rank-1; i > -1; i--) {

        for(  j = i+1; j < rank; j++) {

            C[i][rank] = C[i][rank] - C[i][j] * C[j][rank];
		}

		C[i][rank] = divprotect(C[i][rank],C[i][i]);

	}

for ( i = 0; i < rank; i++) {
		X[i] = C[i][rank];
}

}

kcg_real divprotect( kcg_real n, kcg_real d)
{

	kcg_real q;

	if( d < -1.0e-9 || d > 1.0e-9) {

		q = n / d;

	}else if(d < 0.0) {

		q = n / -1.0e-9;

	}else{

		q = n /  1.0e-9;
	}

	return q;

}
