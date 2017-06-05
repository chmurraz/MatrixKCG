#include "kcg_types.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{
	array_float32_3_3 A;
	array_float32_3_3 ARREF;
	array_float32_3 b;
	array_float32_3 bRREF;
	array_float32_3 x;
	kcg_int32 n = 3;
	kcg_int32 rank = 0;
	kcg_bool consistent = kcg_true;

	/*
	//	Case 1
	A[0][0] = 9.0;
	A[0][1] = 3.0;
	A[0][2] = 4.0;
	A[1][0] = 4.0;
	A[1][1] = 3.0;
	A[1][2] = 4.0;
	A[2][0] = 1.0;
	A[2][1] = 1.0;
	A[2][2] = 1.0;
	b[0] = 7.0;
	b[1] = 8.0;
	b[2] = 3.0;
	//	Solution will be -1/5, 4, -4/5
	*/

	//	Case 2
	A[0][0] = 2.0;
	A[0][1] = 4.0;
	A[0][2] = 8.0;
	A[1][0] = 3.0;
	A[1][1] = 6.0;
	A[1][2] = 5.0;
	A[2][0] = 4.0;
	A[2][1] = 8.0;
	A[2][2] = 2.0;
	b[0] = 16.0;
	b[1] = 3.0;
	b[2] = -10.0;
	//	Solution will be 4, 3, true

	/*
	//	Case 3
	A[0][0] = 5.0;
	A[0][1] = -2.0;
	A[0][2] = 1.0;
	A[1][0] = 4.0;
	A[1][1] = 1.0;
	A[1][2] = -7.0;
	A[2][0] = 1.0;
	A[2][1] = -3.0;
	A[2][2] = 8.0;
	b[0] = 2.0;
	b[1] = 10.0;
	b[2] = -4.0;
	*/

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			ARREF[i][j] = A[i][j];
		bRREF[i] = b[i];
	}

	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;

	Gauss_3(&A, &ARREF, &b, &bRREF, n, &x, &rank, &consistent);

	

	int i, j, r = 3, c = 3;
	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			printf("%f", ARREF[i][j]);
			printf(" ");
			if (j == c - 1)
			{
				printf("%f", x[i]);
				printf("\n");
			}
		}
	}

	return 0;
}