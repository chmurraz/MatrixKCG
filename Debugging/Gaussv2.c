#include "kcg_types.h"

void Gaussv2(kcg_float32 **A, kcg_float32 *x, kcg_float32 *b, int r, int c)
{
	int row;
	int col;
	int pivotRow;
	int pivotCol;
	int forwardElimRow;
	int forwardElimCol;
	int backwardElimRow;
	int backwardElimCol;
	kcg_float32 forwardElimRowLeadingCoeff;
	kcg_float32 backwardElimRowLeadingCoeff;
	//kcg_float32 forwardElimRowLeadingCoeff;
	kcg_float32 pivot;

	//forwardElimRowLeadingCoeff = (kcg_float32 *)malloc(sizeof(kcg_float32));

	for (row = 0; row<3; row++)
	{
		// Search for the pivot in the current row
		for (col = 0; col<3; col++)
		{
			if (A[row][col] == kcg_lit_float32(0))
			{
				// Add code here to handle singular matrix case where a row has no valid pivot
			}
			else
			{
				pivot = A[row][col];
				pivotRow = row;
				pivotCol = col;
				break;
			}
		}

		// Normalize the current row using the pivot
		for (col = pivotCol; col<3; col++)
		{
			A[row][col] = A[row][col] / pivot;
		}

		// Normalize the RHS
		b[row] = b[row] / pivot;

		// Begin forward elimination
		for (forwardElimRow = pivotRow + 1; forwardElimRow<3; forwardElimRow++)
		{
			forwardElimRowLeadingCoeff = A[forwardElimRow][pivotCol];
			//forwardElimRowLeadingCoeff = *A[forwardElimRow][pivotCol];
			for (forwardElimCol = pivotCol; forwardElimCol<3; forwardElimCol++)
			{
				A[forwardElimRow][forwardElimCol] = A[forwardElimRow][forwardElimCol] - forwardElimRowLeadingCoeff * (A[pivotRow][forwardElimCol]);
				//*A[forwardElimRow][forwardElimCol] = *A[forwardElimRow][forwardElimCol] - forwardElimRowLeadingCoeff * (*A[pivotRow][forwardElimCol]);
			}

		}

	}
}