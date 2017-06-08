#include "kcg_types.h"
#include "string.h"
kcg_float32 divprotect(kcg_float32 n, kcg_float32 d);
kcg_bool floatequalzero(kcg_float32 x);

void Gauss_3x3(array_float32_3_3 *A, array_float32_3 *b, kcg_int32 n, array_float32_3 *x, array_float32_3_3 *ARREF, kcg_int32 *rank, kcg_bool *consistent)
{
	kcg_int32 row;
	kcg_int32 col;
	kcg_int32 pivotRow;
	kcg_int32 pivotCol;
	kcg_float32 pivot;
	kcg_int32 elimRow;
	kcg_int32 elimCol;
	kcg_float32 elimRowLeadingCoeff;
	kcg_float32 *Pivot_ij;
	kcg_float32 *Elim_ij;
	kcg_int32 zeroEntryCount;

	/////////////////////////////////////////////////////////
	//
	//	Initialize output variables:  x, rank and consistent
	//
	/////////////////////////////////////////////////////////
	
	//	Note:  The final values for x are to be calculated using row reduction
	memcpy(*x, *b, sizeof(*x));

	//	Note:  ARREF is a copy of the original input matrix.  ARREF will be put in reduced row echelon form.
	//	so it's values will be different from the original matrix A
	memcpy(*ARREF, *A, sizeof(*A));
	
	//	Initialize matrix rank at zero
	*rank = 0;

	//  Initialize the system as consistent
	*consistent = kcg_true;

	/////////////////////////////////////////////
	//
	//	Begin forward elimination
	//
	/////////////////////////////////////////////
	for (row = 0; row<n; row++)
	{
		// Initialize the count of zero entries in this row and search for the pivot in the current row
		zeroEntryCount = 0;
		for (col = 0; col<n; col++)
		{
			Pivot_ij = (**ARREF + row*n + col);
			if (floatequalzero(*Pivot_ij))
			{
				// Add code here to handle non-full rank case where row has no pivot (no non-zero entry)
				zeroEntryCount = zeroEntryCount + 1;
			}
			else
			{
				pivot = *Pivot_ij;
				pivotRow = row;
				pivotCol = col;
				*rank = *rank + 1;
				break;
			}
		}

		//	If a valid pivot exists in this row, proceed with forward elimination
		if (zeroEntryCount < n)
		{
			// Normalize the current row by dividing by the pivot value
			for (col = pivotCol; col < n; col++)
			{
				Pivot_ij = **ARREF + row*n + col;
				*Pivot_ij = divprotect(*Pivot_ij, pivot);
			}

			// Normalize the solution by dividing by the pivot value
			*(*x + row) = divprotect(*(*x + row), pivot);

			// Begin row operations on all rows below the pivot row
			for (elimRow = pivotRow + 1; elimRow < n; elimRow++)
			{
				//	Row operations on each coefficient row below the pivot row
				elimRowLeadingCoeff = *(**ARREF + elimRow*n + pivotCol);
				for (elimCol = pivotCol; elimCol < n; elimCol++)
				{
					Elim_ij = **ARREF + elimRow*n + elimCol;
					Pivot_ij = **ARREF + pivotRow*n + elimCol;
					*Elim_ij = *Elim_ij - elimRowLeadingCoeff * (*Pivot_ij);
				}

				//	Row operation on solution vector
				*(*x + elimRow) = *(*x + elimRow) - elimRowLeadingCoeff * (*(*x + pivotRow));

			}
		}

		// Check for inconsistent system.  If the current row is all zeros and the b-value for this row is non-zero, then inconsistent
		if ((zeroEntryCount == n) && !floatequalzero(*(*x + row)))
			*consistent = kcg_false;
	}

	/////////////////////////////////////////////
	//
	//	Begin reverse elimination
	//
	/////////////////////////////////////////////
	for (row = n-1; row>-1; row--)
	{
		// Initialize the count of zero entries in this row and search for the pivot in the current row
		zeroEntryCount = 0;
		for (col = n-1; col>0; col--)
		{
			Pivot_ij = (**ARREF + row*n + col);
			if (floatequalzero(*Pivot_ij))
			{
				zeroEntryCount = zeroEntryCount + 1;
			}
			else
			{
				pivot = *Pivot_ij;
				pivotRow = row;
				pivotCol = col;
				break;
			}
		}

		//	If a valid pivot exists in this row, proceed with reverse elimination
		if (zeroEntryCount < n)
		{
			// Normalize the current row by dividing by the pivot value
			for (col = pivotCol; col < n; col++)
			{
				Pivot_ij = **ARREF + row*n + col;
				*Pivot_ij = divprotect(*Pivot_ij, pivot);
			}

			// Normalize the solution by dividing by the pivot value
			*(*x + row) = divprotect(*(*x + row), pivot);


			// Begin row operations on all rows above the pivot row
			for (elimRow = pivotRow - 1; elimRow > -1; elimRow--)
			{
				//	Row operations on each coefficient row above the pivot row
				elimRowLeadingCoeff = *(**ARREF + elimRow*n + pivotCol);
				for (elimCol = pivotCol; elimCol < n; elimCol++)
				{
					Elim_ij = **ARREF + elimRow*n + elimCol;
					Pivot_ij = **ARREF + pivotRow*n + elimCol;
					*Elim_ij = *Elim_ij - elimRowLeadingCoeff * (*Pivot_ij);
				}

				//	Row operation on solution vector
				*(*x + elimRow) = *(*x + elimRow) - elimRowLeadingCoeff * (*(*x + pivotRow));

			}
		}
	}
}

void Gauss_4x4(array_float32_4_4 *A, array_float32_4 *b, kcg_int32 n, array_float32_4 *x, array_float32_4_4 *ARREF, kcg_int32 *rank, kcg_bool *consistent)
{
	kcg_int32 row;
	kcg_int32 col;
	kcg_int32 pivotRow;
	kcg_int32 pivotCol;
	kcg_float32 pivot;
	kcg_int32 elimRow;
	kcg_int32 elimCol;
	kcg_float32 elimRowLeadingCoeff;
	kcg_float32 *Pivot_ij;
	kcg_float32 *Elim_ij;
	kcg_int32 zeroEntryCount;

	/////////////////////////////////////////////////////////
	//
	//	Initialize output variables:  x, rank and consistent
	//
	/////////////////////////////////////////////////////////

	//	Note:  The final values for x are to be calculated using row reduction
	memcpy(*x, *b, sizeof(*x));

	//	Note:  ARREF is a copy of the original input matrix.  ARREF will be put in reduced row echelon form.
	//	so it's values will be different from the original matrix A
	memcpy(*ARREF, *A, sizeof(*A));

	//	Initialize matrix rank at zero
	*rank = 0;

	//  Initialize the system as consistent
	*consistent = kcg_true;

	/////////////////////////////////////////////
	//
	//	Begin forward elimination
	//
	/////////////////////////////////////////////
	for (row = 0; row<n; row++)
	{
		// Initialize the count of zero entries in this row and search for the pivot in the current row
		zeroEntryCount = 0;
		for (col = 0; col<n; col++)
		{
			Pivot_ij = (**ARREF + row*n + col);
			if (floatequalzero(*Pivot_ij))
			{
				// Add code here to handle non-full rank case where row has no pivot (no non-zero entry)
				zeroEntryCount = zeroEntryCount + 1;
			}
			else
			{
				pivot = *Pivot_ij;
				pivotRow = row;
				pivotCol = col;
				*rank = *rank + 1;
				break;
			}
		}

		//	If a valid pivot exists in this row, proceed with forward elimination
		if (zeroEntryCount < n)
		{
			// Normalize the current row by dividing by the pivot value
			for (col = pivotCol; col < n; col++)
			{
				Pivot_ij = **ARREF + row*n + col;
				*Pivot_ij = divprotect(*Pivot_ij, pivot);
			}

			// Normalize the solution by dividing by the pivot value
			*(*x + row) = divprotect(*(*x + row), pivot);

			// Begin row operations on all rows below the pivot row
			for (elimRow = pivotRow + 1; elimRow < n; elimRow++)
			{
				//	Row operations on each coefficient row below the pivot row
				elimRowLeadingCoeff = *(**ARREF + elimRow*n + pivotCol);
				for (elimCol = pivotCol; elimCol < n; elimCol++)
				{
					Elim_ij = **ARREF + elimRow*n + elimCol;
					Pivot_ij = **ARREF + pivotRow*n + elimCol;
					*Elim_ij = *Elim_ij - elimRowLeadingCoeff * (*Pivot_ij);
				}

				//	Row operation on solution vector
				*(*x + elimRow) = *(*x + elimRow) - elimRowLeadingCoeff * (*(*x + pivotRow));

			}
		}

		// Check for inconsistent system.  If the current row is all zeros and the b-value for this row is non-zero, then inconsistent
		if ((zeroEntryCount == n) && !floatequalzero(*(*x + row)))
			*consistent = kcg_false;
	}

	/////////////////////////////////////////////
	//
	//	Begin reverse elimination
	//
	/////////////////////////////////////////////
	for (row = n - 1; row>-1; row--)
	{
		// Initialize the count of zero entries in this row and search for the pivot in the current row
		zeroEntryCount = 0;
		for (col = n - 1; col>0; col--)
		{
			Pivot_ij = (**ARREF + row*n + col);
			if (floatequalzero(*Pivot_ij))
			{
				zeroEntryCount = zeroEntryCount + 1;
			}
			else
			{
				pivot = *Pivot_ij;
				pivotRow = row;
				pivotCol = col;
				break;
			}
		}

		//	If a valid pivot exists in this row, proceed with reverse elimination
		if (zeroEntryCount < n)
		{
			// Normalize the current row by dividing by the pivot value
			for (col = pivotCol; col < n; col++)
			{
				Pivot_ij = **ARREF + row*n + col;
				*Pivot_ij = divprotect(*Pivot_ij, pivot);
			}

			// Normalize the solution by dividing by the pivot value
			*(*x + row) = divprotect(*(*x + row), pivot);


			// Begin row operations on all rows above the pivot row
			for (elimRow = pivotRow - 1; elimRow > -1; elimRow--)
			{
				//	Row operations on each coefficient row above the pivot row
				elimRowLeadingCoeff = *(**ARREF + elimRow*n + pivotCol);
				for (elimCol = pivotCol; elimCol < n; elimCol++)
				{
					Elim_ij = **ARREF + elimRow*n + elimCol;
					Pivot_ij = **ARREF + pivotRow*n + elimCol;
					*Elim_ij = *Elim_ij - elimRowLeadingCoeff * (*Pivot_ij);
				}

				//	Row operation on solution vector
				*(*x + elimRow) = *(*x + elimRow) - elimRowLeadingCoeff * (*(*x + pivotRow));

			}
		}
	}
}

void Gauss_6x6(array_float32_6_6 *A, array_float32_6 *b, kcg_int32 n, array_float32_6 *x, array_float32_6_6 *ARREF, kcg_int32 *rank, kcg_bool *consistent)
{
	kcg_int32 row;
	kcg_int32 col;
	kcg_int32 pivotRow;
	kcg_int32 pivotCol;
	kcg_float32 pivot;
	kcg_int32 elimRow;
	kcg_int32 elimCol;
	kcg_float32 elimRowLeadingCoeff;
	kcg_float32 *Pivot_ij;
	kcg_float32 *Elim_ij;
	kcg_int32 zeroEntryCount;

	/////////////////////////////////////////////////////////
	//
	//	Initialize output variables:  x, rank and consistent
	//
	/////////////////////////////////////////////////////////

	//	Note:  The final values for x are to be calculated using row reduction
	memcpy(*x, *b, sizeof(*x));

	//	Note:  ARREF is a copy of the original input matrix.  ARREF will be put in reduced row echelon form.
	//	so it's values will be different from the original matrix A
	memcpy(*ARREF, *A, sizeof(*A));

	//	Initialize matrix rank at zero
	*rank = 0;

	//  Initialize the system as consistent
	*consistent = kcg_true;

	/////////////////////////////////////////////
	//
	//	Begin forward elimination
	//
	/////////////////////////////////////////////
	for (row = 0; row<n; row++)
	{
		// Initialize the count of zero entries in this row and search for the pivot in the current row
		zeroEntryCount = 0;
		for (col = 0; col<n; col++)
		{
			Pivot_ij = (**ARREF + row*n + col);
			if (floatequalzero(*Pivot_ij))
			{
				// Add code here to handle non-full rank case where row has no pivot (no non-zero entry)
				zeroEntryCount = zeroEntryCount + 1;
			}
			else
			{
				pivot = *Pivot_ij;
				pivotRow = row;
				pivotCol = col;
				*rank = *rank + 1;
				break;
			}
		}

		//	If a valid pivot exists in this row, proceed with forward elimination
		if (zeroEntryCount < n)
		{
			// Normalize the current row by dividing by the pivot value
			for (col = pivotCol; col < n; col++)
			{
				Pivot_ij = **ARREF + row*n + col;
				*Pivot_ij = divprotect(*Pivot_ij, pivot);
			}

			// Normalize the solution by dividing by the pivot value
			*(*x + row) = divprotect(*(*x + row), pivot);

			// Begin row operations on all rows below the pivot row
			for (elimRow = pivotRow + 1; elimRow < n; elimRow++)
			{
				//	Row operations on each coefficient row below the pivot row
				elimRowLeadingCoeff = *(**ARREF + elimRow*n + pivotCol);
				for (elimCol = pivotCol; elimCol < n; elimCol++)
				{
					Elim_ij = **ARREF + elimRow*n + elimCol;
					Pivot_ij = **ARREF + pivotRow*n + elimCol;
					*Elim_ij = *Elim_ij - elimRowLeadingCoeff * (*Pivot_ij);
				}

				//	Row operation on solution vector
				*(*x + elimRow) = *(*x + elimRow) - elimRowLeadingCoeff * (*(*x + pivotRow));

			}
		}

		// Check for inconsistent system.  If the current row is all zeros and the b-value for this row is non-zero, then inconsistent
		if ((zeroEntryCount == n) && !floatequalzero(*(*x + row)))
			*consistent = kcg_false;
	}

	/////////////////////////////////////////////
	//
	//	Begin reverse elimination
	//
	/////////////////////////////////////////////
	for (row = n - 1; row>-1; row--)
	{
		// Initialize the count of zero entries in this row and search for the pivot in the current row
		zeroEntryCount = 0;
		for (col = n - 1; col>0; col--)
		{
			Pivot_ij = (**ARREF + row*n + col);
			if (floatequalzero(*Pivot_ij))
			{
				zeroEntryCount = zeroEntryCount + 1;
			}
			else
			{
				pivot = *Pivot_ij;
				pivotRow = row;
				pivotCol = col;
				break;
			}
		}

		//	If a valid pivot exists in this row, proceed with reverse elimination
		if (zeroEntryCount < n)
		{
			// Normalize the current row by dividing by the pivot value
			for (col = pivotCol; col < n; col++)
			{
				Pivot_ij = **ARREF + row*n + col;
				*Pivot_ij = divprotect(*Pivot_ij, pivot);
			}

			// Normalize the solution by dividing by the pivot value
			*(*x + row) = divprotect(*(*x + row), pivot);


			// Begin row operations on all rows above the pivot row
			for (elimRow = pivotRow - 1; elimRow > -1; elimRow--)
			{
				//	Row operations on each coefficient row above the pivot row
				elimRowLeadingCoeff = *(**ARREF + elimRow*n + pivotCol);
				for (elimCol = pivotCol; elimCol < n; elimCol++)
				{
					Elim_ij = **ARREF + elimRow*n + elimCol;
					Pivot_ij = **ARREF + pivotRow*n + elimCol;
					*Elim_ij = *Elim_ij - elimRowLeadingCoeff * (*Pivot_ij);
				}

				//	Row operation on solution vector
				*(*x + elimRow) = *(*x + elimRow) - elimRowLeadingCoeff * (*(*x + pivotRow));

			}
		}
	}
}

kcg_float32 divprotect(kcg_float32 n, kcg_float32 d)
{
	kcg_float32 q;
	if (d < -1.0e-9 || d > 1.0e-9) {
		q = n / d;
	}
	else if (d < 0.0) {
		q = n / -1.0e-9;
	}
	else {
		q = n / 1.0e-9;
	}
	return q;
}

kcg_bool floatequalzero(kcg_float32 x)
{
	if (x < 1.0e-9 && x > -1.0e-9) {
		return kcg_true;
	}
	else {
		return kcg_false;
	}
}