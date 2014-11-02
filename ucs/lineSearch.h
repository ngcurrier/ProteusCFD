#ifndef LINE_SEARCH_H__
#define LINE_SEARCH_H__

#include "general.h"
#include "portFunctions.h"

void lineSearch(int ndv, double* x, double* bounds, void(*function)(FUNC_ARGS), 
		void(*gradient)(GRAD_ARGS), int* uiparm, double* urparm);

#endif
