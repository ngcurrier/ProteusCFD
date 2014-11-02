#ifndef _THREADED_H__
#include "general.h"
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

Int LockVars(std::vector<pthread_mutex_t*>& mutexes);
Int UnlockVars(std::vector<pthread_mutex_t*>& mutexes);

#endif
