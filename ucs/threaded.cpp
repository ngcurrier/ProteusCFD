#include "threaded.h"

#ifdef _OPENMP
Int LockVars(std::vector<pthread_mutex_t*>& mutexes)
{
  Int err;
  err = 0;
  for(std::vector<pthread_mutex_t*>::iterator it = mutexes.begin(); it != mutexes.end(); ++it){
    pthread_mutex_t* mutex = *it;
    err = pthread_mutex_lock(mutex);
    if(err != 0){
      std::cerr << "MUTEX: locking error!" << std::endl;
      return (err);
    }
  }
  return (0);
}

Int UnlockVars(std::vector<pthread_mutex_t*>& mutexes)
{
  Int err;
  err = 0;
  for(std::vector<pthread_mutex_t*>::iterator it = mutexes.begin(); it != mutexes.end(); ++it){
    pthread_mutex_t* mutex = *it;
    err = pthread_mutex_unlock(mutex);
    if(err != 0){
      std::cerr << "MUTEX: unlocking error!" << std::endl;
      return (err);
    }
  }
  return (0);
}

#else
Int LockVars(std::vector<pthread_mutex_t*>& mutexes)
{
  return 0;
}

Int UnlockVars(std::vector<pthread_mutex_t*>& mutexes)
{
  return 0;
}
#endif

