#include "timer.h"
TimerList::TimerList(Int numTimers) :
  list(NULL)
{
  InitList(numTimers);
}

TimerList::TimerList() :
  list(NULL)
{
  //this is a dummy function, the work is done in InitList()
  return;
}

TimerList::~TimerList()
{
  delete [] list;
}

void TimerList::InitList(Int numTimers)
{
  Int i;
  countTimers = 0;
  maxTimers = numTimers;
  list = new Timer[maxTimers];
  for(i = 0; i < numTimers; i++){
    list[i].stopped = false;
    list[i].started = false;
    list[i].accumulated = 0.0;
    list[i].acc_started = false;
  }
}

Int TimerList::CreateTimer(std::string timerName)
{
  if(countTimers >= maxTimers){
    std::cerr << "Timer not created b/c it would result in more than allocated timers" << std::endl;
    std::cerr << "See constructor call to TimerList()" << std::endl;
    return (-999);
  }
  list[countTimers].name = timerName;
  countTimers++;
  //return timer id if access time is an issue
  return(countTimers - 1);
}

Int TimerList::StartTimer(std::string timerName)
{
  Int tnum = GetNumber(timerName);
  return (StartTimer(tnum));
}

Int TimerList::StartTimer(Int number)
{
  if(number < countTimers && number >= 0){
    gettimeofday(&list[number].startTime, NULL);
    list[number].splitTime = list[number].startTime;
    list[number].started = true;
    return (0);
  }
  else{
    std::cerr << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}

Int TimerList::StopTimer(std::string timerName)
{
  Int tnum = GetNumber(timerName);
  return (StopTimer(tnum));
}

Int TimerList::StopTimer(Int number)
{
  if(number < countTimers && number >= 0){
    gettimeofday(&list[number].endTime, NULL);
    list[number].stopped = true;
    return (0);
  }
  else{
    std::cerr << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}

Int TimerList::PrintBreak(std::string timerName, std::ostream& str)
{
  Int tnum = GetNumber(timerName);
  return (PrintBreak(tnum, str));
}

Int TimerList::PrintBreak(Int number, std::ostream& str)
{
  struct timeval breakTime;
  Real deltaTime = 0.0;
  if(number < countTimers && number >= 0){
    if(list[number].started){
      gettimeofday(&breakTime, NULL);
      deltaTime = (Real)(breakTime.tv_sec - list[number].startTime.tv_sec);
      deltaTime += (Real)(breakTime.tv_usec - list[number].startTime.tv_usec)/1.0e+6;
      str << list[number].name << ": " << deltaTime << " sec" << std::endl;
      return (0);
    }
    else{
      str << list[number].name << ": Timer has not been started" << std::endl;
      return (1);
    }
  }
  else{
    std::cerr << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}
Int TimerList::PrintSplit(std::string timerName, std::ostream& str)
{
  Int tnum = GetNumber(timerName);
  return (PrintSplit(tnum, str));
}

Int TimerList::PrintSplit(Int number, std::ostream& str)
{
  struct timeval breakTime;
  Real deltaTime = 0.0;
  if(number < countTimers && number >= 0){
    if(list[number].started){
      gettimeofday(&breakTime, NULL);
      deltaTime = (Real)(breakTime.tv_sec - list[number].splitTime.tv_sec);
      deltaTime += (Real)(breakTime.tv_usec - list[number].splitTime.tv_usec)/1.0e+6;
      //now reset for the next split
      list[number].splitTime = breakTime;
      str << list[number].name << ": " << deltaTime << " sec" << std::endl;
      return (0);
    }
    else{
      str << list[number].name << ": Timer has not been started" << std::endl;
      return (1);
    }
  }
  else{
    std::cerr << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}

Int TimerList::PrintTimer(std::string timerName, std::ostream& str)
{
  Int tnum = GetNumber(timerName);
  return (PrintTimer(tnum, str));
}

Int TimerList::PrintTimer(Int number, std::ostream& str)
{
  Real deltaTime = 0.0;
  struct timeval temp;
  if(number < countTimers && number >= 0){
    if(list[number].started){
      if(!list[number].stopped){
	gettimeofday(&temp, NULL);
      }
      else{
	temp = list[number].endTime;
      }
      deltaTime = (Real)(temp.tv_sec - list[number].startTime.tv_sec);
      deltaTime += (Real)(temp.tv_usec - list[number].startTime.tv_usec)/1.0e+6;
      str << list[number].name << ": " << deltaTime << " sec" << std::endl;
      return (0);
    }
    else{
      str << list[number].name << ": Timer has not been started" << std::endl;
      return (1);
    }
  }
  else{
    std::cerr << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}

void TimerList::PrintAllTimers(std::ostream& str)
{
  Int i;
  Real deltaTime;
  struct timeval temp;
  for(i = 0; i < countTimers; i++){
    if(list[i].started){
      if(!list[i].stopped){
	gettimeofday(&temp, NULL);
      }
      else{
	temp = list[i].endTime;
      }
      deltaTime = 0.0;
      deltaTime = (Real)(temp.tv_sec - list[i].startTime.tv_sec);
      deltaTime += (Real)(temp.tv_usec - list[i].startTime.tv_usec)/1.0e+6;
      str << list[i].name << ": " << deltaTime << " sec" << std::endl;
    }
    //don't print anything, we presumably only care about timers which were running
  }
  return;
}

void TimerList::PrintAllAccumulate(std::ostream& str)
{
  for(Int i = 0; i < countTimers; i++){
    //only print timers which have accumulators started
    if(list[i].acc_started){
      PrintAccumulate(i, str);
    }
  }
}

Int TimerList::ResetAccumulate(std::string timerName)
{
  Int tnum = GetNumber(timerName);
  return (ResetAccumulate(tnum));
}

Int TimerList::ResetAccumulate(Int number){
  if(number < countTimers && number >= 0){
    list[number].accumulated = 0.0;
    return (0);
  }
  else{
    std::cerr << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}

Int TimerList::StartAccumulate(std::string timerName)
{
  Int tnum = GetNumber(timerName);
  return (StartAccumulate(tnum));
}

Int TimerList::StartAccumulate(Int number){
  if(number < countTimers && number >= 0){
    gettimeofday(&list[number].acc_startTime, NULL);
    list[number].acc_started = true;
    return (0);
  }
  else{
    std::cerr << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}

Int TimerList::PauseAccumulate(std::string timerName)
{
  Int tnum = GetNumber(timerName);
  return (PauseAccumulate(tnum));
}

Int TimerList::PauseAccumulate(Int number){
  Real deltaTime;
  struct timeval temp;
  if(number < countTimers && number >= 0){
    gettimeofday(&temp, NULL);
    deltaTime = 0.0;
    deltaTime = (Real)(temp.tv_sec - list[number].acc_startTime.tv_sec);
    deltaTime += (Real)(temp.tv_usec - list[number].acc_startTime.tv_usec)/1.0e+6; 
    list[number].accumulated += deltaTime;
    return (0);
  }
  else{
    std::cerr << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}

Int TimerList::PauseAccumulateAndPrint(std::string timerName, std::ostream& str)
{
  Int tnum = GetNumber(timerName);
  return (PauseAccumulateAndPrint(tnum, str));
}

Int TimerList::PauseAccumulateAndPrint(Int number, std::ostream& str){
  Real deltaTime;
  struct timeval temp;
  if(number < countTimers && number >= 0){
    gettimeofday(&temp, NULL);
    deltaTime = 0.0;
    deltaTime = (Real)(temp.tv_sec - list[number].acc_startTime.tv_sec);
    deltaTime += (Real)(temp.tv_usec - list[number].acc_startTime.tv_usec)/1.0e+6; 
    list[number].accumulated += deltaTime;
    str << list[number].name << ": " << deltaTime << " sec" << std::endl;
    return (0);
  }
  else{
    str << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}

Int TimerList::PrintAccumulate(std::string timerName, std::ostream& str)
{
  Int tnum = GetNumber(timerName);
  return (PrintAccumulate(tnum, str));
}

Int TimerList::PrintAccumulate(Int number, std::ostream& str)
{
  if(number < countTimers && number >= 0){
    str << list[number].name << ": " << list[number].accumulated << " sec" << std::endl;
    return (0);
  }
  else{
    std::cerr << "Timer number " << number << " not found!!" << std::endl;
    return (1);
  }
}

Int TimerList::GetNumber(std::string timerName)
{
  for(Int i = 0; i < countTimers; i++){
    if(timerName == list[i].name){
      return i;
    }
  }
  std::cerr << "Timer named " << timerName << " not found!!" << std::endl;
  return -1;
}
