#ifndef TIMER_H__
#define TIMER_H__

#include "general.h"
#include <sys/time.h>
#include <string>
#include <iostream>
#include <fstream>

class Timer{
  friend class TimerList;

public:
  Timer(){};
  ~Timer(){};

  std::string name;
  bool started;
  bool stopped;
  bool acc_started;
  struct timeval startTime;
  struct timeval endTime;
  struct timeval splitTime;
  Real accumulated;
  struct timeval acc_startTime;
private:
};

class TimerList
{
 public:
  //the preference is to use this constructor
  TimerList(Int numTimers);
  //never use this directly, only use if you need the timerList object
  //inside another class, i.e. you need the pointer
  TimerList();
  //likewise, only use this if you need the timerList object inside a class
  //which will allocate the list later
  void InitList(Int numTimers);
  ~TimerList();

  //creates a timer and returns its id
  //on failure returns negative
  Int CreateTimer(std::string timerName);
  //starts a timer
  Int StartTimer(Int number);
  Int StartTimer(std::string timerName);
  //stops a timer
  Int StopTimer(Int number);
  Int StopTimer(std::string timerName);
  //will print time from start but will not stop timer
  Int PrintBreak(Int number, std::ostream& str = std::cout);
  Int PrintBreak(std::string timerName, std::ostream& str = std::cout);
  //will print time from the last time a split was printed but will not stop timer
  Int PrintSplit(Int number, std::ostream& str = std::cout);
  Int PrintSplit(std::string timerName, std::ostream& str = std::cout);
  Real GetSplit(std::string timerName);
  Real GetSplit(Int number);
  //will only print time from start to stop
  //WARNING: if timer has not been stopped time will be garbage
  Int PrintTimer(Int number, std::ostream& str = std::cout);
  Int PrintTimer(std::string timerName, std::ostream& str = std::cout);
  //Will reset accumulation timer
  Int ResetAccumulate(Int number);
  Int ResetAccumulate(std::string timerName);
  //Will start an accumulation timer without a reset, i.e. start accumulating
  Int StartAccumulate(Int number);
  Int StartAccumulate(std::string timerName);
  //Will stop an accumulation timer without a reset, i.e. stop accumulating
  Int PauseAccumulate(Int number);
  Int PauseAccumulate(std::string timerName);
  //Will stop an accumulation timer without a reset and print the split, i.e. stop accumulating
  Int PauseAccumulateAndPrint(Int number, std::ostream& str = std::cout);
  Int PauseAccumulateAndPrint(std::string timerName, std::ostream& str = std::cout);
  //will print the time accumulated with the correct functions
  Int PrintAccumulate(Int number, std::ostream& str = std::cout);
  Int PrintAccumulate(std::string timerName, std::ostream& str = std::cout);
  //will print all std timers from start to stop
  void PrintAllTimers(std::ostream& str = std::cout);
  //will print all accumulation total from timers which have accumulators started
  void PrintAllAccumulate(std::ostream& str = std::cout);

 private:
  Int GetNumber(std::string timerName);

  Int countTimers;
  Int maxTimers;
  Timer* list;
};


#endif
