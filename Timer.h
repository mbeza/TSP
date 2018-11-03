#ifndef TIMER_H
#define TIMER_H


#include <chrono>
using namespace std::chrono;

class Timer
{
public:
	void timerStart();
	void timerStop();
	__int64 calculateTime();
	high_resolution_clock::time_point startTime;
	high_resolution_clock::time_point stopTime;
	Timer();
	~Timer();
};


#endif // !TIMER_H