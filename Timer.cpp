#include "pch.h"
#include "Timer.h"
#include "time.h"


Timer::Timer()
{
}
__int64 Timer::calculateTime() {

	return duration_cast<nanoseconds>(Timer::stopTime - Timer::startTime).count();

}

void Timer::timerStart() {
	startTime = high_resolution_clock::now();
}

void Timer::timerStop() {
	stopTime = high_resolution_clock::now();
}

Timer::~Timer()
{
}