/*
 * Timer.h
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <cstdlib>
#include <time.h>
#include <windows.h>

class Timer {
public:
	Timer() { m_start = timestamp(); }
	void restart() { m_start = timestamp(); }
	long long elapsed() { return timestamp() - m_start; }
	int gettimeofday(struct timeval *tp, void *tzp)
	{
		time_t clock;
		struct tm tm;
		SYSTEMTIME wtm;
		GetLocalTime(&wtm);
		tm.tm_year = wtm.wYear - 1900;
		tm.tm_mon = wtm.wMonth - 1;
		tm.tm_mday = wtm.wDay;
		tm.tm_hour = wtm.wHour;
		tm.tm_min = wtm.wMinute;
		tm.tm_sec = wtm.wSecond;
		tm.tm_isdst = -1;
		clock = mktime(&tm);
		tp->tv_sec = clock;
		tp->tv_usec = wtm.wMilliseconds * 1000;
		return (0);
	}

private:
	long long m_start;

	// Returns a timestamp ('now') in microseconds
	long long timestamp() {
		struct timeval tp;
		gettimeofday(&tp, nullptr);
		return ((long long)(tp.tv_sec)) * 1000000 + tp.tv_usec;
	}
};

#endif /* TIMER_H_ */
