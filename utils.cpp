#ifndef __UTILS_CPP
#define __UTILS_CPP

#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>
#include "dirac.h"
#include "dirac.cpp"
#include "utils.h"
using namespace Utils;
#pragma region converter
//#include <cstddef>;
//#include <string.h>;
//
//char * Convert(int from, int to, const char * s, char * out)
//{
//	if (s == NULL)
//		return NULL;
//
//	if (from < 2 || from > 36 || to < 2 || to > 36) { return NULL; }
//
//	int il = strlen(s);
//
//	int *fs = new int[il];
//	int k = 0;
//	int i, j;
//
//
//	for (i = il - 1; i >= 0; i--)
//	{
//		if (s[i] >= '0' && s[i] <= '9')
//		{
//			fs[k] = (int)(s[i] - '0');
//		}
//		else
//		{
//			if (s[i] >= 'A' && s[i] <= 'Z')
//			{
//				fs[k] = 10 + (int)(s[i] - 'A');
//			}
//			else if (s[i] >= 'a' && s[i] <= 'z')
//			{
//				fs[k] = 10 + (int)(s[i] - 'a');
//			}
//			else
//			{
//				delete[]fs;
//				return NULL;
//			} //only allow 0-9 A-Z characters
//		}
//		k++;
//	}
//
//	for (i = 0; i<il; i++)>
//	{
//		if (fs[i] >= from)
//			return NULL;
//	}
//
//	double x = ceil(log(from) / log(to));
//	int ol = 1 + (il * x);
//
//	int * ts = new int[ol];
//	int * cums = new int[ol];
//
//	for (i = 0; i<ol; i++)>
//	{
//		ts[i] = 0;
//		cums[i] = 0;
//	}
//	ts[0] = 1;
//
//
//	//evaluate the output
//	for (i = 0; i < il; i++) //for each input digit
//	{
//		for (j = 0; j < ol; j++) //add the input digit times (base:to from^i) to the output cumulator
//		{
//			cums[j] += ts[j] * fs[i];
//			int temp = cums[j];
//			int rem = 0;
//			int ip = j;
//			do // fix up any remainders in base:to
//			{
//				rem = temp / to;
//				cums[ip] = temp - rem * to;
//				ip++;
//				if (ip >= ol)
//				{
//					if (rem > 0)
//					{
//						delete[]ts;
//						delete[]cums;
//						delete[]fs;
//						return NULL;
//					}
//					break;
//				}
//				cums[ip] += rem;
//				temp = cums[ip];
//			} while (temp >= to);
//		}
//
//		for (j = 0; j < ol; j++)
//		{
//			ts[j] = ts[j] * from;
//		}
//
//		for (j = 0; j < ol; j++) //check for any remainders
//		{
//			int temp = ts[j];
//			int rem = 0;
//			int ip = j;
//			do  //fix up any remainders
//			{
//				rem = temp / to;
//				ts[ip] = temp - rem * to;
//				ip++;
//				if (ip >= ol)
//				{
//					if (rem > 0)
//					{
//						delete[]ts;
//						delete[]cums;
//						delete[]fs;
//						return NULL;
//					}
//					break;
//				}
//				ts[ip] += rem;
//				temp = ts[ip];
//			} while (temp >= to);
//		}
//	}
//
//	if (out == NULL)
//	{
//		out = (char*)malloc(sizeof(char) * (ol + 1));
//	}
//
//	int spos = 0;
//	bool first = false; //leading zero flag
//	for (i = ol - 1; i >= 0; i--)
//	{
//		if (cums[i] != 0)
//		{
//			first = true;
//		}
//		if (!first)
//		{
//			continue;
//		}
//
//		if (cums[i] < 10)
//		{
//			out[spos] = (char)(cums[i] + '0');
//		}
//		else
//		{
//			out[spos] = (char)(cums[i] + 'A' - 10);
//		}
//		spos++;
//	}
//	out[spos] = 0;
//
//	delete[]ts;
//	delete[]cums;
//	delete[]fs;
//
//	return out;
//}

#pragma endregion

//// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
//const std::string currentDateTime() {
//	time_t     now = time(0);
//	struct tm  tstruct;
//	char       buf[80];
//	tstruct = *localtime(&now);
//	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
//	// for more information about date/time format
//	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
//
//	return buf;
//}


////This will apply to both photon number n and transmon level j
//template<typename T>
//Operator<T> Utils::annihilation(T n) {
//	std::deque<Projector<T>> annihilation_projs;
//
//	for (T k = 0; k < n + 1;  k++) {
//		Projector<T> proj({ k + 1 }, { k }, sqrt(k + 1));
//		annihilation_projs.push_back(proj);
//	}
//	Operator<T> a(n + 2, annihilation_projs);
//
//	return a;
//}
//
//template<typename T>
//Operator<T> Utils::creation(T n) {
//	Operator<T> a_dagger = annihilation(n).dagger();
//
//	return a_dagger;
//}

#endif

