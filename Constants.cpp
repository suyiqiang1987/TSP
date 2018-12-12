#include "stdafx.h"
#include "Constants.h"


double const Constants::my_PI = 3.14159265358979323846;
double const  Constants::LargeM = 100000;
double const  Constants::INF = 0x3f3f3f3f;
double const   Constants::EarthRaidumMiles = 3959; //miles;
int const  Constants::threeOptPermutation[7][4] = {
	{ 0,1,3,2 },
	{ 1,0,2,3 },
	{ 1,0,3,2 },
	{ 2,3,0,1 },
	{ 2,3,1,0 },
	{ 3,2,0,1 },
	{ 3,2,1,0 }
};
