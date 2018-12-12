#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "LatLong.h"
#include <math.h>
using namespace std;

class Constants {
public:
	static double const my_PI;
	static double const  LargeM;
	static double const  INF;
	static double const   EarthRaidumMiles; //miles;
	static int const  threeOptPermutation[7][4];
	
	static double DegreeToRadian(double angle) { return my_PI * angle / 180.0; }
	static double haversine(const double lat1, const double long1, const double lat2, const double long2) {
		double latRad1 = DegreeToRadian(lat1);
		double latRad2 = DegreeToRadian(lat2);
		double lonRad1 = DegreeToRadian(long1);
		double lonRad2 = DegreeToRadian(long2);

		double diffLa = latRad2 - latRad1;
		double doffLo = lonRad2 - lonRad1;

		double computation = asin(sqrt(sin(diffLa / 2) * sin(diffLa / 2) + cos(latRad1) * cos(latRad2) * sin(doffLo / 2) * sin(doffLo / 2)));
		return 2 * EarthRaidumMiles * computation;
	}
	static double haversine(LatLong l1, LatLong l2) 
	{ return haversine(l1.Latitude(), l1.Longitude(), l2.Latitude(), l2.Longitude()); }
};
#endif
