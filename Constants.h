#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <math.h>
# define M_PI 3.14159265358979323846 
using namespace std;
const double LargeM = 100000;
const double INF = 0x3f3f3f3f;
const static double EarthRaidumMiles = 3959; //miles;
const int threeOptPermutation[7][4] = {
	{0,1,3,2},
	{1,0,2,3},
	{1,0,3,2},
	{2,3,0,1},
	{2,3,1,0},
	{3,2,0,1},
	{3,2,1,0}
};

class LatLong {
public:
	double Latitude;
	double Longitude;
	LatLong(double rLatitude, double rLongitude) {
		this->Latitude = rLatitude;
		this->Longitude = rLongitude;
	}
};

inline double DegreeToRadian(double angle){return M_PI * angle / 180.0;}
inline double haversine(const double lat1, const double long1, const double lat2, const double long2) {
	double dist = 0;
	double latRad1 = DegreeToRadian(lat1);
	double latRad2 = DegreeToRadian(lat2);
	double lonRad1 = DegreeToRadian(long1);
	double lonRad2 = DegreeToRadian(long2);

	double diffLa = latRad2 - latRad1;
	double doffLo = lonRad2 - lonRad1;

	double computation = asin(sqrt(sin(diffLa / 2) * sin(diffLa / 2) + cos(latRad1) * cos(latRad2) * sin(doffLo / 2) * sin(doffLo / 2)));
	return 2 * EarthRaidumMiles * computation;

	return dist;
}
inline double haversine(LatLong l1, LatLong l2) { return haversine(l1.Latitude,l1.Longitude,l2.Latitude,l2.Longitude);}
#endif