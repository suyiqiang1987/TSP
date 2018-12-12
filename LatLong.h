#pragma once
#ifndef LATLONG_H
#define LATLONG_H
class LatLong {
private:
	double latitude;
	double longitude;
public:
	LatLong(double rLatitude, double rLongitude) {
		this->latitude = rLatitude;
		this->longitude = rLongitude;
	}
	double Latitude() { return this->latitude; }
	double Longitude() { return this->longitude; }
};
#endif
