//
//  Points.hpp
//  simplemps
//
//
//// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000,
//// there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
//// Export of this program may require a license from the United States Government.
////
////  Created by Scott Alan Mitchell on 3/27/15.
//  Copyright (c) 2015 Sandia National Labs darts team. All rights reserved.
//

#ifndef simplemps_Points_hpp
#define simplemps_Points_hpp


#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <cmath>
#include <assert.h>


// ============= functions

void print_point(double *p, std::ostream &out = std::cout); // coordinates

// Point methods
// consider using " double p[2]; " instead of new and delete
double *new_point();  // allocate memory
void delete_point(double *&p); // deallocate
void point_assign(double *p, const double *a); // p = a
void point_minus(double *p, const double *a, const double *b);  // p = a - b
void point_add(double *p, const double *a, const double *b); // p = a + b
void point_add( double *p, double scalar); // add scalar to each of p's coordinates
void point_scale( double *p, double scale );    // p *= scale;
void point_midpoint( double *midpoint, const double *a, const double *b);
double point_dot( const double *a, const double *b ); // d = a * b
double norm_squared( const double *a ); // || a ||^2
double norm( const double *a );
double point_determinant(const double *a, const double *b);

double distance_squared( const double *a, const double *b ); // || a - b ||^2

// distance from p to the edge ab
double distance_to_edge_squared( double *p, double *a, double *b );

double cw_angle( const double *a, const double *b, const double *c );

// ========================
// definition of an Arc

class Arc
{
public:
  // arc has center of circle c, with radius r, and extends clockwise from thetaA to thetaB
  double _center[2];
  double _r, _thetaA, _thetaB;
  
  // distance from point p to the closest point of the arc
  double distance_to_arc( const double *p ) const;

  // convert angle theta to point coordinates p
  void angle_to_point( double *p, double theta ) const;

  // convert angle to [0,2pi)
  double zero_to_2pi( double angle ) const;


};

// =================== implementations



inline
double *new_point()
{
  return new double[2];
}
inline
void delete_point(double *&p)
{
  delete [] p;
  p = 0;
}
inline
void point_assign(double *p, const double *q)
{
  p[0] = q[0];
  p[1] = q[1];
}
inline
void point_minus(double *p, const double *a, const double *b)
{
  p[0] = a[0] - b[0];
  p[1] = a[1] - b[1];
}
inline
void point_add(double *p, const double *a, const double *b)
{
  p[0] = a[0] + b[0];
  p[1] = a[1] + b[1];
}
inline
void point_add( double *p, double scalar )
{
  p[0] += scalar;
  p[1] += scalar;
}
inline
void point_midpoint( double *midpoint, const double *a, const double *b)
{
  midpoint[0] = (a[0] + b[0]) / 2.;
  midpoint[1] = (a[1] + b[1]) / 2.;
}
inline
void point_scale( double *p, double scale )
{
  p[0] *= scale;
  p[1] *= scale;
}
inline
double point_dot( const double *a, const double *b )
{
  return a[0] * b[0] + a[1] * b[1];
}
inline
double norm_squared( const double *a )
{
  return a[0] * a[0] + a[1] * a[1];
}
inline
double norm( const double *a )
{
  return sqrt( point_dot(a,a) );
}
inline
double point_determinant(const double *a, const double *b)
{
  return a[0]*b[1] - a[1]*b[0];
}
inline
double distance_squared( const double *a, const double *b )
{
  return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]);
}

#endif
