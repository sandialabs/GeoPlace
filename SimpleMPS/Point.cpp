//
//  Points.cpp
//  simplemps
//
//
//// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000,
//// there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
//// Export of this program may require a license from the United States Government.
////
////  Created by Scott Alan Mitchell on 3/27/15.
//  Copyright (c) 2015 darts team. All rights reserved.
//

#include "Point.hpp"

// clockwise angle at a between ab and ac,
// in [-pi, pi]
double cw_angle( const double *a, const double *b, const double *c )
{
  double ab[2];
  double ac[2];
  point_minus(ab, b,a);
  point_minus(ac, c,a);
  
  double dot = point_dot(ab, ac); // dot product, proportional to cos
  double det = point_determinant(ab, ac); // determinant, proportional to sin
  // atan2(y, x) or atan2(sin, cos)
  double angle = (dot == 0. && det == 0.) ? M_PI : atan2(det, dot);
  
  assert( angle >= -M_PI );
  assert( angle <=  M_PI );
  
  return angle;
}

double Arc::distance_to_arc( const double *p ) const
{
  // angle p makes at c
  double x[2];
  x[0] = _center[0];
  x[1] = _center[1] + 1.;
  double angle_p = cw_angle( _center, x, p );
  
  // convert to 0 to 2pi range
  angle_p = zero_to_2pi( angle_p );
  const double theta_A = zero_to_2pi( _thetaA );
  const double theta_B = zero_to_2pi( _thetaB );
  
  const bool in_arc = (theta_B > theta_A) ? ( angle_p <= theta_A || angle_p >= theta_B ) : (angle_p <= theta_A && angle_p >= theta_B );
  
  if (in_arc)
  {
    const double pc2 = distance_squared( _center, p );
    const double d = sqrt( pc2 ) - _r;
    return d;
  }
//  else
  {
    double a[2], b[2];
    angle_to_point( a, _thetaA );
    angle_to_point( b, _thetaA );
    const double pa2 = distance_squared(p,a);
    const double pb2 = distance_squared(p,b);
    const double d = sqrt( ( pa2 < pb2 ) ? pa2 : pb2 );
    return d;
  }
}

void Arc::angle_to_point( double *p, double theta ) const
{
  point_assign(p, _center);
  p[0] += _r * cos(theta);
  p[1] += _r * sin(theta);
}

double Arc::zero_to_2pi( double angle ) const
{
  const double two_pi = 2 * M_PI;
  while ( angle >= two_pi )
    angle -= two_pi;
  while ( angle < 0 )
    angle += two_pi;
  assert( angle >= 0. );
  assert( angle < two_pi );
  return angle;
}
