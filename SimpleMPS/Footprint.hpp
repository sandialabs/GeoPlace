//
//  Footprint.hpp
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

#ifndef simplemps_Footprint_hpp
#define simplemps_Footprint_hpp

#include <fstream>

// abstract class defining the interface
class Footprint
{
protected:
  double _scale;
public:
  // the maximum distance from the toe to the farthest covered point
  virtual double max_radius() = 0;
  // is point p covered by a footprint with the given toe point
  // assumes both are in the same coordinate system, as the scaling of the footprint
  virtual bool is_covered(const double *toe, const double *p) = 0;
  
  virtual void plot( std::fstream &file, double graphics_scale, double toe_x, double toe_y ) = 0;
  
  // scale the footprint by the scale factor
  virtual void set_scale( double scale ) { _scale = scale; }
  
  // >= 0, the amount a point has to be inside the footprint to be considered covered
  double _offset;
  
  Footprint() : _scale(0.), _offset(0.) {}
  virtual ~Footprint() {}
  
  static Footprint *read_footprint( std::string fname );

protected:
  virtual void parse_footprint( std::ifstream &file ) = 0;

};

// footprint types
class FootprintCircle : public Footprint
{
private:
  // scaled radius, not original radius
  double _radius;
  double _radius_squared;
public:
  FootprintCircle( double radius = 1. ) : _radius(radius), _radius_squared( radius * radius ) {}
  
  // the maximum distance from the toe to the farthest covered point
  virtual double max_radius()
  { return _radius; }

  // is point p covered by a footprint with the given toe point
  bool is_covered(const double *toe, const double *p);
  
  virtual void set_scale( double scale )
  {
    Footprint::set_scale(scale);
    _radius *= scale;
    _radius_squared = _radius * _radius;
  }

  virtual void plot( std::fstream &file, double graphics_scale, double toe_x, double toe_y );

protected:
  virtual void parse_footprint( std::ifstream &file );

};

#endif
