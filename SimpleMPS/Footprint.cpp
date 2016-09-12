//
//  Footprint.cpp
//  simplemps
//
//
//// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000,
//// there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
//// Export of this program may require a license from the United States Government.
////
//  Created by Scott Alan Mitchell on 3/27/15.
//  Copyright (c) 2015 darts team. All rights reserved.
//

#include "Footprint.hpp"
#include "Point.hpp"
#include "IO.hpp"

void FootprintCircle::plot( std::fstream &file, double graphics_scale, double toe_x, double toe_y )
{
  // outer radius
  set_color(file, green);
  plot_circle( file, graphics_scale, toe_x, toe_y, _radius, false /*true*/ );
  set_color(file, dark_blue);
  // inner radius
  plot_circle( file, graphics_scale, toe_x, toe_y, _radius - _offset, false );
  // center point
  plot_circle( file, graphics_scale, toe_x, toe_y, 0.01, false );
  set_color(file);
}

bool FootprintCircle::is_covered(const double *toe, const double *p)
{
  const double d2 = distance_squared( toe, p );
  return ( sqrt(d2) < _radius - _offset );
}


Footprint *Footprint::read_footprint( std::string fname )
{
  std::ifstream file( fname.c_str(), std::ios::in );
  std::string footprint_name;
  getline( file, footprint_name );
  
  // create
  Footprint * foot(0);
  if ( footprint_name.compare("circle") == 0 )
  {
    foot = new FootprintCircle();
  }
//  else if ( footprint_name.compare("name") )
//  {
//    
//  }
  else
  {
    std::cerr << "unknown footprint type " << footprint_name << std::endl;
    return 0;
  }
  
  // finish parsing it in the appropriate virtual function
  foot->parse_footprint(file);
  
  return foot;
}


void FootprintCircle::parse_footprint( std::ifstream &file )
{
  double scale;
  file >> scale;
  set_scale( scale );
}
