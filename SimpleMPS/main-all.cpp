//
//  main-all.cpp
//  run all the steps, sans the external optimization problem
//  mainly for development and testing
//
//
//// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000,
//// there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
//// Export of this program may require a license from the United States Government.
////
////  Created by Scott Alan Mitchell on 4/09/15.
//  Copyright (c) 2015 Sandia National Labs darts team. All rights reserved.
//


// method submitted to Euro Graphics
// Author: Mohamed S. Ebeida
// Modified 03/01/2015, for remote sensing disk cover purposes
// Author: Scott A. Mitchell


// done for diskcover project:
// added a polygonal boundary
// define polygon boundary
//    read in polygon from file
//    draw in plot2d so we can debug and understand
//
// add "has boundary" flag to top level cells, list of edges of the boundary passing through this cell
//
// implement checking if a sample point is inside a parent cell that "has boundary"
// do this by project the query point in +y direction, count crossings. Define a segment as closed on its endpoint with the smaller x value, and open on the other endpoint. For vertical segments, always count them as two == zero crossings.
// implement checking if a sample point is inside a refined square that "has boundary" - check parent cell,
//
// set initial squares as inside/outside
//   pick a boundary vertex and find the cell it is in. Follow its edges to march around the polygon boundary, crossing from cell to cell.
//   These are marked as "has boundary"
// floodfill would speed stuff up but is not needed in this context =  mark each of the four cell vertices as inside or outside - for each cell edge that is inside-inside and outside-outside edge, flood to the edge-adjacent neighbors
//
// define footprint : polygon? arc-gon? implement "distance from point to segment" and "distance from point to arc"
// center footprint on sample points, determine which other samples are covered by (footprint-r) (closest point on line?)

// =========
// in separate stand-alone codes, do the individual steps:
//  sampling
//  determining coverage
//  displaying the solution

#include "diskcover.hpp"

// =========================

int main(int argc, char *argv[])
{
#pragma region Main:
  std::cout << "Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for use of this work by or on behalf of the U.S. Government. Export of this program may require a license from the United States Government." << std::endl;

  clock_t start_time, end_time; double cpu_time;
  start_time = clock();
  
  std::cout << std::endl << "Diskcover-all: creating domain sample points for coverage and placement, finding coverage, and plotting solution." << std::endl;

  // initialize RNG
  init_RNG(1);
  
  // sampling densities
  // _epsilon is the sampling density for the coverage samples
  // _delta is the sampling density for the center samples
  // these are absolute values, as is the domain coordinates
  // these all are scaled internally
  //  double _epsilon = 0.0585;
  //  double _epsilon = 0.00585;
  double _epsilon = 0.04; // 0.07;
  double _delta = _epsilon * 2.; // _epsilon * 4.; * 1.5 etc
  // assert( _delta >= _epsilon );
  const bool reuse_all_samples = false; // not implemented yet
  const bool reuse_nonconflict_samples = true; // true in production
  // if true, then the _epsilon samples for representation are used for footprint placement
  
  // define the footprint and its size
  // define the footprint and its size
  Footprint *_footprint = Footprint::read_footprint( "footA.txt" );
  assert(_footprint);
  // Footprint * _footprint = new FootprintCircle(0.25);
  // this is used for the frame and the offset
  
  Domain _domain;
  Grid _grid;
  _grid._fname_prefix = "points_for_coverage";
  
  
  _domain._frame_size = _footprint->max_radius() * 1.01; // use this in production
  // const double frame_size = _delta * 1.1 + 0.1;
  _domain.read( "polygon.txt" );
  
  // ============ representation sampling
  
  _domain._offset = 0.;
  // _domain._offset = 0.3; // debug
  _grid.setup( &_domain, _epsilon );
  _grid.MPS();
  
  //std::cout<< "\n*** Testing Output ... " << std::endl;
  //test_output();
  
  //plot();
  
  _grid.save_point_cloud("points_for_coverage.txt");
  _grid.save_point_cloud("points_for_coverage-debug.txt", true);
  
  
  // ============ placement sampling
  Grid _grid_placement;
  _grid_placement._fname_prefix = "points_for_placement";
  
  // farthest a footprint can be and still cover part of the domain
  // the offset is the second part where the footprint that affects the sampling, so just be conservative over all possible footprints
  _domain._offset = _footprint->max_radius() * _domain._scale;
  
  _grid_placement.setup( &_domain, _delta );
  
  // const double scaled_delta = _delta * _domain._scale;
  
  std::cout << "adding coverage samples conditionally " << std::endl;
  if (reuse_all_samples)
    _grid_placement.add_samples( &_grid );
  else if (reuse_nonconflict_samples)
    _grid_placement.add_samples_conditionally( &_grid );
  
  std::cout << "adding placement samples randomly " << std::endl;
  _grid_placement.MPS();
  
  _grid_placement.save_point_cloud("points_for_placement.txt");
  _grid_placement.save_point_cloud("points_for_placement-debug.txt", true);
  
  delete _footprint;
  
  // ============= which points are covered by which placements?
  find_footprint_coverage( "footA.txt" );
  
  
  // ======= run the LP
  
  
  // =========== display the final solution
  find_solution_coverage( "solution.txt" );
  
  // ============= finish
  
  std::cout<< "\n*** Mission Accomplished ***." << std::endl;
  
  // performance stats
  end_time = clock();
  cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
  std::cout << "Total runtime for all algorithms = " << cpu_time;
  
  // wait for keypress for debugging
//  char dummy_cc;
//  std::cout<< "\n*** enter q to quit ...";
//  std::cin >> dummy_cc;
  
  std::cout << std::endl;
  return 0;
  
#pragma endregion
};

