// diskcover.cpp

//
//// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000,
//// there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
//// Export of this program may require a license from the United States Government.
////
////  Copyright (c) 2015 Sandia National Labs darts team. All rights reserved.
//
// method submitted to Euro Graphics
// Author: Mohamed S. Ebeida
// Modified 03/01/2015, for remote sensing disk cover purposes
// modification
// Author: Scott A. Mitchell

// local

#include "diskcover.hpp"
#include <numeric>

// ================ random number generator ============
class RandomNumberGenerator
{
private:
  size_t _MY_RAND_MAX;
public:
  RandomNumberGenerator()
  {
    unsigned int seed = 1425430940; // fixed
    // unsigned int seed = (unsigned int) size_t(time(0)); // random
    // std::cout << "default random seed = " << seed << std::endl;
    srand(seed);
    
    // Constants
    _MY_RAND_MAX = RAND_MAX; // 15 bits
    _MY_RAND_MAX = _MY_RAND_MAX << 15;
    _MY_RAND_MAX += RAND_MAX; // it's now 30 bits
  }
  
  // used a lot
  double generate_a_random_number()
  {
    // samitch: I think we should be using KISS64 instead of "rand"
    
    // put two 15 bit random numbers together into one 30 bit random number
    size_t num(rand());
    num = num << 15;
    num += rand(); // another random number
    return static_cast<double>(num * 1.0 / _MY_RAND_MAX);
  }
};

RandomNumberGenerator THE_RNG;

void init_RNG(const unsigned int seed)
{
  std::cout << "setting random seed = " << seed << std::endl;
  srand(seed);
}

// =====================================================

double * SamplingIterator::next()
{
  const size_t max_index = _sampling->cell_points_size();
  if (_celli >= max_index )
    return 0;
  while (_sampling->_cell_points[_celli++] == 0)
  {
    if (_celli >= max_index)
      return 0;
  }
  return _sampling->_cell_points[_celli - 1];
}

void Sampling::new_cell_points(size_t cell_points_size)
{
  // get rid of old memory, if any
  if (_cell_points_size)
    delete_cell_points();
  
  _cell_points_size = cell_points_size;
  _cell_points = new double*[_cell_points_size];
  for (size_t i = 0; i < _cell_points_size; ++i)
  {
    _cell_points[i] = 0;
  }
}

void Sampling::delete_cell_points()
{
  for (size_t i = 0; i < _cell_points_size; ++i)
  {
    delete [] _cell_points[i];
  }
  delete [] _cell_points;
  _cell_points_size = 0;
  _cell_points = 0;
  _num_inserted_points = 0;
}

Grid::~Grid()
{
  // clear memory
  
  // active pool
  
  // active children
  for (size_t ichild = 0; ichild < _num_max_children; ichild++)
  {
    delete[] _active_children[ichild];
  }
  delete[] _active_children;
  
  delete[] _i;
  delete[] _j;
  delete[] _u;
  delete[] _x;
}



// ================= implementation

Domain::~Domain()
{
  for (size_t v = 0; v < _vertices.size(); ++v)
    delete [] _vertices[v];
  _vertices.clear();
  for (size_t p = 0; p < _loops.size(); ++p)
  {
    Polygon *poly = _loops[p];
    delete poly;
  }
  _loops.clear();
}
double Polygon::calculate_angle_sum( const double * p, size_t first_edge )
{
  double angle_sum = 0;
  for (size_t e = first_edge; e < _edges.size(); ++e )
  {
    double angle = cw_angle(p, v0(e), v1(e));
    // if (angle > M_PI) angle -= 2. * M_PI;
    angle_sum += angle;
  }
  return angle_sum;
}

bool Polygon::is_point_geometrically_inside( const double *p )
{
  // winding number test
  // some abiguity of what it will return for very close to the boundary, but robust that it will return an answer, consistently
  // vs. projection tests requiring special case handling
  
  double angle_sum = calculate_angle_sum(p);
  
  // 0 is outside, 2pi is inside
  bool geometrically_inside = fabs( angle_sum ) > M_PI ;
  
  return geometrically_inside;
}

bool Polygon::is_point_inside( const double *p )
{
  bool geometrically_inside = is_point_geometrically_inside(p);
  
  // xor, with whether the polygon defines and inner or outer boundary
  bool logically_inside = ( geometrically_inside == is_outer_boundary() );
  
  return logically_inside;
}

bool Polygon::is_outer_boundary()
{
  if (_edges.empty())
    return false;
  
  // calculate and store the first time
  if (_is_outer_boundary == 0) // is unset
  {
    double midpoint[2];
    point_midpoint( midpoint, v0(0), v1(0));
    
    double angle_sum = calculate_angle_sum( midpoint, 1 );
    
    // -PI for an outer boundary, +PI for an inner boundary
    _is_outer_boundary = angle_sum < 0. ? 1 : -1; // set
  }
  return (_is_outer_boundary > 0);
}

bool Domain::is_point_inside( const double *p )
{
  if (_loops.empty())
    return true;
  
  bool in_outer = false;
  bool out_inner = true;
  
  for (size_t el = 0; el < _loops.size(); ++el )
  {
    Polygon *poly = _loops[el];
    if (poly->is_outer_boundary())
    {
      // if we finished a set of loops defining a connected part of the domain,
      // and were inside it, then return true
      if (in_outer && out_inner)
        return true;
      in_outer = poly->is_point_geometrically_inside(p);
      out_inner = true; // so far
    }
    else
    {
      out_inner = out_inner && !poly->is_point_geometrically_inside(p);
    }
  }
  return (in_outer && out_inner);
}

bool Grid::is_cell_outside( const size_t* i )
{
  //  // debug
  //  if (i[0] == 12)
  //    std::cout << " x = 12 " << std::endl;
  
  double lower[2];
  get_childcell_lower(i, lower);
  return _domain->is_cell_outside(lower, current_diagonal() );
}

bool Domain::is_cell_outside( double* ll, double ss )
{
  // ll == lower left coordinates of cell, ur == upper right
  // ss == side length
  
  // if the cell contains a domain vertex
  // then it is not outside
  for (size_t v = 0; v < _vertices.size(); ++v)
  {
    // domain vertex coordinates
    double *p = _vertices[v];
    
    if (! ( p[0] < ll[0] || p[0] > ll[0] + ss ||
            p[1] < ll[1] || p[1] > ll[1] + ss  ) )
    {
      // p is in rectangular cell
      return false;
    }
  }
  
  // now cell doesn't contain a vertex, but could be pierced by two edges
  //
  // check the distance from the center to the domain.
  // if that minus the cell diagonal is bigger than offset, then the cell is outside for sure.
  double center[2];
  point_assign(center, ll);
  const double cell_radius = ss / 2.;
  point_add( center, cell_radius );
  
  const double cd = signed_distance_to_boundary(center);
  const bool too_far = ( cd > _offset + sqrt(2) * cell_radius );
  
  if (DEBUG_BOUNDARY)
  {
    std::cout << "cell center (" << center[0] << " , " << center[1] << ") at distance " << cd << " vs. " << _offset + cell_radius << " is ";
    if ( too_far )
      std::cout <<" too far from the domain boundary to contain an interior point " << std::endl;
    else
      std::cout <<" AMBIGOUSLY CLOSE to the boundary" << std::endl;
  }
  
  if (too_far)
    return true;
  
  // inside convex hull?
  if (_convex_hull.is_point_geometrically_inside(center))
    return false;
  
  // far outside convex hull?
  if (_convex_hull.distance_to_boundary_squared(center) > ss * ss * 0.5)
    return true;
  
  return false; // indeterminant
  
}

// 100 is big since the box is scaled to 0-1
const double infinite_distance = 100;

// these distances are positive
double distance_to_edge_squared( double *p, double *a, double *b )
{
  // distance from p to line ab
  // if the projection of p to ab is beyond a and b;
  // then return a big number since we check distances to vertices elsewhere
  
  double ab[2];
  point_minus(ab, b, a); //ab = b-a
  const double nab2 = norm_squared(ab); // || ab ||^2
  double ap[2];
  point_minus(ap, p, a); //ap = p-a
  const double ab_dot_ap = point_dot(ab, ap);
  // if closest point is a, return infinity
  if (ab_dot_ap < 0.)
    return infinite_distance;
  const double d2 = ab_dot_ap * ab_dot_ap / nab2; // d^2
  // if closest point is b, return infinity
  if ( d2 > nab2 )
    return infinite_distance;
  
  const double nap2 = norm_squared(ap); // || ap ||^2
  const double e2 = (nap2 - d2 ); //  e^2 = f^2 - d^2
  
  // we could get the signed distance to the edge if we wanted to, but
  // that doesn't help in the case that the closest distance is a vertex
  return e2;
}

// these distances are positive
double Polygon::distance_to_boundary_squared(double *p)
{
  // distance is min distance to an edg of the polygon
  double min_distance2 = infinite_distance;
  for (size_t e = 0; e < _edges.size(); ++e)
  {
    const double d2 = distance_to_edge_squared( p, v0(e), v1(e) );
    if ( d2 < min_distance2 )
      min_distance2 = d2;
  }
  return min_distance2;
}

void Polygon::print()
{
  std::cout << "Polygon:\n  Edges:" << std::endl;
  for (size_t e = 0; e < _edges.size(); ++e)
    std::cout << "    " << _edges[e].first << "  -  " << _edges[e].second << std::endl;
  std::cout << "  Vertices:" << std::endl;
  for (size_t e = 0; e < _edges.size(); ++e)
    print_point( (*_vertices)[ _edges[e].first ] );
  std::cout << std::endl;
}

void Polygon::convex_hull()
{
  _is_outer_boundary = 1; // normal polygon, not a hole
  
  // find vertex with minimum-y ( lex min x )
  unsigned minp = 0;
  for (unsigned p=1; p<_vertices->size(); ++p)
  {
    const double *pp    = (*_vertices)[p];
    const double *pminp = (*_vertices)[minp];
    if (pp[1] < pminp[1] || (pp[1] == pminp[1] && pp[0] < pminp[0]) )
      minp = p;
  }
  
  
  // sort by angle
  // vector of indices into the vertices
  std::vector<unsigned> v(_vertices->size());
  std::iota( v.begin(), v.end(), 0); // sequential integers 0..size-1
  std::sort( v.begin(), v.end(), less_than_angle(_vertices, minp) );
  
  // search for minimum y-vertex in sorted list
  int start = 0;
  while (v[start] != minp )
    ++start;
  assert( v[start] == minp );
  
  // start is now the index of minp

  // hull vertices in order
  std::vector<unsigned> h;
  h.reserve(_vertices->size()+1);
  {
    int i = start;
    // ccw traversal, keeping only convex steps
    do
    {
      unsigned next = v[i];
      bool addnext = false;
      if ( h.size() < 2 )
      {
        addnext = true;
      }
      else
      {
        const unsigned curr = h[ h.size() -1 ];
        const unsigned prev = h[ h.size() -2 ];
        // check angle
        const double cw = cw_angle( (*_vertices)[curr], (*_vertices)[prev], (*_vertices)[next] );
        addnext = (cw >= 0. );
      }
      
      if (addnext)
      {
        h.push_back(next);
        
        // increment i
        if (--i < 0) i = (int)v.size() - 1;
      }
      else
      {
        // pop current, try again without incrementing i
        h.pop_back();
      }
    }  while ( i != start );
  }

  // convert h into edges
  _edges.clear();
  _edges.reserve(h.size());
  // duplicage last entry of h to simplify loop
  h.push_back( h.front() );
  for (size_t i = 0; i < h.size()-1; ++i)
    _edges.push_back( Edge(h[i], h[i+1]) );
  
  // debug
  std::cout << "Convex Hull" << std::endl;
  print();
}


// these distances are positive
double Domain::distance_to_boundary(double *p)
{
  // distance is min distance to a polygon
  // if there are vertices with no edges, then we must also check the distance to individual points
  // we check that first, as an optimization, anyway so the lower level routines do not need to
  // compute them.
  
  // empty domain?
  if (_vertices.empty())
    return infinite_distance;
  
  // min distance to a vertex
  double min_distance2 = infinite_distance;
  for (size_t v = 0; v < _vertices.size(); ++v)
  {
    double *q = _vertices[v];
    const double d2 = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]);
    if (d2 < min_distance2)
      min_distance2 = d2;
  }
  
  // min distance to an edge
  for (size_t el = 0; el < _loops.size(); ++el )
  {
    const double d2 = _loops[el]->distance_to_boundary_squared(p);
    if (d2 < min_distance2)
    {
      min_distance2 = d2;
    }
  }
  
  return sqrt(min_distance2);
}

double Domain::signed_distance_to_boundary(double *p)
{
  // - if inside, + if outside, 0 if on
  
  double d = distance_to_boundary(p);
  if (d == 0.)
    return 0.;
  if (is_point_inside(p))
    return -d;
  return d;
}

// ============================================================

void Grid::report_stats(double cpu_time)
{
  std::cout<< "==================================================================" << std::endl;
  std::cout<< "Number of inserted points = " << _num_inserted_points << std::endl;
  std::cout<< "Number of thrown darts = " << _num_darts << std::endl;
  std::cout<< "Execution Time = " << cpu_time << " seconds." << std::endl;
  std::cout<< "Remaining Cells Ratio = " << _num_valid * 1.0 / _no << std::endl;
  std::cout<< "Point Insertion Rate = " << _num_inserted_points * 1.0 / cpu_time << " points / sec" << std::endl;
  std::cout<< "Point Insertion Rate = " << _num_inserted_points * 1.0 / _no << " points / cell" << std::endl;
  std::cout<< "==================================================================" << std::endl;
}

// convert the indices for each dimension (first argument) into  a linear cell index (second argument)
inline void Grid::get_cell_center(const size_t* i, double* center)
{
  center[0] = _xo + ( (double) i[0] + 0.5 ) * _s;
  center[1] = _xo + ( (double) i[1] + 0.5 ) * _s;
}
inline void Grid::get_childcell_center(const size_t* i, double* center)
{
  center[0] = _xo + ( (double) i[0] + 0.5 ) * _ss;
  center[1] = _xo + ( (double) i[1] + 0.5 ) * _ss;
}
inline void Grid::get_childcell_lower(const size_t* i, double* lower)
{
  lower[0] = _xo + ( (double) i[0] ) * _ss;
  lower[1] = _xo + ( (double) i[1] ) * _ss;
}
inline void Grid::get_childcell_upper(const size_t* i, double* upper)
{
  upper[0] = _xo + ( (double) (i[0] + 1) ) * _ss;
  upper[1] = _xo + ( (double) (i[1] + 1) ) * _ss;
}
inline void Grid::get_cell_lower(const size_t* i, double* lower)
{
  lower[0] = _xo + ( (double) i[0] ) * _s;
  lower[1] = _xo + ( (double) i[1] ) * _s;
}
inline void Grid::get_cell_upper(const size_t* i, double* upper)
{
  upper[0] = _xo + ( (double) (i[0] + 1) ) * _s;
  upper[1] = _xo + ( (double) (i[1] + 1) ) * _s;
}

inline void Grid::get_cell_index(const size_t* i, size_t &icell)
{
  icell = i[1] + _n_pow[0] *  i[0];
};

inline void Grid::get_cell_indices(size_t icell, size_t* i)
{
  i[0] = icell / _n_pow[0];
  icell -= i[0] * _n_pow[0];
  i[1] = icell;
};

inline void Grid::get_cell_indicies( const double *p, size_t *i )
{
  // _s is side length
  {
    i[0] = floor( (p[0] - _xo) / _s );
    i[1] = floor( (p[1] - _xo) / _s );
  }
}

// mark cells by the  domain boundary
void Grid::initialize_cells_by_domain( Domain *domain )
{
  // uses global data:
  //   _domain, _s
  // sets global data:
  //   _boundary_cells, _invalid_cell
  
  // Domain can be empty of polygons, but can not be NULL
  assert(domain);
  _domain = domain;
  
  _domain->set_convex_hull();
  
  size_t i[2]; // cell indices
  size_t icell; // linear cell index
  
  // vertices, mark containing cells
  for (size_t v = 0; v < _domain->_vertices.size(); ++v)
  {
    // domain vertex coordinates
    double *p = _domain->_vertices[v];
    
    // find top level cell containing vert
    get_cell_indicies( p, i );
    get_cell_index( i, icell);
    
    // mark it as containing the boundary
    _boundary_cells[icell] = true;
    
    // debug
    if (DEBUG_BOUNDARY)
    {
      print_cell(i);
      std::cout << " contains ";
      _domain->print_vertex(v);
      std::cout << std::endl;
    }
  }
  
  // loops, walk edges
  for (size_t p = 0; p < _domain->_loops.size(); ++p)
  {
    Polygon *poly = _domain->_loops[p];
    for (size_t e = 0; e < poly->_edges.size(); ++e)
    {
      double *v0 = poly->v0(e);
      double *v1 = poly->v1(e);
      
      // find top level cell containing v0
      get_cell_indicies( v0, i );
      get_cell_index( i, icell);
      
      if (DEBUG_BOUNDARY)
      {
        std::cout << "Walking edge " << e << " between ";
        _domain->print_vertex(poly->_edges[e].first);
        std::cout << " in ";
        print_cell(i);
        std::cout << " and ";
        _domain->print_vertex(poly->_edges[e].second);
        std::cout << std::endl;
      }
      
      // t is fraction of edge from v0 to v1
      double t = 0.;
      const double dx = v1[0] - v0[0];
      const double dy = v1[1] - v0[1];
      // is x,y increasing or decreasing?
      const int dxi = dx > 0. ? 1 : -1;
      const int dyi = dy > 0. ? 1 : -1;
      // is the next crossing on the lower or upper edge of the next square?
      const int side_x = dx > 0. ? 0 : +1;
      const int side_y = dy > 0. ? 0 : +1;
      
      int current_xi = (int) i[0];
      int current_yi = (int) i[1];
      //      double current_x = current_xi * _s;
      //      double current_y = current_yi * _s;
      const double big_t = 999;
      
      // increase t, crossing grid edges, until the far vertex is reached
      while (t <= 1.)
      {
        // find next value of t that crosses an x grid edge, and a y grid edge.
        double next_xi = current_xi + dxi;
        double next_yi = current_yi + dyi;
        double next_x = _xo + (next_xi + side_x) * _s;
        double next_y = _xo + (next_yi + side_y) * _s;
        
        double next_xt = fabs(dx) > 0. ? (next_x - v0[0]) / dx : big_t;
        double next_yt = fabs(dy) > 0. ? (next_y - v0[1]) / dy : big_t;
        assert( next_xt >= 0. );
        assert( next_yt >= 0. );
        
        // increment either the x or y cell index
        if (next_xt < next_yt)
        {
          t = next_xt;
          current_xi = next_xi;
        }
        else
        {
          t = next_yt;
          current_yi = next_yi;
        }
        
        // mark that cell as having the boundary
        if (t <= 1.)
        {
          size_t next_i[2];
          next_i[0] = current_xi;
          next_i[1] = current_yi;
          get_cell_index( next_i, icell);
          _boundary_cells[ icell ] = true;
          
          if (DEBUG_BOUNDARY)
          {
            std::cout << "Edge cuts through ";
            print_cell( next_i );
            std::cout << " fraction " << t;
            std::cout << std::endl;
          }
        }
        else if (DEBUG_BOUNDARY)
          std::cout << "end of edge, fraction " << t << std::endl;
      }
    }
  }
  
  // future optimization: flood fill
  
  // set all non-boundary cells as being inside or outside the boundary,
  // depending on whether their center point lies inside or outside the boundary
  double center[2];
  for (size_t icell = 0; icell < _n; icell++)
  {
    if( _boundary_cells[ icell ] )
    {
      if (DEBUG_BOUNDARY)
      {
        get_cell_indices(icell, i);
        get_cell_center(i, center);
        std::cout << "bdy      ";
        std::cout << icell << " ";
        print_cell( i );
        std::cout << std::endl;
      }
      
      continue;
    }
    
    if (! _invalid_cells[icell] )
    {
      get_cell_indices(icell, i);
      
      // interior cell?
      double center[2];
      get_cell_center(i, center);
      assert( domain->_offset >= 0.); // below is not correct for negative offsets
      const bool is_interior = domain->is_point_inside( center );
      
      bool is_outside = !is_interior;
      if ( is_outside )
      {
        // if there is no offset, then testing whether the center is inside or outside is sufficient,
        // because we already set the boundary cells in the prior loop
        if ( domain->_offset > 0. )
        {
          // mark it as outside if it is outside the convex hull, regardless of offset, or outside the domain boundary+offset
          is_outside = is_cell_outside( i );
        }
      }
      
      // invalidate exterior cells
      if ( is_outside )
      {
        _invalid_cells[icell] = true;
        _num_valid--;
      }
      
      // mark ambigous cells as boundary
      if ( !is_interior && !is_outside )
      {
        _boundary_cells[ icell ] = true;
      }
      
      // do nothing for inside cells
      
      if (DEBUG_BOUNDARY)
      {
        if ( is_outside )
          std::cout << "EXTERIOR ";
        else if (is_interior)
          std::cout << "interior ";
        else
          std::cout << "ambigous marked as boundary ";
        // else silent
        std::cout << icell << " ";
        print_cell( i );
        std::cout << std::endl;
      }
      
    }
  }
}


void Grid::allocate_working_memory()
{
  // indices of two cells, i is usually a bottom-level cell, and j its containing top level cell
  _i = new size_t[2];
  _j = new size_t[2];
  
  // array of random numbers, a random location
  _u = new double[2];
  _x = new double[2];
  _xc = new double[2];
}

void Grid::generate_top_grid()
{
  // have ghost layers out to distance equal to _dm from the bounding box, beyond 0 and 1
  _num_ghost_layers = ceil(sqrt(double(2))) + 1;
  _num_cells_in_row = _num_ghost_layers * 2 + 1;
  
  // _s = top level grid cell side length
  _s = _dm / sqrt(2.0);
  _ss = _s;
  
  // assign _n_pow
  // n = 1/_s = number of top level grid cells in a row of the domain along each dimension, not including ghosts
  // _n_pow[k] = n^(k+1). That is, _n_pow[1] = n^2
  _n_pow = new size_t[2];
  _n_pow[0] = size_t(floor(1.0 / _s));
  _s = 1.0 / _n_pow[0];
  _xo = - double(_num_ghost_layers) * _s;
  size_t n_non_ghosts = _n_pow[0];
  _n_pow[0] += (_num_ghost_layers * 2);
  _n_pow[1] = _n_pow[0] * _n_pow[0];
  
  // adjusting distance function to match grid size
  // that is, _dm is shrunk to be exactly equal to a top level grid square diagonal
  _dm = _s * sqrt(2.0);
  _dm_squared = _dm * _dm;
  
  // _n = number of top level grid cells
  _n = _n_pow[1];
  std::cout<< "*** Number of top level grid cells = " << floor( pow( n_non_ghosts, 2) + 0.5) << ", including ghosts = " << _n << std::endl;
  std::cout<< "*** Cells per dimension = " << n_non_ghosts << ", including ghosts = " << _n_pow[0] << std::endl;
  std::cout<< "*** Cell side length = " << _s << std::endl;
  std::cout<< "*** Modified distance function, cell diagonal = " << _dm << std::endl;
  std::cout<< "======================" << std::endl << std::endl;
  
  
  // initialize all top level grid cells to be (boolean) empty, (boolean) active, and (double*) contain no sample point coordinates
  // if a cell is outside the domain, e.g. a ghost cell, then it is marked inactive
  new_cell_points(_n);
  _empty_cells = new bool[_n];
  _boundary_cells = new bool[_n];
  _num_valid = 0;
  _invalid_cells = new bool[_n];
  for (size_t icell = 0; icell < _n; icell++)
  {
    _empty_cells[icell] = true;
    _boundary_cells[icell] = false;
    
    // _invalid_cells
    get_cell_indices(icell, _i);
    if (valid_cell(_i))
    {
      _invalid_cells[icell] = false;
      _num_valid++;
    }
    else
    {
      _invalid_cells[icell] = true;
    }
  }
}


void Grid::initialize_children()
{
  // _active_children = an array of index-array of child cells
  _num_max_children = 4;
  _active_children = new size_t*[_num_max_children];
  for (size_t ichild = 0; ichild < _num_max_children; ichild++)
  {
    _active_children[ichild] = new size_t[2];
  }
}

void Grid::initialize_neighbors()
{
  // jcell_o is the upper-right ghost cell, the max cell?
  size_t jcell, jcell_o;
  _j[0] = _j[1] = _num_ghost_layers;
  get_cell_index(_j, jcell_o);
  _j[0]= _j[1] = 0;
  
  // samitch: I don't know what this block does
  // it may have something to do with the next block, to
  // compute the template of neighbor indices
  _num_neighbors = 0;
  size_t k_dim(2);
  while (true)
  {
    while (_j[k_dim] < _num_cells_in_row)
    {
      _j[k_dim]++; _num_neighbors++;
    }
    size_t kk_dim(k_dim - 1);
    
    bool done(false);
    while (true)
    {
      _j[kk_dim]++;
      if (_j[kk_dim] == _num_cells_in_row)
      {
        _j[kk_dim] = 0;
        if (kk_dim == 0)
        {
          done = true;
          break;
        }
        kk_dim--;
      }
      else break;
    }
    if (done) break;
    _j[k_dim] = 0;
  }
  
  // collect neighbors
  // that is, build the template of neighbor indices that we check for disks covering a sample point.
  _j[0] = _j[1] = 0;
  
  _neighbors = new int[_num_neighbors];
  double* neighbors_dmin = new double [_num_neighbors];
  double* neighbors_dmax = new double [_num_neighbors];
  _num_neighbors = 0; double mid_point((_num_ghost_layers + 0.5) * _s), min_point(_num_ghost_layers * _s), max_point((_num_ghost_layers + 1.0) * _s);
  while (true)
  {
    while (_j[k_dim] < _num_cells_in_row)
    {
      get_cell_index(_j, jcell);
      if (jcell != jcell_o)
      {
        _neighbors[_num_neighbors] = -int(jcell_o) + int(jcell);
        neighbors_dmin[_num_neighbors] = 0.0;
        for (size_t dim = 0; dim < 2; dim++)
        {
          double corner = _j[dim] * _s;
          if (fabs(corner - mid_point) > fabs(corner + _s - mid_point)) corner += _s;
          double dst = (corner - min_point);
          if (fabs(corner - max_point) < fabs(dst)) dst = corner - max_point;
          neighbors_dmin[_num_neighbors] += dst * dst;
        }
        neighbors_dmax[_num_neighbors] = 0.0;
        for (size_t dim = 0; dim < 2; dim++)
        {
          double corner = _j[dim] * _s;
          if (fabs(corner - mid_point) < fabs(corner + _s - mid_point)) corner += _s;
          double dst = (corner - min_point);
          if (fabs(corner - max_point) < fabs(dst)) dst = corner - max_point;
          neighbors_dmax[_num_neighbors] += dst * dst;
        }
        _num_neighbors++;
      }
      
      _j[k_dim]++;
    }
    size_t kk_dim(k_dim - 1);
    
    bool done(false);
    while (true)
    {
      _j[kk_dim]++;
      if (_j[kk_dim] == _num_cells_in_row)
      {
        _j[kk_dim] = 0;
        if (kk_dim == 0)
        {
          done = true;
          break;
        }
        kk_dim--;
      }
      else break;
    }
    if (done) break;
    _j[k_dim] = 0;
  }
  
  // arrange dimension based on their fathest corner from the center of jcell_o
  
  for (size_t ii = 0; ii < _num_neighbors; ii++)
  {
    for (size_t jj = ii + 1; jj < _num_neighbors; jj++)
    {
      if (neighbors_dmax[jj] < neighbors_dmax[ii])
      {
        int tmp = _neighbors[ii]; _neighbors[ii] = _neighbors[jj]; _neighbors[jj] = tmp;
        double tmp_d = neighbors_dmin[ii]; neighbors_dmin[ii] = neighbors_dmin[jj]; neighbors_dmin[jj] = tmp_d;
        tmp_d = neighbors_dmax[ii]; neighbors_dmax[ii] = neighbors_dmax[jj]; neighbors_dmax[jj] = tmp_d;
      }
    }
  }
  // remove neighbors that are too far away from the center
  for (int ii = int(_num_neighbors) - 1; ii >= 0; ii--)
  {
    if (neighbors_dmin[ii] > _dm_squared - 1.0E-10)
    {
      _num_neighbors--;
    }
  }
}


void Grid::setup(Domain *_domain, double &sampling_radius)
{
  _start_time = clock();
  
  const double desired_radius = sampling_radius * _domain->_scale;
  
  initialize(desired_radius);
  
  // update cells with the polygon
  initialize_cells_by_domain( _domain );
  
  initialize_top_pool();
  
  sampling_radius = _dm / _domain->_scale;
}

void Grid::MPS()
{
  // Simple_MPS,
  // Throw darts, refine cells, etc.
  srand(1);
  size_t ref_levels(30);
  for (size_t iref = 0; iref <= ref_levels; iref++)
  {
    if ( use_grid_level(iref) )
      break;
    subdivide_cells();
  }
  
  _end_time = clock();
  const double cpu_time = ((double) (_end_time - _start_time)) / CLOCKS_PER_SEC;
  
  report_stats(cpu_time);
  
}

void Grid::add_samples(Sampling *sampling)
{
  assert(0); // not implemented properly
  // force them all in
  size_t i[2], icell;
  SamplingIterator si ( sampling );
  double *p = si.first();
  while (p)
  {
    get_cell_indicies( p, i );
    get_cell_index( i, icell );
    _cell_points[icell] = p; // convert to some other assignment
    
    p = si.next();
  }
}

void Grid::add_samples_conditionally(Sampling *sampling)
{
  // conditionally add points, checking conflict for each
  // todo: go in random order
  // force them all in
  SamplingIterator si ( sampling );
  double *p = si.first();
  while (p)
  {
    // old grid coordinates are the same as the new grid coordinates, same domain and frame
    add_sample_conditionally(p);
    
    p = si.next();
  }
}

void Grid::initialize(double sampling_radius)
{
  
  _dm = sampling_radius;
		
  allocate_working_memory();
  
  // generate background grid
  generate_top_grid();
  
  initialize_children();
  
  initialize_neighbors();
}

void write_matrix( Matrix &matrix, size_t column_size, std::string fname )
{
  std::cout<< "writing sparse matrix, with 0-1 values " << fname << std::endl;
  std::cout<< "indices starting at 1, rows and columns in first line" << std::endl;
  std::fstream file(fname.c_str(), std::ios::out);
  file << matrix.size() << " " << column_size << std::endl;  // number of lines
  for (size_t j = 0; j < matrix.size(); ++j )
  {
    file << j+1 << "      ";
    for (size_t i = 0; i < matrix[j].size(); ++i )
    {
      file << " " << matrix[j][i] + 1;
      assert( matrix[j][i] < column_size );
    }
    file << std::endl;
  }
}

void read_matrix( Matrix &matrix, size_t &column_size, std::string fname )
{
  std::cout<< "reading sparse matrix, with 0-1 entries " << fname << std::endl;
  std::cout<< "indices starting at 1, rows and columns in first line" << std::endl;
  std::ifstream file(fname.c_str(), std::ios::in);
  size_t row_size(0);
  file >> row_size >> column_size;
  matrix.resize(row_size);
  for (size_t j = 0; j < matrix.size(); ++j )
  {
    size_t row;
    file >> row;
    assert( row == j+1 );
    std::string line;
    getline( file, line );
    std::istringstream fin(line);
    size_t col;
    while (fin >> col)
    {
      assert( col > 0 && col <= column_size);
      matrix[j].push_back(col-1);
    }
  }
}

void write_matrix_entries( MatrixEntries &matrix, size_t column_size, size_t max_entry, std::string fname )
{
  std::cout<< "writing sparse matrix, with floating point values " << fname << std::endl;
  std::cout<< "indices starting at 1, rows and columns in first line" << std::endl;
  std::fstream file(fname.c_str(), std::ios::out);
  file << matrix.size() << " " << column_size << std::endl;  // number of lines
  const double max_entry_d = (double) max_entry;
  for (size_t j = 0; j < matrix.size(); ++j )
  {
    file << j+1 << "        ";
    for (size_t i = 0; i < matrix[j].size(); ++i )
    {
      file << "   " << matrix[j][i].first + 1 << " " << matrix[j][i].second / max_entry_d ;
      assert( matrix[j][i].first < column_size );
    }
    file << std::endl;
  }

}


// inverse of save_point_cloud
void read_point_coordinates( PointVector &points, double &r, std::string fname )
{
  std::ifstream file( fname.c_str(), std::ios::in );
  
  size_t num_points(0);
  file >> num_points;
  points.resize(num_points);
  
  file >> r;
  
  for ( size_t i = 0; i < num_points; ++i )
  {
    file >> points[i][0] >> points[i][1];
  }
}

void find_footprint_coverage( std::string footprint_fname )
{
  Footprint *footprint = Footprint::read_footprint( footprint_fname );
  assert(footprint);
  // footprint->set_scale( 1. );
  
  // in the integrated demo we used the grid coordinates, but in the stand-alone we will use absolute coordinates,
  // that is, epsilon instead of scaled_epsilon,
  PointVector coverage_points, placement_points;
  double epsilon, delta;
  read_point_coordinates( coverage_points, epsilon, "points_for_coverage.txt");
  read_point_coordinates( placement_points, delta, "points_for_placement.txt");
  
  Matrix placement_covers; // for placement i, the list of j coverage points covered
  Matrix covered_by_placements; // for coverage point i, the list of placemetns that cover it
  placement_covers.resize(placement_points.size());
  covered_by_placements.resize(coverage_points.size());
  
  // define the offset, how much we have to shrink the footprint to ensure the original one covers the whole domain
  footprint->_offset = epsilon;
  
  // calculate coverage matrices
  for (size_t j = 0; j < coverage_points.size(); ++j)
  {
    double *p = coverage_points[j].data();
    
    // todo optimization:
    // use the grid to efficiently iterate over just the coverage points within max_radius;
    // if we have less than 10,000 points in each, this is still instantaneous
    // virtual double Footprint::max_radius()
    
    for (size_t i = 0; i < placement_points.size(); ++i)
    {
      double *c = placement_points[i].data();
      
      if ( footprint->is_covered( p, c ) )
      {
        placement_covers[i].push_back(j);
        covered_by_placements[j].push_back(i);
      }
    }
  }
  
  // write coordinates of points we need to cover, 1 to n
  // write coordinates of placements, 1 to m
  
  // for each placement, write which samples it covers
  write_matrix( placement_covers, covered_by_placements.size(), "footprint_covers.txt");
  write_matrix( covered_by_placements, placement_covers.size(), "covered_by_footprints.txt");
  
  // for each placement, determine the amount of overlap with all other placements
  size_t max_overlap = 0;
  // each row is a placement index, each entry is another placement index, and the amount of overlap
  MatrixEntries overlapping_pairs(placement_covers.size());
  for (size_t i = 0; i < placement_covers.size(); ++i)
  {
    // placement_covers[i] overlaps
    for (size_t j = 0; j < placement_covers.size(); ++j)
    {
      // skip writing self-overlaps
      if ( i != j )
      {
        // placement_covers[j] overlaps
        size_t overlap = num_common_entries( placement_covers[i], placement_covers[j] );
        if (overlap)
        {
          overlapping_pairs[i].push_back( std::make_pair(j, overlap ) );
          if (overlap>max_overlap)
            max_overlap = overlap;
        }
      }
    }
  }
  write_matrix_entries( overlapping_pairs, placement_covers.size(), max_overlap, "footprint_overlaps.txt" );
  
  // plot coverage, for each placement, in a separate file
  {
    std::vector< Footprint *> footprints(1);
    footprints[0] = footprint;
    std::vector<size_t> toes(1);
    for (size_t i = 0; i < placement_points.size(); ++i)
    {
      toes[0] = i;
      
      // make is_covered vector, set flags for this particular placement
      std::vector<char> is_covered( coverage_points.size() );
      for ( size_t j = 0; j < placement_covers[i].size(); ++j )
      {
        const size_t covered_point = placement_covers[i][j];
        is_covered[covered_point] = '\1'; // true;
      }
      
      std::ostringstream oss;
      // oss << footprint_fname << "_footprint_covers_" << i + 1 << ".ps";
      oss << "footprint_covers_" << i + 1 << ".ps";
      plot_coverage( footprints, toes, placement_points, coverage_points, is_covered, oss.str() );
    }
  }
  
  delete footprint;
}
size_t num_common_entries( std::vector<size_t> & rowA, std::vector<size_t> &rowB )
{
  size_t a(0), b(0), c(0);
  while ( a < rowA.size() && b < rowB.size() )
  {
    if ( rowA[a] == rowB[b] )
    {
      ++a;
      ++b;
      ++c;
    }
    else if ( rowA[a] < rowB[b] )
    {
      ++a;
    }
    else
    {
      assert( rowA[a] > rowB[b] );
      ++b;
    }
  }
  return c;
}



void read_placements( std::vector<Footprint*> &footprints,  std::vector<size_t> &toes, std::string solution_fname )
{
  std::ifstream file( solution_fname.c_str(), std::ios::in );
  
  size_t num_toes(0);
  file >> num_toes;
  footprints.resize(num_toes);
  toes.resize(num_toes);
  
  for ( size_t i = 0; i < num_toes; ++i )
  {
    std::string footprint_filename;
    size_t toeid;
    file >> toeid >> footprint_filename;
    assert( toeid > 0 );
    // assert( toeid < num_placements );
    toes[i] = toeid - 1;
    footprints[i] = Footprint::read_footprint(footprint_filename);
  }
}

void find_solution_coverage( std::string solution_fname )
{
  std::vector<Footprint*> footprints;
  std::vector<size_t> toes;
  read_placements( footprints, toes, solution_fname );
  
  PointVector coverage_points, placement_points;
  double epsilon, delta;
  read_point_coordinates( coverage_points, epsilon, "points_for_coverage.txt");
  read_point_coordinates( placement_points, delta, "points_for_placement.txt");
  
  // define the offset, how much we have to shrink the footprint to ensure the original one covers the whole domain
  for (size_t i = 0; i < footprints.size(); ++i)
    footprints[i]->_offset = epsilon;
  
  Matrix placement_covers; // for placement i, the list of j coverage points covered
  Matrix covered_by_placements; // for coverage point i, the list of placemetns that cover it
  size_t coverage_size(0), placement_size(0);
  read_matrix( placement_covers, coverage_size, "footprint_covers.txt");
  read_matrix( covered_by_placements, placement_size, "covered_by_footprints.txt");
  assert( coverage_size == covered_by_placements.size() );
  assert( placement_size == placement_covers.size() );
  
  plot_solution( footprints, toes, placement_covers, coverage_points, placement_points );
  
  // cleanup
  for (size_t i = 0; i < footprints.size(); ++i)
    delete footprints[i];
}

void plot_solution( std::vector<Footprint*> &footprints, std::vector<size_t> &toes, Matrix &placement_covers, PointVector &coverage_points, PointVector &placement_points )
{
  assert( footprints.size() == toes.size() );
  
  // plot coverage, for chosen placements
  std::vector<char> is_covered( coverage_points.size() );
  {
    for (size_t j = 0; j < toes.size(); ++j )
    {
      const size_t k = toes[j]; // index of placement
      for (size_t i = 0; i < placement_covers[k].size(); ++i)
      {
        is_covered[ placement_covers[k][i] ] = '\1'; // true
      }
    }
  }
  
  std::ostringstream oss;
  // oss << footprint_fname << "_footprint_covers_" << i + 1 << ".ps";
  oss << "footprint_covers" << ".ps";
  plot_coverage( footprints, toes, placement_points, coverage_points, is_covered, oss.str() );
}

void Grid::initialize_top_pool()
{
  _pool_size = 32 * _num_valid;
  
  _active_pool_i = new size_t[_pool_size * 2];
  
  // fill the active pool with the indices of the (valid in-domain) top level grid cells
  size_t tmp(0);
  size_t no(0);
  for (size_t icell = 0; icell < _n; icell++)
  {
    if (!_invalid_cells[icell])
    {
      get_cell_indices(icell, _i);
      for (size_t idim = 0; idim < 2; idim++)
      {
        _active_pool_i[tmp] = _i[idim];
        tmp++;
      }
      no++;
    }
  }
  
  // initialize the rest of the active poool with zeros
  for (size_t ii = no; ii < _pool_size; ii++)
  {
    for (size_t idim = 0; idim < 2; idim++)
    {
      _active_pool_i[tmp] = 0;
      tmp++;
    }
  }
  
  _no = no;
  _TWO_POW_REF_LEVEL = 1;
  _ss = _s;
  
  std::ostringstream oss;
  oss << _fname_prefix << "_0domain.ps";
  plot(false, false, oss.str() );
}

void Grid::subdivide_cells()
{
  // size_t new_two_pow_ref_lef  = _TWO_POW_REF_LEVEL * 2;
  
  size_t num_active_children;
  size_t ii(_pool_size - 1), jj((ii + 1) *2 - 1), kk(_num_valid - 1);
  // size_t prevLevel(refLevel - 1);
  // kk is the bottom of the queue of active this-level children
  // jj is the start of the d-indices of the subdivided next-level children,
  //   placed at the bottom of the array upwards.
  while (true)
  {
    // icell is the top level cell
    // _i is the d-indices of an active child cell, taken from the bottom of the queue
    // _j is the d-indices of its top level cell
    // icell is the linear index of its top level cell
    size_t kkst = kk * 2;
    for (size_t idim = 0; idim <2; idim++)
    {
      _i[idim] = _active_pool_i[kkst + idim];
      _j[idim] = _i[idim] / _TWO_POW_REF_LEVEL; // indices of icell
    }
    size_t icell;
    get_cell_index(_j, icell);
    
    // if top level cell is not covered
    if (!_invalid_cells[icell])
    {
      // subdivide current level cell into children
      get_active_children(_i, icell, _active_children, num_active_children);
      assert( num_active_children <= _max_children );
      
      // put subdivided children into the pool,
      for (size_t ichild = 0; ichild < num_active_children; ichild++)
      {
        for (size_t idim = 0; idim < 2; idim++)
        {
          _active_pool_i[jj] = _active_children[ichild][2 - idim - 1];
          jj--;
        }
        if (ii < kk)
        {
          // what does this error indicate?
          std::cout << "out of memory when creating children" << std::endl;
          std::cout << "ERROR in filling active pool" << std::endl;
        }
        assert (ii >= kk); // else we will start overwriting the pool with bad indices
        ii--;
      }
    }
    if (kk == 0) break;
    
    kk--;
  }
  ii++;
  
  // move child cells to the front of the array?
  // alternatively, one could reverse the index directions and avoid all this data movement
  _num_valid = _pool_size - ii;
  // no = _num_valid;
  for (size_t ivalid = 0; ivalid < _num_valid; ivalid++)
  {
    for (size_t idim = 0; idim < 2; idim++) _active_pool_i[ivalid * 2 + idim] = _active_pool_i[(ivalid + ii) * 2 + idim];
  }
  _TWO_POW_REF_LEVEL *= 2;
  _ss *= 0.5;
}

bool Grid::use_grid_level(size_t refLevel)
{
#pragma region Grid Level refLevel:
  
  clock_t start_time, end_time; double cpu_time;
  start_time = clock();
  
  size_t icell; // jcell;
  size_t num_inserted(_num_inserted_points);
  
  // quit if no children are uncovered (valid)
  if (_num_valid == 0)
  {
    end_time = clock();
    cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    
    std::cout<< "Refinement Level : (" << refLevel << ")" << std::endl;
    std::cout<< "spacing : (" << _ss << ")" << std::endl;
    std::cout<< "Pool is empty .. maximal condition is achieved!" << std::endl;
    std::cout<< "Execution Time = " << cpu_time << " seconds." << std::endl;
    std::cout<< std::endl;
    return true;
  }
  
  // plot the current cells
  {
    std::ostringstream oss;
    oss << _fname_prefix << "_" << (refLevel + 1) << ".ps";
    plot( true, false, oss.str() );
  }
  
  double urand;
  size_t irand;
  size_t num_invalidated(0), num_invalidated_invalid(0);
  
  // throw darts, nit times
  const size_t no(_num_valid);
  
  size_t nit(0);
  size_t nit_max(size_t(_RF * no));
  if (refLevel == 0)
    nit_max *= 2;
  while (nit < nit_max)
  {
    _num_darts++;
    
    // pick a random cell i from active pool.
    // its top level grid cell is j
    urand = THE_RNG.generate_a_random_number();
    irand = size_t(urand * (_num_valid - 1));
    
    // kst is the index into the active pool of the cell we picked to sample from
    size_t kst = irand * 2;
    for (size_t idim = 0; idim < 2; idim++)
    {
      _i[idim] = _active_pool_i[kst + idim];
      _j[idim] =  _i[idim] / _TWO_POW_REF_LEVEL;
    }
    get_cell_index(_j, icell);
    
    // Check if top level cell is invalid.
    // This happens when some other cell at this level receives a sample, which covers the top level cell
    if (_invalid_cells[icell])
    {
      // swap the invalid cell at this level with the end of the active list
      _num_valid--;
      size_t tmp(0);
      size_t kkst = _num_valid * 2;
      for (size_t idim = 0; idim < 2; idim++)
      {
        size_t i1 = kst + idim;
        size_t i2 = kkst + idim;
        tmp =  _active_pool_i[i1];
        _active_pool_i[i1] = _active_pool_i[i2];
        _active_pool_i[i2] = tmp;
      }
      num_invalidated_invalid++;
      nit++;
      if (_num_valid == 0) break;
      continue;
    }
    
    // pick a random point
    for (size_t idim = 0; idim < 2; idim++)
    {
      _u[idim] = THE_RNG.generate_a_random_number();
      _x[idim] = _xo +  (_i[idim] + _u[idim]) * _ss;
    }
    
    if ( add_sample_conditionally( _x, icell ) )
    {
      // deactivate that cell
      num_invalidated++;
      _num_valid--;
      size_t tmp(0);
      size_t kkst = _num_valid * 2;
      for (size_t idim = 0; idim < 2; idim++)
      {
        size_t i1 = kst + idim;
        size_t i2 = kkst + idim;
        tmp =  _active_pool_i[i1];
        _active_pool_i[i1] = _active_pool_i[i2];
        _active_pool_i[i2] = tmp;
      }
    }
    
    if (_num_valid == 0) break;
    
    nit++;
  }
  
  end_time = clock();
  cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
  
  std::cout<< "Refinement Level : (" << refLevel << ")" << std::endl;
  std::cout<< "Active cells initially " << no << std::endl
  << "    x/max_cells_memory = " << no * 1.0 / _pool_size << std::endl;
  //  std::cout<< "Maximum number of cells that can fit in the allocated memory " << _pool_size << std::endl;
  std::cout<< "Number of inserted points this level = " << (_num_inserted_points - num_inserted) << std::endl;
  std::cout<< std::endl;
  std::cout<< "Iterations, number of times pick a square, pick a point = " << nit << std::endl
  << "    x/initial = " << nit * 1.0 / no << std::endl;
  std::cout<< "Remaining cells = " << _num_valid << std::endl << "    x/initial = " << _num_valid * 1.0 / no << std::endl;
  std::cout<< "Invalidated cells due to point insertion in this cell " <<  num_invalidated
  << std::endl << "    x/initial = " << num_invalidated * 1.0 / no << std::endl;
  std::cout<< "Invalidated cells due to a nearby point insertion = " << num_invalidated_invalid << std::endl
  << "    x/initial = " << num_invalidated_invalid * 1.0 / no << std::endl;
  std::cout<< "Insertion Rate = " << (_num_inserted_points - num_inserted) * 1.0 / cpu_time << " points/sec." << std::endl;
  std::cout<< "Execution Time this level = " << cpu_time << " seconds." << std::endl;
  std::cout<< std::endl;
  std::cout<< "--------------------------------------------------------------------" << std::endl;
  std::cout<< std::endl;
  
  // final output plot, with disks
  {
    std::ostringstream oss;
    oss << _fname_prefix << ".ps";
    plot( true, true, oss.str() );
  }
  
  // done if no more valid cells
  return (_num_valid == 0);
  
#pragma endregion
};

bool Grid::add_sample_conditionally( double * x, size_t icell )
{
  // invalid cell? e.g. already a point in icell?
  if ( _invalid_cells[icell] )
    return false;
  
  // check for conflicts with prior disks
  bool valid_point(true);
  for (size_t index = 0; index < _num_neighbors; index++)
  {
    size_t jcell = icell + _neighbors[index];
    if (_empty_cells[jcell]) continue;
    
    double dd(0.0);
    for (size_t idim = 0; idim < 2; idim++)
    {
      double dx = x[idim] - _cell_points[jcell][idim];
      dd += dx * dx;
    }
    
    if (dd < _dm_squared)
    {
      valid_point = false;
      break;
    }
  }
  
  // check if point is inside the domain
  // taking the offset into account
  
  // valid_point = valid_point && (!_boundary_cells[icell] || _domain->is_point_inside(x) "up to offset" );
  if (valid_point)
  {
    if ( _boundary_cells[icell] ) // including cells closer to the boundary than offset
    {
      if (_domain->_offset == 0.)
        valid_point = _domain->is_point_inside(x);
      else
      {
        valid_point = (_domain->signed_distance_to_boundary( x ) < _domain->_offset) &&
                       _domain->_convex_hull.is_point_geometrically_inside( x );
      }
    }
  }
  
  // sample is good
  // save it in the top level cell
  // remove its current-level cell from the active list
  if (valid_point)
  {
    // Create a new points and assign it to that cell
    _invalid_cells[icell] = true;
    _empty_cells[icell] = false;
    
    _cell_points[icell] = new double[2];
    for (size_t idim = 0; idim < 2; idim++) _cell_points[icell][idim] = x[idim];
    _num_inserted_points++;
    
    return true;
  }
  return false;
}

inline void Grid::get_active_children(size_t* i, size_t icell, size_t** active_children, size_t &num_active_children)
{
#pragma region Retrieve Active Children:
  // a child cell is "active" if it is not covered by a disk, and not outside the domain
  
  // i is the cell we are subdividing into children
  // icell is the top level grid cell containing subcell i
  
  num_active_children = 0;
  
  // _x = coordinates of cell i's lower left (smallest) corner
  for (size_t idim = 0; idim < 2; idim++) _x[idim] = _xo + i[idim] * _ss;
  
  // ss = size of children
  double ss = 0.5 * _ss;
  
  // _j = (0,1), indices of a child cell of the parent
  for (size_t idim = 0; idim < 2; idim++) _j[idim] = 0;
  
  // a nested loop over children
  // at each loop we add one to _j, taken as a binary number: 0000, 0001, 0010, 0011, 0100, etc
  size_t k_dim(1);
  while (true)
  {
    while (_j[k_dim] < 2)
    {
      // coordinates of lower left corner of a child
      for (size_t idim = 0; idim < 2; idim++)
      {
        _xc[idim] = _x[idim] + _j[idim] * ss;
      }
      
      // check if this child is covered by a disk
      bool covered(false);
      // check nearby top level cells
      for (size_t index = 0; index < _num_neighbors;  index++)
      {
        size_t jcell = icell + _neighbors[index];
        if (_empty_cells[jcell]) continue;
        
        double dd(0.0);
        for (size_t idim = 0; idim < 2; idim++)
        {
          double xx = _cell_points[jcell][idim];
          
          double xf(_xc[idim]), xfp(xf + ss);
          if (fabs(xf - xx) < fabs(xfp - xx)) xf = xfp;
          double dx(xf - xx);
          dd += dx * dx;
        }
        if (dd < _dm_squared) {covered = true; break;}
      }
      
      // check if child is outside the domain
      // here _boundary_cells includes those ambigously close to an offset boundary
      if (!covered && _boundary_cells[icell])
      {
        bool is_exterior = _domain->is_cell_outside(_xc, ss);
        if (is_exterior)
          covered = true;
      }
      
      if (!covered)
      {
        for (size_t idim = 0; idim < 2; idim++)
        {
          active_children[num_active_children][idim] = 2 * i[idim] + _j[idim];
        }
        num_active_children++;
      }
      
      _j[k_dim]++; // move to the next child
    }
    size_t kk_dim(k_dim - 1);
    
    bool done(false);
    while (true)
    {
      _j[kk_dim]++;
      if (_j[kk_dim] == 2)
      {
        _j[kk_dim] = 0;
        if (kk_dim == 0)
        {
          done = true;
          break;
        }
        kk_dim--;
      }
      else break;
    }
    if (done) break;
    _j[k_dim] = 0;
  }
#pragma endregion
};


inline bool Grid::valid_cell(size_t* i)
{
#pragma region Check if cell is valid with regard to boundaries:
  for (size_t dim = 0; dim < 2; dim++)
  {
    double x = _xo + i[dim] * _s;
    if (x > 1.0 - 1E-10 || x + _s < 1E-10) return false;
  }
  return true;
#pragma endregion
};

void Grid::save_point_cloud( std::string file_name, bool debug )
{
#pragma region Save Point Cloud:
  if (debug) std::cout << "\ndebug file format ";
  std::cout<< "\n*** Saving point cloud in " << file_name << " ... " << std::endl;
  
  std::fstream file(file_name.c_str(), std::ios::out);
  
  // # Points
  file << _num_inserted_points;           // number of points
  if (debug) file << " number of points ";
  file << std::endl;
  
  // # sampling density r
  if (debug)
  {
    file << "sampling density r = " << _dm << " grid-abs    ";
  }
  file << _dm / _domain->_scale << std::endl;
  
  file.precision(15); file << std::fixed;
  double a[2];
  for (size_t icell = 0; icell < _n; icell++)
  {
    if (_empty_cells[icell]) continue;
    
    bool ghost_cell(false);
    get_cell_indices(icell, _i); // index of ghost cell
    for (size_t idim = 0; idim < 2; idim++)
    {
      if (_i[idim] < _num_ghost_layers)
      {ghost_cell = true; break;}
      
      if (_i[idim] + _num_ghost_layers >= _n_pow[0])
      {ghost_cell = true; break;}
    }
    if (ghost_cell) continue;
    
    // grid coordinates in (0,1)
    if (debug)
    {
      file << _cell_points[icell][0];
      for (size_t idim = 1; idim < 2; idim++)
        file << " " << _cell_points[icell][idim];
      file << " grid-abs   ";
    }
    
    // absolute coordinates
    _domain->grid_coordinates_to_absolute( _cell_points[icell], a );
    file << a[0] << " " << a[1] << std::endl;
    
  }
#pragma endregion
}

void Grid::test_output()
{
#pragma region Test Output:
  double dmin = 100 * _dm_squared;
  
  for (size_t icell = 0; icell < _n; icell++)
  {
    if (_empty_cells[icell]) continue;
    for (size_t jcell = icell + 1; jcell < _n; jcell++)
    {
      if (_empty_cells[jcell]) continue;
      
      double dstsq(0.0);
      for (size_t idim = 0; idim < 2; idim++)
      {
        double dx = _cell_points[icell][idim] - _cell_points[jcell][idim];
        dstsq += dx * dx;
      }
      if (dstsq < dmin) dmin = dstsq;
    }
  }
  
  std::cout<<"    ^ Minimum distance between any two points is " << sqrt(dmin) << std::endl;
  
#pragma endregion
};

void Grid::print_cell(size_t *i, std::ostream &out)
{
  out << "cell (" << i[0];
  for (size_t d = 1; d < 2; ++d)
    out << "," << i[d];
  out << ")";
  out << " coordinates ";
  double p[2];
  get_cell_lower(i, p);
  print_point(p, out);
  out << " to ";
  get_cell_upper(i, p);
  print_point(p, out);
}

void Domain::print_vertex(size_t v, std::ostream &out)
{
  out << "vertex " << v << " ";
  print_point( _vertices[v], out );
}

void print_point(double *p, std::ostream &out)
{
  out << "(" << p[0];
  for (size_t d = 1; d < 2; ++d)
    out << "," << p[d];
  out << ")";
}

void Domain::plot_ps( std::fstream &file, double scale)
{
  for (size_t p = 0; p < _loops.size(); ++p)
  {
    Polygon *poly = _loops[p];
    for (size_t e = 0; e < poly->_edges.size(); ++e)
    {
      double *v0 = poly->v0(e);
      double *v1 = poly->v1(e);
      plot_seg( file, scale, v0[0], v0[1], v1[0], v1[1] );
    }
  }
}

void plot_coverage( std::vector< Footprint *> &footprints, std::vector<size_t> &toes, PointVector &placement_points, PointVector &coverage_points, std::vector<char> &is_covered, std::string fname )
{
  // preamble
  std::fstream file;
  open_ps_file( file, fname );
  double scale = define_plot_box( file, 0, 0, 1, 1 );
  define_shapes( file );
  const double s = 0.01;
  
  
  // plot footprints - solid colors?
  for (size_t j=0; j < footprints.size(); ++j)
  {
    std::array<double,2> p = placement_points[ toes[j] ];
    footprints[j]->plot(file, scale, p[0], p[1]);
  }
  
  // plot placement points +
  for (size_t i=0; i < placement_points.size(); ++i)
  {
    plot_plus( file, scale, placement_points[i][0], placement_points[i][1], s );
  }
  
  // plot toe placement points
  set_color(file, dark_blue);
  for (size_t j=0; j < toes.size(); ++j)
  {
    size_t i = toes[j];
    //  plot_plus( file, scale, placement_points[i][0], placement_points[i][1], s );
    plot_circle( file, scale, placement_points[i][0], placement_points[i][1], s, false );
  }
  set_color(file);
  
  
  // plot coverage points x
  assert( is_covered.empty() || is_covered.size() == coverage_points.size() );
  bool last_covered = false;
  set_color(file);
  for (size_t i=0; i < coverage_points.size(); ++i)
  {
    // switch colors, if need be
    bool covered = is_covered.empty() ? false : is_covered[i] == '\1';
    if ( covered && ! last_covered)
      set_color(file, red);
    else if (!covered && last_covered)
      set_color(file);
    
    plot_ex( file, scale, coverage_points[i][0], coverage_points[i][1], s );
    
    last_covered = covered;
  }
  set_color(file);

  // plotting polygon 
  // currently doesn't work correctly, the scaling is off

  // Domain domain;
  // domain.read( "polygon.txt" );
  // domain.plot_ps(file, scale);

}

void Grid::plot_grid_cells( std::fstream &file, double scale )
{
  size_t* i = new size_t [2];
  for (size_t icell = 0; icell < _n; icell++)
  {
    get_cell_indices(icell, i);
    double xo = _xo + i[0] * _s;
    double yo = _xo + i[1] * _s;
    double xn = xo + _s;
    double yn = yo + _s;
    
    // plot cell
    if (xo > -1E-10 && yo > -1E-10 && xn < 1.0 + 1E-10 && yn < 1.0 + 1E-10)
    {
      // boundary = medium
      // exterior = dark
      // interior = light
      std::string quad_weight = _boundary_cells[icell] ? "medium" : ((_invalid_cells[icell]) ? "dark" : "light" );
      
      plot_quad( file, scale, xo, yo, xn, yn, quad_weight );
      
    }
  }
}

void Grid::plot_active_cells( std::fstream &file, double scale)
{
  for (size_t icell = 0; icell < _pool_size; icell++)
  {
    _i[0] = _active_pool_i[icell * 2];
    _i[1] = _active_pool_i[icell * 2 + 1];
    _j[0] = _i[0] / _TWO_POW_REF_LEVEL;
    _j[1] = _i[1] / _TWO_POW_REF_LEVEL;
    
    size_t jcell;
    get_cell_index(_j, jcell);
    
    if (_invalid_cells[jcell]) continue;
    
    if (icell < _num_valid)
    {
      // light grey
      double xo = _xo + _i[0] * _ss;
      double yo = _xo + _i[1] * _ss;
      double xn = xo + _ss;
      double yn = yo + _ss;
      plot_quad( file, scale, xo, yo, xn, yn, "light" );
    }
  }
}

void Grid::plot_samples( std::fstream &file, double scale, bool disks)
{
  // plot disks
  if (disks)
  {
    // plot filled circles
    for (size_t icell = 0; icell < _n; icell++)
    {
      if (!_empty_cells[icell])
        plot_circle( file, scale, _cell_points[icell][0], _cell_points[icell][1], _dm, true );
    }
    
    // plot discs boundaries
    for (size_t icell = 0; icell < _n; icell++)
    {
      if (!_empty_cells[icell])
        plot_circle( file, scale, _cell_points[icell][0], _cell_points[icell][1], _dm, false );
    }
  }
  
  
  // plot center points of sample disks
  double s(_dm * 0.12);
  for (size_t icell = 0; icell < _n; icell++)
  {
    if (!_empty_cells[icell])
      plot_ex( file, scale, _cell_points[icell][0], _cell_points[icell][1], s );
  }
}

void Grid::plot( bool active_cells, bool disks, std::string fname ) // true, "points.ps"
{
  // file
  std::fstream file;
  open_ps_file( file, fname );
  
  // bounding box
  double xmin(-2.0 * double(_num_ghost_layers) * _s);
  double ymin(-2.0 * double(_num_ghost_layers) * _s);
  double Lx(1.0 + 4.0 * double(_num_ghost_layers) * _s);
  double Ly(1.0 + 4.0 * double(_num_ghost_layers) * _s);
  
  double scale = define_plot_box( file, xmin, ymin, Lx, Ly );

  // shapes
  define_shapes( file );
  
  // plot background grid
  if (!active_cells)
  {
    plot_grid_cells(file, scale);
  }
  
  // plot active pool, those that are not invalid
  if (/* DISABLES CODE */ (active_cells))
  {
    plot_active_cells(file, scale);
  }
  
  // samples, and their disks
  plot_samples( file, scale, disks );
  
  // plot Domain polygons
  if (1)
  {
    _domain->plot_ps( file, scale );
  }
  
  plot_bbox( file, scale );
  
  file << "showpage" << std::endl;
  
#pragma endregion
};


void Domain::read( std::string file_name )
{
  std::ifstream file_in( file_name.c_str(), std::ios::in );
  
  // vertices
  size_t vertices_size(0);
  file_in >> vertices_size;
  _vertices.resize(vertices_size);
  double x, y;
  for (size_t v = 0; v < vertices_size; ++v)
  {
    file_in >> x >> y;
    _vertices[v] = new double[2];
    _vertices[v][0] = x;
    _vertices[v][1] = y;
  }
  
  // rescale vertices to unit box
  // could move to separate function
  if (vertices_size)
  {
    // get extent of domain
    double max_x, max_y;
    max_x = _min_x = _vertices[0][0];
    max_y = _min_y = _vertices[0][1];
    for (size_t v = 1; v < vertices_size; ++v)
    {
      const double x = _vertices[v][0];
      const double y = _vertices[v][1];
      if (x > max_x) max_x = x;
      if (x < _min_x) _min_x = x;
      if (y > max_y) max_y = y;
      if (y < _min_y) _min_y = y;
    }
    
    // add the frame around the domain
    _min_x -= _frame_size;
    max_x += _frame_size;
    _min_y -= _frame_size;
    max_y += _frame_size;
    
    // vertices will be placed within epsilon of the domain boundary
    const double length_x = max_x - _min_x;
    const double length_y = max_y - _min_y;
    _length = (length_x > length_y) ? length_x : length_y;
    // the "1." in 1./length arises because the grid bounds are from 0 to 1
    _scale = (_length > 0.) ? 1. /_length : 1.;
    
    // do the scaling
    for (size_t v = 0; v < vertices_size; ++v)
    {
      absolute_coordinates_to_grid(_vertices[v], _vertices[v]);
    }
  }
  
  
  // read edges
  size_t edges_size(0);
  while (1)
  {
    file_in >> edges_size;
    if (file_in.eof() || edges_size == 0)
      break;
    Polygon *poly = new Polygon;
    _loops.push_back(poly);
    poly->_edges.resize(edges_size);
    poly->_vertices = &_vertices;
    size_t e1, e2;
    for (size_t e = 0; e < edges_size; ++e)
    {
      file_in >> e1 >> e2;
      poly->_edges[e] = Edge(e1, e2);
    }
  }
}


