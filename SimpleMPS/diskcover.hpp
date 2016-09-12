// diskcover.hpp
//
//  Copyright (c) 2015 Sandia National Labs darts team. All rights reserved.
//
// method submitted to Euro Graphics
// Author: Mohamed S. Ebeida
// Modified 03/01/2015, for remote sensing disk cover purposes
// modification
// Author: Scott A. Mitchell

#ifndef simplemps_DISKCOVER_HPP
#define simplemps_DISKCOVER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

#include <time.h>

#include <cmath>
#include <vector>
#include <array>
#include <list>
#include <utility>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <assert.h>

#include "Point.hpp"
#include "Footprint.hpp"
#include "IO.hpp"

// =========================
// initialize the RNG that is used internally
void init_RNG(const unsigned int seed);

class Domain;
class Polygon;

const double _RF = 5.;     // How many dart throws per initial cell at each refinement level. E.g. 5

// print diagnostics about the domain boundary
// #define DEBUG_BOUNDARY 1
#define DEBUG_BOUNDARY 0

#define DEBUG_OUT_CENTERS 1

// dimension of the problem, hardcoded to 2 for domain issues
const size_t _max_children(4);

class Sampling
{
public:

  size_t _num_darts;

  size_t _num_inserted_points;

  // to do
  // provide an iterator to go over the created points

  void new_cell_points(size_t cell_points_size);
  size_t cell_points_size() const { return _cell_points_size; }
  
  Sampling() : _num_darts(0), _num_inserted_points(0), _cell_points(0), _cell_points_size(0) {}
  
  virtual ~Sampling();

protected:
  friend class SamplingIterator;
  // coordinates of sample points, by top level cell
  // note many will be null pointers, for empty cells
  // note the sampling is stored according to the Grid, a dependency
  size_t _cell_points_size;
  double** _cell_points;
private:
  void delete_cell_points();
};

class SamplingIterator
{
private:
  size_t _celli;
  Sampling *_sampling;
public:
  double *first();
  double *next();
  
  SamplingIterator(Sampling *sampling) : _sampling(sampling), _celli(0) {}
};

// Idiom
//SamplingIterator si ( sampling );
//double *p = si.first();
//while (p)
//{
//  // do something with p
//  
//  p = si.next();
//}

inline
double * SamplingIterator::first()
{
  _celli = 0;
  return next();
}

inline
Sampling::~Sampling()
{
  delete_cell_points();
}


class Grid : public Sampling
{
public:

  Grid() : _domain(0), _TWO_POW_REF_LEVEL(1), _fname_prefix("points") {}
  virtual ~Grid();

  // call this first, set up grid
  // note sampling_radius is adjusted slightly to be 1/n for some integer n
  void setup(Domain *domain, double &sampling_radius);

  // optional, initialize the grid with the given sample points
  // not implemented because the data structure for points only allows one point per cell
  void add_samples( Sampling *sampling );
  
  // optional, initialize the grid with a maximal subset of the given sample points, using the conflict radius _dm
  void add_samples_conditionally( Sampling *sampling );
  
  // call this, create the sample points
  void MPS();
  
  
  //-------------------
  //-- i/o, debug
  //-------------------
  void report_stats(double cpu_time);

  void plot( bool active_cells = true, bool disks = true, std::string fname = "points.ps" );
  
  void print_cell(size_t *i, std::ostream &out = std::cout); // indices, coordinates of corners

  // output sample points to a file
  // if debug flag is set, then the scaled-grid coordinates and explanatory text is written
  void save_point_cloud( std::string fname = "mps_points.txt", bool debug = false );
  
  // ensure no two points are too close together, etc
  void test_output();
  
  std::string _fname_prefix;
  
private:
  Domain *_domain;
  
  size_t _num_neighbors;
  int* _neighbors; // Array of indices of the neighbors in two layers

  // _n is the number of top level cells, or ??
  size_t _n, _no, _pool_size;

  size_t* _n_pow; // number of cells in dimesnions 1, 2

  // _dm = disk radius = square diagonal
  // todo: separate radius and diagonal (radius >= diagonal ) so the disk radius does not have to 1/k for k integer.
  // _s = side length of top level grid, _ss = side length of current grid squares
  double _dm, _dm_squared, _s, _ss;
  double current_diagonal() const { return _dm / ( (double) _TWO_POW_REF_LEVEL ); }
  double current_diagonal_squared() const { return _dm_squared / ( (double) (_TWO_POW_REF_LEVEL * _TWO_POW_REF_LEVEL) ); }
  
  size_t _TWO_POW_REF_LEVEL; // 2^current_refinement_level

  double _xo, _Lx; // lower left corner (same value for both x and y

  // cells outside the domain (e.g. in the ghost layer) are invalid
  // cells covered by a disk are invalid
  bool* _invalid_cells;

  bool* _empty_cells;   // array of top level grid cells. if true, then that cell does not contain a sample point (disk center)

  // cells cut by a boundary edge or containing a boundary point.
  // future optimization: store which portions of the boundary pass through this cell
  bool* _boundary_cells;

  // flat quadtree cells that we generate sample candidates in
  size_t* _active_pool_i; // 2 dimension array, of cell indices

  // refined child cells
  size_t _num_max_children; // array length
  // array of index-array of child cells that make it to the next refinement level
  size_t** _active_children;

  size_t _num_valid;

  size_t _num_ghost_layers; // number of layers of squares outside the actual domain

  // indices of a cell, 2 long
  size_t* _i; // cell_i
  size_t* _j; // cell_j

  // coordinates of a point, 2 long
  double* _u;
  double* _x;
  double* _xc;

  clock_t _start_time, _end_time;

  // ========================= methods

  // set up the grid, ready for sampling
  void initialize(double sampling_radius); // = _dm
  
  // updates invalid_cells, boundary_cells, based on the domain
  void initialize_cells_by_domain(Domain *domain);
  
  // run simple mps on the remaining cells at the current refinement level
  // returns true if the algorithm is done
  void initialize_top_pool();
  void subdivide_cells();
  bool use_grid_level(size_t refLevel);
  
  bool add_sample_conditionally( double * x )
  {
    size_t i[2];
    get_cell_indicies( x, i );
    size_t icell;
    get_cell_index(i, icell);
    return add_sample_conditionally( x, icell );
  }
  bool add_sample_conditionally( double * x, size_t icell );

  
  // convert the indices for each dimension (first argument) into  a linear cell index (second argument)
  void get_cell_index(const size_t* i, size_t &icell);
  // inverse of prior function
  // convert a linear cell index (first argument) into the indices for each dimension (second argument)
  void get_cell_indices(size_t icell, size_t* i);
  // convert coordinates of point p into the indices of the top level cell that contains it
  void get_cell_indicies( const double *p, size_t *i );
  
  // get the coordinates of the upper or lower corners, or the center, of the cell
  void get_cell_center(const size_t* i, double* center);
  void get_cell_upper(const size_t* i, double* upper);
  void get_cell_lower(const size_t* i, double* lower);
  void get_childcell_center(const size_t* i, double* center);
  void get_childcell_upper(const size_t* i, double* upper);
  void get_childcell_lower(const size_t* i, double* lower);
  
  // if true, then cell is definately outside
  // if false, the cell may or may not be inside
  // assumes current level cell, not top
  bool is_cell_outside( const size_t* i );
  
  inline bool valid_cell(size_t* i);
  
  inline void get_active_children(size_t* i, size_t icell, size_t** jj, size_t &num_active_children);
  
protected:
  
  // plotting ps files
  void plot_grid_cells( std::fstream &file, double scale );
  void plot_active_cells( std::fstream &file, double scale);
  void plot_samples( std::fstream &file, double scale, bool disks);

  
  size_t _num_cells_in_row;
  
  void allocate_working_memory();
  void generate_top_grid();
  void initialize_children();
  void initialize_neighbors();

};

typedef double* Vertex; // length 2
//typedef std::pair<double,double> Vertex;
typedef std::vector< Vertex > Vertices;
typedef std::pair<unsigned,unsigned> Edge;
typedef std::vector<Edge> Edges;


// sort vertices by angle around the origin
// usage
// std::sort(vec.begin(), vec.end(), less_than_angle(_vertices, reference_vertex_index));
struct less_than_angle
{
  Vertices &_vertices;
  int _v0;
  less_than_angle(Vertices *vertices, int v0) : _vertices( *vertices ), _v0(v0) {}

  
  double cross(double p1x, double p1y, double p2x, double p2y) {return p1x * p2y - p1y * p2x;}
  // return true if p1 < p2 in terms of angle around _v0
  inline bool operator() (const unsigned& p1, const unsigned& p2)
  {
    // _v0 is always the smallest
    if ( p1 == _v0 ) return true;
    if ( p2 == _v0 ) return false;
    
    // vector to p1 or p2 from v0
    const double p1x = _vertices[p1][0] - _vertices[_v0][0];
    const double p1y = _vertices[p1][1] - _vertices[_v0][1];
    const double p2x = _vertices[p2][0] - _vertices[_v0][0];
    const double p2y = _vertices[p2][1] - _vertices[_v0][1];
    
//    if( p1y==0. && p1x>0.) return true; //angle of p1 is 0, thus p2>=p1
//    if( p2y==0. && p2x>0.) return false; //angle of p2 is 0 , thus p1>=p2
//    if( p1y>0. &&  p2y<0.) return true; //p1 is between 0 and 180, p2 between 180 and 360
//    if( p1y<0. &&  p2y>0.) return false;
    return cross(p1x,p1y,p2x,p2y)>0.; //return true if p1 is clockwise from p2
  }
};


class Polygon // can have holes
{
private:
  // 1 = a normal polygon, -1 = a hole in a polygon, 0 = unset
  int _is_outer_boundary;
  
  // calculate the sum of the signed cw angles (-pi,pi) between p and the edges, from the first edge to the last
  double calculate_angle_sum( const double * p, size_t first_edge = 0);

public:
  Vertices *_vertices;
  Edges _edges;
  Polygon() : _vertices(0), _edges(0), _is_outer_boundary(0) {}
  ~Polygon() {_vertices = 0; _edges.clear(); }
  
  // shorthand to get the first and second vertices of an edge
  double * v0(size_t e) { return (*_vertices)[ _edges[ e ].first ]; };
  double * v1(size_t e) { return (*_vertices)[ _edges[ e ].second ]; };
  
  bool is_point_geometrically_inside( const double *p );
  bool is_point_inside( const double *p ); // takes the outer-boundary orientation of the polygon into account
  bool is_outer_boundary();

  double distance_to_boundary_squared(double *p);
  
  // sets polygon edges to be the convex hull of its Vertices
  void convex_hull();
  
  void print();
};

class Domain
{
public:
  Domain() : _offset(0.), _scale(1.), _length(1.), _min_x(0.), _min_y(0.), _frame_size(0.) {}
  virtual ~Domain();

  Vertices _vertices;
  std::vector<Polygon*> _loops;
  Polygon _convex_hull;
  void set_convex_hull()
  {
    _convex_hull._vertices = &_vertices;
    _convex_hull.convex_hull();
  }
  
  void read( std::string file_name );

  void print_vertex(size_t v, std::ostream &out = std::cout); // index and coordinates

  
  void plot_ps( std::fstream &file, double scale );

  // +/-, 0.
  // Boundary is move outward by _offset
  // Currently only +,0 offsets are supported for creating sample points.
  // But - is supported for testing if points are inside an offset polygon
  // Sample points are required to be within _offset of the boundary
  // The is_point_inside and is_cell_outside take this into account
  double _offset;
  
  bool is_point_inside( const double *p );

  // use this for conservatively excluding/discarding cells
  // if true, then the cell is definitively outside the boundary+offset
  // if false, then the cell may or may not be inside the boundary+offset
  // is_cell_outside(
  bool is_cell_outside( double* ll, double ss );
  // ll == lower left coordinates of cell
  // ss == side length
  
  double distance_to_boundary(double *p); // >= 0
  double signed_distance_to_boundary(double *p); // - if inside, + if outside, 0 if on
  
  // the input domain extends from _min_x to _min_x + _length, _min_y to _min_y + _length
  // it is scaled bye _scale and offset(-) by _min_* so that it fits in a 0,1 square.
  // the output sample points should be inversely scaled and offset
  double _scale, _length, _min_x, _min_y;
  
  // the frame is an buffer space outside the domain but inside the grid, in absolute coordinates
  double _frame_size;

  void absolute_coordinates_to_grid( const double *a , double *g )
  {
    g[0] = (a[0] - _min_x) * _scale;
    g[1] = (a[1] - _min_y) * _scale;
  }

  void grid_coordinates_to_absolute( const double *g, double *a )
  {
    a[0] = g[0] / _scale + _min_x;
    a[1] = g[1] / _scale + _min_y;
  }
  

};


typedef  std::vector< std::array< double, 2 > > PointVector;
void plot_coverage( std::vector< Footprint *> &footprints, std::vector<size_t> &toes, PointVector &placement_points, PointVector &coverage_points, std::vector<char> &is_covered, std::string fname );
void find_footprint_coverage( std::string footprint_fname );

// 0-1 matrix. For each row, list of columns such that (row,column) = 1
typedef std::vector< std::vector<size_t> > Matrix;
void write_matrix( Matrix &matrix, size_t column_size, std::string fname );
void read_matrix( Matrix &matrix, size_t &column_size, std::string fname );
// count the number of common entries to both rows
// relies on each row being in ascending order
size_t num_common_entries( std::vector<size_t> & rowA, std::vector<size_t> &rowB );


// matrix of size_t values. For each row, list the column and entry, such that (row,column) = entry.
// e.g. 1   4 3   7 2  ==> (1,4) = 3, and (1,7) = 2
typedef std::vector< std::vector< std::pair<size_t,size_t> > > MatrixEntries;
// write out entries, but scaled values as floats, with values between 0 and 1.
void write_matrix_entries( MatrixEntries &matrix, size_t column_size, size_t max_entry, std::string fname );


void plot_solution( std::vector<Footprint*> &footprints, std::vector<size_t> &toes, Matrix &placement_covers, PointVector &coverage_points, PointVector &placement_points );
void find_solution_coverage( std::string solution_fname );

// read in the points that were written out by Grid::save_point_cloud
// but don't put them in a grid
void read_point_coordinates( PointVector &points, double &r, std::string fname );

// read the solution from the integer program, the list of placements (toes) and type of footprints
void read_placements( std::vector<Footprint*> &footprints,  std::vector<size_t> &toes, std::string solution_fname = "solution.txt");


#endif
