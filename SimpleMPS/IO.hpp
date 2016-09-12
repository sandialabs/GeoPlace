//
//  IO.hpp
//  simplemps
//
//
//// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000,
//// there is a non-exclusive license for use of this work by or on behalf of the U.S. Government.
//// Export of this program may require a license from the United States Government.
////
//
////  Created by Scott Alan Mitchell on 3/27/15.
//  Copyright (c) 2015 darts team. All rights reserved.
//

#ifndef simplemps_IO_hpp
#define simplemps_IO_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>

// ======== ps file

void open_ps_file(std::fstream &file, std::string fname );

// return the scaling factor
// inputs are the lower left corner, and its distances to the upper right corner
double define_plot_box( std::fstream &file, double xmin, double ymin, double Lx, double Ly );

// note the last-drawn shape will be on top and visible. Lower objects may be covered

// defines shapes:  /circ  /fcirc  /seg  /seg2  /quad_  /quad_bold  /quad_dark  /quad_medium  /quad_light  /quad_white
void define_shapes( std::fstream &file );

// quad_weight is one of "", "bold", "dark", "medium", "light",  "white"
void plot_quad( std::fstream &file, double scale, double xo, double yo, double xn, double yn, std::string quad_weight );

// filled makes a solid shaded circle. false just draws the perimieter.
void plot_circle( std::fstream &file, double scale, double x, double y, double r, bool filled );

enum plot_colors {red, green, dark_blue, light_blue, black, white};
void set_color( std::fstream &file, plot_colors color = black );

// plot the letter "x" centered at x,y with size s
void plot_ex( std::fstream &file, double scale, double x, double y, double s );
// plot the plus sign "+" centered at x,y with size s
void plot_plus( std::fstream &file, double scale, double x, double y, double s );

void plot_seg( std::fstream &file, double scale, double Ax, double Ay, double Bx, double By );

// plot the 0-1 bounding box surrounding the domain
void plot_bbox(std::fstream &file, double scale);

#endif
