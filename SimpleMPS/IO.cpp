//
//  IO.cpp
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

#include "IO.hpp"

void open_ps_file(std::fstream &file, std::string fname )
{
  // std::fstream file();
  file.open( fname.c_str(), std::ios::out);
  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;
}

std::vector< std::string > color_strings;

// defines shapes:  /circ  /fcirc  /seg  /seg2  /quad_  /quad_bold  /quad_dark  /quad_medium  /quad_light  /quad_white
void define_shapes( std::fstream &file )
{
  
  color_strings.resize( 6 );
  color_strings[red]        = "1 0 0";
  color_strings[green]      = "0 1 0";
  color_strings[dark_blue]  = "0 0 1";
  color_strings[light_blue] = "0 1 1";
  color_strings[black]      = "0 0 0";
  color_strings[white]      = "1 1 1";
  
  file << "/circ    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
  
  file << "/fcirc    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.6 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
  
  file << "/seg      % stack: x1 y1 x2 y2" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.02 setlinewidth" << std::endl; // was 0.05 for extra thick
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
  
  file << "/seg2      % stack: x1 y1 x2 y2" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.008 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
  
  // the quad will be filled with whatever color was used last
  file << "/quad_      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.0 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
  
  file << "/quad_bold      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.01 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
  
  file << "/quad_dark      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.2 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
  
  file << "/quad_medium      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.5 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
  
  file << "/quad_light      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl; // 0.9
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
  
  file << "/quad_white      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 1.0 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << "} def" << std::endl;
#pragma endregion
}

double define_plot_box( std::fstream &file, double xmin, double ymin, double Lx, double Ly )
{
  
  double scale_x, scale_y, scale;
  double shift_x, shift_y;
  
  scale_x = 6.5 / Lx;
  scale_y = 9.0 / Ly;
  
  if (scale_x < scale_y)
  {
    scale = scale_x;
    shift_x = 1.0 - xmin * scale;
    shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
  }
  else
  {
    scale = scale_y;
    shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
    shift_y = 1.0 - ymin * scale;
  }
  file << shift_x << " " << shift_y << " translate" << std::endl;
  return scale;
}


void plot_quad( std::fstream &file, double scale, double xo, double yo, double xn, double yn, std::string quad_weight )
{
  file << xo * scale << "  " << yo * scale << "  ";
  file << xn * scale << "  " << yo * scale << "  ";
  file << xn * scale << "  " << yn * scale << "  ";
  file << xo * scale << "  " << yn * scale << "  ";
  file << "quad_" << quad_weight << std::endl;
}

void plot_circle( std::fstream &file, double scale, double x, double y, double r, bool filled )
{
  file << x * scale << "  " << y * scale << "  " << r * scale << "  ";
  if (filled)
    file << "f";
  file << "circ"     << std::endl;
}

void set_color(std::fstream &file, plot_colors color )
{
  file << color_strings[color] << " setrgbcolor" << std::endl;
}

void plot_plus( std::fstream &file, double scale, double x, double y, double s )
{
  // horizontal
  file << (x - s) * scale << "  " << (y) * scale << "  ";
  file << (x + s) * scale << "  " << (y) * scale << "  ";
  file << "seg2"     << std::endl;
  
  // vertical
  file << (x) * scale << "  " << (y - s) * scale << "  ";
  file << (x) * scale << "  " << (y + s) * scale << "  ";
  file << "seg2"     << std::endl;
}
  
void plot_ex( std::fstream &file, double scale, double x, double y, double s )
{
  // ll to ur
  file << (x - s) * scale << "  " << (y - s) * scale << "  ";
  file << (x + s) * scale << "  " << (y + s) * scale << "  ";
  file << "seg2"     << std::endl;

  // ul to lr
  file << (x + s) * scale << "  " << (y - s) * scale << "  ";
  file << (x - s) * scale << "  " << (y + s) * scale << "  ";
  file << "seg2"     << std::endl;
}

void plot_seg( std::fstream &file, double scale, double Ax, double Ay, double Bx, double By )
{
  file << Ax * scale << " " << Ay * scale << "  " << Bx * scale << " " << By * scale << "  ";
  file << "seg" << std::endl;
}


void plot_bbox(std::fstream &file, double scale)
{
  // plot a big white (blank) square around the bounding box
  if (true)
  {
    plot_quad( file, scale, -1, 2, -1, 0, "white" );
    plot_quad( file, scale, -1, 2,  1, 2, "white" );
    plot_quad( file, scale,  1, 2, -1, 2, "white" );
    plot_quad( file, scale, -1, 0, -1, 2, "white" );
  }
  
  // plot 0-1 box boundary
  {
    plot_quad( file, scale, 0, 0, 1, 1, "bold" );
  }
}
