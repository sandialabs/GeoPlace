-- quick start
If pyomo, gurobi or cplex, etc., is all in place, you may

% make all sample coverage
% ./runfootprint.sh
You should eventually see an image pop up.
Additional output will be in the run directory

-- runfootprint.sh does the following

./make-all
cd run
../diskcover-sample
../diskcover-coverage
#  run your optimization code, to overwrite solution.txt
../diskcover-plot-solution


-- description

diskcover suite of tools for remote sensor footprint placement to cover a region of interest

diskcover-all performs all of the operations below, except the optimization code. Experimental.

Normal workflow follows: 
  diskcover-sample: creating domain sample points for coverage and placement - only.
  diskcover-coverage: finding coverage of sample points by footprints.
  run the external optimization code to produce the solution, where to place footprints
  diskcover-plot-solution: plotting solution, the footprints and their placements over the domain.


Run the following scripts to make the different tools:
make-all
make-sample
make-coverage
make-plot-solution	

Get the following executables:
diskcover-all
diskcover-sample
diskcover-coverage
diskcover-plot-solution


Run them in a directory with data files

input files -> tool -> output files

all: diskcover-all

sample 
  input: polygon.txt
  visualization output: points_for_coverage*.ps points_for_placement*.ps
  output: points_for_coverage.txt points_for_placement.txt
coverage 
  input: footA.txt points_for_coverage.txt points_for_placement.txt
  visualization output: footprint_covers_*.ps
  output: footprint_covers.txt covered_by_footprints.txt footprint_overlaps.txt
optimization code (external)
  input: covered_by_footprints.txt footprint_overlaps.txt
  output: solution.txt
plot: 
  input: solution.txt footA.txt points_for_coverage.txt points_for_placement.txt
  visualization output: footprint_covers.ps
  output: 

see the *-format.txt files for descriptions of the contents and layout of these files

epsilon (sample point spacing) is specified in the C++ file main-all.cpp
delta (coverage point spacing) is 0.04
footprint radius is assumed to be 0.25

these all should be replaced by user input values

only circular footprints are implemented

