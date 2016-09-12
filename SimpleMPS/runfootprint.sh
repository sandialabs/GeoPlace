cd run
./clean.sh
../diskcover-all
cd ../../geo_coverage
./createplacements.sh
cd ../SimpleMPS_EG2012/run
../diskcover-plot-solution
display footprint_covers.ps