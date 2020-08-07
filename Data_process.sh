#!/bin/bash
echo Interpolating part 0 
./Orbit_interp.out -i Y34D0522_0610_SC1 -o SC1_part_0 -s 0 -n 1001

echo Interpolating part 1
./Orbit_interp.out -i Y34D0522_0610_SC1 -o SC1_part_1 -s 1000 -n 1001

echo Interpolating part 2
./Orbit_interp.out -i Y34D0522_0610_SC1 -o SC1_part_2 -s 2000 -n 1001

cat SC1_part_0 >> GEI_SC1
cat SC1_part_1 >> GEI_SC1
cat SC1_part_2 >> GEI_SC1
rm -f SC1_part_*

echo Converting to GCS
./getGCS_TT.out -i GEI_SC1 -o GCS_SC1

mv GEI_SC1 ./Processed/
mv GCS_SC1 ./Processed/

echo All Compelete!
