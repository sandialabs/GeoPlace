This software project makes use of and extends the Simple MPS algorithm developed by Ebeida, et al. and reported here: http://dl.acm.org/citation.cfm?id=2322144.

Some results of applying GeoPlace software to place multiple different footprint shapes to a variety of different coverage regions are reported here: ...


#Format file for subfootprint problem definition
#Line 1: row dimension (m), column dimension (n)
#Line 2: Maximum potential priority
#Lines 3-infinity: m x n dimensional array of chips with each entry denoting chip priority. Entries are comma delimited. Line 3 starts at the mth row and 0th column.

#Example
4,3
10

1,0,0
0,3,0
0,0,2
0,4,6
