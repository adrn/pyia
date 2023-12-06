Gaia DR2: Re-normalised Unit Weight Error (RUWE) - tables of u0(G,C) and u0(G)

L. Lindegren (V1, 2018 Aug 10)

The files table_u0_*.txt contain lookup tables for the functions u0(G,C) and u0(G)
defined in GAIA-C3-TN-LU-LL-124-01. u0(G,C) is the reference value of the UWE as
function of magnitude (G) and BP-RP colour (C) such that RUWE = UWE/u0(G,C).
For sources without a colour the function u0(G) may be used instead.

Detailed file descriptions:

All files contain comma separated values (CSV) including a single header line 
with unique column names. The files can be read e.g. in TOPCAT using format 
specification CSV.
 
table_u0_g_col.txt (5375015 bytes) is a lookup table for u0(G,C). 
It has 3 columns (g_mag, bp_rp, u0) and 193251 data lines (+ header) 
for G = 3.60(0.01)21.00 and C = -1.0(0.1)10.0.

table_u0_g.txt (32502 bytes) is a lookup table for u0(G). 
It has 2 columns (g_mag, u0) and 1741 data lines (+ header) 
for G = 3.60(0.01)21.00.

table_u0_2D.txt (1782073 bytes) contains exactly the same data as the previous
two tables, only arranged in a way that may be more convenient in some cases. 
It has 113 columns (g_mag, u0g, u0m010, u0m009, ... u0p100) and 1741 data lines 
(+ header) for G = 3.60(0.01)21.00. The columns are:
g_mag = G magnitude
u0g = u0(G)
u0m010 = u0(G,-1.0) (i.e. for C = -1.0)
u0m009 = u0(G,-0.9)
...
u0p100 = u0(G,10.0)
Thus columns i = 3 to 113 contain u0(G,C) for C = (i-13)/10.

The values of G and C always refer to the centre of the bin.
