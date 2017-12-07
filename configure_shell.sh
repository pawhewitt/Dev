
# Disable Everything 
#./configure prefix=/usr/local/SU2_CAD  --disable-CFD --disable-GEO --disable-IDE --disable-MSH --disable-SOL --disable-DEF

# Enable Everything
#./configure prefix=/usr/local/SU2_CAD

# Enable CFD

./configure prefix=/usr/local/SU2_CAD --disable-GEO --disable-IDE --disable-MSH --disable-SOL --disable-DEF --disable-DOT
