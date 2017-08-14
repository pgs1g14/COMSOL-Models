Data of impedances named: Impedances_<frequency>Hz.mat
contains the vector Z with the following lines:
1 - FFT(Liner)
2 - FFT(Plane 1) at 70 mm
3 - FFT(Plane 2) at 40*h_h = 25.4 mm
4 - FFT (Mic 1) at 119mm
5 - FFT (Mic 2) at 99 mm
6 - IRW (Inlet) at 208.2 mm
7 - FC (Mic 2)
8 - ASTM
9 - FC (Liner)

The results data file named: Results_<frequency>Hz.mat 
contains the time domain results on the 3D Matrix (time, value, parameter in an position). The parameter in an positions (3rd position) are:
1 - Incident pressure on the inlet
2 - Reflected pressure on the inlet
3 - Pressure on the plane 1, multiphysics/acoustics domain
4 - Velocity on the plane 1
5 - Pressure on the plane 2, acoustics domain
6 - Velocift on the plane 2
7 - Pressure on the Mic 1
8 - Pressure on the Mic 2
9 - Aceleration on the Mic 1
10 - Aceleration on the Mic 2
11 - Average Pressure on the liner surface
12 - Average velocity on the liner surface


