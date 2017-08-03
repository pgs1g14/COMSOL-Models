%
% LinerAxi_PT_155dB_500Hz.m
%
% Model exported on Jun 30 2017, 14:23 by COMSOL 5.3.0.260.

%% Constants Header = comments, name of file, path, ambient, positions, geometry, sound properties
tic;
%Model physical constants
freq=500; %Frequency to excite and generate the mesh
%------------------------------------------------
w=2*pi()*freq; %angular frequency  
TC=20; %Temperature of the measurements to run the simulation �C
T= 273.15+TC; %Temperature of the fluid [K]
ATMp= 101450; %Atmosferic pressure [Pa]
R=287.058; %Specific air constant [J/kg.K]
rho=ATMp/(R*T); %Specific gas density [kg/m^3]
mi= 1.8205e-5;%dynamic viscosity [kg/m.s]
ni=mi/rho; %kinematic viscosity [m^2/s]
c0= 331.29*sqrt(1+((T-273.15)/273.15)); %Fluid sound velocity [m/s]
x1=0.099+0.02; %Mic 1 position
x2=0.099; %Mic 2 position
SPL=155; %Sound pressure level dB
k0=w/c0; %wave number
terms=1; %Number of terms related to tone components
j=1; %incremental variable;
%Liner geometry
w_tube=0.00517; %impedance tube diameter [m]
h_r=0.019; %Liner cavity depth - resonator height [m]
h_slit=0.000635; %Liner hole thickness [m]
POA=0.0366; %percentage of open area of the sample in the holder
%POA=0.0518; %percentage of open area of the flanged liner sample
h_in=0.2082; %impedance tube length [m]
d_h=0.00099; %Hole diameter 0,99 mm - on the holder
yPlane1=0.07;
yPlane2=40*h_slit;
elementHole=2;
elementNear=15;
elementFar=24;
dt_sol=1*10^-4;
Tfim=0.01;

comment=['Liner Axisimetrical PT thickness ' num2str(h_slit) ' mm ' num2str(SPL) 'dB ' num2str(freq) 'Hz. Impedance tube model with a sample in a holder with POA 3.66%, d=' num2str(d_h) ' mm, h=19.0mm nominal values.'];
file=['LinerAxi_PT_155dB_' num2str(freq) 'Hz'];
%Comsol Instructions Windows
pasta='C:\Local\Pablo\COMSOL\Circular Axi model\LinerAxi_PT_155dB_1500Hz\';
%Comsol instructions to Linux
%pasta='~/Models/Results/';
%addpath('/local/software/comsol/5.2a/mli')
%addpath('C:\Program Files\COMSOL\COMSOL53\Multiphysics\mli') on windows
%mphstart;


%% Start the model
import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.showProgress(true); %Show the progress bar
ModelUtil.setServerBusyHandler(ServerBusyHandler(10));

model = ModelUtil.create('Model');

model.modelPath(['C:\Local\Pablo\COMSOL\Circular Axi model\LinerAxi_PT_155dB_' num2str(freq) 'Hz']);

model.label([file '.mph']);

model.comments(comment);

mphsave(model,['../' file '.mph']);

%% Parameters
model.param.set('w_tube', [num2str(w_tube) '[m]'], 'Tube width');
model.param.set('h_r', [num2str(h_r) '[m]'], 'Resonator height');
model.param.set('h_slit', [num2str(h_slit) '[m]'], 'Slit thickness');
model.param.set('w_slit', [num2str(d_h) '[m]'], 'width');
model.param.set('h_in', [num2str(h_in) '[m]'], 'Inlet height');
model.param.set('y0', [num2str(h_in+h_slit/2) '[m]'], 'Inlet location y-coordinate');
model.param.set('h_tube', [num2str(h_r+h_in+h_slit) '[m]'], 'Total tube length');
model.param.set('rho0', [num2str(rho) '[kg/m^3]'], 'Density');
model.param.set('c0', [num2str(c0) '[m/s]'], 'Speed of sound');
model.param.set('beta', '1.2', 'Nonlinear coefficent');
model.param.set('L0', num2str(SPL), 'Incident wave amplitude dB');
model.param.set('p0', [num2str(10^(SPL/20)*20e-6) '[Pa]'], 'Incident wave amplitude Pa');
model.param.set('f0', [num2str(freq) '[Hz]'], 'Driving frequency');
model.param.set( 'Tstart' , '7/f0', 'Stat time for post processing 7 cycles');
model.param.set('lambda0', [num2str(c0/freq) '[m]'], 'Wavelength at f0');
model.param.set('dvisc', '220[um]*sqrt(100[Hz]/f0)', 'Viscous boundary layer thickness at f0');
model.param.set('k0', '2*pi/lambda0', 'Wave number at f0');
model.param.set('omega0', w, 'Angular frequency');
model.param.set('T0', [num2str(1/freq) '[s]'], 'Period');
model.param.set('dt_sol', num2str(dt_sol), 'Solver time step (resolve f0 with T0/30)');
model.param.set('lmesh', 'c0*dt_sol', 'Minimum element size acoustic mesh');
model.param.set('Tend',['Tstart+' num2str(Tfim)], 'End time for post processing');


model.component.create('comp1', false);

model.component('comp1').geom.create('geom1', 2);

model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');
model.result.table.create('tbl4', 'Table');
model.result.table.create('tbl5', 'Table');
model.result.table.create('tbl6', 'Table');
model.result.table.create('tbl9', 'Table');
model.result.table.create('tbl10', 'Table');
model.result.table.create('tbl11', 'Table');
model.result.table.create('tbl12', 'Table');
model.result.table.create('tbl13', 'Table');
model.result.table.create('tbl14', 'Table');
model.result.table.create('evl2', 'Table');
model.result.table.create('tbl15', 'Table');
model.result.table.create('tbl16', 'Table');
%% Geometry
model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', {'0' 'h_tube/2-(h_r+h_slit/2)'});
model.component('comp1').geom('geom1').feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('r1').set('layerbottom', false);
model.component('comp1').geom('geom1').feature('r1').set('layertop', true);
model.component('comp1').geom('geom1').feature('r1').set('size', {'w_tube' 'h_tube'});
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('base', 'center');
model.component('comp1').geom('geom1').feature('r2').set('layername', {'Layer 1'});
model.component('comp1').geom('geom1').feature('r2').setIndex('layer', '(w_tube-w_slit)/2', 0);
model.component('comp1').geom('geom1').feature('r2').set('layerleft', true);
model.component('comp1').geom('geom1').feature('r2').set('layerright', true);
model.component('comp1').geom('geom1').feature('r2').set('layerbottom', false);
model.component('comp1').geom('geom1').feature('r2').set('size', {'w_tube' 'h_slit'});
model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('r3').set('base', 'center');
model.component('comp1').geom('geom1').feature('r3').set('size', {'10*w_slit' '40*h_slit'});
model.component('comp1').geom('geom1').create('uni1', 'Union');
model.component('comp1').geom('geom1').feature('uni1').selection('input').set({'r1' 'r2' 'r3'});
model.component('comp1').geom('geom1').create('del1', 'Delete');
model.component('comp1').geom('geom1').feature('del1').selection('input').init(2);
model.component('comp1').geom('geom1').feature('del1').selection('input').set('uni1(1)', [2 5 8 9]);
model.component('comp1').geom('geom1').create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('pol1').set('type', 'open');
model.component('comp1').geom('geom1').feature('pol1').set('x', '-w_tube/2 w_tube/2');
model.component('comp1').geom('geom1').feature('pol1').set('y', '40*h_slit 40*h_slit');
model.component('comp1').geom('geom1').create('pol2', 'Polygon');
model.component('comp1').geom('geom1').feature('pol2').set('type', 'open');
model.component('comp1').geom('geom1').feature('pol2').set('x', '-w_tube/2 w_tube/2');
model.component('comp1').geom('geom1').feature('pol2').set('y', '0.07 0.07');
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').feature('pt1').setIndex('p', 'w_tube/2', 0, 0);
model.component('comp1').geom('geom1').feature('pt1').setIndex('p', '0.119', 1, 0);
model.component('comp1').geom('geom1').create('pt2', 'Point');
model.component('comp1').geom('geom1').feature('pt2').setIndex('p', 'w_tube/2', 0, 0);
model.component('comp1').geom('geom1').feature('pt2').setIndex('p', '0.099', 1, 0);
model.component('comp1').geom('geom1').run;
%% Variables
model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('Pin', 'p0*sin(omega0*t+k0*(y-y0))', 'Incident pressure');
model.component('comp1').variable('var1').set('rhoL', 'rho0+p/c0^2', 'Linear density relation for air');
model.component('comp1').variable('var1').set('rhoNL', 'rho0+p/c0^2-1/(rho0*c0^4)*(beta-1)*p^2', 'Nonlinear density relation for air');

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat1').propertyGroup.create('idealGas', 'Ideal gas');

model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl.create('aveop2', 'Average');
model.component('comp1').cpl.create('aveop3', 'Average');
model.component('comp1').cpl('intop1').selection.geom('geom1', 1);
model.component('comp1').cpl('intop1').selection.set([10]);
model.component('comp1').cpl('aveop1').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop1').selection.set([5 15 19 22 26]);
model.component('comp1').cpl('aveop2').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop2').selection.set([7]);
model.component('comp1').cpl('aveop3').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop3').selection.set([9]);
%% Physics
model.component('comp1').physics.create('actd', 'TransientPressureAcoustics', 'geom1');
model.component('comp1').physics('actd').field('pressure').field('p2');
model.component('comp1').physics('actd').selection.set([3 4]);
model.component('comp1').physics('actd').create('pwr1', 'PlaneWaveRadiation', 1);
model.component('comp1').physics('actd').feature('pwr1').selection.set([10]);
model.component('comp1').physics('actd').feature('pwr1').create('ipf1', 'IncidentPressureField', 1);
model.component('comp1').physics('actd').create('nvel1', 'NormalVelocity', 1);
model.component('comp1').physics('actd').feature('nvel1').selection.set([7]);
model.component('comp1').physics.create('spf', 'LaminarFlow', 'geom1');
model.component('comp1').physics('spf').selection.set([1 2 5 6 7]);
model.component('comp1').physics('spf').create('bs1', 'BoundaryStress', 1);
model.component('comp1').physics('spf').feature('bs1').selection.set([7]);
model.component('comp1').physics('spf').create('wall2', 'Wall', 1);
model.component('comp1').physics('spf').feature('wall2').selection.set([1 2 3 4 5 24 26 27 28]);
%% Mesh
model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('map1').selection.set([3 4]);
model.component('comp1').mesh('mesh1').feature('map1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').create('size2', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').create('size3', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').create('cr1', 'CornerRefinement');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').selection.set([1 2]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').selection.set([5 6]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').selection.set([7]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('cr1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('cr1').selection.set([5 6 7]);

model.component('comp1').mesh('mesh1').feature('map1').set('adjustedgdistr', true);
model.component('comp1').mesh('mesh1').feature('map1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('map1').feature('size1').set('hmax', 'lmesh');
model.component('comp1').mesh('mesh1').feature('map1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('map1').feature('size1').set('hcurve', 1);
model.component('comp1').mesh('mesh1').feature('map1').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('map1').feature('size1').set('hcurveactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmax', ['dvisc*' num2str(elementFar)]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmin', ['dvisc*' num2str(elementNear)]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hgrad', 1.15);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hcurveactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hnarrowactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hgradactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('hauto', 4);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('hmax', ['dvisc*' num2str(elementNear)]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('hmin', ['dvisc*' num2str(elementHole)]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('hgrad', 1.13);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('hcurveactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('hnarrowactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size2').set('hgradactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('hmax', ['dvisc*' num2str(elementHole)]);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('hmin', 'dvisc');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('hgrad', 1.05);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('hcurveactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('hnarrowactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size3').set('hgradactive', false);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('cr1').selection('boundary').set([13 15 17 20 21 22]);
model.component('comp1').mesh('mesh1').run;
%% Tables
model.result.table('tbl1').comments('Point Evaluation 1 (Pin)');
model.result.table('tbl1').label('Point Evaluation 1 (Pin)');
model.result.table('tbl2').comments('Point Evaluation 1 (Pout)');
model.result.table('tbl2').label('Point Evaluation 1 (Pout)');
model.result.table('tbl3').label('Pressure Liner');
model.result.table('tbl3').comments('(aveop1(p))');
model.result.table('tbl4').label('Velocity Liner');
model.result.table('tbl4').comments('(aveop1(v))');
model.result.table('tbl5').label('Ave P Plane Multiphysics');
model.result.table('tbl6').label('Ave v Plane Mutiphysics');
model.result.table('tbl9').label('Pressure Mic 2');
model.result.table('tbl10').label('Aceleration Mic 2');
model.result.table('tbl11').label('Pressure Mic 1');
model.result.table('tbl12').comments('actd.ay');
model.result.table('tbl13').label('Pressure incident Inlet');
model.result.table('tbl13').comments('Pin');
model.result.table('tbl14').label('Pressure reflected Inlet');
model.result.table('tbl14').comments('p2-Pin');
model.result.table('evl2').label('Evaluation 2D');
model.result.table('evl2').comments('Interactive 2D values');
model.result.table('tbl15').label('Ave P Plane Acoustic');
model.result.table('tbl15').comments('P plane ac (aveop3(p2))');
model.result.table('tbl16').label('Av vel Plane Acoustic');
model.result.table('tbl16').comments('V plane ac (aveop3(actd.a_inst/(2*pi*f0)))');

model.component('comp1').view('view1').axis.set('xmin', -0.012346153147518635);
model.component('comp1').view('view1').axis.set('xmax', 0.01584230735898018);
model.component('comp1').view('view1').axis.set('ymin', -0.007656451314687729);
model.component('comp1').view('view1').axis.set('ymax', 0.023603826761245728);
model.component('comp1').view('view1').axis.set('abstractviewlratio', 0.0742705687880516);
model.component('comp1').view('view1').axis.set('abstractviewrratio', 0.08861720561981201);
model.component('comp1').view('view1').axis.set('abstractviewbratio', 0.06006757169961929);
model.component('comp1').view('view1').axis.set('abstractviewtratio', -0.8190760612487793);
model.component('comp1').view('view1').axis.set('abstractviewxscale', 5.2716935897478834E-5);
model.component('comp1').view('view1').axis.set('abstractviewyscale', 5.271693953545764E-5);

model.component('comp1').material('mat1').label('Air');
model.component('comp1').material('mat1').set('family', 'air');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/8.314/T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('dermethod', 'manual');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argders', {'pA' 'd(pA*0.02897/8.314/T,pA)'; 'T' 'd(pA*0.02897/8.314/T,T)'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA' '0' '1'; 'T' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*287*T)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('dermethod', 'manual');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argders', {'T' 'd(sqrt(1.4*287*T),T)'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T[1/K])[Pa*s]');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T[1/K])[J/(kg*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(pA[1/Pa],T[1/K])[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]'});
model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T[1/K])[m/s]');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', '0');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('def').addInput('pressure');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.component('comp1').material('mat1').propertyGroup('idealGas').set('Rs', '287[J/kg/K]');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('pressure');

model.component('comp1').cpl('aveop1').label('Average Liner');
model.component('comp1').cpl('aveop2').label('Average Plane Interface');
model.component('comp1').cpl('aveop3').label('Average Plane Acoustic');

model.component('comp1').physics('actd').feature('tpam1').set('rho', 'rho0');
model.component('comp1').physics('actd').feature('tpam1').set('c', 'c0');
model.component('comp1').physics('actd').feature('pwr1').feature('ipf1').set('p', 'Pin');
model.component('comp1').physics('actd').feature('nvel1').set('Type', 'vel');
model.component('comp1').physics('actd').feature('nvel1').set('vel', {'u'; 'v'; '0'});
model.component('comp1').physics('spf').prop('PhysicalModelProperty').set('Compressibility', 'CompressibleMALT03');
model.component('comp1').physics('spf').feature('fp1').set('rho', 'rhoL');
model.component('comp1').physics('spf').feature('bs1').set('BoundaryCondition', 'NormalStress');
model.component('comp1').physics('spf').feature('bs1').set('f0', 'p2');
model.component('comp1').physics('spf').feature('wall2').set('BoundaryCondition', 'Slip');

model.frame('material1').sorder(1);

model.component('comp1').physics('actd').feature('tpam1').set('rho_mat', 'userdef');
model.component('comp1').physics('actd').feature('tpam1').set('c_mat', 'userdef');
model.component('comp1').physics('spf').feature('fp1').set('rho_mat', 'userdef');

model.result.numerical.create('pev1', 'EvalPoint');
model.result.numerical.create('pev2', 'EvalPoint');
model.result.numerical.create('gev3', 'EvalGlobal');
model.result.numerical.create('gev4', 'EvalGlobal');
model.result.numerical.create('gev5', 'EvalGlobal');
model.result.numerical.create('gev6', 'EvalGlobal');
model.result.numerical.create('pev3', 'EvalPoint');
model.result.numerical.create('pev4', 'EvalPoint');
model.result.numerical.create('pev5', 'EvalPoint');
model.result.numerical.create('pev6', 'EvalPoint');
model.result.numerical.create('pev7', 'EvalPoint');
model.result.numerical.create('pev8', 'EvalPoint');
model.result.numerical.create('gev7', 'EvalGlobal');
model.result.numerical.create('gev8', 'EvalGlobal');
%% Study - Solution
model.study.create('std1');
model.study('std1').create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.study('std1').label('Study 1 - ACTD and SPF');
model.study('std1').feature('time').set('tlist', 'range(0,dt_sol,Tend)');
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', '1e-3');
model.study('std1').feature('time').set('plot', true);
model.study('std1').feature('time').set('probesel', 'none');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('v1').set('clist', {'range(0,T0/30,Tend)'});
model.sol('sol1').feature('t1').set('control', 'user');
model.sol('sol1').feature('t1').set('tlist', 'range(0,T0/30,Tend)');
model.sol('sol1').feature('t1').set('rtol', 0.001);
model.sol('sol1').feature('t1').set('plot', true);
model.sol('sol1').feature('t1').set('probesel', 'none');
model.sol('sol1').feature('t1').feature('fc1').set('dtech', 'auto');

mphsave(model,['C:\Local\Pablo\COMSOL\Perforated 5.4% Laminar flow PT\Multiphysics model\Result files\' file '\MATLAB\' file '.m']);
mphsave(model,['C:\Local\Pablo\COMSOL\Perforated 5.4% Laminar flow PT\Multiphysics model\COMSOL files\' file '.mph']);

model.sol('sol1').runAll;
%% Set Results

model.result.numerical('pev1').set('data', 'dset1');
model.result.numerical('pev2').set('data', 'dset1');
model.result.numerical('gev3').set('data', 'dset1');
model.result.numerical('gev4').set('data', 'dset1');
model.result.numerical('gev5').set('data', 'dset1');
model.result.numerical('gev6').set('data', 'dset1');
model.result.numerical('pev3').set('data', 'dset1');
model.result.numerical('pev4').set('data', 'dset1');
model.result.numerical('pev5').set('data', 'dset1');
model.result.numerical('pev6').set('data', 'dset1');
model.result.numerical('pev7').set('data', 'dset1');
model.result.numerical('pev8').set('data', 'dset1');
model.result.numerical('gev7').set('data', 'dset1');
model.result.numerical('gev8').set('data', 'dset1');

model.result.numerical('pev1').label('Point Evaluation: RMS Incident');
model.result.numerical('pev1').set('table', 'tbl1');
model.result.numerical('pev1').set('expr', {'Pin'});
model.result.numerical('pev1').set('unit', {'Pa'});
model.result.numerical('pev1').set('descr', {'Incident pressure'});
model.result.numerical('pev1').set('dataseries', 'rms');
model.result.numerical('pev2').label('Point Evaluation: RMS Reflected');
model.result.numerical('pev2').set('table', 'tbl2');
model.result.numerical('pev2').set('expr', {'p2-Pin'});
model.result.numerical('pev2').set('unit', {'Pa'});
model.result.numerical('pev2').set('descr', {''});
model.result.numerical('pev2').set('dataseries', 'rms');
model.result.numerical('gev3').label('P Liner');
model.result.numerical('gev3').set('table', 'tbl3');
model.result.numerical('gev3').set('expr', {'aveop1(p)'});
model.result.numerical('gev3').set('unit', {'Pa'});
model.result.numerical('gev3').set('descr', {''});
model.result.numerical('gev4').label('V Liner');
model.result.numerical('gev4').set('table', 'tbl4');
model.result.numerical('gev4').set('expr', {'aveop1(v)'});
model.result.numerical('gev4').set('unit', {'m/s'});
model.result.numerical('gev4').set('descr', {''});
model.result.numerical('gev5').label('P plane 1');
model.result.numerical('gev5').set('table', 'tbl5');
model.result.numerical('gev5').set('expr', {'aveop2(p)'});
model.result.numerical('gev5').set('unit', {'Pa'});
model.result.numerical('gev5').set('descr', {''});
model.result.numerical('gev6').label('V plane 1');
model.result.numerical('gev6').set('table', 'tbl6');
model.result.numerical('gev6').set('expr', {'aveop2(v)'});
model.result.numerical('gev6').set('unit', {'m/s'});
model.result.numerical('gev6').set('descr', {''});
model.result.numerical('pev3').label('P Mic 2');
model.result.numerical('pev3').selection.set([24]);
model.result.numerical('pev3').set('table', 'tbl9');
model.result.numerical('pev3').set('expr', {'p2'});
model.result.numerical('pev3').set('unit', {'Pa'});
model.result.numerical('pev3').set('descr', {'Pressure'});
model.result.numerical('pev4').label('V Mic 2');
model.result.numerical('pev4').selection.set([24]);
model.result.numerical('pev4').set('table', 'tbl10');
model.result.numerical('pev4').set('expr', {'actd.ay'});
model.result.numerical('pev4').set('unit', {'m/s^2'});
model.result.numerical('pev4').set('descr', {'Local acceleration, y component'});
model.result.numerical('pev5').label('P Mic 1');
model.result.numerical('pev5').set('table', 'tbl11');
model.result.numerical('pev5').set('expr', {'p2'});
model.result.numerical('pev5').set('unit', {'Pa'});
model.result.numerical('pev5').selection.set([25]);
model.result.numerical('pev5').set('descr', {'Pressure'});
model.result.numerical('pev6').label('V Mic 1');
model.result.numerical('pev6').set('table', 'tbl12');
model.result.numerical('pev6').set('expr', {'actd.ay'});
model.result.numerical('pev6').set('unit', {'m/s^2'});
model.result.numerical('pev6').selection.set([25]);
model.result.numerical('pev6').set('descr', {'Local acceleration, y component'});
model.result.numerical('pev7').label('Pressure incident Inlet');
model.result.numerical('pev7').set('table', 'tbl13');
model.result.numerical('pev7').set('expr', {'Pin'});
model.result.numerical('pev7').set('unit', {'Pa'});
model.result.numerical('pev7').selection.set([26]);
model.result.numerical('pev7').set('descr', {'Incident pressure'});
model.result.numerical('pev8').label('Pressure reflected Inlet');
model.result.numerical('pev8').selection.set([26]);
model.result.numerical('pev8').set('table', 'tbl14');
model.result.numerical('pev8').set('expr', {'p2-Pin'});
model.result.numerical('pev8').set('unit', {'Pa'});
model.result.numerical('pev8').set('descr', {''});
model.result.numerical('gev7').label('P plane ac');
model.result.numerical('gev7').set('table', 'tbl15');
model.result.numerical('gev7').set('expr', {'aveop3(p2)'});
model.result.numerical('gev7').set('unit', {'Pa'});
model.result.numerical('gev7').set('descr', {'Average Pressure Plane Acoustic'});
model.result.numerical('gev8').label('V plane ac');
model.result.numerical('gev8').set('table', 'tbl16');
model.result.numerical('gev8').set('expr', {'aveop3(actd.ay)'});
model.result.numerical('gev8').set('unit', {'m/s'});
model.result.numerical('gev8').set('descr', {'Average Velocity Plane Acoustic'});


model.result.numerical('pev1').setResult;
model.result.numerical('pev2').setResult;
model.result.numerical('gev3').setResult;
model.result.numerical('gev4').setResult;
model.result.numerical('gev5').setResult;
model.result.numerical('gev6').setResult;
model.result.numerical('pev3').setResult;
model.result.numerical('pev4').setResult;
model.result.numerical('pev5').setResult;
model.result.numerical('pev6').setResult;
model.result.numerical('pev7').setResult;
model.result.numerical('pev8').setResult;
model.result.numerical('gev7').setResult;
model.result.numerical('gev8').setResult;

%% Set figures
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup1D');
model.result.create('pg3', 'PlotGroup2D');
model.result.create('pg5', 'PlotGroup2D');
model.result.create('pg6', 'PlotGroup1D');
model.result.create('pg7', 'PlotGroup1D');
model.result.create('pg8', 'PlotGroup2D');
model.result.create('pg9', 'PlotGroup1D');
model.result.create('pg10', 'PlotGroup2D');
model.result.create('pg11', 'PlotGroup1D');
model.result.create('pg13', 'PlotGroup1D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').create('surf2', 'Surface');
model.result('pg2').create('ptgr2', 'PointGraph');
model.result('pg2').feature('ptgr2').selection.set([13 14 15 18 19 22]);
model.result('pg3').create('surf1', 'Surface');
model.result('pg5').create('surf1', 'Surface');
model.result('pg6').create('ptgr1', 'PointGraph');
model.result('pg6').feature('ptgr1').selection.set([24 25 26]);
model.result('pg7').create('ptgr1', 'PointGraph');
model.result('pg7').create('ptgr2', 'PointGraph');
model.result('pg7').feature('ptgr1').selection.set([26]);
model.result('pg7').feature('ptgr2').selection.set([26]);
model.result('pg8').create('surf1', 'Surface');
model.result('pg9').create('ptgr1', 'PointGraph');
model.result('pg9').feature('ptgr1').selection.set([26]);
model.result('pg10').create('con1', 'Contour');
model.result('pg10').create('con2', 'Contour');
model.result('pg11').create('tblp1', 'Table');
model.result('pg11').create('tblp2', 'Table');
model.result('pg11').create('tblp3', 'Table');
model.result('pg13').create('tblp1', 'Table');
model.result('pg13').create('tblp2', 'Table');
model.result('pg13').create('ptgr1', 'PointGraph');
model.result('pg13').create('tblp3', 'Table');
model.result('pg13').feature('ptgr1').selection.set([24 25 26]);
model.result.export.create('anim1', 'Animation');
model.result.export.create('anim3', 'Animation');
model.result.export.create('data1', 'Data');
model.result.export.create('img1', 'Image2D');
model.result.export.create('plot1', 'Plot');
model.result.export.create('plot2', 'Plot');
model.result.export.create('plot3', 'Plot');

model.result('pg1').label('Presure');
model.result('pg1').feature('surf1').set('expr', 'p');
model.result('pg1').feature('surf1').set('descr', 'Pressure');
model.result('pg1').feature('surf1').set('rangecoloractive', true);
model.result('pg1').feature('surf1').set('rangecolormin', -3000);
model.result('pg1').feature('surf1').set('rangecolormax', 3000);
model.result('pg1').feature('surf1').set('rangedataactive', true);
model.result('pg1').feature('surf1').set('rangedatamax', 3000);
model.result('pg1').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('surf2').set('expr', 'p2');
model.result('pg1').feature('surf2').set('descr', 'Pressure');
model.result('pg1').feature('surf2').set('rangecoloractive', true);
model.result('pg1').feature('surf2').set('rangecolormin', -3000);
model.result('pg1').feature('surf2').set('rangecolormax', 3000);
model.result('pg1').feature('surf2').set('resolution', 'normal');
model.result('pg1').feature('surf2').set('rangedataactive', true);
model.result('pg1').feature('surf2').set('rangedatamax', 3000);
model.result('pg1').feature('surf2').set('rangedatamin', -3000);
model.result('pg2').label('Pressure at some positions');
model.result('pg2').set('xlabel', 'Time (s)');
model.result('pg2').set('ylabel', 'Pressure (Pa)');
model.result('pg2').set('xlabelactive', false);
model.result('pg2').set('ylabelactive', false);
model.result('pg2').feature('ptgr2').set('expr', 'p');
model.result('pg2').feature('ptgr2').set('descr', 'Pressure');
model.result('pg2').feature('ptgr2').set('legend', true);
model.result('pg2').feature('ptgr2').set('legendmethod', 'manual');
model.result('pg2').feature('ptgr2').set('legends', {'Hole inside edge' 'Hole outside edge' 'Middle cavity' 'Middle tube' 'Back cavity' 'Multiphysics boundary'});
model.result('pg3').label('Velocity');
model.result('pg3').set('looplevel', [313]);
model.result('pg3').set('axisactive', true);
model.result('pg3').set('axisnotation', 'scientific');
model.result('pg3').set('axiscommonexp', false);
model.result('pg3').set('axistrailingzeros', true);
model.result('pg3').set('axisprecision', 2);
model.result('pg3').feature('surf1').set('expr', 'spf.U');
model.result('pg3').feature('surf1').set('unit', 'm/s');
model.result('pg3').feature('surf1').set('descr', 'Velocity magnitude');
model.result('pg3').feature('surf1').set('rangecoloractive', true);
model.result('pg3').feature('surf1').set('rangecolormax', 55);
model.result('pg3').feature('surf1').set('rangedataactive', true);
model.result('pg3').feature('surf1').set('rangedatamin', -55);
model.result('pg3').feature('surf1').set('rangedatamax', 55);
model.result('pg3').feature('surf1').set('colortable', 'Wave');
model.result('pg3').feature('surf1').set('resolution', 'normal');
model.result('pg5').label('Vorticity');
model.result('pg5').set('looplevel', [316]);
model.result('pg5').feature('surf1').set('expr', '10*log(spf.vort_magn)');
model.result('pg5').feature('surf1').set('unit', '');
model.result('pg5').feature('surf1').set('descr', '10*log(spf.vort_magn)');
model.result('pg5').feature('surf1').set('resolution', 'normal');
model.result('pg6').label('Pressure at Slit FFT');
model.result('pg6').set('xlabel', 'Time (s)');
model.result('pg6').set('xlabelactive', true);
model.result('pg6').set('ylabel', 'Pressure (Pa)');
model.result('pg6').set('ylabelactive', true);
model.result('pg6').feature('ptgr1').set('expr', 'p2');
model.result('pg6').feature('ptgr1').set('descr', 'Pressure');
model.result('pg6').feature('ptgr1').set('xdata', 'spectrum');
model.result('pg6').feature('ptgr1').set('nfreqs', 322);
model.result('pg6').feature('ptgr1').set('freqrangeactive', true);
model.result('pg6').feature('ptgr1').set('freqmin', 100);
model.result('pg6').feature('ptgr1').set('freqmax', 5000);
model.result('pg6').feature('ptgr1').set('scale', true);
model.result('pg6').feature('ptgr1').set('linewidth', 2);
model.result('pg6').feature('ptgr1').set('legend', true);
model.result('pg6').feature('ptgr1').set('legendmethod', 'manual');
model.result('pg6').feature('ptgr1').set('legends', {'Tap - near' 'Tap - far' 'Inlet'});
model.result('pg7').label('Incident and Reflected');
model.result('pg7').set('xlabel', 'Time (s)');
model.result('pg7').set('legendpos', 'lowerleft');
model.result('pg7').set('xlabelactive', false);
model.result('pg7').feature('ptgr1').set('expr', 'Pin');
model.result('pg7').feature('ptgr1').set('descr', 'Incident pressure');
model.result('pg7').feature('ptgr1').set('linewidth', 2);
model.result('pg7').feature('ptgr1').set('legend', true);
model.result('pg7').feature('ptgr1').set('legendmethod', 'manual');
model.result('pg7').feature('ptgr1').set('legends', {'Incident' 'Pressure'});
model.result('pg7').feature('ptgr2').set('expr', 'p2-Pin');
model.result('pg7').feature('ptgr2').set('descr', 'p2-Pin');
model.result('pg7').feature('ptgr2').set('linewidth', 2);
model.result('pg7').feature('ptgr2').set('legend', true);
model.result('pg7').feature('ptgr2').set('legendmethod', 'manual');
model.result('pg7').feature('ptgr2').set('legends', {'Reflected' 'Pressure'});
model.result('pg8').label('Linearity for Density Approximation');
model.result('pg8').set('looplevel', [300]);
model.result('pg8').feature('surf1').set('expr', 'rhoL/rho0');
model.result('pg8').feature('surf1').set('unit', '1');
model.result('pg8').feature('surf1').set('descr', 'rhoL/rho0');
model.result('pg8').feature('surf1').set('resolution', 'normal');
model.result('pg9').label('Reflections Coefficient');
model.result('pg9').set('titletype', 'manual');
model.result('pg9').set('title', 'Reflection Coefficient');
model.result('pg9').set('xlabel', 'f (Hz)');
model.result('pg9').set('xlabelactive', true);
model.result('pg9').set('ylabel', 'abs(R)');
model.result('pg9').set('ylabelactive', true);
model.result('pg9').set('axislimits', true);
model.result('pg9').set('xmin', 100);
model.result('pg9').set('xmax', 3100);
model.result('pg9').set('ymin', 0);
model.result('pg9').feature('ptgr1').set('expr', 'sqrt(timeint(Tstart,Tend-T0,(p2-Pin)^2))/sqrt(timeint(Tstart,Tend-T0,Pin^2))');
model.result('pg9').feature('ptgr1').set('unit', '');
model.result('pg9').feature('ptgr1').set('descr', 'sqrt(timeint(Tstart,Tend-T0,(p2-Pin)^2))/sqrt(timeint(Tstart,Tend-T0,Pin^2))');
model.result('pg9').feature('ptgr1').set('xdata', 'expr');
model.result('pg9').feature('ptgr1').set('xdataexpr', 'f0');
model.result('pg9').feature('ptgr1').set('xdataunit', 'Hz');
model.result('pg9').feature('ptgr1').set('xdatadescr', '');
model.result('pg9').feature('ptgr1').set('linestyle', 'none');
model.result('pg9').feature('ptgr1').set('linewidth', 3);
model.result('pg9').feature('ptgr1').set('linemarker', 'circle');
model.result('pg9').feature('ptgr1').set('markerpos', 'datapoints');
model.result('pg10').label('Pressure Contours');
model.result('pg10').set('looplevel', [241]);
model.result('pg10').set('showlegendsmaxmin', true);
model.result('pg10').feature('con1').set('expr', 'p2');
model.result('pg10').feature('con1').set('descr', 'Pressure');
model.result('pg10').feature('con1').set('levelmethod', 'levels');
model.result('pg10').feature('con1').set('levels', 'range(-2000,100,2000)');
model.result('pg10').feature('con1').set('resolution', 'normal');
model.result('pg10').feature('con2').set('expr', 'p');
model.result('pg10').feature('con2').set('descr', 'Pressure');
model.result('pg10').feature('con2').set('levelmethod', 'levels');
model.result('pg10').feature('con2').set('levels', 'range(-2000,100,2000)');
model.result('pg10').feature('con2').set('resolution', 'normal');
model.result('pg11').label('v time');
model.result('pg11').set('xlabel', 'Time (s)');
model.result('pg11').set('xlabelactive', false);
model.result('pg11').feature('tblp1').set('table', 'tbl4');
model.result('pg11').feature('tblp1').set('legend', true);
model.result('pg11').feature('tblp1').set('legendmethod', 'manual');
model.result('pg11').feature('tblp1').set('legends', {'Average_v_Plane_Liner_(m/s)'});
model.result('pg11').feature('tblp2').set('table', 'tbl6');
model.result('pg11').feature('tblp2').set('legend', true);
model.result('pg11').feature('tblp2').set('legendmethod', 'manual');
model.result('pg11').feature('tblp2').set('legends', {'Average_v_Plane_Multiphysics_Interface_(m/s)'});
model.result('pg11').feature('tblp3').set('table', 'tbl16');
model.result('pg11').feature('tblp3').set('legend', true);
model.result('pg11').feature('tblp3').set('legendmethod', 'manual');
model.result('pg11').feature('tblp3').set('legends', {'Average_v_Plane_Acoustics_Domain_(m/s)'});
model.result('pg13').label('p time');
model.result('pg13').set('xlabel', 'Time (s)');
model.result('pg13').set('xlabelactive', false);
model.result('pg13').feature('tblp1').set('table', 'tbl3');
model.result('pg13').feature('tblp1').set('legend', true);
model.result('pg13').feature('tblp1').set('legendmethod', 'manual');
model.result('pg13').feature('tblp1').set('legends', {'Av_P_Plane_Liner_(Pa)'});
model.result('pg13').feature('tblp2').set('table', 'tbl5');
%model.result('pg13').feature('tblp2').set('xaxisdata', 1);
model.result('pg13').feature('tblp2').set('legend', true);
model.result('pg13').feature('tblp2').set('legendmethod', 'manual');
model.result('pg13').feature('tblp2').set('legends', {'Av_P_Plane_Multi_(Pa)'});
model.result('pg13').feature('ptgr1').set('expr', 'p2');
model.result('pg13').feature('ptgr1').set('descr', 'Pressure');
model.result('pg13').feature('ptgr1').set('legend', true);
model.result('pg13').feature('ptgr1').set('legendmethod', 'manual');
model.result('pg13').feature('ptgr1').set('legends', {'P_Tap_near' 'P_Tap_far' 'P_Inlet'});
model.result('pg13').feature('tblp3').set('table', 'tbl15');
model.result('pg13').feature('tblp3').set('legend', true);
model.result('pg13').feature('tblp3').set('legendmethod', 'manual');
model.result('pg13').feature('tblp3').set('legends', {'Av_P_Plane_Acoustic_(Pa)'});
model.result.export('anim1').set('plotgroup', 'pg3');
model.result.export('anim1').set('target', 'player');
model.result.export('anim1').set('framesel', 'all');
model.result.export('anim1').set('showframe', 171);
model.result.export('anim1').set('shownparameter', '0.011333');
model.result.export('anim1').set('frametime', 0.15);
model.result.export('anim1').set('title', 'on');
model.result.export('anim1').set('legend', 'on');
model.result.export('anim1').set('logo', 'on');
model.result.export('anim1').set('options', 'off');
model.result.export('anim1').set('fontsize', '9');
model.result.export('anim1').set('customcolor', [1 1 1]);
model.result.export('anim1').set('background', 'color');
model.result.export('anim1').set('axisorientation', 'on');
model.result.export('anim1').set('grid', 'on');
model.result.export('anim1').set('axes', 'on');
model.result.export('anim1').set('showgrid', 'on');
model.result.export('anim3').set('plotgroup', 'pg5');
model.result.export('anim3').set('target', 'player');
model.result.export('anim3').set('looplevelinput', 'manual');
model.result.export('anim3').set('looplevel', [157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323]);
model.result.export('anim3').set('framesel', 'all');
model.result.export('anim3').set('showframe', 154);
model.result.export('anim3').set('shownparameter', '0.0206');
model.result.export('anim3').set('frametime', '0.20');
model.result.export('anim3').set('title', 'on');
model.result.export('anim3').set('legend', 'on');
model.result.export('anim3').set('logo', 'on');
model.result.export('anim3').set('options', 'off');
model.result.export('anim3').set('fontsize', '9');
model.result.export('anim3').set('customcolor', [1 1 1]);
model.result.export('anim3').set('background', 'color');
model.result.export('anim3').set('axisorientation', 'on');
model.result.export('anim3').set('grid', 'on');
model.result.export('anim3').set('axes', 'on');
model.result.export('anim3').set('showgrid', 'on');
model.result.export('data1').set('looplevelinput', {'first'});
model.result.export('data1').set('expr', {'v'});
model.result.export('data1').set('unit', {'m/s'});
model.result.export('data1').set('descr', {'Velocity field, y component'});
model.result.export('data1').set('filename', 'Velocidades2.txt');
model.result.export('data1').set('struct', 'sectionwise');
model.result.export('img1').set('plotgroup', 'pg3');
model.result.export('img1').set('pngfilename', 'VelMag_t0.0021863s');
model.result.export('img1').set('printunit', 'mm');
model.result.export('img1').set('webunit', 'px');
model.result.export('img1').set('printheight', '90');
model.result.export('img1').set('webheight', '600');
model.result.export('img1').set('printwidth', '120');
model.result.export('img1').set('webwidth', '800');
model.result.export('img1').set('printlockratio', 'off');
model.result.export('img1').set('weblockratio', 'off');
model.result.export('img1').set('printresolution', '300');
model.result.export('img1').set('webresolution', '96');
model.result.export('img1').set('size', 'current');
model.result.export('img1').set('antialias', 'on');
model.result.export('img1').set('zoomextents', 'off');
model.result.export('img1').set('title', 'on');
model.result.export('img1').set('legend', 'on');
model.result.export('img1').set('logo', 'on');
model.result.export('img1').set('options', 'off');
model.result.export('img1').set('fontsize', '9');
model.result.export('img1').set('customcolor', [1 1 1]);
model.result.export('img1').set('background', 'color');
model.result.export('img1').set('axes', 'on');
model.result.export('img1').set('qualitylevel', '92');
model.result.export('img1').set('qualityactive', 'off');
model.result.export('img1').set('imagetype', 'png');
model.result.export('plot1').set('plotgroup', 'pg3');
model.result.export('plot1').set('filename', 'Velocity');
model.result.export('plot1').set('header', false);
model.result.export('plot1').set('fullprec', false);
model.result.export('plot2').set('filename', 'PressureField.txt');
model.result.export('plot3').set('plotgroup', 'pg3');
model.result.export('plot3').set('exporttype', 'vtu');
model.result.export('plot3').set('filename', 'VelocityTestVTK.vtu');

%% Save
mphsave(model,['../' file '.mph']);

model.result().table('tbl3').save([pasta file '\Data\PLine.txt']);
model.result().table('tbl4').save([pasta file '\Data\VLine.txt']);
model.result().table('tbl5').save([pasta file '\Data\PrePlane1.txt']);
model.result().table('tbl6').save([pasta file '\Data\VelPlane1.txt']);
model.result().table('tbl9').save([pasta file '\Data\PreMic2.txt']);
model.result().table('tbl10').save([pasta file '\Data\AcMic2.txt']);
model.result().table('tbl11').save([pasta file '\Data\PreMic1.txt']);
model.result().table('tbl12').save([pasta file '\Data\AcMic1.txt']);
model.result().table('tbl13').save([pasta file '\Data\Pin.txt']);
model.result().table('tbl14').save([pasta file '\Data\Pout.txt']);
model.result().table('tbl15').save([pasta file '\Data\PrePlane2.txt']);
model.result().table('tbl16').save([pasta file '\Data\VelPlane2.txt']);

toc