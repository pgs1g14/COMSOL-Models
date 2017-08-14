function LinerAxi_PT_155dB_parametric_v5_2(d_h,POA,h_h,h_r,cycles,freq)
%
% LinerAxi_PT_155dB_parametric_v5_2.m
%
% Model exported on Aug 9 2017, 17:45 by COMSOL 5.2.1.229.
% Edited to include parametric analysis on MATLAB on 11/08/2017
%   freq=600;
%   cycles=5;
%% Constants
tic;
%Model physical constants
%------------------------------------------------
w=2*pi()*freq; %Angular frequency
TC=20; %Temperature of the measurements to run the simulation ºC
T= 273.15+TC; %Temperature of the fluid [K]
ATMp= 101450; %Atmosferic pressure [Pa]
R=287.058; %Specific air constant [J/kg.K]
rho=ATMp/(R*T); %Specific gas density [kg/m^3]
mi= 1.8205e-5;%dynamic viscosity [kg/m.s]
ni=mi/rho; %kinematic viscosity [m^2/s]
c0= 331.29*sqrt(1+((T-273.15)/273.15)); %Fluid sound velocity [m/s]
SPL=155; %Sound pressure level dB
k0=w/c0; %wave number

%Liner geometry
% d_h=0.00099; %Hole diameter 1,0614 mm - on the holder, Hole diameter 1,5022 mm - flanged
% POA=0.0366; %percentage of open area of the sample in the holder
% h_h=0.000635; %Liner slit thickness [m]
% h_r=0.019; %Liner cavity depth - resonator height [m]

h_in=0.2082; %impedance tube length [m]
yPlane1=0.07; %Acoustics
yPlane2=40*h_h; %Multiphysics
x1=0.099+0.02; %Mic 1 position
x2=0.099; %Mic 2 position
elementHole=2;
elementNear=15;
elementFar=24;

%-----
ppl=10; % Points per lambda to have good description on time domain
%-----
%% Start Model
import com.comsol.model.*
import com.comsol.model.util.*

    model = ModelUtil.create('Model');

    ModelUtil.showProgress(true); %Show the progress bar
    ModelUtil.setServerBusyHandler(ServerBusyHandler(10));

    model.modelPath('C:\Local\Pablo\COMSOL\Circular Axi model\LinerAxi_PT_155dB_1500Hz');

    model.label('LinerAxi_PT_155dB_parametric_v5.2.mph');

    model.comments(['LinerAxi PT 155dB 500Hz v5.2\n\nLinerAxi PT 155dB 600Hz v5.2\n\nLiner Axisimetrical PT thickness 0.000635 m 155dB 500Hz. Impedance tube model with a sample in a holder with POA 3.66%, d=0.00099 m, h=0.019m nominal values.']);
%% Parameters
    model.param.set('h_r', [num2str(h_r) '[m]'], 'Resonator height');
    model.param.set('h_h', [num2str(h_h) '[m]'], 'Hole thickness');
    model.param.set('d_h', [num2str(d_h) '[m]'], 'Hole diameter');
    model.param.set('h_in', [num2str(h_in) '[m]'], 'Inlet height');
    model.param.set('y0', [num2str(h_r+h_h+h_in) '[m]'], 'Inlet location y-coordinate');
    model.param.set('rho0', [num2str(rho) '[kg/m^3]'], 'Density');
    model.param.set('c0', [num2str(c0) '[m/s]'], 'Speed of sound');
    model.param.set('yPlane1', [num2str(yPlane1) '[m]'], 'Distance to plane 1');
    model.param.set('yPlane2', [num2str(yPlane2) '[m]'], 'Distance to plane 2');
    model.param.set('POA', num2str(POA) , 'Percentage of open area');model.param.set('beta', '1.2', 'Nonlinear coefficent');
    model.param.set('L0', num2str(SPL), 'Incident wave amplitude dB');
    model.param.set('f0', [num2str(freq) '[Hz]'], 'Driving frequency');
    model.param.set('cycles', num2str(cycles), 'Cycles to extract from the beggining');
    model.param.set('k0', [num2str(k0) '[1/m]'], 'Wave number at f0');
    model.param.set('ppl', num2str(ppl), 'Points per lambda');
    model.param.set('d_tube', 'sqrt(d_h^2/POA) [m]', 'Tube diameter considering POA=3.66%');
    model.param.set('p0', '10^(L0/20)*2*10^-5[Pa]', 'Incident wave amplitude Pa');
    model.param.set('Tstart', 'cycles/f0', 'Start time for a number of cycles');
    model.param.set('lambda0', 'c0/f0[m]', 'Wavelength at f0');
    model.param.set('dvisc', '0.125[m]*sqrt(4*pi*1.51*10^-5[Hz]/f0)', 'Viscous boundary layer thickness at f0');
    model.param.set('omega0', '2*pi*f0', 'Angular frequency');
    model.param.set('T0', '1/f0', 'Period');
    model.param.set('dt_sol', 'T0/ppl', 'Solver time step (resolve f0 with T0/ppl)');
    model.param.set('lmesh', 'lambda0[m]/ppl', 'Maximum element size acoustic mesh');
    model.param.set('Tend', 'Tstart+100*dt_sol', 'End time for post processing');
    model.param.set('d_r', 'sqrt(d_h^2/POA) [m]', 'Resonator diameter');
   
%% Geometry
    model.modelNode.create('comp1');

    model.geom.create('geom1', 2);

    model.func.create('rect1', 'Rectangle');
    model.func.create('an1', 'Analytic');
    model.func('rect1').model('comp1');
    model.func('rect1').set('lower', '-1.5');
    model.func('rect1').set('upper', '1.5');
    model.func('an1').model('comp1');
    model.func('an1').set('expr', 'rect(p)');
    model.func('an1').set('args', {'p'});
    model.func('an1').set('plotargs', {'p' '0' '0.005'});
    model.func('an1').set('funcname', 'rect');

    model.geom('geom1').axisymmetric(true);

    model.mesh.create('mesh1', 'geom1');

    model.geom('geom1').create('r1', 'Rectangle');
    model.geom('geom1').feature('r1').set('size', {'d_r/2' 'h_r'});
    model.geom('geom1').create('r2', 'Rectangle');
    model.geom('geom1').feature('r2').set('size', {'d_h/2' 'h_h'});
    model.geom('geom1').feature('r2').set('pos', {'0' 'h_r'});
    model.geom('geom1').create('r3', 'Rectangle');
    model.geom('geom1').feature('r3').set('size', {'d_r/2' 'h_in'});
    model.geom('geom1').feature('r3').set('pos', {'0' 'h_r+h_h'});
    model.geom('geom1').create('pol1', 'Polygon');
    model.geom('geom1').feature('pol1').set('x', '0 d_r/2');
    model.geom('geom1').feature('pol1').set('y', 'yPlane1 yPlane1');
    model.geom('geom1').create('pol2', 'Polygon');
    model.geom('geom1').feature('pol2').set('x', '0 d_r/2');
    model.geom('geom1').feature('pol2').set('y', 'yPlane2+h_r+h_h yPlane2+h_r+h_h');
    model.geom('geom1').create('pt1', 'Point');
    model.geom('geom1').feature('pt1').setIndex('p', 'd_r/2', 0, 0);
    model.geom('geom1').feature('pt1').setIndex('p', '0.119+h_h+h_r', 1, 0);
    model.geom('geom1').create('pt2', 'Point');
    model.geom('geom1').feature('pt2').setIndex('p', 'd_r/2', 0, 0);
    model.geom('geom1').feature('pt2').setIndex('p', '0.099+h_h+h_r', 1, 0);
    model.geom('geom1').run;
    model.geom('geom1').run('fin');
%% Variables
    model.variable.create('var1');
    model.variable('var1').model('comp1');
    model.variable('var1').set('Pin', 'p0*sin(omega0*t)', 'Incident pressure');
    model.variable('var1').set('rhoL', 'rho0+p/c0^2', 'Linear density relation for air');
    model.variable('var1').set('rhoNL', 'rho0+p/c0^2-1/(rho0*c0^4)*(beta-1)*p^2', 'Nonlinear density relation for air');

    model.view('view1').tag('view2');
    model.view.create('view3', 3);
    model.view.create('view4', 2);
    model.view.create('view5', 2);
    model.view.create('view6', 2);
    model.view.create('view7', 2);

    model.material.create('mat1', 'Common', 'comp1');
    model.material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
    model.material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
    model.material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
    model.material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
    model.material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
    model.material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
    model.material('mat1').propertyGroup.create('idealGas', 'Ideal gas');

    model.cpl.create('aveop1', 'Average', 'geom1');
    model.cpl.create('aveop2', 'Average', 'geom1');
    model.cpl.create('aveop3', 'Average', 'geom1');
    model.cpl('aveop1').selection.geom('geom1', 1);
    model.cpl('aveop1').selection.set([6 14]);
    model.cpl('aveop2').selection.geom('geom1', 1);
    model.cpl('aveop2').selection.set([8]);
    model.cpl('aveop3').selection.geom('geom1', 1);
    model.cpl('aveop3').selection.set([10]);
%% Physics
    model.physics.create('actd', 'TransientPressureAcoustics', 'geom1');
    model.physics('actd').field('pressure').field('p2');
    model.physics('actd').selection.set([4 5]);
    model.physics('actd').create('nvel1', 'NormalVelocity', 1);
    model.physics('actd').feature('nvel1').selection.set([8]);
    model.physics('actd').create('pwr1', 'PlaneWaveRadiation', 1);
    model.physics('actd').feature('pwr1').selection.set([11]);
    model.physics('actd').feature('pwr1').create('ipf1', 'IncidentPressureField', 1);
    model.physics('actd').create('ssb1', 'SoundSoft', 1);
    model.physics('actd').feature('ssb1').selection.set([17 18 19 20]);
    model.physics('actd').create('pr1', 'Pressure', 1);
    model.physics('actd').feature('pr1').selection.set([8]);
    model.physics.create('spf', 'LaminarFlow', 'geom1');
    model.physics('spf').selection.set([1 2 3]);
    model.physics('spf').create('bs1', 'BoundaryStress', 1);
    model.physics('spf').feature('bs1').selection.set([8]);
    model.physics('spf').create('wall2', 'Wall', 1);
    model.physics('spf').feature('wall2').selection.set([15 16]);
        model.physics('actd').prop('EquationForm').set('form', 'Transient');
    model.physics('actd').feature('tpam1').set('FluidModel', 'IdealGas');
    model.physics('actd').feature('nvel1').set('nvel', 'w');
    model.physics('actd').feature('pwr1').feature('ipf1').set('p', 'Pin');
    model.physics('actd').feature('ssb1').active(false);
    model.physics('actd').feature('pr1').set('p0', 'p');
    model.physics('actd').feature('pr1').active(false);
    model.physics('spf').prop('PhysicalModelProperty').set('Compressibility', 'CompressibleMALT03');
    model.physics('spf').feature('bs1').set('BoundaryCondition', 'NormalStress');
    model.physics('spf').feature('bs1').set('f0', 'p2');
    model.physics('spf').feature('wall2').set('BoundaryCondition', 'Slip');
%% Tables
    model.result.table.create('tbl1', 'Table');
    model.result.table.create('tbl2', 'Table');
    model.result.table.create('tbl3', 'Table');
    model.result.table.create('tbl17', 'Table');
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
    model.result.table.create('tbl18', 'Table');
    model.result.table.create('tbl19', 'Table');

    model.view('view2').label('View 2');
    model.view('view2').axis.set('abstractviewxscale', '3.0793496989645064E-5');
    model.view('view2').axis.set('abstractviewtratio', '-0.8695948719978333');
    model.view('view2').axis.set('abstractviewlratio', '-1.7451910972595215');
    model.view('view2').axis.set('abstractviewyscale', '3.0793496989645064E-5');
    model.view('view2').axis.set('abstractviewrratio', '1.5749908685684204');
    model.view('view2').axis.set('abstractviewbratio', '0.03957964852452278');
    model.view('view2').axis.set('ymax', '0.029710859060287476');
    model.view('view2').axis.set('xmax', '0.006662531290203333');
    model.view('view2').axis.set('ymin', '0.00901762954890728');
    model.view('view2').axis.set('xmin', '-0.004515507258474827');
    model.view('view4').axis.set('abstractviewxscale', '3.696437051985413E-4');
    model.view('view4').axis.set('abstractviewtratio', '0.04999997466802597');
    model.view('view4').axis.set('abstractviewlratio', '-24.643888473510742');
    model.view('view4').axis.set('abstractviewyscale', '3.6964379251003265E-4');
    model.view('view4').axis.set('abstractviewrratio', '24.64388656616211');
    model.view('view4').axis.set('abstractviewbratio', '-0.04999997466802597');
    model.view('view4').axis.set('ymax', '0.23922674357891083');
    model.view('view4').axis.set('xmax', '0.06635098904371262');
    model.view('view4').axis.set('ymin', '-0.011391744017601013');
    model.view('view4').axis.set('xmin', '-0.06376359611749649');
    model.view('view5').axis.set('ymax', '1.4776785373687744');
    model.view('view5').axis.set('ymin', '-1.4776785373687744');
    model.view('view6').axis.set('ymax', '810.0000610351562');
    model.view('view6').axis.set('xmax', '7.7674736976623535');
    model.view('view6').axis.set('ymin', '790');
    model.view('view6').axis.set('xmin', '-5.7674736976623535');
    model.view('view7').axis.set('abstractviewxscale', '3.718374064192176E-4');
    model.view('view7').axis.set('abstractviewtratio', '0.04999997466802597');
    model.view('view7').axis.set('abstractviewlratio', '-19.475807189941406');
    model.view('view7').axis.set('abstractviewyscale', '3.718375228345394E-4');
    model.view('view7').axis.set('abstractviewrratio', '19.475807189941406');
    model.view('view7').axis.set('abstractviewbratio', '-0.04999997466802597');
    model.view('view7').axis.set('ymax', '0.23922671377658844');
    model.view('view7').axis.set('xmax', '0.05297910049557686');
    model.view('view7').axis.set('ymin', '-0.011391714215278625');
    model.view('view7').axis.set('xmin', '-0.050391700118780136');

    model.material('mat1').label('Air');
    model.material('mat1').set('family', 'air');
    model.material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
    model.material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
    model.material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
    model.material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
    model.material('mat1').propertyGroup('def').func('rho').set('dermethod', 'manual');
    model.material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/8.314/T');
    model.material('mat1').propertyGroup('def').func('rho').set('argders', {'pA' 'd(pA*0.02897/8.314/T,pA)'; 'T' 'd(pA*0.02897/8.314/T,T)'});
    model.material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
    model.material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA' '0' '1'; 'T' '0' '1'});
    model.material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
    model.material('mat1').propertyGroup('def').func('k').set('arg', 'T');
    model.material('mat1').propertyGroup('def').func('cs').set('dermethod', 'manual');
    model.material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*287*T)');
    model.material('mat1').propertyGroup('def').func('cs').set('argders', {'T' 'd(sqrt(1.4*287*T),T)'});
    model.material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
    model.material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T' '0' '1'});
    model.material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T[1/K])[Pa*s]');
    model.material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
    model.material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
    model.material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T[1/K])[J/(kg*K)]');
    model.material('mat1').propertyGroup('def').set('density', 'rhoNL');
    model.material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]'});
    model.material('mat1').propertyGroup('def').set('soundspeed', 'c0');
    model.material('mat1').propertyGroup('def').addInput('temperature');
    model.material('mat1').propertyGroup('def').addInput('pressure');
    model.material('mat1').propertyGroup('RefractiveIndex').set('n', '');
    model.material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
    model.material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.material('mat1').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
    model.material('mat1').propertyGroup('idealGas').set('Rs', '287');
    model.material('mat1').propertyGroup('idealGas').addInput('temperature');
    model.material('mat1').propertyGroup('idealGas').addInput('pressure');

    model.cpl('aveop1').label('Average 1 - Liner');
    model.cpl('aveop1').set('axisym', true);
    model.cpl('aveop2').label('Average 2 - Plane multiphysics');
    model.cpl('aveop2').set('axisym', true);
    model.cpl('aveop3').label('Average 3 - Plane acoustics');
    model.cpl('aveop3').set('axisym', true);    
%% Mesh
    model.mesh('mesh1').create('map1', 'Map');
    model.mesh('mesh1').create('ftri1', 'FreeTri');
    model.mesh('mesh1').feature('map1').selection.geom('geom1', 2);
    model.mesh('mesh1').feature('map1').selection.set([4 5]);
    model.mesh('mesh1').feature('map1').create('size1', 'Size');
    model.mesh('mesh1').feature('ftri1').create('size1', 'Size');
    model.mesh('mesh1').feature('ftri1').create('size2', 'Size');
    model.mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 2);
    model.mesh('mesh1').feature('ftri1').feature('size1').selection.set([1 3]);
    model.mesh('mesh1').feature('ftri1').feature('size2').selection.geom('geom1', 2);
    model.mesh('mesh1').feature('ftri1').feature('size2').selection.set([2]);

    model.mesh('mesh1').feature('map1').set('adjustedgdistr', true);
    model.mesh('mesh1').feature('map1').feature('size1').set('hauto', 2);
    model.mesh('mesh1').feature('map1').feature('size1').set('custom', 'on');
    model.mesh('mesh1').feature('map1').feature('size1').set('hmaxactive', true);
    model.mesh('mesh1').feature('map1').feature('size1').set('hmax', 'lmesh');
    model.mesh('mesh1').feature('map1').create('dis1', 'Distribution');
    model.mesh('mesh1').feature('map1').feature('dis1').selection.set([8]);
    model.mesh('mesh1').feature('map1').feature('dis1').set('numelem', 3);
  
    model.mesh('mesh1').feature('ftri1').feature('size1').set('custom', 'on');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hminactive', true);
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', true);
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hmin', 'dvisc');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', 'dvisc*24');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hgradactive', true);
    model.mesh('mesh1').feature('ftri1').feature('size1').set('table', 'cfd');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hgrad', '1.3');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hnarrowactive', false);
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hmin', 'dvisc');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', 'dvisc*24');
    
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hcurveactive', false);
    model.mesh('mesh1').feature('ftri1').feature('size2').set('custom', 'on');
    model.mesh('mesh1').feature('ftri1').feature('size2').set('hmaxactive', true);
    model.mesh('mesh1').feature('ftri1').feature('size2').set('hgrad', '1.05');
    model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc');
    model.mesh('mesh1').feature('ftri1').feature('size2').set('hgradactive', true);
    model.mesh('mesh1').run;

    model.frame('material1').sorder(1);

    model.result.table('tbl1').label('Point Evaluation 1 (Pin)');
    model.result.table('tbl1').comments('Point Evaluation 1 (Pin)');
    model.result.table('tbl2').label('Point Evaluation 1 (Pout)');
    model.result.table('tbl2').comments('Point Evaluation 1 (Pout)');
    model.result.table('tbl3').label('Pressure Liner');
    model.result.table('tbl3').comments('(aveop1(p))');
    model.result.table('tbl17').label('Velocity Liner');
    model.result.table('tbl17').comments('V Liner (aveop1(w))');
    model.result.table('tbl5').label('Ave P Plane Multiphysics');
    model.result.table('tbl6').label('Ave v Plane Mutiphysics');
    model.result.table('tbl9').label('Pressure Mic 2');
    model.result.table('tbl10').label('Aceleration Mic 2');
    model.result.table('tbl11').label('Pressure Mic 1');
    model.result.table('tbl12').label('Aceleration Mic 1');
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
    model.result.table('tbl18').label('Integral of Pressure');
    model.result.table('tbl18').comments('Line Integration 1 (p)');
    model.result.table('tbl19').label('Integral of Velocity');
    model.result.table('tbl19').comments('Integrate Velocity Liner ()');

    model.study.create('std1');
    model.study('std1').create('time', 'Transient');

    model.sol.create('sol1');
    model.sol('sol1').study('std1');
    model.sol('sol1').attach('std1');
    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').create('t1', 'Time');
    model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('t1').create('d1', 'Direct');
    model.sol('sol1').feature('t1').feature.remove('fcDef');

    model.result.dataset.create('rev1', 'Revolve2D');
    model.result.dataset.create('rev2', 'Revolve2D');
    model.result.dataset.create('rev3', 'Revolve1D');
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
    model.result.numerical.create('int1', 'IntLine');
    model.result.numerical.create('int2', 'IntLine');
    model.result.numerical('pev1').selection.set([16]);
    model.result.numerical('pev1').set('probetag', 'none');
    model.result.numerical('pev2').selection.set([16]);
    model.result.numerical('pev2').set('probetag', 'none');
    model.result.numerical('gev3').set('probetag', 'none');
    model.result.numerical('gev4').set('probetag', 'none');
    model.result.numerical('gev5').set('probetag', 'none');
    model.result.numerical('gev6').set('probetag', 'none');
    model.result.numerical('pev3').selection.set([14]);
    model.result.numerical('pev3').set('probetag', 'none');
    model.result.numerical('pev4').selection.set([14]);
    model.result.numerical('pev4').set('probetag', 'none');
    model.result.numerical('pev5').selection.set([15]);
    model.result.numerical('pev5').set('probetag', 'none');
    model.result.numerical('pev6').selection.set([15]);
    model.result.numerical('pev6').set('probetag', 'none');
    model.result.numerical('pev7').selection.set([16]);
    model.result.numerical('pev7').set('probetag', 'none');
    model.result.numerical('pev8').selection.set([16]);
    model.result.numerical('pev8').set('probetag', 'none');
    model.result.numerical('gev7').set('probetag', 'none');
    model.result.numerical('gev8').set('probetag', 'none');
    model.result.numerical('int1').selection.set([6 14]);
    model.result.numerical('int1').set('probetag', 'none');
    model.result.numerical('int2').selection.set([6 14]);
    model.result.numerical('int2').set('probetag', 'none');
    
    model.result.numerical('pev1').label('Point Evaluation: RMS Incident');
    model.result.numerical('pev1').set('descr', {'Incident pressure'});
    model.result.numerical('pev1').set('table', 'tbl1');
    model.result.numerical('pev1').set('unit', {''});
    model.result.numerical('pev1').set('dataseries', 'rms');
    model.result.numerical('pev1').set('expr', {'Pin'});
    model.result.numerical('pev1').set('unit', {''});
    model.result.numerical('pev2').label('Point Evaluation: RMS Reflected');
    model.result.numerical('pev2').set('descr', {''});
    model.result.numerical('pev2').set('table', 'tbl2');
    model.result.numerical('pev2').set('unit', {'Pa'});
    model.result.numerical('pev2').set('dataseries', 'rms');
    model.result.numerical('pev2').set('expr', {'Pin-p2'});
    model.result.numerical('gev3').label('P Liner');
    model.result.numerical('gev3').set('descr', {''});
    model.result.numerical('gev3').set('table', 'tbl3');
    model.result.numerical('gev3').set('unit', {'Pa'});
    model.result.numerical('gev3').set('expr', {'aveop1(p)'});
    model.result.numerical('gev4').label('V Liner');
    model.result.numerical('gev4').set('descr', {''});
    model.result.numerical('gev4').set('table', 'tbl17');
    model.result.numerical('gev4').set('unit', {'m/s'});
    model.result.numerical('gev4').set('expr', {'aveop1(w)'});
    model.result.numerical('gev5').label('P plane 1');
    model.result.numerical('gev5').set('descr', {''});
    model.result.numerical('gev5').set('table', 'tbl5');
    model.result.numerical('gev5').set('unit', {'Pa'});
    model.result.numerical('gev5').set('expr', {'aveop2(p)'});
    model.result.numerical('gev6').label('V plane 1');
    model.result.numerical('gev6').set('descr', {''});
    model.result.numerical('gev6').set('table', 'tbl6');
    model.result.numerical('gev6').set('unit', {'m/s'});
    model.result.numerical('gev6').set('expr', {'aveop2(w)'});
    model.result.numerical('pev3').label('P Mic 2');
    model.result.numerical('pev3').set('descr', {'Pressure'});
    model.result.numerical('pev3').set('table', 'tbl9');
    model.result.numerical('pev3').set('unit', {'Pa'});
    model.result.numerical('pev3').set('expr', {'p2'});
    model.result.numerical('pev4').label('V Mic 2');
    model.result.numerical('pev4').set('descr', {'Local acceleration, z component'});
    model.result.numerical('pev4').set('table', 'tbl10');
    model.result.numerical('pev4').set('unit', {'m/s^2'});
    model.result.numerical('pev4').set('expr', {'actd.az'});
    model.result.numerical('pev5').label('P Mic 1');
    model.result.numerical('pev5').set('descr', {'Pressure'});
    model.result.numerical('pev5').set('table', 'tbl11');
    model.result.numerical('pev5').set('unit', {'Pa'});
    model.result.numerical('pev5').set('expr', {'p2'});
    model.result.numerical('pev6').label('V Mic 1');
    model.result.numerical('pev6').set('descr', {'Local acceleration, z component'});
    model.result.numerical('pev6').set('table', 'tbl12');
    model.result.numerical('pev6').set('unit', {'m/s^2'});
    model.result.numerical('pev6').set('expr', {'actd.az'});
    model.result.numerical('pev7').label('Pressure incident Inlet');
    model.result.numerical('pev7').set('descr', {'Incident pressure'});
    model.result.numerical('pev7').set('table', 'tbl13');
    model.result.numerical('pev7').set('unit', {''});
    model.result.numerical('pev7').set('expr', {'Pin'});
    model.result.numerical('pev7').set('unit', {''});
    model.result.numerical('pev8').label('Pressure reflected Inlet');
    model.result.numerical('pev8').set('descr', {''});
    model.result.numerical('pev8').set('table', 'tbl14');
    model.result.numerical('pev8').set('unit', {'Pa'});
    model.result.numerical('pev8').set('expr', {'Pin-p2'});
    model.result.numerical('gev7').label('P plane ac');
    model.result.numerical('gev7').set('descr', {'Average Pressure Plane Acoustic'});
    model.result.numerical('gev7').set('table', 'tbl15');
    model.result.numerical('gev7').set('unit', {'Pa'});
    model.result.numerical('gev7').set('expr', {'aveop3(p2)'});
    model.result.numerical('gev8').label('V plane ac');
    model.result.numerical('gev8').set('descr', {'Average Velocity Plane Acoustic'});
    model.result.numerical('gev8').set('table', 'tbl16');
    model.result.numerical('gev8').set('unit', {'m/s^2'});
    model.result.numerical('gev8').set('expr', {'aveop3(actd.az)'});
    model.result.numerical('int1').label('Integrate Pressure Liner');
    model.result.numerical('int1').set('descr', {'Pressure'});
    model.result.numerical('int1').set('table', 'tbl18');
    model.result.numerical('int1').set('unit', {''});
    model.result.numerical('int1').set('intsurface', true);
    model.result.numerical('int1').set('expr', {'exp(-i*omega0*t)*p'});
    model.result.numerical('int1').set('unit', {''});
    model.result.numerical('int2').label('Integrate Velocity Liner');
    model.result.numerical('int2').set('descr', {''});
    model.result.numerical('int2').set('table', 'tbl19');
    model.result.numerical('int2').set('unit', {''});
    model.result.numerical('int2').set('intsurface', true);
    model.result.numerical('int2').set('expr', {'exp(-i*omega0*t)*w'});
    model.result.numerical('int2').set('unit', {''});
   
    model.result.create('pg1', 'PlotGroup2D');
    model.result.create('pg2', 'PlotGroup1D');
    model.result.create('pg3', 'PlotGroup2D');
    model.result.create('pg6', 'PlotGroup1D');
    model.result.create('pg8', 'PlotGroup2D');
    model.result.create('pg14', 'PlotGroup3D');
    model.result.create('pg16', 'PlotGroup3D');
    model.result.create('pg15', 'PlotGroup1D');
    model.result.create('pg17', 'PlotGroup1D');
    model.result('pg1').create('surf1', 'Surface');
    model.result('pg1').create('surf2', 'Surface');
    model.result('pg2').create('ptgr2', 'PointGraph');
    model.result('pg17').create('ptgr3', 'PointGraph');
    model.result('pg2').feature('ptgr2').set('data', 'dset1');
    model.result('pg2').feature('ptgr2').selection.set([1 2 3 4]);
    model.result('pg17').feature('ptgr3').set('data', 'dset1');
    model.result('pg17').feature('ptgr3').selection.set([5 6 14 15]);
    model.result('pg3').create('surf1', 'Surface');
    model.result('pg3').create('surf2', 'Surface');
    model.result('pg6').create('ptgr1', 'PointGraph');
    model.result('pg6').feature('ptgr1').selection.set([2 3]);
    model.result('pg8').create('surf1', 'Surface');
    model.result('pg14').create('vol1', 'Volume');
    model.result('pg14').create('vol2', 'Volume');
    model.result('pg16').create('vol1', 'Volume');
    model.result('pg15').create('glob1', 'Global');
    model.result('pg15').create('glob2', 'Global');
    model.result.export.create('anim1', 'Animation');
    model.result.export.create('anim2', 'Animation');
    model.result.export.create('anim3', 'Animation');

    model.study('std1').feature('time').set('probefreq', 'tout');
    model.study('std1').feature('time').set('tlist', 'range(0,dt_sol,Tend)');

    model.sol('sol1').attach('std1');
    model.sol('sol1').feature('v1').set('clist', {'range(0,dt_sol,Tend)'});
    model.sol('sol1').feature('t1').set('probefreq', 'tout');
    model.sol('sol1').feature('t1').set('estrat', 'exclude');
    model.sol('sol1').feature('t1').set('tlist', 'range(0,dt_sol,Tend)');
    model.sol('sol1').feature('t1').set('timemethod', 'genalpha');
    model.sol('sol1').feature('t1').set('atolglobal', '5.0E-4');
    model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
    model.sol('sol1').feature('t1').feature('fc1').set('maxiter', '6');
%% Ready to run

    model.sol('sol1').runAll;
%% Set Results
    model.result.dataset('rev1').label('Revolution 2D');
    model.result.dataset('rev1').set('revangle', '225');
    model.result.dataset('rev1').set('startangle', '-90');
    model.result.dataset('rev1').set('genpoints', {'0' '0'; '0' '1'});
    model.result.dataset('rev2').label('Revolution 2D 1');
    model.result.dataset('rev2').set('revangle', '225');
    model.result.dataset('rev2').set('startangle', '-90');
    model.result.dataset('rev2').set('genpoints', {'0' '0'; '0' '1'});
  
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
    model.result.numerical('int1').setResult;
    model.result.numerical('int2').setResult;
    model.result('pg1').label('Presure');
    model.result('pg1').set('solrepresentation', 'solnum');
    model.result('pg1').set('view', 'view7');
    model.result('pg1').feature('surf1').set('rangedatamax', '2000');
    model.result('pg1').feature('surf1').set('rangedataactive', 'on');
    model.result('pg1').feature('surf1').set('descr', 'Pressure');
    model.result('pg1').feature('surf1').set('rangecolormax', '2000');
    model.result('pg1').feature('surf1').set('rangedatamin', '-2000');
    model.result('pg1').feature('surf1').set('rangecolormin', '-2000');
    model.result('pg1').feature('surf1').set('rangecoloractive', 'on');
    model.result('pg1').feature('surf1').set('expr', 'p2');
    model.result('pg1').feature('surf1').set('resolution', 'normal');
    model.result('pg1').feature('surf2').set('rangedatamax', '2000');
    model.result('pg1').feature('surf2').set('rangedataactive', 'on');
    model.result('pg1').feature('surf2').set('descr', 'Pressure');
    model.result('pg1').feature('surf2').set('rangecolormax', '2000');
    model.result('pg1').feature('surf2').set('rangedatamin', '-2000');
    model.result('pg1').feature('surf2').set('rangecolormin', '-2000');
    model.result('pg1').feature('surf2').set('rangecoloractive', 'on');
    model.result('pg1').feature('surf2').set('expr', 'p');
    model.result('pg1').feature('surf2').set('resolution', 'normal');
    model.result('pg2').label('Pressure at some positions');
    model.result('pg2').set('data', 'none');
    model.result('pg2').set('xlabel', 'Time (s)');
    model.result('pg2').set('xlabelactive', false);
    model.result('pg2').feature('ptgr2').label('Laminar Flow domain');
    model.result('pg2').feature('ptgr2').set('descr', 'Pressure');
    model.result('pg2').feature('ptgr2').set('legends', {'Back cavity' 'Hole inside edge' 'Hole outside edge' 'Multiphysics boundary'});
    model.result('pg2').feature('ptgr2').set('legend', true);
    model.result('pg2').feature('ptgr2').set('legendmethod', 'manual');
    model.result('pg2').feature('ptgr2').set('expr', 'p');
    
    model.result('pg17').label('Pressure Acoustics domain');
    model.result('pg17').set('data', 'none');
    model.result('pg17').set('xlabel', 'Time (s)');
    model.result('pg17').set('xlabelactive', false);
    model.result('pg17').feature('ptgr3').label('Pressure acoustics domain');
    model.result('pg17').feature('ptgr3').set('descr', 'Pressure');
    model.result('pg17').feature('ptgr3').set('legends', {'Acoustics Plane' 'Inlet' 'Mic 2 near' 'Mic 1 far' ''});
    model.result('pg17').feature('ptgr3').set('legend', true);
    model.result('pg17').feature('ptgr3').set('linewidth', '2');
    model.result('pg17').feature('ptgr3').set('legendmethod', 'manual');
    model.result('pg17').feature('ptgr3').set('expr', 'comp1.actd.p_t');
    model.result('pg3').label('Velocity');
    model.result('pg3').set('axisprecision', '2');
    model.result('pg3').set('axisactive', 'on');
    model.result('pg3').set('axistrailingzeros', true);
    model.result('pg3').set('axisnotation', 'scientific');
    model.result('pg3').set('axiscommonexp', false);
    model.result('pg3').set('looplevel', {'11'});
    model.result('pg3').feature('surf1').set('descractive', true);
    model.result('pg3').feature('surf1').set('rangedatamax', '30');
    model.result('pg3').feature('surf1').set('rangedataactive', 'on');
    model.result('pg3').feature('surf1').set('descr', 'Velocity (CFD)');
    model.result('pg3').feature('surf1').set('rangecolormax', '30');
    model.result('pg3').feature('surf1').set('unit', 'm/s');
    model.result('pg3').feature('surf1').set('rangecoloractive', 'on');
    model.result('pg3').feature('surf1').set('expr', 'abs(w)');
    model.result('pg3').feature('surf1').set('resolution', 'normal');
    model.result('pg3').feature('surf2').set('descractive', true);
    model.result('pg3').feature('surf2').set('rangedatamax', '30');
    model.result('pg3').feature('surf2').set('rangedataactive', 'on');
    model.result('pg3').feature('surf2').set('descr', 'Velocity (Acoustics)');
    model.result('pg3').feature('surf2').set('rangecolormax', '30');
    model.result('pg3').feature('surf2').set('unit', 'm/s^2');
    model.result('pg3').feature('surf2').set('rangecoloractive', 'on');
    model.result('pg3').feature('surf2').set('expr', 'actd.a_inst/2/pi/omega0');
    model.result('pg3').feature('surf2').set('unit', 'm/s^2');
    model.result('pg3').feature('surf2').set('resolution', 'normal');
    model.result('pg6').label('Velocity on the resonator');
    model.result('pg6').set('ylabel', 'Velocity field, z component (m/s)');
    model.result('pg6').set('xlabel', 'Time (s)');
    model.result('pg6').set('ylabelactive', false);
    model.result('pg6').set('xlabelactive', false);
    model.result('pg6').feature('ptgr1').label('Points on the hole');
    model.result('pg6').feature('ptgr1').set('descr', 'Velocity field, z component');
    model.result('pg6').feature('ptgr1').set('legends', {'Back plate' 'Hole inner' 'Hole outer'});
    model.result('pg6').feature('ptgr1').set('unit', 'm/s');
    model.result('pg6').feature('ptgr1').set('linemarker', 'cycle');
    model.result('pg6').feature('ptgr1').set('markerpos', 'datapoints');
    model.result('pg6').feature('ptgr1').set('legend', true);
    model.result('pg6').feature('ptgr1').set('legendmethod', 'manual');
    model.result('pg6').feature('ptgr1').set('expr', 'w');
    model.result('pg8').label('Vorticity');
    model.result('pg8').set('looplevel', {'11'});
    model.result('pg8').feature('surf1').label('Laminar flow domain');
    model.result('pg8').feature('surf1').set('rangedatamax', '60');
    model.result('pg8').feature('surf1').set('rangedataactive', 'on');
    model.result('pg8').feature('surf1').set('descr', '10*log10(spf.vort_magn)');
    model.result('pg8').feature('surf1').set('rangecolormax', '60');
    model.result('pg8').feature('surf1').set('unit', '');
    model.result('pg8').feature('surf1').set('rangecoloractive', 'on');
    model.result('pg8').feature('surf1').set('expr', '10*log10(spf.vort_magn)');
    model.result('pg8').feature('surf1').set('resolution', 'normal');
    model.result('pg14').label('Pressure Axi');
    model.result('pg14').set('looplevel', {'49'});
    model.result('pg14').feature('vol1').set('rangedatamax', '2000');
    model.result('pg14').feature('vol1').set('rangedataactive', 'on');
    model.result('pg14').feature('vol1').set('rangecolormax', '2000');
    model.result('pg14').feature('vol1').set('rangedatamin', '-2000');
    model.result('pg14').feature('vol1').set('rangecolormin', '-2000');
    model.result('pg14').feature('vol1').set('rangecoloractive', 'on');
    model.result('pg14').feature('vol1').set('resolution', 'normal');
    model.result('pg14').feature('vol2').set('rangedatamax', '2000');
    model.result('pg14').feature('vol2').set('rangedataactive', 'on');
    model.result('pg14').feature('vol2').set('descr', 'Pressure');
    model.result('pg14').feature('vol2').set('rangecolormax', '2000');
    model.result('pg14').feature('vol2').set('rangedatamin', '-2000');
    model.result('pg14').feature('vol2').set('rangecolormin', '-2000');
    model.result('pg14').feature('vol2').set('rangecoloractive', 'on');
    model.result('pg14').feature('vol2').set('expr', 'p');
    model.result('pg14').feature('vol2').set('resolution', 'normal');
    model.result('pg16').label('Velocity Axi');
    model.result('pg16').set('looplevel', {'93'});
    model.result('pg16').set('edges', 'off');
    model.result('pg16').feature('vol1').set('descr', 'Velocity field, z component');
    model.result('pg16').feature('vol1').set('unit', 'm/s');
    model.result('pg16').feature('vol1').set('expr', 'w');
    model.result('pg16').feature('vol1').set('resolution', 'normal');
    model.result('pg15').label('Impedance (time)');
    model.result('pg15').set('xlabel', 'Time (s)');
    model.result('pg15').set('xlabelactive', false);
    model.result('pg15').feature('glob1').label('Liner');
    model.result('pg15').feature('glob1').set('descr', {''});
    model.result('pg15').feature('glob1').set('unit', {'Pa*s/m'});
    model.result('pg15').feature('glob1').set('expr', {'real(aveop1(-p))/real(aveop1(w))'});
    model.result('pg15').feature('glob2').active(false);
    model.result('pg15').feature('glob2').label('Plane 1');
    model.result('pg15').feature('glob2').set('descr', {''});
    model.result('pg15').feature('glob2').set('unit', {'Pa*s/m'});
    model.result('pg15').feature('glob2').set('expr', {'aveop2(p)/aveop2(w)'});
    model.result.export('anim1').label('Pressure Animation');
    model.result.export('anim1').set('showframe', '101');
    model.result.export('anim1').set('shownparameter', '0.01');
    model.result.export('anim1').set('framesel', 'all');
    model.result.export('anim1').set('target', 'player');
    model.result.export('anim1').set('frametime', '0.2');
    model.result.export('anim1').set('solnumtype', 'inner');
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
    model.result.export('anim2').label('Velocity magnitude animation');
    model.result.export('anim2').set('showframe', '101');
    model.result.export('anim2').set('shownparameter', '0.01');
    model.result.export('anim2').set('framesel', 'all');
    model.result.export('anim2').set('plotgroup', 'pg3');
    model.result.export('anim2').set('target', 'player');
    model.result.export('anim2').set('frametime', '0.2');
    model.result.export('anim2').set('title', 'on');
    model.result.export('anim2').set('legend', 'on');
    model.result.export('anim2').set('logo', 'on');
    model.result.export('anim2').set('options', 'off');
    model.result.export('anim2').set('fontsize', '9');
    model.result.export('anim2').set('customcolor', [1 1 1]);
    model.result.export('anim2').set('background', 'color');
    model.result.export('anim2').set('axisorientation', 'on');
    model.result.export('anim2').set('grid', 'on');
    model.result.export('anim2').set('axes', 'on');
    model.result.export('anim2').set('showgrid', 'on');
    model.result.export('anim3').label('Vorticity magnitude');
    model.result.export('anim3').set('showframe', '101');
    model.result.export('anim3').set('shownparameter', '0.01');
    model.result.export('anim3').set('framesel', 'all');
    model.result.export('anim3').set('plotgroup', 'pg8');
    model.result.export('anim3').set('target', 'player');
    model.result.export('anim3').set('frametime', '0.2');
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
%% Exporting the results
folder= [num2str(freq) 'Hz']
pathD='..\Data\';
mkdir(pathD,folder)

    model.result().table('tbl3').save([pathD num2str(freq) 'Hz\PLine.txt']);
    model.result().table('tbl17').save([pathD num2str(freq) 'Hz\VLine.txt']);
    model.result().table('tbl5').save([pathD num2str(freq) 'Hz\PrePlane1.txt']);
    model.result().table('tbl6').save([pathD num2str(freq) 'Hz\VelPlane1.txt']);
    model.result().table('tbl9').save([pathD num2str(freq) 'Hz\PreMic2.txt']);
    model.result().table('tbl10').save([pathD num2str(freq) 'Hz\AcMic2.txt']);
    model.result().table('tbl11').save([pathD num2str(freq) 'Hz\PreMic1.txt']);
    model.result().table('tbl12').save([pathD num2str(freq) 'Hz\AcMic1.txt']);
    model.result().table('tbl13').save([pathD num2str(freq) 'Hz\Pin.txt']);
    model.result().table('tbl14').save([pathD num2str(freq) 'Hz\Pout.txt']);
    model.result().table('tbl15').save([pathD num2str(freq) 'Hz\PrePlane2.txt']);
    model.result().table('tbl16').save([pathD num2str(freq) 'Hz\VelPlane2.txt']);
    Results=zeros(length(load([pathD num2str(freq) 'Hz\Pin.txt'])),2,12);
    Results(:,:,1)=load([pathD num2str(freq) 'Hz\Pin.txt']); %incident pressure inlet
    Results(:,:,2)=load([pathD num2str(freq) 'Hz\Pout.txt']); %reflected pressure inlet
    Results(:,:,3)=load([pathD num2str(freq) 'Hz\PrePlane1.txt']); %pressure plane 1 near (multiphysics)
    Results(:,:,4)=load([pathD num2str(freq) 'Hz\VelPlane1.txt']); %velocity plane 1
    Results(:,:,5)=load([pathD num2str(freq) 'Hz\PrePlane2.txt']); %pressure plane 2 away
    Results(:,:,6)=load([pathD num2str(freq) 'Hz\VelPlane2.txt']); %velocity plane 2
    Results(:,:,7)=load([pathD num2str(freq) 'Hz\PreMic1.txt']); %pressure mic1
    Results(:,:,8)=load([pathD num2str(freq) 'Hz\PreMic2.txt']); %pressure mic2
    Results(:,:,9)=load([pathD num2str(freq) 'Hz\AcMic1.txt']); %velocity mic1 (far)
    Results(:,:,10)=load([pathD num2str(freq) 'Hz\AcMic2.txt']); %velocity mic2
    Results(:,:,11)=load([pathD num2str(freq) 'Hz\PLine.txt']); %Pressure on the liner surface
    Results(:,:,12)=load([pathD num2str(freq) 'Hz\VLine.txt']); %Velocity on the liner surface
str=[pathD num2str(freq) 'Hz\Results_' num2str(freq) 'Hz']
    save(str,'Results')
  mphsave(model,['..\Data\' num2str(freq) 'Hz\Solved_' num2str(freq) 'Hz'])  
end
