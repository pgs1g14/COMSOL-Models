function out = model
%
% LinerAxi_PT_155dB_500Hz_nonsolved.m
%
% Model exported on Aug 5 2017, 17:53 by COMSOL 5.3.0.260.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Local\Pablo\COMSOL\Circular Axi model\LinerAxi_PT_155dB_500Hz');

model.label('LinerAxi_PT_155dB_500Hz.mph');

model.comments('Liner Axisimetrical PT thickness 0.000635 mm 155dB 500Hz. Impedance tube model with a sample in a holder with POA 3.66%, d=0.00099 mm, h=19.0mm nominal values.');

model.param.set('w_tube', '0.00517[m]', 'Tube width');
model.param.set('h_r', '0.019[m]', 'Resonator height');
model.param.set('h_slit', '0.000635[m]', 'Slit thickness');
model.param.set('w_slit', '0.00099[m]', 'width');
model.param.set('h_in', '0.2082[m]', 'Inlet height');
model.param.set('y0', '0.20852[m]', 'Inlet location y-coordinate');
model.param.set('h_tube', '0.22783[m]', 'Total tube length');
model.param.set('rho0', '1.2056[kg/m^3]', 'Density');
model.param.set('c0', '343.2043[m/s]', 'Speed of sound');
model.param.set('beta', '1.2', 'Nonlinear coefficent');
model.param.set('L0', '155', 'Incident wave amplitude dB');
model.param.set('p0', '1124.6827[Pa]', 'Incident wave amplitude Pa');
model.param.set('f0', '500[Hz]', 'Driving frequency');
model.param.set('Tstart', '7/f0', 'Stat time for post processing 7 cycles');
model.param.set('lambda0', '0.68641[m]', 'Wavelength at f0');
model.param.set('dvisc', '220[um]*sqrt(100[Hz]/f0)', 'Viscous boundary layer thickness at f0');
model.param.set('k0', '2*pi/lambda0', 'Wave number at f0');
model.param.set('omega0', 3141.592653589793, 'Angular frequency');
model.param.set('T0', '0.002[s]', 'Period');
model.param.set('dt_sol', '0.0001', 'Solver time step (resolve f0 with T0/30)');
model.param.set('lmesh', 'c0*dt_sol', 'Minimum element size acoustic mesh');
model.param.set('Tend', 'Tstart+0.01', 'End time for post processing');

model.comments(['LinerAxi PT 155dB 500Hz\n\nLiner Axisimetrical PT thickness 0.000635 mm 155dB 500Hz. Impedance tube model with a sample in a holder with POA 3.66%, d=0.00099 mm, h=19.0mm nominal values.']);
model.comments('Liner Axisimetrical PT thickness 0.000635 m 155dB 500Hz. Impedance tube model with a sample in a holder with POA 3.66%, d=0.00099 mm, h=19.0mm nominal values.');
model.comments('Liner Axisimetrical PT thickness 0.000635 m 155dB 500Hz. Impedance tube model with a sample in a holder with POA 3.66%, d=0.00099 mm, h=19.0mm nominal values.');

model.geom.create('part1', 'Part', 2);

model.comments(['LinerAxi PT 155dB 500Hz\n\nLiner Axisimetrical PT thickness 0.000635 m 155dB 500Hz. Impedance tube model with a sample in a holder with POA 3.66%, d=0.00099 mm, h=19.0mm nominal values.']);

model.modelNode.create('comp1', true);

model.geom.create('geom1', 2);
model.geom('geom1').model('comp1');
model.geom('geom1').axisymmetric(true);

model.mesh.create('mesh1', 'geom1');

model.geom.remove('part1');

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

model.geom('geom1').create('r1', 'Rectangle');
model.geom('geom1').feature('r1').set('pos', {'0' 'h_tube/2-(h_r+h_slit/2)'});

model.param.set('d_h', '0.00099[m]', 'width');

model.geom('geom1').feature('r1').setIndex('pos', 'w_tube/2', 0);
model.geom('geom1').feature('r1').set('pos', {'w_tube/2' 'h_r/2'});
model.geom('geom1').run('r1');
model.geom('geom1').feature('r1').set('pos', {'' 'h_r/2'});
model.geom('geom1').feature('r1').set('size', {'w_tube/2' '1'});
model.geom('geom1').feature('r1').set('pos', {'' ''});
model.geom('geom1').feature('r1').set('size', {'w_tube/2' 'h_r/2'});
model.geom('geom1').feature('r1').set('pos', [0 0]);
model.geom('geom1').run('r1');
model.geom('geom1').run('r1');
model.geom('geom1').create('r2', 'Rectangle');
model.geom('geom1').feature('r2').set('size', {'1' 'h_r/2'});
model.geom('geom1').feature('r2').set('pos', {'0' 'h_r/2'});

model.param.set('h_h', '0.000635[m]', 'Slit thickness');
model.param.set('h_h', '0.000635[m]', 'Slit thickness');

model.geom('geom1').feature('r2').set('size', {'1' 'h_h'});
model.geom('geom1').feature('r1').set('size', {'w_tube/2' 'h_r'});
model.geom('geom1').run('r1');
model.geom('geom1').feature('r2').set('size', {'d_h/2' 'h_h'});
model.geom('geom1').run('r2');
model.geom('geom1').feature('r2').set('pos', {'0' 'h_r'});
model.geom('geom1').run('r2');
model.geom('geom1').run('r2');
model.geom('geom1').create('sq1', 'Square');
model.geom('geom1').feature.remove('sq1');
model.geom('geom1').run('r2');
model.geom('geom1').create('r3', 'Rectangle');
model.geom('geom1').feature('r3').set('pos', {'0' 'h_r+h_h'});
model.geom('geom1').feature('r3').set('size', {'w_tube/2' 'h_in'});
model.geom('geom1').run('r3');

model.param.set('w_tube', '0.00517[m]', 'Tube diameter considering POA=3.66%');
model.param.set('d_tube', '0.00517[m]', 'Tube diameter considering POA=3.66%');

model.geom('geom1').feature('r3').set('size', {'d_tube/2' 'h_in'});
model.geom('geom1').run('r3');
model.geom('geom1').feature('r1').set('size', {'d_tube/2' 'h_r'});
model.geom('geom1').run('r1');
model.geom('geom1').run('r2');
model.geom('geom1').run('r3');
model.geom('geom1').create('pol1', 'Polygon');
model.geom('geom1').feature('pol1').set('type', 'open');
model.geom('geom1').feature('pol1').set('x', '0 d_tube/2');
model.geom('geom1').feature('pol1').set('y', '40*h_h 40*h_h');
model.geom('geom1').run('pol1');
model.geom('geom1').create('pol2', 'Polygon');
model.geom('geom1').feature('pol2').set('type', 'open');
model.geom('geom1').feature('pol2').set('x', '0 d_tube/2');
model.geom('geom1').feature('pol2').set('y', '0.07+h_h+h_r 0.07+h_h+h_r');
model.geom('geom1').run('pol2');
model.geom('geom1').run('pol2');

model.param.set('h_r', '0.019[m]', 'Resonator height');
model.param.set('h_h', '0.000635[m]', 'Slit thickness');
model.param.remove('h_slit');
model.param.remove('w_slit');
model.param.remove('w_tube');

model.geom('geom1').run('pol1');
model.geom('geom1').run('pol2');
model.geom('geom1').feature('pol1').set('y', '40*(h_h+h_r) 40*(h_h+h_r)');
model.geom('geom1').run('pol1');
model.geom('geom1').feature('pol1').set('y', '40*d_h 40*d_h');
model.geom('geom1').run('pol1');
model.geom('geom1').run('pol2');
model.geom('geom1').run('fin');

model.variable.create('var1');
model.variable.var1.model('comp1');
model.variable('var1').set('Pin', 'p0*sin(omega0*t+k0*(y-y0))', 'Incident pressure');
model.variable('var1').set('rhoL', 'rho0+p/c0^2', 'Linear density relation for air');
model.variable('var1').set('rhoNL', 'rho0+p/c0^2-1/(rho0*c0^4)*(beta-1)*p^2', 'Nonlinear density relation for air');

model.material.create('mat1', 'Common', 'comp1');
model.material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
model.material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
model.material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat1').propertyGroup.create('idealGas', 'Ideal gas');

model.cpl.create('intop1', 'Integration', 'geom1');
model.cpl.create('aveop1', 'Average', 'geom1');
model.cpl.create('aveop2', 'Average', 'geom1');
model.cpl.create('aveop3', 'Average', 'geom1');
model.cpl('intop1').selection.geom('geom1', 1);
model.cpl('intop1').selection.set([11]);
model.cpl('aveop1').selection.geom('geom1', 1);
model.cpl('aveop1').selection.set([6 14]);
model.cpl('aveop2').selection.geom('geom1', 1);
model.cpl('aveop2').selection.set([8]);
model.cpl('aveop3').selection.geom('geom1', 1);
model.cpl('aveop3').selection.set([10]);
model.cpl('aveop3').label('Acoustic Plane 3');
model.cpl('aveop2').label('Multiphysics Plane 2');
model.cpl('aveop1').label('Liner Plane 1');
model.cpl('intop1').label('Inlet 1');

model.physics.create('actd', 'TransientPressureAcoustics', 'geom1');
model.physics('actd').model('comp1');

model.material('mat1').set('family', 'air');

model.physics('actd').field('pressure').field('p2');
model.physics('actd').selection.set([4 5]);
model.physics('actd').create('pwr1', 'PlaneWaveRadiation', 1);
model.physics('actd').feature('pwr1').selection.set([11]);
model.physics('actd').feature('pwr1').create('ipf1', 'IncidentPressureField', 1);
model.physics('actd').create('nvel1', 'NormalVelocity', 1);
model.physics('actd').feature('nvel1').selection.set([8]);
model.physics('actd').feature('nvel1').set('Type', 'nvel');
model.physics('actd').feature('nvel1').set('nvel', 'v');
model.physics('actd').feature('nvel1').selection.set([8]);
model.physics.create('spf', 'LaminarFlow', 'geom1');
model.physics('spf').model('comp1');
model.physics('spf').selection.set([1 2 3]);
model.physics('spf').prop('PhysicalModelProperty').set('Compressibility', 'CompressibleMALT03');
model.physics('spf').create('bs1', 'BoundaryStress', 1);
model.physics('spf').feature('bs1').selection.set([8]);
model.physics('spf').feature('bs1').set('BoundaryCondition', 'NormalStress');
model.physics('actd').feature('nvel1').set('Type', 'vel');
model.physics('actd').feature('nvel1').set('vel', {'0' '0' 'v'});
model.physics('spf').create('wall2', 'Wall', 1);
model.physics('spf').feature.remove('wall2');
model.physics('actd').prop('TransientSettings').set('TimeStepping', 'Automatic');
model.physics('actd').prop('TransientSettings').set('fmax', '6400[Hz]');
model.physics('spf').prop('PhysicalModelProperty').set('StokesFlowProp', false);
model.physics('actd').feature('nvel1').set('vel', {'u' '0' 'v'});
model.physics('spf').feature('bs1').set('f0', 'p2');
model.physics('spf').create('wallbc2', 'WallBC', 1);
model.physics('spf').feature('wallbc2').set('BoundaryCondition', 'Slip');
model.physics('spf').feature('wallbc2').selection.set([15 16]);

model.material('mat1').set('family', 'air');
model.material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
model.material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/8.314/T');
model.material('mat1').propertyGroup('def').func('rho').set('dermethod', 'manual');
model.material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA' '0' '1'; 'T' '0' '1'});
model.material('mat1').propertyGroup('def').func('rho').set('argders', {'pA' 'd(pA*0.02897/8.314/T,pA)'; 'T' 'd(pA*0.02897/8.314/T,T)'});
model.material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('cs').set('args', 'T');
model.material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*287*T)');
model.material('mat1').propertyGroup('def').func('cs').set('dermethod', 'manual');
model.material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T' '0' '1'});
model.material('mat1').propertyGroup('def').func('cs').set('argders', {'T' 'd(sqrt(1.4*287*T),T)'});
model.material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T[1/K])[Pa*s]');
model.material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
model.material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T[1/K])[J/(kg*K)]');
model.material('mat1').propertyGroup('def').set('density', 'rho(pA[1/Pa],T[1/K])[kg/m^3]');
model.material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]'});
model.material('mat1').propertyGroup('def').set('soundspeed', 'cs(T[1/K])[m/s]');
model.material('mat1').propertyGroup('def').set('bulkviscosity', '0');
model.material('mat1').propertyGroup('def').addInput('temperature');
model.material('mat1').propertyGroup('def').addInput('pressure');
model.material('mat1').propertyGroup('RefractiveIndex').set('n', []);
model.material('mat1').propertyGroup('RefractiveIndex').set('ki', []);
model.material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat1').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.material('mat1').propertyGroup('idealGas').set('Rs', '287[J/kg/K]');
model.material('mat1').propertyGroup('idealGas').addInput('temperature');
model.material('mat1').propertyGroup('idealGas').addInput('pressure');

model.mesh('mesh1').automatic(false);
model.mesh('mesh1').feature('size').set('custom', true);
model.mesh('mesh1').create('map1', 'Map');
model.mesh('mesh1').feature('map1').create('size1', 'Size');
model.mesh('mesh1').feature('map1').feature('size1').selection.geom('geom1', 1);
model.mesh('mesh1').feature('map1').feature('size1').selection.set([11]);
model.mesh('mesh1').feature('map1').feature('size1').set('custom', true);
model.mesh('mesh1').feature('map1').feature('size1').set('hmaxactive', true);
model.mesh('mesh1').feature('map1').feature('size1').set('hmax', '6/f0');
model.mesh('mesh1').run('map1');
model.mesh('mesh1').feature.remove('bl1');
model.mesh('mesh1').feature.remove('cr1');
model.mesh('mesh1').feature.remove('size2');
model.mesh('mesh1').feature('size').set('table', 'cfd');
model.mesh('mesh1').feature.remove('size1');
model.mesh('mesh1').feature.remove('ftri1');
model.mesh('mesh1').run('size');
model.mesh('mesh1').run;
model.mesh('mesh1').create('ftri1', 'FreeTri');
model.mesh('mesh1').feature('size').set('custom', true);
model.mesh('mesh1').feature('size').set('hmax', 'lambda0/6');
model.mesh('mesh1').feature('size').set('hmin', 'dvisc/3');
model.mesh('mesh1').feature('size').set('hgrad', 1.1);
model.mesh('mesh1').feature('size').set('hnarrow', 3);
model.mesh('mesh1').feature('map1').feature('size1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('map1').feature('size1').selection.geom('geom1', 1);
model.mesh('mesh1').feature('map1').feature('size1').selection.set([12]);
model.mesh('mesh1').feature('map1').feature('size1').set('hmax', 'dvisc');
model.mesh('mesh1').feature('map1').feature('size1').set('hminactive', true);
model.mesh('mesh1').feature('map1').feature('size1').set('hmin', 'dvisc/2');

model.geom('geom1').run('pol2');
model.geom('geom1').create('r4', 'Rectangle');
model.geom('geom1').run;
model.geom('geom1').feature.remove('r4');
model.geom('geom1').run;

model.mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.mesh('mesh1').feature('map1').feature('size1').selection.set([11]);
model.mesh('mesh1').feature('map1').feature('size1').set('hminactive', false);
model.mesh('mesh1').feature('map1').feature('size1').set('hmax', 'd_tube/2/6');
model.mesh('mesh1').feature('map1').feature('size1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('map1').feature('size1').selection.set([5]);
model.mesh('mesh1').feature('map1').feature('size1').set('hmax', 'lmesh');
model.mesh('mesh1').feature('ftri1').feature('size1').set('custom', true);
model.mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 1);
model.mesh('mesh1').feature('ftri1').feature('size1').selection.set([12]);
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', true);
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', 'dvisc/2');
model.mesh('mesh1').feature('ftri1').create('size2', 'Size');
model.mesh('mesh1').feature('ftri1').feature('size2').selection.geom('geom1', 2);
model.mesh('mesh1').feature('ftri1').feature('size2').selection.set([1 2 3]);
model.mesh('mesh1').feature('ftri1').feature('size2').set('custom', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmaxactive', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc*elementNear');
model.mesh('mesh1').feature('ftri1').feature('size2').set('hminactive', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmin', 'dvisc*elementNear');

model.param.set('elementNear', '15');
model.param.descr('elementNear', 'Times that the viscous boundary layer is considered in the domain');
model.param.set('elementFar', '24');
model.param.descr('elementFar', 'Times that the viscous boundary layer is considered in the domain far from the hole');
model.param.descr('elementNear', 'Times that the viscous boundary layer is considered in the domain near the hole');

model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc*elementFar');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('map1').feature('size1').set('hmax', 'd_tube/8');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size1').set('table', 'cfd');
model.mesh('mesh1').feature('ftri1').feature('size1').set('custom', true);
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', true);
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', 'dvisc/2');
model.mesh('mesh1').feature('ftri1').feature('size2').set('table', 'cfd');
model.mesh('mesh1').feature('ftri1').feature('size2').set('custom', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmaxactive', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc*elementFar');
model.mesh('mesh1').feature('ftri1').feature('size2').set('hminactive', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmin', 'dvisc*elementNear');
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', false);
model.mesh('mesh1').feature('ftri1').feature('size1').set('custom', false);
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('custom', false);
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature.remove('size1');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hauto', 1);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('map1').feature('size1').selection.set([4 5]);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').selection.set([1 2 3]);
model.mesh('mesh1').feature('ftri1').feature('size2').set('custom', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmaxactive', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc/2');
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri1').active(false);
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature('size').set('custom', false);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('size').set('custom', true);
model.mesh('mesh1').feature('size').set('hmax', 'lambda/6');
model.mesh('mesh1').run('size');
model.mesh('mesh1').feature('size').set('hmax', 'lambda0/6');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('size').set('custom', false);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('size').set('hauto', 9);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('map1').active(false);
model.mesh('mesh1').feature('ftri1').active(true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hminactive', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmin', 'dvisc/2');
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc*elementFar');
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').active(false);
model.mesh('mesh1').feature('size').set('custom', false);
model.mesh('mesh1').feature('size').set('table', 'default');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('size').set('hauto', 5);
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run('size');
model.mesh('mesh1').feature('size').set('custom', false);
model.mesh('mesh1').feature('size').set('table', 'cfd');
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh.remove('mesh1');
model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').automatic(false);
model.mesh('mesh1').automatic(true);
model.mesh('mesh1').run;
model.mesh('mesh1').autoMeshSize(6);
model.mesh('mesh1').run;

model.physics('spf').feature('wallbc2').selection.set([15 16]);

model.mesh('mesh1').run;
model.mesh('mesh1').autoMeshSize(5);
model.mesh('mesh1').run;
model.mesh('mesh1').autoMeshSize(4);
model.mesh('mesh1').run;
model.mesh('mesh1').autoMeshSize(1);
model.mesh('mesh1').run;
model.mesh('mesh1').autoMeshSize(3);
model.mesh('mesh1').run;
model.mesh('mesh1').autoMeshSize(5);
model.mesh('mesh1').run;

model.result.table('tbl1').clearTableData;
model.result.table('tbl2').clearTableData;
model.result.table('tbl3').clearTableData;
model.result.table('tbl4').clearTableData;
model.result.table('tbl5').clearTableData;
model.result.table('tbl6').clearTableData;
model.result.table('tbl9').clearTableData;
model.result.table('tbl10').clearTableData;
model.result.table('tbl11').clearTableData;
model.result.table('tbl12').clearTableData;
model.result.table('tbl13').clearTableData;
model.result.table('tbl14').clearTableData;
model.result.table('evl2').clearTableData;
model.result.table('tbl15').clearTableData;
model.result.table('tbl16').clearTableData;
model.result.table.clear;

model.cpl('aveop1').label('Average Liner');
model.cpl('aveop2').label('Average Plane Interface');
model.cpl('aveop3').label('Average Plane Acoustic');

model.physics('actd').feature('tpam1').set('c', 'c0');
model.physics('spf').prop('PhysicalModelProperty').set('EnablePorousMediaDomains', false);
model.physics('actd').feature('tpam1').set('rho', 'rho0');

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

model.physics('actd').prop('ReferencePressure').set('ReferenceType', 'ReferencePressureAir');

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

model.geom('geom1').create('pt1', 'Point');
model.geom('geom1').feature('pt1').setIndex('p', 'd_tube/2', 0, 0);
model.geom('geom1').feature('pt1').setIndex('p', '0.119+h_h+h_r', 1, 0);
model.geom('geom1').runPre('fin');
model.geom('geom1').create('pt2', 'Point');
model.geom('geom1').feature('pt2').setIndex('p', 'd_tube/2', 0, 0);
model.geom('geom1').feature('pt2').setIndex('p', '0.099+h_h+h_r', 1, 0);
model.geom('geom1').run;

model.mesh('mesh1').run;

model.label('LinerAxi_PT_155dB_500Hz.mph');

model.physics('actd').feature('pwr1').feature('ipf1').set('p', 'Pin');

model.variable('var1').set('Pin', 'p0*sin(omega0*t)');

model.param.set('Tstart', '5/f0');

model.mesh('mesh1').create('map1', 'Map');
model.mesh('mesh1').feature('map1').set('adjustedgdistr', true);
model.mesh('mesh1').feature('map1').create('size1', 'Size');
model.mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('map1').selection.set([4 5]);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('map1').feature('size1').set('table', 'cfd');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('map1').feature('size1').set('table', 'default');
model.mesh('mesh1').feature('map1').feature('size1').set('hauto', 3);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('map1').feature('size1').set('hauto', 1);
model.mesh('mesh1').run;
model.mesh('mesh1').create('ftri1', 'FreeTri');
model.mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 1);
model.mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('ftri1').feature('size1').selection.set([1 2 3]);
model.mesh('mesh1').feature('ftri1').feature('size1').set('table', 'cfd');
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 7);
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri1').create('cr1', 'CornerRefinement');
model.mesh('mesh1').feature('ftri1').feature('cr1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('ftri1').feature.remove('cr1');
model.mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 1);
model.mesh('mesh1').feature('ftri1').feature('size1').selection.set([12]);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('ftri1').feature('size1').selection.set([1 2 3]);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').create('size2', 'Size');
model.mesh('mesh1').feature('ftri1').feature('size2').selection.geom('geom1', 1);
model.mesh('mesh1').feature('ftri1').feature('size2').selection.set([12]);
model.mesh('mesh1').run('ftri1');
model.mesh('mesh1').feature('ftri1').feature('size2').set('hauto', 3);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('table', 'cfd');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hauto', 2);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hauto', 3);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('custom', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmaxactive', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc*2');
model.mesh('mesh1').feature('ftri1').feature('size2').set('hminactive', true);
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmin', 'dvisc/4');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmin', 'dvisc/8');
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 8);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmin', 'dvisc/16');
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc/2');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 9);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc/4');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hgradactive', true);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmin', 'dvisc/8');
model.mesh('mesh1').run;
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc/2');
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hgrad', 1.05);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hcurveactive', true);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hcurveactive', false);
model.mesh('mesh1').feature('ftri1').feature('size2').selection.set([4 6 12]);
model.mesh('mesh1').run;
model.mesh('mesh1').feature('ftri1').feature('size2').set('hmax', 'dvisc');
model.mesh('mesh1').run;

model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').label('Acoustic Pressure (actd)');
model.result('pg1').set('data', 'dset1');
model.result('pg1').feature.create('surf1', 'Surface');
model.result('pg1').feature('surf1').label('Surface');
model.result('pg1').feature('surf1').set('data', 'parent');
model.result.dataset.create('rev1', 'Revolve2D');
model.result.dataset('rev1').label('Revolution 2D');
model.result.dataset('rev1').set('startangle', -90);
model.result.dataset('rev1').set('revangle', 225);
model.result.dataset('rev1').set('data', 'dset1');
model.result.create('pg2', 'PlotGroup3D');
model.result('pg2').label('Acoustic Pressure, 3D (actd)');
model.result('pg2').set('data', 'rev1');
model.result('pg2').feature.create('surf1', 'Surface');
model.result('pg2').feature('surf1').label('Surface');
model.result('pg2').feature('surf1').set('data', 'parent');
model.result.create('pg3', 'PlotGroup2D');
model.result('pg3').label('Velocity (spf)');
model.result('pg3').set('frametype', 'spatial');
model.result('pg3').set('data', 'dset1');
model.result('pg3').feature.create('surf1', 'Surface');
model.result('pg3').feature('surf1').label('Surface');
model.result('pg3').feature('surf1').set('expr', 'spf.U');
model.result('pg3').feature('surf1').set('data', 'parent');
model.result.create('pg4', 'PlotGroup2D');
model.result('pg4').label('Pressure (spf)');
model.result('pg4').set('frametype', 'spatial');
model.result('pg4').set('data', 'dset1');
model.result('pg4').feature.create('con1', 'Contour');
model.result('pg4').feature('con1').label('Contour');
model.result('pg4').feature('con1').set('expr', 'p');
model.result('pg4').feature('con1').set('number', 40);
model.result('pg4').feature('con1').set('data', 'parent');
model.result.dataset('rev1').set('data', 'dset1');
model.result.create('pg5', 'PlotGroup3D');
model.result('pg5').label('Velocity (spf) 1');
model.result('pg5').set('frametype', 'spatial');
model.result('pg5').set('data', 'rev1');
model.result('pg5').feature.create('surf1', 'Surface');
model.result('pg5').feature('surf1').label('Surface');
model.result('pg5').feature('surf1').set('expr', 'spf.U');
model.result('pg5').feature('surf1').set('data', 'parent');
model.result.remove('pg2');
model.result.remove('pg4');
model.result.remove('pg3');
model.result.remove('pg5');

model.study('std1').feature('time').set('plot', false);

out = model;
