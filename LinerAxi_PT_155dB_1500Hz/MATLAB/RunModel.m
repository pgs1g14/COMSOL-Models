%% Run FEM Model
% The parametric model choosing the start stop and frequency step

clear all
clc
tic

freq=600:500:1100;
cycles=5;

hi = waitbar(0,'Calculating FEM model, parametric on frequency.');

for i=2:1:length(freq)
    waitbar(i/length(freq),hi)
    freq(i)
    %Calls the numerical model using COMSOL
    LinerAxi_PT_155dB_parametric_v5_2(cycles,freq(i)) %Run the model and extract results
    % Load the Time domain results
    DataPath=['..\Data\' num2str(freq(i)) 'Hz\Results_' num2str(freq(i)) 'Hz.mat'];
    eval(['data=load(DataPath)'])
    %Calculate the impedance and save the figures to analysis
    ImpedanceCalc(j,freq(i),load(DataPath))
    ImpPath=['..\Data\' num2str(freq(i)) 'Hz\Impedances_' num2str(freq(i)) 'Hz_8_Cycles.mat'];
    imp=load(ImpPath)
    % Create a vector with the FFT(Liner) impedance result
    Num(i,1)=freq(i);
    Num(i,2)=imp.Z(1);
    % Plot the intermediate liner results
    figure(30)
    plot(Num(i,1),real(Num(:,2)),'-ob')
    hold on;
    plot(Num(i,1),imag(Num(:,2)),'-xb')
end
close(hi)

%% Plot Numerical x Experimental results
%Export Experimental Result
Exp=load('C:\Local\Pablo\PIM measurements\Post Process\Sample Holder\Measurement 14.08.2015\PureTone Process\PT_IMP.mat')
h=figure(1)
hold on;
plot(Exp.PT_IMP(:,1),Exp.PT_IMP(:,12),'-or')
plot(Num(:,1),real(Num(:,2)),'-ob')
hold on;
plot(Exp.PT_IMP(:,1),Exp.PT_IMP(:,13),'-xr')
plot(Num(:,1),imag(Num(:,2)),'-xb')
legend('\theta - Exper. Sample Holder','\theta - Numerical','\chi - Exper. Sample Holder','\chi - Numerical')
   string = ['..\Data\Validation2ndTry']
                saveas(h,[string '.eps'],'epsc2')
                saveas(h,[string '.fig'])
                saveas(h,[string '.png'])
toc

