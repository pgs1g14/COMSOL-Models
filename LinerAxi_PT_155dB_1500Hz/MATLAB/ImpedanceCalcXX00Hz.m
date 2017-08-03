%function ImpedanceCalc3000Hz(num)  % Extracting average pressures on the planes and microfones from txt
   % file exported from COMSOL
%% Constants
    num=1;
    
for num=2:9

    freq=3000;

TC=20; %Temperature of the measurements to run the simulation ºC
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
terms=1;
%Liner geometry
w_tube=0.029; %impedance tube diameter [m]
h_r=0.0191; %Liner cavity depth - resonator height [m]
h_slit=0.00083; %Liner slit thickness [m]
POA_holder=0.0366; %percentage of open area of the sample in the holder
POA_liner=0.0518; %percentage of open area of the flanged liner sample
h_in=0.2082; %impedance tube length [m]
d_h=POA_holder*w_tube; %Hole diameter 1,0614 mm
d_l=POA_liner*w_tube; %Hole diameter 1,5022 mm    
yPlane2=0.07-h_slit/2; %Y plane 2 acoustics
yPlane1=40*h_slit; %Y plane 1 multiphysics
y0=h_in+h_slit/2;   

%num is the number of cycles taken of the solution at the beginning to
%avoid the transient solution
Tstart=3.5*y0/c0+num/freq;
j=1;
 w=2*pi()*freq(j);   
 k0=w/c0;
 
DataPath='../Data/'
            Pin=load([DataPath 'Pin.txt']); %incident pressure inlet
             Pout=load([DataPath 'Pout.txt']); %reflected pressure inlet
             PrePlane1=load([DataPath 'PrePlane1.txt']); %pressure plane 1 near (multiphysics)
             VelPlane1=load([DataPath 'VelPlane1.txt']); %velocity plane 1
             PrePlane2=load([DataPath 'PrePlane2.txt']); %pressure plane 2 away
             VelPlane2=load([DataPath 'VelPlane2.txt']); %velocity plane 2
             PreMic1=load([DataPath 'PreMic1.txt']); %pressure mic1
             PreMic2=load([DataPath 'PreMic2.txt']); %pressure mic2
             VelMic1=load([DataPath 'AcMic1.txt']); %velocity mic1 (far)
             VelMic2=load([DataPath 'AcMic2.txt']); %velocity mic2
            PLine=load([DataPath 'PLine.txt']); %Pressure on the liner surface
            VLine=load([DataPath 'VLine.txt']); %Velocity on the liner surface
            %% Extract the data of the struct
             PIN=Pin;
             POUT=Pout;
             P1=PrePlane1;
             V1=VelPlane1;
             P2=PrePlane2;
             V2=VelPlane2;
             PM1=PreMic1;
             PM2=PreMic2;
             VM1=VelMic1;
             VM2=VelMic2;
            Pl=PLine; 
            Vl=VLine;        
   % Preparing the signal 
     % take out 3000 micro seconds of the signal to keep only the steadystate
      %[Tstart,u,def,d] = Tstart;
            samp=round(Tstart/(Pl(2,1)-Pl(1,1)));
            pinc=PIN(samp:length(PIN),2); %permanent solution of incident pressure in the inlet
            pref=POUT(samp:length(POUT),2); %permanent solution of reflected pressure in the inlet
            p1=P1(samp:length(P1),2); %permanent solution of pressure plane 1
            v1=V1(samp:length(V1),2); %permanent solution of velocity plane 1
            p2=P2(samp:length(P2),2); %permanent solution of pressure plane 1
            v2=V2(samp:length(V2),2); %permanent solution of velocity plane 1
            pm1=PM1(samp:length(PM1),2); %extract just the permanent solution of pressure mic1
            vm1=VM1(samp:length(VM1),2); %extract just the permanent solution of velocity mic1
            pm2=PM2(samp:length(PM2),2); %extract just the permanent solution of pressure mic2
            vm2=VM2(samp:length(VM2),2); %extract just the permanent solution of velocity mic2
            pliner=Pl(samp:length(Pl),2); %extract just the permanent solution of pressure
            vliner=Vl(samp:length(Vl),2); %extract just the permanent solution of velocity
            time=Pl(samp:length(Pl),1);          
            %timehole=PreH(samp:length(PreH),1)*10^(-6); %extract time response
            timeliner=Pl(samp:length(Pl),1); %extract time response
            % Frequency domain
            step=(Pl(2,1)-Pl(1,1)); %time step       
            Dt=step;
            maxtime=(time(length(time))-time(1));
            tempo=0:step:maxtime;   %time vector  
            Fs=1/Dt;
            ls=length(time); %Length of the time signal
            Fmax=Fs/2; % Maximum Frequency analysed 5kHz to 100ms steps
            deltaf=Fs/ls; %delta frequency
            fs=-Fmax:deltaf:Fmax-Fs/ls; %frequency vector   
   % Impedance P/V at the liner surface
            Pliner_f = fft(-pliner);
            Vliner_f = fft(vliner);
            HP=fftshift(Pliner_f);
            HV=fftshift(Vliner_f);         
            zetaLiner(j)=interp1(fs,HP,freq(j))./interp1(fs,HV,freq(j))/rho/c0 %Numerical value calculated at the surface
            zetaLiner1(j)=interp1(fs,HP./HV,freq(j))/rho/c0 %Numerical value calculated at the surface
            erreLiner(j)=-(1-zetaLiner(j))./(1+zetaLiner(j))
            abs(erreLiner(j))
            
            PPP1=fit(time,-pliner,['fourier' num2str(terms)]); %Furthest
            coeffs=coeffvalues(PPP1);
            coeffs(terms*2+2)=w;
            PP1=fit(time,-pliner,['fourier' num2str(terms)],'StartPoint',coeffs);
            coeffs=coeffvalues(PP1);
            A0=coeffs(1);
            A1=coeffs(2)
            B1=coeffs(3)
            omega=coeffs(terms*2+2);
        
         
           PPP2=fit(time,vliner,['fourier' num2str(terms)]); %Nearer
            coeffs=coeffvalues(PPP2);
            coeffs(terms*2+2)=w;
            PP2=fit(time,vliner,['fourier' num2str(terms)],'StartPoint',coeffs);
            coeffs=coeffvalues(PP2);
            A02=coeffs(1);
            A2=coeffs(2)
            B2=coeffs(3)
            omega2=coeffs(terms*2+2);
         theta=atan(B1/A1)
         theta2=atan(B2/A2)         
            zzz=(A1-1i*B1)/(A2-1i*B2)/rho/c0;
            
    % Impedance P/V at the Plane 1
            P1_f = fft(-p1);
            V1_f = fft(v1);
            HP1=fftshift(P1_f);
            HV1=fftshift(V1_f);         
            zetaP1(j)=interp1(fs,HP1,freq(j))./interp1(fs,HV1,freq(j))/rho/c0; %Numerical value calculated at the surface
            erreP1(j)=-(1-zetaP1(j))./(1+zetaP1(j))
            abs(erreP1(j)) 
            xp1=yPlane1;
            zLinerP1(j)=(zetaP1(j)-1i*tan(k0*xp1))./(1-zetaP1(j)*1i*tan(k0*xp1))
      % Impedance P/V at the Plane 2
            P2_f = fft(-p2);
            V2_f = fft(v2);
            HP2=fftshift(P2_f);
            HV2=fftshift(V2_f);         
            zetaP2(j)=interp1(fs,HP2,freq(j))./interp1(fs,HV2,freq(j))/rho/c0*1i*w; %Numerical value calculated at the surface
            erreP2(j)=-(1-zetaP2(j))./(1+zetaP2(j))
            abs(erreP2(j)) 
            xp2=yPlane2;
            zLinerP2(j)=(zetaP2(j)-1i*tan(k0*xp2))./(1-zetaP2(j)*1i*tan(k0*xp2))
      % Impedance P/V at the Mic 1
            PM1_f = fft(-pm1);
            VM1_f = fft(vm1);
            HPM1=fftshift(PM1_f);
            HVM1=fftshift(VM1_f);         
            zetaPM1(j)=interp1(fs,HPM1,freq(j))./interp1(fs,HVM1,freq(j))/rho/c0*1i*w; %Numerical value calculated at the surface
            errePM1(j)=(zetaPM1(j)-1)./(1+zetaPM1(j))
            abs(errePM1(j))
            zLinerPM1(j)=(zetaPM1(j)-1i*tan(k0*x1))./(1-zetaPM1(j)*1i*tan(k0*x1))
      % Impedance P/V at the Mic 2
            PM2_f = fft(-pm2);
            VM2_f = fft(vm2);
            HPM2=fftshift(PM2_f);
            HVM2=fftshift(VM2_f);         
            zetaPM2(j)=interp1(fs,HPM2,freq(j))./interp1(fs,HVM2,freq(j))/rho/c0*1i*w; %Numerical value calculated at the surface
            errePM2(j)=-(1-zetaPM2(j))./(1+zetaPM2(j))
            abs(errePM2(j)) 
            zLinerPM2(j)=(zetaPM2(j)-1i*tan(k0*x2))./(1-zetaPM2(j)*1i*tan(k0*x2))
      % ASTM transfer function method    
      H21(j)=interp1(fs,HPM2,freq(j))./interp1(fs,HPM1,freq(j));
      erreASTM(j)=(H21(j)-exp(-1i*k0*(x1-x2)))./(exp(1i*k0*(x1-x2))-H21(j)).*exp(1i*2*k0*x1)
      abs(erreASTM(j))
      zLinerASTM(j)=(1+erreASTM(j))./(1-erreASTM(j))
      % Impedance at the Inlet
            PINC_f = fft(pinc);
            PREF_f = fft(pref);
            PINC=fftshift(PINC_f);
            PREF=fftshift(PREF_f);         
            erreInlet(j)=interp1(fs,PREF,freq(j))./interp1(fs,PINC,freq(j)) %Numerical value calculated at the inlet
            zetaInlet(j)=(1+erreInlet(j))./(1-erreInlet(j))
            abs(erreInlet(j)) 
            L=h_in;
            zLinerInlet(j)=(zetaInlet(j)-1i*tan(k0*L))./(1-zetaInlet(j)*1i*tan(k0*L))
      % Calculating the Impedance based on the A and B at the 2 planes
    %  Fitting curve
            PP1=fit(time,pm1,['fourier' num2str(terms)]); %Furthest
            coeffs=coeffvalues(PP1);
            coeffs(terms*2+2)=w;
            PP1=fit(time,pm1,['fourier' num2str(terms)],'StartPoint',coeffs);
            coeffs=coeffvalues(PP1);
            A0=coeffs(1);
            A2=coeffs(2)
            B2=coeffs(3)
            omega=coeffs(terms*2+2);
            
            PP2=fit(time,pm2,['fourier' num2str(terms)]); %Nearer
            coeffs=coeffvalues(PP2);
            coeffs(terms*2+2)=w;
            PP2=fit(time,pm2,['fourier' num2str(terms)],'StartPoint',coeffs);
            coeffs=coeffvalues(PP2);
            A02=coeffs(1);
            A1=coeffs(2)
            B1=coeffs(3)
            omega2=coeffs(terms*2+2);
    % Calculating a and b
            C1=A1-1i*B1;
            C2=A2-1i*B2;
            d=x1-x2;
            h1=x2;
            a(j)=1i*(C1*exp(-1i*k0*d)-C2)/sin(k0*d)/2
            b(j)=-1i*(C1*exp(1i*k0*d)-C2)/sin(k0*d)/2
   %  calculating impedance
            zP1(j)=(a(j)+b(j))/(a(j)-b(j));
            R(j)=-(1-zP1(j))/(1+zP1(j));
            abs(R(j))
  %   putting the impedance on the liner position
            zSurf1(j)=(exp(-1i*2*k0*h1)+R(j))/(exp(-2*1i*k0*h1)-R(j));
            zSurf1(j)=(a(j)*exp(-1i*k0*h1)+b(j)*exp(1i*k0*h1))/(a(j)*exp(-1i*k0*h1)-b(j)*exp(1i*k0*h1))           
  %% Cheching signals
     h=figure(2)
     plotTimeFreq(freq,j,time,vliner,pliner,fs,HV,HP,'liner')
     hold on;
     plotTimeFreq(freq,j,time,v1,p1,fs,HV1,HP1,'plane1')
     plotTimeFreq(freq,j,time,v2/w,p2,fs,HV2/1i/w,HP2,'plane2')
     plotTimeFreq(freq,j,time,vm1/w,pm1,fs,HVM1/1i/w,HPM1,'Mic1')
     plotTimeFreq(freq,j,time,vm2/w,pm2,fs,HVM2/1i/w,HPM2,'Mic2')
     legend('Liner Surface','Plane Multiphysics','Plane Acoustics','Mic 1','Mic 2','Location','Best')

   %   Saving signals
    set(gcf, 'Position', get(0,'ScreenSize')); 
    string = ['..\Figures\Signals\Signals_' num2str(SPL) 'dB_' num2str(freq(j)) 'Hz']
    saveas(h,[string '.eps'],'epsc2')
    saveas(h,[string '.fig'])
    saveas(h,[string '.png'])
close(h)
    %% Impedances using several methods
    h=figure(800)
    subplot(2,1,1)
       title('Impedances at the planes')
        hold on
         plot(freq(j),real(zetaLiner(j)),'<k')
         plot(freq(j),real(zLinerP2(j)),'*k')
         plot(freq(j),real(zLinerP1(j)),'>k')
         plot(freq(j),real(zLinerPM2(j)),'ok')
         plot(freq(j),real(zLinerPM1(j)),'^k')
         plot(freq(j),real(zLinerInlet(j)),'dk')
         plot(freq(j),real(zSurf1(j)),'vk')
         plot(freq(j),real(zLinerASTM(j)),'xk')
         grid on;
         grid minor;
         legend('FFT(Liner)','FFT(Plane1)','FFT(Plane2)','FFT(Mic1)','FFT(Mic2)','IRW(Inlet)','FC(Mic2)','ASTM','Location','Best')
         xlabel('Frequency [Hz]')
         ylabel('Normalized Resistance, \theta')
         axis([0 3500 0 3])
    subplot(2,1,2)
        hold on
         plot(freq(j),imag(zetaLiner(j)),'<k')
         plot(freq(j),imag(zLinerP1(j)),'>k')
         plot(freq(j),imag(zLinerP2(j)),'*k')
         plot(freq(j),imag(zLinerPM2(j)),'ok')
         plot(freq(j),imag(zLinerPM1(j)),'^k')
         plot(freq(j),imag(zLinerInlet(j)),'dk')
         plot(freq(j),imag(zSurf1(j)),'vk')
         plot(freq(j),imag(zLinerASTM(j)),'xk')
         grid on;
         grid minor;
        xlabel('Frequency [Hz]','FontSize',12)
         ylabel('Normalized Rectance, \chi','FontSize',12)
       axis([0 3500 -10 10])
     % Saving Impedance
     %set(gcf, 'Position', get(0,'ScreenSize')); 
                string = ['..\Figures\Impedance\Impedance_' num2str(SPL) 'dB_' num2str(freq(j)) 'Hz']
                saveas(h,[string '.eps'],'epsc2')
                saveas(h,[string '.fig'])
                saveas(h,[string '.png'])
                close(h)
    %% Ploting Impedance spectrum and harmonics          
        h=figure(1100)
        plot(fs,real(HP./HV/rho/c0),'o',fs,imag(HP./HV/rho/c0),'o')
        hold on;
        xlabel('Frequency [Hz]','FontSize',12)
        ylabel('Normalized Impedance','FontSize',12)
        legend('Resistance','Reactance')
        plot([freq freq],get(gca,'ylim'),'-k')
        axis([0 6400 -10 10])
        txt =['\leftarrow Fundamental Resi. = ' num2str(real(zetaLiner1(j)))];
        txt2 =['\leftarrow Fundamental Reac. = ' num2str(imag(zetaLiner1(j)))];
        text(freq,real(zetaLiner1(j)),txt)
        text(freq,imag(zetaLiner1(j)),txt2)

        z3h=interp1(fs,HP,3*freq(j))./interp1(fs,HV,3*freq(j))/rho/c0;
        terms=2
        for i=3:2:terms*3
        z3h1(i)=interp1(fs,HP./HV,i*freq)/rho/c0;

        plot([i*freq i*freq],get(gca,'ylim'),'-k')
        txt =['\leftarrow harmonic Resi. = ' num2str(real(z3h1(i)))];
        txt2 =['\leftarrow harmonic Reac. = ' num2str(imag(z3h1(i)))];
        text(i*freq,real(z3h1(i)),txt)
        text(i*freq,imag(z3h1(i)),txt2)
        end        
        % Saving Impedance of the harmonics
        set(gcf, 'Position', get(0,'ScreenSize')); 
                string = ['..\Figures\Impedance\ImpedanceSpectrum_' num2str(SPL) 'dB_' num2str(freq(j)) 'Hz']
                saveas(h,[string '.eps'],'epsc2')
                saveas(h,[string '.fig'])
                saveas(h,[string '.png'])
    %% Reflection coefficients
   close(h)
   h=figure(750)
    subplot(2,1,1)
        hold on
        plot(freq(j),abs(erreLiner(j)),'<k')
        plot(freq(j),abs(erreP1(j)),'>k')
        plot(freq(j),abs(erreP2(j)),'*k')
        plot(freq(j),abs(errePM2(j)),'ok')
        plot(freq(j),abs(errePM1(j)),'^k')
        plot(freq(j),abs(erreInlet(j)),'dk')
        plot(freq(j),abs(R(j)),'vk')
        plot(freq(j),abs(erreASTM(j)),'xk')
        grid on;
        grid minor;
        legend('FFT(Liner)','FFT(Plane1)','FFT(Plane2)','FFT(Mic1)','FFT(Mic2)','IRW(Inlet)','FC(Mic2)','ASTM','Location','Best')
        xlabel('Frequency [Hz]','FontSize',12)
        ylabel('\Gamma Magnitude','FontSize',12)
        axis([0 3500 0 1])
    subplot(2,1,2)
        hold on
        plot(freq(j),abs(atan(imag(erreLiner(j))/real(erreLiner(j)))*180/2/pi),'<k')
        plot(freq(j),abs(atan(imag(erreP1(j))/real(erreP1(j)))*180/2/pi),'>k')
        plot(freq(j),abs(atan(imag(erreP2(j))/real(erreP2(j)))*180/2/pi),'*k')
        plot(freq(j),abs(atan(imag(errePM2(j))/real(errePM2(j)))*180/2/pi),'ok')
        plot(freq(j),abs(atan(imag(errePM1(j))/real(errePM1(j)))*180/2/pi),'^k')
        plot(freq(j),abs(atan(imag(erreInlet(j))/real(erreInlet(j)))*180/2/pi),'dk')
        plot(freq(j),abs(atan(imag(R(j))/real(R(j)))*180/2/pi),'vk')
        plot(freq(j),abs(atan(imag(erreASTM(j))/real(erreASTM(j)))*180/2/pi),'xk')
        grid on;
        grid minor;
        xlabel('Frequency [Hz]','FontSize',12)
        ylabel('\Gamma Phase [degrees]','FontSize',12)   
        axis([0 3500 -180 180])
       %Save Reflection coefficient
       %set(gcf, 'Position', get(0,'ScreenSize')); 
                string = ['..\Figures\Impedance\ReflectionCoefficient_' num2str(SPL) 'dB_' num2str(freq(j)) 'Hz']
                saveas(h,[string '.eps'],'epsc2')
                saveas(h,[string '.fig'])
                saveas(h,[string '.png'])
close(h)
   
%          wl=c0/freq;
%         phaseP1=yPlane1/wl
%         phaseP2=yPlane2/wl
%         phasePM2=x2/wl
%         phasePM1=x1/wl
%         phaseIn=h_in/wl
%         
%         plot(freq(j),atan(imag(erreLiner(j))/real(erreLiner(j)))*180/pi,'<k')
%         plot(freq(j),(-phaseP1+atan(imag(erreP1(j))/real(erreP1(j))))*180/pi,'>k')
%         plot(freq(j),(-phaseP2+atan(imag(erreP2(j))/real(erreP2(j))))*180/pi,'*k')
%         plot(freq(j),(-phasePM2+atan(imag(errePM2(j))/real(errePM2(j))))*180/pi,'ok')
%         plot(freq(j),(-phasePM1+atan(imag(errePM1(j))/real(errePM1(j))))*180/pi,'^k')
%         plot(freq(j),(-phaseIn+atan(imag(erreInlet(j))/real(erreInlet(j))))*180/pi,'dk')
%         plot(freq(j),(-phasePM1+atan(imag(R(j))/real(R(j))))*180/pi,'vk')
%         plot(freq(j),(-phasePM1+atan(imag(erreASTM(j))/real(erreASTM(j))))*180/pi,'xk')
%         
%               
  %% Cycle analysis     
        
    h=figure(900)
    subplot(2,1,1)
       title(['Total time convergence analysis for pure tone ' num2str(freq) ' Hz'])
        hold on
         plot(10-num,real(zetaLiner(j)),'<k')
         plot(10-num,real(zLinerP2(j)),'*k')
         plot(10-num,real(zLinerP1(j)),'>k')
         plot(10-num,real(zLinerPM2(j)),'ok')
         plot(10-num,real(zLinerPM1(j)),'^k')
         plot(10-num,real(zLinerInlet(j)),'dk')
         plot(10-num,real(zSurf1(j)),'vk')
         plot(10-num,real(zLinerASTM(j)),'xk')
         
         X(1,10-num)=zetaLiner(j);
         X(2,10-num)=zLinerP1(j);
         X(3,10-num)=zLinerP2(j);
         X(4,10-num)=zLinerPM2(j);
         X(5,10-num)=zLinerPM1(j);
         X(6,10-num)=zLinerInlet(j);
         X(7,10-num)=zSurf1(j);
         X(8,10-num)=zLinerASTM(j);
         X(9,10-num)=zzz;
          
         x=num;
         y=mean(real(X(:,10-num)));
         err=std(real(X(:,10-num)));
         errorbar(10-x,y,err);
         grid on;
         grid minor;
         xlabel('Number of Cycles to perform the FFT')
         ylabel('Normalized Resistance, \theta','FontSize',12)
         axis([0.5 9.5 -3 5])
         
    subplot(2,1,2)
        hold on
         plot(10-num,imag(zetaLiner(j)),'<k')
         plot(10-num,imag(zLinerP1(j)),'>k')
         plot(10-num,imag(zLinerP2(j)),'*k')
         plot(10-num,imag(zLinerPM2(j)),'ok')
         plot(10-num,imag(zLinerPM1(j)),'^k')
         plot(10-num,imag(zLinerInlet(j)),'dk')
         plot(10-num,imag(zSurf1(j)),'vk')
         plot(10-num,imag(zLinerASTM(j)),'xk')
         legend('FFT(Liner)','FFT(Plane1)','FFT(Plane2)','FFT(Mic1)','FFT(Mic2)','IRW(Inlet)','FC(Mic2)','ASTM','Location','Best')
         
          x=num;
         y=mean(imag(X(:,10-num)));
         err=std(imag(X(:,10-num)));
         errorbar(10-x,y,err);
                  grid on;
         grid minor;
        xlabel('Number of Cycles to perform the FFT','FontSize',12)
         ylabel('Normalized Rectance, \chi','FontSize',12)
       axis([0.5 9.5 3 6])
end   
 %% Saving Impedance
     set(gcf, 'Position', get(0,'ScreenSize')); 
                 string = ['..\Figures\Impedance\ImpedanceConvergence_' num2str(SPL) 'dB_' num2str(freq(j)) 'Hz']
                saveas(h,[string '.eps'],'epsc2')
                saveas(h,[string '.fig'])
                saveas(h,[string '.png'])
                
                save Impedances.mat X
                
%    %% Save Struct
%                 
%                 field1='Liner';
%                 value1={'Slit liner properties';
%                    
%                 %Struct liner
%                 field1='Hole_diameter';
%                 value1=d_h;
%                 field2='Face_sheet_thickness'
%                 value2=h_slit;
%                 field3='Cavity_depth'
%                 value3=h_r;
%                 field4='POA'
%                 value4=POA_holder;
%                 liner_description=struct(field1,value1,field2,value2,field3,value3,field4,value4)
%                 
%                 %Struct constants
%                 field1='Temperature_celsius';
%                 value1=T;
%                 field2='Density'
%                 value2=rho;
%                 field3='Speed_of_sound'
%                 value3=c0;
%                 field4='Mic_1_position'
%                 value4=x1;
%                 field5='Mic_2_position'
%                 value5=x2;
%                 field6='Plane_1_position'
%                 value6=yPlane1;
%                 field7='Plane_2_position'
%                 value7=yPlane2;
%                 constants=struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7)
%                 
%                 %Struct frequency
%                 field1='Excitation_frequency'
%                 value1=freq;
%                 field2='Array_frequency_responses'
%                 value2=fs;
%                 
% 
%                 
%   ['SPL=' num2str(SPL) 'dB'];

%        psave(:,1)=time
%        psave(:,2)=pliner
%        
%        vsave(:,1)=time
%        vsave(:,2)=vliner
%        
%        pfsave(:,1)=fs
%        pfsave(:,2)=real(HP)
%        pfsave(:,3)=imag(HP)
%        
%        vfsave(:,1)=fs
%        vfsave(:,2)=real(HV)
%        vfsave(:,3)=imag(HV)
%        
%        save pressure_time.txt psave -ascii
%        save velocity_time.txt vsave -ascii
%        
%        save pressure_fft.txt pfsave -ascii
%        save velocity_fft.txt vfsave -ascii
% end