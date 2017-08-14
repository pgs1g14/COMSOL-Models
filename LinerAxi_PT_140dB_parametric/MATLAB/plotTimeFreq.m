function plotTimeFreq(freq,j,time,vliner,pliner,fs,HV,HP,position)
subplot(2,2,1)
        plot(time,vliner)
         hold on
         title(['Pure tone f=' num2str(freq(j)) ' Hz'])
         grid on
         grid minor
         xlabel('Time [s]')
         ylabel('Normal Average Velocity, v(t) [m/s]')
         xlim([time(length(time))-5/freq time(length(time))]) 
     subplot(2,2,2)
         semilogy(fs,abs(HV))
         hold on
         grid on
         grid minor
         xlabel('Frequency [Hz]')
         ylabel('Normal Average Velocity, abs(v(\omega)) [-]')
         axis([0 6400 0 max(abs(HV))+100])
     subplot(2,2,3)
          plot(time,pliner)
         hold on
         grid on
         grid minor
         xlabel('Time [s]')
         ylabel('Average Pressure, p(t) [Pa]') 
         xlim([time(length(time))-5/freq time(length(time))])  
     subplot(2,2,4)
         semilogy(fs,abs(HP))
         hold on
         grid on
         grid minor
         xlabel('Frequency [Hz]')
         ylabel('Average Pressure, abs(p(\omega)) [-]')
         axis([0 6400 0 max(abs(HP))+100000])
end