%% Run the parametric model choosing the start stop and frequency step

 freq=2000
cycles=5
    LinerAxi_PT_155dB_parametric_v5_2(cycles,freq) %Run the model and extract results
    DataPath=['..\Data\Results_' num2str(freq) 'Hz.mat'];
    eval(['data=load(DataPath)'])
    
    for cycles=34:38
        ImpedanceCalc(cycles,freq,load(DataPath))
    end
