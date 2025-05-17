% Script to evaluate the " 18 m/s Hurdles" for the LAC Summer Games 2025.
% Task:
% Get the best possible wind preview (2s) by improving the baseline lidar 
% data processing (LDP_v3)!
% Results 4BeamPulsed:
% v3: Cost for Summer Games 2025 ("18 m/s hurdles"):  0.516059 m/s 
% Results CircularCW:
% v3: Cost for Summer Games 2025 ("18 m/s hurdles"):  0.463949 m/s 

%% Setup
clearvars;close all;clc;
addpath(genpath('..\WetiMatlabFunctions'))

% select LidarType
LidarType           = '4BeamPulsed'; % [4BeamPulsed/CircularCW]

% Seeds (can be adjusted, but will provide different results)
nSeed               = 6;                        % [-]	    number of stochastic turbulence field samples
Seed_vec            = [1:nSeed]+18*100;         % [-]  	    vector of seeds

% Parameters postprocessing (can be adjusted, but will provide different results)
t_start             = 60;                       % [s] 	    ignore data before for STD and spectra
TMax                = 659.75;                   % [s]       total run time (longest time in input files)
DT                  = 0.01;                     % [s]       time step
time                = [0:DT:TMax]';             % [s]       time vector
R                   = 120;                      % [m]  	    rotor radius to calculate REWS

% Parameter for Cost (Summer Games 2025)
tau                 = 2;                        % [s]       time to overcome pitch actuator, from Example 1: tau = T_Taylor - T_buffer, since there T_filter = T_scan = 0

switch LidarType
    case '4BeamPulsed'
        % configuration from LDP_v1_4BeamPulsed.IN and FFP_v1_4BeamPulsed.IN
        LDP.NumberOfBeams       = 4;            % [-]       Number of beams measuring at different directions               
        LDP.AngleToCenterline   = 19.176;       % [deg]     Angle around centerline
        LDP.IndexGate           = 6;            % [-]       IndexGate
        LDP.FlagLPF             = 1;            % [0/1]     Enable low-pass filter (flag)
        LDP.omega_cutoff        = 0.13;         % [rad/s]   Corner frequency (-3dB) of the low-pass filter
        LDP.T_buffer            = 0.2;          % [s]       Buffer time for filtered REWS signal
    case 'CircularCW'
        % configuration from LDP_v1_CircularCW.IN and FFP_v1_CircularCW.IN
        LDP.NumberOfBeams       = 50;           % [-]       Number of beams measuring at different directions               
        LDP.AngleToCenterline   = 15;           % [deg]     Angle around centerline
        LDP.IndexGate           = 1;            % [-]       IndexGate
        LDP.FlagLPF             = 1;            % [0/1]     Enable low-pass filter (flag)
        LDP.omega_cutoff        = 0.20;         % [rad/s]   Corner frequency (-3dB) of the low-pass filter
        LDP.T_buffer            = 4.2;          % [s]       Buffer time for filtered REWS signal        
end

% Files (should not be be changed)
SimulationFolderLAC     = 'solis_lidar_data';

%% Postprocessing: evaluate data

% Allocation
MAE     = NaN(1,nSeed); % mean absolute error [m/s]

% Loop over all seeds
for iSeed = 1:nSeed    

    % Load data
    Seed                                = Seed_vec(iSeed);
	WindFileName                        = ['URef_18_Seed_',num2str(Seed,'%02d')];
    SolisResultFile                     = fullfile(SimulationFolderLAC,[WindFileName,'_lidar_data_',LidarType,'.csv']);
    Data                                = readtable(SolisResultFile);    
    beamID                              = interp1(Data.time,Data.beamID,time,'previous');
    isValid                             = interp1(Data.time,Data.("isValid"+LDP.IndexGate),time,'previous');
    lineOfSightWindSpeed                = interp1(Data.time,Data.("lineOfSightWindSpeed"+LDP.IndexGate),time,'previous');

    % Get REWS from the wind field 
    TurbSimResultFile                 	= ['TurbulentWind\URef_18_Seed_',num2str(Seed,'%02d'),'.wnd'];   
    [REWS_WindField,Time_WindField]  	= CalculateREWSfromWindField(TurbSimResultFile,R,2);
               
    % Calculate REWS
    clear LDP_v3 % clearing all persistent variables from previous call
    [REWS,REWS_f,REWS_b]                = LDP_v3(time,isValid,beamID,lineOfSightWindSpeed,DT,LDP);

    % Calculate mean absolute error
    REWS_WindField_Fs                   = interp1(Time_WindField,REWS_WindField,time);
    REWS_WindField_Fs_shifted           = interp1(Time_WindField-tau,REWS_WindField,time);     % shift the REWS from wind field by tau (intented prediction time) into the future (lower times)
    Error                               = REWS_WindField_Fs_shifted-REWS_b;                    % error is  REWS from wind field shifted minus REWS from lidar filtered and buffered.
    MAE(iSeed)                          = mean(abs(detrend(Error(time>=t_start),'constant'))); % only consider error after t_start 
             
    % Plot REWS for absolute error             
    figure('Name',['REWS seed ',num2str(Seed)])
    subplot(311)
    hold on; grid on; box on
    plot(time,  REWS_WindField_Fs);
    plot(time,  REWS);
    ylabel('REWS [m/s]');
    legend('wind field','lidar estimate')
    subplot(312)
    hold on; grid on; box on
    plot(time,  REWS_WindField_Fs_shifted);
    plot(time,  REWS_b);
    ylabel('REWS [m/s]');
    legend('wind field shifted','lidar estimate filtered and buffered')
    subplot(313)
    hold on; grid on; box on
    plot(time,  Error);
    ylabel('error [m/s]');
    xlabel('time [s]')

end

%% Calculation of Cost for Summer Games 2025
Cost                        = mean(MAE);
fprintf('Cost for Summer Games 2025 ("18 m/s hurdles"):  %f \n',Cost);