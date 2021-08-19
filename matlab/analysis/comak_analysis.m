%% COMAK Analysis

h5_file = 'walking_meniscus.h5';
osim_file = '../../models/healthy/smith2019/smith2019.osim';
jam = jam_analysis(osim_file,h5_file);

%% Plot Muscles
% - Activation
% - Force

jam.forceset.plot_muscle_param({'addbrev_r'},'activation');
jam.forceset.plot_muscle_param_groups(jam.model,'activation');
%% Plot Ligaments
% - Strain
% - Force

%% Plot Reserve Actuators
% - Activation
% - Force

%% Plot Damping Forces

%% Secondary Coordinates
% - Values

%% Primary Coordinates
% - Values
% - Moments

%% Ground Reaction Forces