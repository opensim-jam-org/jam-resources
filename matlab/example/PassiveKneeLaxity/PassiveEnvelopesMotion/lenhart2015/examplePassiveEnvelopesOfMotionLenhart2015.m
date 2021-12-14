%% Passive Envelopes of Motion
%==========================================================================
% Author: Colin Smith 
%
%
%==========================================================================

%% Setup Environment and Folders
clear; close all;
import org.opensim.modeling.*
Logger.setLevelString('info');

useVisualizer = true;
model_file = '../../../../models/knee_healthy/smith2019/smith2019.osim';

secondary_coordinates = {...
    'knee_add_r',...
    'knee_rot_r',...;
    'knee_tx_r',...;
    'knee_ty_r',...;
    'knee_tz_r',...;
    'pf_flex_r',...;
    'pf_rot_r',...;
    'pf_tilt_r',...;
    'pf_tx_r',...;
    'pf_ty_r',...;
    'pf_tz_r',...;
    };

numSecondaryCoordinates = length(secondary_coordinates);
secondary_joints ={...
    'knee_r',...
    'knee_r',...;
    'knee_r',...;
    'knee_r',...;
    'knee_r',...;
    'pf_r',...;
    'pf_r',...;
    'pf_r',...;
    'pf_r',...;
    'pf_r',...;
    'pf_r',...;
    };

laxity_coordinates = {...
    'knee_add_r',...
    'knee_rot_r',...;
    'knee_tx_r',...;
    'knee_ty_r',...;
    'knee_tz_r'};

laxity_max_load =
%%
