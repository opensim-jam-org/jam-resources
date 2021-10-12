%% Setup Environment and Folders
clear; close all;
import org.opensim.modeling.*
Logger.setLevelString('debug');

useVisualizer = true;

model_file = '../../../models/knee_healthy/lenhart2015/lenhart2015.osim';
results_basename = 'passive_flexion';
forsim_result_dir = './results/forsim';
jnt_mech_result_dir = './results/joint-mechanics';

if(exist('./inputs','dir')~=7)
    mkdir('./inputs')
end
if(exist('./results','dir')~=7)
    mkdir('./results')
end
if(exist('./results/graphics','dir')~=7)
    mkdir('./results/graphics')
end


%% Create Input Files
prescribed_coord_file = './inputs/prescribed_coordinates.sto';

% Simulation Time
% Simulation consists of two phases:
%   settle : allow unconstraind knee DOFs to settle into equilbrium 
%   flex : prescribe the tibiofemoral flexion
% All time units are in seconds 

time_step = 0.01;

settle_duration = 0.5; 
flex_duration = 2.0; 

settle_time = 0 : time_step : settle_duration;
flex_time = settle_duration + time_step : time_step : flex_duration + settle_duration;

time = [settle_time, flex_time];

time_points = [0,settle_duration,settle_duration + flex_duration];

nSettleSteps = length(settle_time);
nFlexSteps = length(flex_time);

% Prescribe Knee Flexion
max_knee_flex = 90;

knee_flex = [0,0,max_knee_flex];
smooth_knee_flex = interp1(time_points, knee_flex, time,'pchip');

prescribed_coord_data.knee_flex_r = smooth_knee_flex';
prescribed_coord_data.time = time;
prescribed_coord_data.pelvis_tilt = ones(nSettleSteps+nFlexSteps,1)*90;

prescribed_coord_table = osimTableFromStruct(prescribed_coord_data);

sto = STOFileAdapter();
sto.write(prescribed_coord_table,prescribed_coord_file);

%% Plot Simulation Inputs 
coord_fig = figure('name','prescribed_coordinates','Position',  [100, 100, 667, 300]);

subplot(1,2,1);
plot(time,prescribed_coord_data.pelvis_tilt,'LineWidth',2)
ylim([0.0 100])
xlabel('Time [s]')
ylabel('Angle [^o]')
title('Pelvis Tilt (pelvis\_tilt\_r)')
box off

subplot(1,2,2);
plot(time,prescribed_coord_data.knee_flex_r,'LineWidth',2)
ylim([0.0 100])
xlabel('Time [s]')
ylabel('Angle [^o]')
title('Knee Flexion (knee\_flex\_r)')
box off

saveas(coord_fig,'./results/graphics/prescribed_coordinates.png')

%% Perform Simulation with ForsimTool
forsim = ForsimTool();
forsim.set_model_file(model_file);
forsim.set_results_directory(forsim_result_dir);
forsim.set_results_file_basename(results_basename);
forsim.set_start_time(-1);
forsim.set_stop_time(-1);
forsim.set_integrator_accuracy(1e-2); %Note this should be 1e-6 for research
forsim.set_constant_muscle_control(0.02); %Set all muscles to 2% activation to represent passive state
forsim.set_ignore_activation_dynamics(false);
forsim.set_ignore_tendon_compliance(false);
forsim.set_unconstrained_coordinates(0,'/jointset/knee_r/knee_add_r');
forsim.set_unconstrained_coordinates(1,'/jointset/knee_r/knee_rot_r');
forsim.set_unconstrained_coordinates(2,'/jointset/knee_r/knee_tx_r');
forsim.set_unconstrained_coordinates(3,'/jointset/knee_r/knee_ty_r');
forsim.set_unconstrained_coordinates(4,'/jointset/knee_r/knee_tz_r');
forsim.set_unconstrained_coordinates(5,'/jointset/pf_r/pf_flex_r');
forsim.set_unconstrained_coordinates(6,'/jointset/pf_r/pf_rot_r');
forsim.set_unconstrained_coordinates(7,'/jointset/pf_r/pf_tilt_r');
forsim.set_unconstrained_coordinates(8,'/jointset/pf_r/pf_tx_r');
forsim.set_unconstrained_coordinates(9,'/jointset/pf_r/pf_ty_r');
forsim.set_unconstrained_coordinates(10,'/jointset/pf_r/pf_tz_r');
forsim.set_prescribed_coordinates_file(prescribed_coord_file);
forsim.set_use_visualizer(useVisualizer);
forsim.print('./inputs/forsim_settings.xml');

disp('Running Forsim Tool...')
forsim.run();
%% Perform Analysis with JointMechanicsTool

jnt_mech = JointMechanicsTool();
jnt_mech.set_model_file(model_file);
jnt_mech.set_input_states_file([forsim_result_dir '/' results_basename '_states.sto']);
jnt_mech.set_results_file_basename(results_basename);
jnt_mech.set_results_directory(jnt_mech_result_dir);
jnt_mech.set_start_time(-1);
jnt_mech.set_stop_time(-1);
jnt_mech.set_normalize_to_cycle(false);
jnt_mech.set_contacts(0,'all');
jnt_mech.set_ligaments(0,'all');
jnt_mech.set_muscles(0,'none');
%jmt.set_muscle_outputs(0,'none');
jnt_mech.set_attached_geometry_bodies(0,'/bodyset/femur_distal_r');
jnt_mech.set_attached_geometry_bodies(1,'/bodyset/tibia_proximal_r');
jnt_mech.set_attached_geometry_bodies(2,'/bodyset/patella_r');
jnt_mech.set_output_orientation_frame('ground');
jnt_mech.set_output_position_frame('ground');
jnt_mech.set_write_vtp_files(true);
jnt_mech.set_write_h5_file(true);
jnt_mech.set_h5_kinematics_data(true);
jnt_mech.set_h5_states_data(true);
jnt_mech.set_use_visualizer(useVisualizer);

jnt_mech.print('./inputs/joint_mechanics_settings.xml');

disp('Running JointMechanicsTool...');
% jnt_mech.run();

%% Analyze Simulation Results
results = jam_analysis('lenhart2015',{[jnt_mech_result_dir '/' results_basename '.h5']});

% Secondary Kinematics
sec_coords = {...
    'knee_flex_r','knee_add_r','knee_rot_r',...
    'knee_tx_r','knee_ty_r','knee_tz_r',...
    'pf_flex_r','pf_rot_r','pf_tilt_r',...
    'pf_tx_r','pf_ty_r','pf_tz_r'
    };

figure('name','Secondary Coordinates')
for i = 1:length(sec_coords)
    subplot(4,3,i)
    hold on
    plot(results.coordinateset.(sec_coords{i}).value);    
end

% Ligament Forces
ligament_names = {...
    'MCLd','MCLs'...
    'LCL',...
    'ACLam','ACLpl',...
    'PCLal','PCLpm',...
    'PT',...
    'ITB',...
    
    };

figure('name','Ligament Forces')
for i = 1:length(ligament_names)
    subplot(3,3,i);hold on;
    
    fiber_names = fieldnames(results.forceset.Blankevoort1991Ligament);
    
    fibers = fiber_names(contains(fiber_names,ligament_names{i}));
    
    data = 0;
    for k = 1:length(fibers)
        data = data + results.forceset.Blankevoort1991Ligament.(fibers{k}).total_force;

    end
    plot(data)
    title([ligament_names{i} ' force'])
    
end

% Muscle Forces
% muscle_names = {...
%     'vasint_r',...
%     'vaslat_r',...
%     'vasmed_r',...
%     'recfem_r',...
%     'gaslat_r',...
%     'gasmed_r',...
%     'bflh_r',...
%     'bfsh_r',...
%     'semimem_r',...
%     'semiten_r',...
%     'sart_r',...
%     'tfl_r'
%     };
% 
% for i = 1:length(muscle_names)
%     subplot(3,4,i);hold on;
%     
%     muscle_names = fieldnames(results.forceset.Muscle);
%     
%     msls = fiber_names(contains(muscle_names,ligament_names{i}));
%     
%     data = 0;
%     for k = 1:length(msls)
%         data = data + results.forceset.Blankevoort1991Ligament.(msls{k}).actuation;
% 
%     end
%     plot(data)
%     title([ligament_names{i} ' force'])
%     
% end