%% exampleValidationHealthyKnee
%==========================================================================
% Compare measured tibiofemoral and patellofemoral kinematics from 
% dynamic MRI against simulation predictions for passive and two 
% types of loading (inertia & elastic). 
%

%==========================================================================
%% Setup Environment and Folders
clear; close all;
import org.opensim.modeling.*
Logger.setLevelString('info');

useVisualizer = true;

model_file = '../../../../models/knee_healthy/smith2019/smith2019.osim';
% experimental_data_dir = ...
%     '../../../../models/knee_tka/grand_challenge/DM/experimental_data/motion_analysis/';

trial_names = {'DM_ngait_tmf_slow1','DM_ngait_tmf_slow2'};
num_trials = length(trial_names);

inputs.min_knee_flex = 5;
inputs.max_knee_flex = 40;

if(exist('./inputs','dir')~=7); mkdir('./inputs'); end
if(exist('./results','dir')~=7); mkdir('./results'); end
if(exist('./results/passive','dir')~=7); mkdir('./results/passive'); end
if(exist('./results/inertia','dir')~=7); mkdir('./results/inertia'); end
if(exist('./results/elastic','dir')~=7); mkdir('./results/elastic'); end
if(exist(['./results/graphics'],'dir')~=7)
    mkdir(['./results/graphics'])
end

%% Time
% 0.5 Hz, 2 seconds for 1 cycle
time.step = 0.01;
time.settle_duration = 1.0;
time.cycle_duration = 2.0;
time.duration = time.settle_duration + time.cycle_duration + time.cycle_duration;
time.settle_time = 0 : time.step : time.settle_duration;
time.cycle1_time = time.settle_duration + time.step : time.step : time.settle_duration + time.cycle_duration;
time.cycle2_time = time.settle_duration + time.cycle_duration + time.step : time.step : time.settle_duration + time.cycle_duration + time.cycle_duration;
time.values = [time.settle_time, time.cycle1_time, time.cycle2_time];
time.values_pad = [time.settle_time, time.cycle1_time, time.cycle2_time,time.settle_time];

time_points = [0,time.settle_duration,...
    time.settle_duration + time.cycle_duration,...
    time.settle_duration + time.cycle_duration + time.cycle_duration];

time.num_settle_steps = length(time.settle_time);
time.num_cycle1_steps = length(time.cycle1_time);
time.num_cycle2_steps = length(time.cycle2_time);
time.num_steps = length(time.values);


%% Prescribed Coordinates File
prescribed_coordinates_file = './inputs/prescribed_coordinates.sto';

% inputs.knee_flex = [inputs.min_knee_flex,inputs.min_knee_flex,...
%     inputs.max_knee_flex ,inputs.min_knee_flex,inputs.max_knee_flex,...
%     inputs.min_knee_flex,inputs.min_knee_flex]; 
% 
% time_points2 = linspace(0,time.duration+time.settle_duration,length(inputs.knee_flex));
% inputs.smooth_knee_flex = interp1(time_points2, inputs.knee_flex, time.values_pad,'pchip');
% inputs.smooth_knee_flex = inputs.smooth_knee_flex(1:end-length(time.settle_time));
coord_data.pelvis_ty = ones(time.num_steps,1) * 1;
coord_data.pelvis_tilt = ones(time.num_steps,1) * 90;
coord_data.hip_flex_r = ones(time.num_steps,1) * 20;
% coord_data.knee_flex_r = inputs.smooth_knee_flex';
f=0.5;
tp=[0.0:.01:4]';
coord_data.knee_flex_r = [linspace(0,inputs.min_knee_flex,50)'; inputs.min_knee_flex*ones(50,1); inputs.min_knee_flex+inputs.max_knee_flex/2*(1+cos(f*2*pi*tp+pi))];
coord_data.ankle_flex_r = ones(time.num_steps,1) * -30;
coord_data.time = time.values;

coord_table = osimTableFromStruct(coord_data); %% Function distributed in OpenSim 4.0 resources
coord_table.addTableMetaDataString('inDegrees','yes');
STOFileAdapter.write(coord_table,prescribed_coordinates_file);

%% Primary Coordinates (inertia, elastic)

%% External Loads



%**********************************************
% Elastic and Inertial Loadings
 % frequency in Hz
% t=[0.0:.01:5]';
% f=0.5;
% tp=[0.0:.01:4]';
% kf=[5*ones(100,1); 5+20*(1+cos(f*2*pi*tp+pi))];
% n=length(t);
% z=zeros(n,1);
% o=ones(n,1);
% af=-20*o;
% hdr={'time';'knee_flex_r';'ankle_flex_r';'tibia_r_force_px';'tibia_r_force_py';'tibia_r_force_pz';'tibia_r_force_vx';'tibia_r_force_vy';'tibia_r_force_vz';'tibia_r_torque_x';'tibia_r_torque_y';'tibia_r_torque_z'};
% 
% % 
% % % No-load flexion-extension at same frequency
% F=[zeros(100,1); 0.*cos(f*2*pi*tp)/.35];
% data=[t kf af z -D*o Zoffset*o F z z z z z];
% done=write_motion('kfec_halfhz.mot',data,hdr,'cyclic knee flexion');
% 
% % Elastic load flexion-extension

% data=[t kf af z -D*o Zoffset*o F z z z z z];
% done=write_motion('kfec_elastic.mot',data,hdr,'cyclic knee flexion');
% 
% % Inertial load flexion-extension

% data=[t kf af z -D*o Zoffset*o F z z z z z];
% done=write_motion('kfec_inertia.mot',data,hdr,'cyclic knee flexion');


%% Create External Loads (inertia and elastic)
inputs.elastic_magnitude = 17;
inputs.inertia_magnitude = 15;
inputs.cycle_frequency = 0.5; %Hz

EMAG=inputs.elastic_magnitude;
IMAG=inputs.inertia_magnitude;
freq=0.5;
tp=[0.0:.01:4]';

inputs.elastic_force = [zeros(100,1); -EMAG*cos(freq*2*pi*tp)/.35-3/.35];
inputs.inertia_force = [zeros(100,1); IMAG*cos(freq*2*pi*tp)/.35-3/.35];

elastic_load_data.time = time.values';
elastic_load_data.external_force_vx = inputs.elastic_force;
elastic_load_data.external_force_vy = zeros(time.num_steps,1);
elastic_load_data.external_force_vz = zeros(time.num_steps,1);
elastic_load_data.external_force_px = zeros(time.num_steps,1);
elastic_load_data.external_force_py = -0.35 * ones(time.num_steps,1);
elastic_load_data.external_force_pz = zeros(time.num_steps,1);
elastic_load_data.external_torque_x = zeros(time.num_steps,1);
elastic_load_data.external_torque_y = zeros(time.num_steps,1);
elastic_load_data.external_torque_z = zeros(time.num_steps,1);

elastic_load_table = osimTableFromStruct(elastic_load_data); % Function distributed in OpenSim 4.0 resources
elastic_load_sto_file = './inputs/elastic_load.sto';
STOFileAdapter.write(elastic_load_table,elastic_load_sto_file);

elastic_external_load = ExternalLoads();
elastic_ext_force = ExternalForce();
elastic_ext_force.setName('elastic_force');
elastic_ext_force.set_applied_to_body('tibia_r');
elastic_ext_force.set_force_expressed_in_body('tibia_r');
elastic_ext_force.set_point_expressed_in_body('tibia_r');
elastic_ext_force.set_force_identifier('external_force_v');
elastic_ext_force.set_point_identifier('external_force_p');
elastic_ext_force.set_torque_identifier('external_torque_');
elastic_external_load.adoptAndAppend(elastic_ext_force);
elastic_external_load.setDataFileName('elastic_load.sto');
elastic_external_loads_xml_file = './inputs/elastic_external_load.xml';
elastic_external_load.print(elastic_external_loads_xml_file);



inertia_load_data.time = time.values';
inertia_load_data.external_force_vx = inputs.inertia_force;
inertia_load_data.external_force_vy = zeros(time.num_steps,1);
inertia_load_data.external_force_vz = zeros(time.num_steps,1);
inertia_load_data.external_force_px = zeros(time.num_steps,1);
inertia_load_data.external_force_py = -0.35 * ones(time.num_steps,1);
inertia_load_data.external_force_pz = zeros(time.num_steps,1);
inertia_load_data.external_torque_x = zeros(time.num_steps,1);
inertia_load_data.external_torque_y = zeros(time.num_steps,1);
inertia_load_data.external_torque_z = zeros(time.num_steps,1);

inertia_load_table = osimTableFromStruct(inertia_load_data); % Function distributed in OpenSim 4.0 resources
inertia_load_sto_file = './inputs/inertia_load.sto';
STOFileAdapter.write(inertia_load_table,inertia_load_sto_file);

inertia_external_load = ExternalLoads();
inertia_ext_force = ExternalForce();
inertia_ext_force.setName('inertia_force');
inertia_ext_force.set_applied_to_body('tibia_r');
inertia_ext_force.set_force_expressed_in_body('tibia_r');
inertia_ext_force.set_point_expressed_in_body('tibia_r');
inertia_ext_force.set_force_identifier('external_force_v');
inertia_ext_force.set_point_identifier('external_force_p');
inertia_ext_force.set_torque_identifier('external_torque_');
inertia_external_load.adoptAndAppend(inertia_ext_force);
inertia_external_load.setDataFileName('inertia_load.sto');
inertia_external_loads_xml_file = './inputs/inertia_external_load.xml';
inertia_external_load.print(inertia_external_loads_xml_file);

%% Plot Inputs
figure('name','Dynamic MRI Validation Inputs')
subplot(1,3,1); hold on;
plot(time.values,coord_data.knee_flex_r)
title('Knee Flexion')
xlabel('Time [s]') 
ylabel('Angle [^o]')

subplot(1,3,2); hold on;
plot(time.values,inputs.inertia_force)
title('Inertia Load')
xlabel('Time [s]') 
ylabel('Force [N]')

subplot(1,3,3); hold on;
plot(time.values,inputs.elastic_force)
title('Elastic Load')
xlabel('Time [s]') 
ylabel('Force [N]')

%% Passive 
useActivationDynamics = true;
useTendonCompliance = true;
useMusclePhysiology = true;
integratorAccuracy = 1e-2;

passive_forsim_result_dir = './results/passive/forsim';
passive_joint_mechanics_result_dir = './results/passive/joint-mechanics';
basename = 'passive';

passive_forsim = ForsimTool();
passive_forsim.set_model_file(model_file);
passive_forsim.set_results_directory(passive_forsim_result_dir);
passive_forsim.set_results_file_basename(basename);
passive_forsim.set_start_time(-1);
passive_forsim.set_stop_time(-1);
passive_forsim.set_integrator_accuracy(integratorAccuracy); %Note this should be 1e-6 for research
passive_forsim.set_constant_muscle_control(0.02); %Set all muscles to 2% activation to represent passive state
passive_forsim.set_use_activation_dynamics(useActivationDynamics);
passive_forsim.set_use_tendon_compliance(useTendonCompliance);
passive_forsim.set_use_muscle_physiology(useMusclePhysiology)
passive_forsim.set_unconstrained_coordinates(0,'/jointset/knee_r/knee_add_r');
passive_forsim.set_unconstrained_coordinates(1,'/jointset/knee_r/knee_rot_r');
passive_forsim.set_unconstrained_coordinates(2,'/jointset/knee_r/knee_tx_r');
passive_forsim.set_unconstrained_coordinates(3,'/jointset/knee_r/knee_ty_r');
passive_forsim.set_unconstrained_coordinates(4,'/jointset/knee_r/knee_tz_r');
passive_forsim.set_unconstrained_coordinates(5,'/jointset/pf_r/pf_flex_r');
passive_forsim.set_unconstrained_coordinates(6,'/jointset/pf_r/pf_rot_r');
passive_forsim.set_unconstrained_coordinates(7,'/jointset/pf_r/pf_tilt_r');
passive_forsim.set_unconstrained_coordinates(8,'/jointset/pf_r/pf_tx_r');
passive_forsim.set_unconstrained_coordinates(9,'/jointset/pf_r/pf_ty_r');
passive_forsim.set_unconstrained_coordinates(10,'/jointset/pf_r/pf_tz_r');
passive_forsim.set_prescribed_coordinates_file(prescribed_coordinates_file);
passive_forsim.set_use_visualizer(useVisualizer);
passive_forsim.print('./inputs/forsim_settings.xml');

disp('Running Forsim Tool...')
% passive_forsim.run();

% Joint Mechanics
passive_jnt_mech = JointMechanicsTool();
passive_jnt_mech.set_model_file(model_file);
passive_jnt_mech.set_input_states_file([passive_forsim_result_dir '/' basename '_states.sto']);
passive_jnt_mech.set_use_muscle_physiology(useMusclePhysiology);
passive_jnt_mech.set_results_file_basename(basename);
passive_jnt_mech.set_results_directory(passive_joint_mechanics_result_dir);
passive_jnt_mech.set_start_time(-1);
passive_jnt_mech.set_stop_time(-1);
passive_jnt_mech.set_resample_step_size(-1);
passive_jnt_mech.set_normalize_to_cycle(true);
passive_jnt_mech.set_lowpass_filter_frequency(-1);
passive_jnt_mech.set_contacts(0,'all');
passive_jnt_mech.set_contact_outputs(0,'all');
passive_jnt_mech.set_contact_mesh_properties(0,'none');
passive_jnt_mech.set_ligaments(0,'all');
passive_jnt_mech.set_ligament_outputs(0,'all');
passive_jnt_mech.set_muscles(0,'all');
passive_jnt_mech.set_muscle_outputs(0,'all');
passive_jnt_mech.set_attached_geometry_bodies(0,'all');
passive_jnt_mech.set_output_orientation_frame('ground');
passive_jnt_mech.set_output_position_frame('ground');
passive_jnt_mech.set_write_vtp_files(false);
passive_jnt_mech.set_vtp_file_format('binary');
passive_jnt_mech.set_write_h5_file(true);
passive_jnt_mech.set_h5_kinematics_data(true);
passive_jnt_mech.set_h5_states_data(false);
passive_jnt_mech.set_use_visualizer(useVisualizer);

passive_jnt_mech.print('./inputs/joint_mechanics_settings.xml');

disp('Running JointMechanicsTool...');
% passive_jnt_mech.run();

%% Setup Muscle Driven Simulations
% All coordinates not listed in prescribed_coordinates file or
% below are prescribed to their default values

prescribed_coordinates = {...
    '/jointset/gnd_pelvis/pelvis_tilt',...
    '/jointset/hip_r/hip_flex_r'};

primary_coordinates = {...
    '/jointset/knee_r/knee_flex_r',...
    '/jointset/ankle_r/ankle_flex_r'};

secondary_coordinates(1).name = 'knee_add_r';
secondary_coordinates(1).max_change = 0.005;
secondary_coordinates(1).coordinate = '/jointset/knee_r/knee_add_r';

secondary_coordinates(2).name = 'knee_rot_r';
secondary_coordinates(2).max_change = 0.005;
secondary_coordinates(2).coordinate = '/jointset/knee_r/knee_rot_r';

secondary_coordinates(3).name = 'knee_tx_r';
secondary_coordinates(3).max_change = 0.001;
secondary_coordinates(3).coordinate = '/jointset/knee_r/knee_tx_r';

secondary_coordinates(4).name = 'knee_ty_r';
secondary_coordinates(4).max_change = 0.001;
secondary_coordinates(4).coordinate = '/jointset/knee_r/knee_ty_r';

secondary_coordinates(5).name = 'knee_tz_r';
secondary_coordinates(5).max_change = 0.001;
secondary_coordinates(5).coordinate = '/jointset/knee_r/knee_tz_r';

secondary_coordinates(6).name = 'pf_flex_r';
secondary_coordinates(6).max_change = 0.005;
secondary_coordinates(6).coordinate = '/jointset/pf_r/pf_flex_r';

secondary_coordinates(7).name = 'pf_rot_r';
secondary_coordinates(7).max_change = 0.005;
secondary_coordinates(7).coordinate = '/jointset/pf_r/pf_rot_r';

secondary_coordinates(8).name = 'pf_tilt_r';
secondary_coordinates(8).max_change = 0.005;
secondary_coordinates(8).coordinate = '/jointset/pf_r/pf_tilt_r';

secondary_coordinates(9).name = 'pf_tx_r';
secondary_coordinates(9).max_change = 0.001;
secondary_coordinates(9).coordinate = '/jointset/pf_r/pf_tx_r';

secondary_coordinates(10).name = 'pf_ty_r';
secondary_coordinates(10).max_change = 0.001;
secondary_coordinates(10).coordinate = '/jointset/pf_r/pf_ty_r';

secondary_coordinates(11).name = 'pf_tz_r';
secondary_coordinates(11).max_change = 0.001;
secondary_coordinates(11).coordinate = '/jointset/pf_r/pf_tz_r';

%% Elastic
useMusclePhysiology = false;

elastic_comak_result_dir = './results/elastic/comak';
elastic_joint_mechanics_result_dir = './results/elastic/joint-mechanics';
basename = 'elastic';

elastic_comak = COMAKTool();

elastic_comak.set_use_muscle_physiology(useMusclePhysiology);
elastic_comak.set_model_file(model_file);
elastic_comak.set_coordinates_file(prescribed_coordinates_file);
elastic_comak.set_external_loads_file(elastic_external_loads_xml_file);
elastic_comak.set_results_directory(elastic_comak_result_dir);
elastic_comak.set_results_prefix(basename);
elastic_comak.set_replace_force_set(false);
elastic_comak.set_force_set_file('../../../../models/knee_healthy/smith2019/smith2019_reserve_actuators.xml');
elastic_comak.set_start_time(-1);
elastic_comak.set_stop_time(-1);
elastic_comak.set_time_step(0.01);

for j = 1:length(prescribed_coordinates)
    elastic_comak.set_prescribed_coordinates(j-1,prescribed_coordinates{j});
end

for j = 1:length(primary_coordinates)
    elastic_comak.set_primary_coordinates(j-1,primary_coordinates{j});
end

secondary_coord_set = COMAKSecondaryCoordinateSet();        
for j = 1:length(secondary_coordinates)
    comak_secondary_coord = COMAKSecondaryCoordinate();
    comak_secondary_coord.setName(secondary_coordinates(j).name);
    comak_secondary_coord.set_max_change(secondary_coordinates(j).max_change);
    comak_secondary_coord.set_coordinate(secondary_coordinates(j).coordinate);
    secondary_coord_set.cloneAndAppend(comak_secondary_coord);
end
elastic_comak.set_COMAKSecondaryCoordinateSet(secondary_coord_set);

elastic_comak.set_settle_secondary_coordinates_at_start(true);
elastic_comak.set_settle_threshold(1e-3);
elastic_comak.set_settle_accuracy(1e-2);
elastic_comak.set_settle_internal_step_limit(10000);
elastic_comak.set_print_settle_sim_results(true);
elastic_comak.set_settle_sim_results_directory(elastic_comak_result_dir);
elastic_comak.set_settle_sim_results_prefix('elastic_settle_sim');
elastic_comak.set_max_iterations(25);
elastic_comak.set_udot_tolerance(1);
elastic_comak.set_udot_worse_case_tolerance(50);
elastic_comak.set_unit_udot_epsilon(1e-6);
elastic_comak.set_optimization_scale_delta_coord(1);
elastic_comak.set_ipopt_diagnostics_level(3);
elastic_comak.set_ipopt_max_iterations(500);
elastic_comak.set_ipopt_convergence_tolerance(1e-4);
elastic_comak.set_ipopt_constraint_tolerance(1e-4);
elastic_comak.set_ipopt_limited_memory_history(200);
elastic_comak.set_ipopt_nlp_scaling_max_gradient(10000);
elastic_comak.set_ipopt_nlp_scaling_min_value(1e-8);
elastic_comak.set_ipopt_obj_scaling_factor(1);
elastic_comak.set_activation_exponent(2);
elastic_comak.set_contact_energy_weight(0);
elastic_comak.set_non_muscle_actuator_weight(1000);
elastic_comak.set_model_assembly_accuracy(1e-12);
elastic_comak.set_use_visualizer(useVisualizer);
elastic_comak.print('./inputs/comak_settings.xml');
disp('Running COMAK Tool...')
% elastic_comak.run();

% Joint Mechanics
elastic_jnt_mech = JointMechanicsTool();
elastic_jnt_mech.set_model_file(model_file);
elastic_jnt_mech.set_input_states_file([elastic_comak_result_dir '/' basename '_states.sto']);
elastic_jnt_mech.set_use_muscle_physiology(useMusclePhysiology);
elastic_jnt_mech.set_results_file_basename(basename);
elastic_jnt_mech.set_results_directory(elastic_joint_mechanics_result_dir);
elastic_jnt_mech.set_start_time(-1);
elastic_jnt_mech.set_stop_time(-1);
elastic_jnt_mech.set_resample_step_size(-1);
elastic_jnt_mech.set_normalize_to_cycle(true);
elastic_jnt_mech.set_lowpass_filter_frequency(-1);
elastic_jnt_mech.set_contacts(0,'all');
elastic_jnt_mech.set_contact_outputs(0,'all');
elastic_jnt_mech.set_contact_mesh_properties(0,'none');
elastic_jnt_mech.set_ligaments(0,'all');
elastic_jnt_mech.set_ligament_outputs(0,'all');
elastic_jnt_mech.set_muscles(0,'all');
elastic_jnt_mech.set_muscle_outputs(0,'all');
elastic_jnt_mech.set_attached_geometry_bodies(0,'all');
elastic_jnt_mech.set_output_orientation_frame('ground');
elastic_jnt_mech.set_output_position_frame('ground');
elastic_jnt_mech.set_write_vtp_files(false);
elastic_jnt_mech.set_vtp_file_format('binary');
elastic_jnt_mech.set_write_h5_file(true);
elastic_jnt_mech.set_h5_kinematics_data(true);
elastic_jnt_mech.set_h5_states_data(false);
elastic_jnt_mech.set_use_visualizer(false);

elastic_jnt_mech.print('./inputs/joint_mechanics_settings.xml');

disp('Running JointMechanicsTool...');
% elastic_jnt_mech.run();


%% Inertia Simulations
inertia_comak_result_dir = './results/inertia/comak';
inertia_joint_mechanics_result_dir = './results/inertia/joint-mechanics';
basename = 'inertia';

inertia_comak = elastic_comak;

inertia_comak.set_external_loads_file(inertia_external_loads_xml_file);
inertia_comak.set_results_directory(inertia_comak_result_dir);
inertia_comak.set_results_prefix(basename);
disp('Running COMAK Tool...')
% inertia_comak.run();

inertia_jnt_mech = elastic_jnt_mech;
inertia_jnt_mech.set_input_states_file([inertia_comak_result_dir '/' basename '_states.sto']);
inertia_jnt_mech.set_use_muscle_physiology(useMusclePhysiology);
inertia_jnt_mech.set_results_file_basename(basename);
inertia_jnt_mech.set_results_directory(inertia_joint_mechanics_result_dir);
disp('Running JointMechanicsTool...');
% inertia_jnt_mech.run();

%% Plot Results
h5_file_list = {...
    [passive_joint_mechanics_result_dir '/passive.h5'],...
    [elastic_joint_mechanics_result_dir '/elastic.h5'],...
    [inertia_joint_mechanics_result_dir '/inertia.h5']};
sim_names = {'passive','elastic','inertia'};

jam = jam_analysis(h5_file_list);

line_width = 2;

% Plot TF Kinematics
tf_coords = {'knee_flex_r','knee_add_r','knee_rot_r','knee_tx_r','knee_ty_r','knee_tz_r'};

figure('name','Tibiofemoral Kinematics')
for i = 1:length(tf_coords)
    subplot(2,3,i);hold on;
    for n = 1:jam.num_files
        plot(jam.time,jam.coordinateset.(tf_coords{i}).value(:,n),...
            'LineWidth',line_width)
    end
    xlabel('Time [s]')
    if (i < 4)
        ylabel('Angle [^o]')
    else
        ylabel('Translation [mm]')
    end  
    title(tf_coords{i})
end
legend(sim_names)

% Plot PF Kinematics
pf_coords = {'pf_flex_r','pf_rot_r','pf_tilt_r','pf_tx_r','pf_ty_r','pf_tz_r'};

figure('name','Patellofemoral Kinematics')
for i = 1:length(pf_coords)
    subplot(2,3,i);hold on;
    for n = 1:jam.num_files
        plot(jam.time,jam.coordinateset.(pf_coords{i}).value(:,n),...
            'LineWidth',line_width)
    end
    xlabel('Time [s]')
    if (i < 4)
        ylabel('Angle [^o]')
    else
        ylabel('Translation [mm]')
    end  
    title(pf_coords{i})
end
legend(sim_names)