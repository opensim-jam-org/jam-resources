%% exampleWalkingTKA
%==========================================================================

%% Setup Environment and Folders
clear;
import org.opensim.modeling.*
Logger.setLevelString('info');

useVisualizer = true;

model_file = '../../../models/knee_tka/grand_challenge/DM/DM.osim';
%trial = 'DM_smooth1';
trial = 'DM_ngait_og9';

if(strcmp(trial,'DM_smooth1'))
    start_time = 2.53;
    stop_time = 3.775;
elseif(strcmp(trial,'DM_bouncy5'))
    start_time = 1.876;
    stop_time = 3.075; 
elseif(strcmp(trial,'DM_ngait_og1'))
    start_time = 4.928;
    stop_time = 6.283; 
elseif(strcmp(trial,'DM_ngait_og3'))
    start_time = 2.478;
    stop_time = 3.737; 
elseif(strcmp(trial,'DM_ngait_og4'))
    start_time = 3.081;
    stop_time = 4.288; 
elseif(strcmp(trial,'DM_ngait_og5'))
    start_time = 2.188;
    stop_time = 3.428; 
elseif(strcmp(trial,'DM_ngait_og6'))
    start_time = 2.976;
    stop_time = 4.183; 
elseif(strcmp(trial,'DM_ngait_og7'))
    start_time = 1.136;
    stop_time = 2.337;
elseif(strcmp(trial,'DM_ngait_og9'))
    start_time = 3.214;
    stop_time = 4.45;
end

results_basename = trial;
ik_result_dir = ['./results/' trial '/comak-inverse-kinematics'];
comak_result_dir = ['./results/' trial '/comak'];
jnt_mech_result_dir = ['./results/' trial '/joint-mechanics'];
experimental_data_dir = ...
    '../../../models/knee_tka/grand_challenge/DM/experimental_data/motion_analysis';

marker_file = ...
    [experimental_data_dir '/' trial '_new.trc'];

ik_output_motion_file = [trial '_ik.sto'];

external_loads_file = ...
    [experimental_data_dir '/' trial '_external_loads.xml'];

if(exist('./inputs','dir')~=7)
    mkdir('./inputs')
end
if(exist(['./inputs/' trial],'dir')~=7)
    mkdir(['./inputs/' trial])
end
if(exist('./results','dir')~=7)
    mkdir('./results')
end
if(exist(['./results/' trial],'dir')~=7)
    mkdir(['./results/' trial])
end
if(exist(['./results/' trial '/graphics'],'dir')~=7)
    mkdir(['./results/' trial '/graphics'])
end

%% Perform Inverse Kinematics
comak_ik = COMAKInverseKinematicsTool();
comak_ik.set_model_file(model_file);
comak_ik.set_results_directory(ik_result_dir);
comak_ik.set_results_prefix(results_basename);
comak_ik.set_perform_secondary_constraint_sim(false);
comak_ik.set_secondary_coordinates(0,'/jointset/knee_r/knee_add_r');
comak_ik.set_secondary_coordinates(1,'/jointset/knee_r/knee_rot_r');
comak_ik.set_secondary_coordinates(2,'/jointset/knee_r/knee_tx_r');
comak_ik.set_secondary_coordinates(3,'/jointset/knee_r/knee_ty_r');
comak_ik.set_secondary_coordinates(4,'/jointset/knee_r/knee_tz_r');
comak_ik.set_secondary_coordinates(5,'/jointset/pf_r/pf_flex_r');
comak_ik.set_secondary_coordinates(6,'/jointset/pf_r/pf_rot_r');
comak_ik.set_secondary_coordinates(7,'/jointset/pf_r/pf_tilt_r');
comak_ik.set_secondary_coordinates(8,'/jointset/pf_r/pf_tx_r');
comak_ik.set_secondary_coordinates(9,'/jointset/pf_r/pf_ty_r');
comak_ik.set_secondary_coordinates(10,'/jointset/pf_r/pf_tz_r');
comak_ik.set_secondary_coupled_coordinate('/jointset/knee_r/knee_flex_r');
comak_ik.set_secondary_constraint_sim_settle_threshold(1e-4);
comak_ik.set_secondary_constraint_sim_sweep_time(3.0);
comak_ik.set_secondary_coupled_coordinate_start_value(0);
comak_ik.set_secondary_coupled_coordinate_stop_value(100);
comak_ik.set_secondary_constraint_sim_integrator_accuracy(1e-3);
comak_ik.set_secondary_constraint_sim_internal_step_limit(10000);
secondary_constraint_function_file = ...
    './results/comak-inverse-kinematics/secondary_coordinate_constraint_functions.xml';
comak_ik.set_secondary_constraint_function_file(secondary_constraint_function_file);
comak_ik.set_constraint_function_num_interpolation_points(20);
comak_ik.set_print_secondary_constraint_sim_results(true);
comak_ik.set_constrained_model_file('./results/comak-inverse-kinematics/ik_constrained_model.osim');
comak_ik.set_perform_inverse_kinematics(true);
comak_ik.set_marker_file(marker_file);

comak_ik.set_output_motion_file(ik_output_motion_file);
comak_ik.set_time_range(0, start_time - 0.2);
comak_ik.set_time_range(1, stop_time);
comak_ik.set_report_errors(true);
comak_ik.set_report_marker_locations(false);
comak_ik.set_ik_constraint_weight(100);
comak_ik.set_ik_accuracy(1e-5);
comak_ik.set_use_visualizer(useVisualizer);

tasks(1).Name = 'Sacral'; tasks(1).Weight = 20;    
tasks(2).Name = 'R_Asis'; tasks(2).Weight = 30;    
tasks(3).Name = 'L_Asis'; tasks(3).Weight = 30;    
tasks(4).Name = 'R_Psis'; tasks(4).Weight = 20;    
tasks(5).Name = 'L_Psis'; tasks(5).Weight = 20;    
tasks(6).Name = 'R_Shoulder'; tasks(6).Weight = 1;    
tasks(7).Name = 'L_Shoulder'; tasks(7).Weight = 1;    
tasks(8).Name = 'R_Elbow'; tasks(8).Weight = 1;    
tasks(9).Name = 'R_Wrist'; tasks(9).Weight = 1;    
tasks(10).Name = 'L_Elbow'; tasks(10).Weight = 1;    
tasks(11).Name = 'L_Wrist'; tasks(11).Weight = 1;    
tasks(12).Name = 'R_Patella'; tasks(12).Weight = 1;    
tasks(13).Name = 'R_Thigh_Superior'; tasks(13).Weight = 1;    
tasks(14).Name = 'R_Thigh_Inferior'; tasks(14).Weight = 1;        
tasks(15).Name = 'R_Thigh_Lateral'; tasks(15).Weight = 1;    
tasks(16).Name = 'R_Shank_Superior'; tasks(16).Weight = 1;    
tasks(17).Name = 'R_Shank_Inferior'; tasks(17).Weight = 1;    
tasks(18).Name = 'R_Shank_Lateral'; tasks(18).Weight = 1;    
tasks(19).Name = 'R_Midfoot_Medial'; tasks(19).Weight = 1;    
tasks(20).Name = 'R_Midfoot_Lateral'; tasks(20).Weight = 1;    
tasks(21).Name = 'R_Hindfoot'; tasks(21).Weight = 1;    
tasks(22).Name = 'R_Midfoot_Superior'; tasks(22).Weight = 1;    
tasks(23).Name = 'R_ToeMedial'; tasks(23).Weight = 1;    
tasks(24).Name = 'R_ToeLateral'; tasks(24).Weight = 1;    
tasks(25).Name = 'R_Toe'; tasks(25).Weight = 1;    
tasks(26).Name = 'R_Heel'; tasks(26).Weight = 1;    
tasks(27).Name = 'L_Patella'; tasks(27).Weight = 1;    
tasks(28).Name = 'L_Thigh_Superior'; tasks(28).Weight = 1;    
tasks(29).Name = 'L_Thigh_Inferior'; tasks(29).Weight = 1;    
tasks(30).Name = 'L_Thigh_Lateral'; tasks(30).Weight = 1;    
tasks(31).Name = 'L_Shank_Superior'; tasks(31).Weight = 1;    
tasks(32).Name = 'L_Shank_Inferior'; tasks(32).Weight = 1;    
tasks(33).Name = 'L_Shank_Lateral'; tasks(33).Weight = 1;    
tasks(34).Name = 'L_Midfoot_Medial'; tasks(34).Weight = 1;    
tasks(35).Name = 'L_Midfoot_Lateral'; tasks(35).Weight = 1;    
tasks(36).Name = 'L_Hindfoot'; tasks(36).Weight = 1;    
tasks(37).Name = 'L_Midfoot_Superior'; tasks(37).Weight = 1;    
tasks(38).Name = 'L_ToeMedial'; tasks(38).Weight = 1;    
tasks(39).Name = 'L_ToeLateral'; tasks(39).Weight = 1;    
tasks(40).Name = 'L_Toe'; tasks(40).Weight = 1;        
tasks(41).Name = 'L_Heel'; tasks(41).Weight = 1;

ik_task_set = IKTaskSet();
ik_task = IKMarkerTask();

for j = 1:length(tasks)
    ik_task.setName(tasks(j).Name);
    ik_task.setWeight(tasks(j).Weight);
    ik_task_set.cloneAndAppend(ik_task);
end

comak_ik.set_IKTaskSet(ik_task_set);
comak_ik.print(['./inputs/' trial '/comak_inverse_kinematics_settings.xml']);

disp('Running COMAKInverseKinematicsTool...')
% comak_ik.run();
%% Perform COMAK Simulation

comak = COMAKTool();
comak.set_model_file(model_file);
comak.set_coordinates_file([ik_result_dir '/' ik_output_motion_file]);
comak.set_external_loads_file(external_loads_file),
comak.set_results_directory(comak_result_dir);
comak.set_results_prefix(results_basename);
comak.set_replace_force_set(false);
comak.set_force_set_file('../../../models/knee_tka/grand_challenge/DM/DM_reserve_actuators.xml');
comak.set_start_time(start_time - 0.1);
comak.set_stop_time(stop_time);
comak.set_time_step(0.01);
comak.set_lowpass_filter_frequency(6);
comak.set_print_processed_input_kinematics(false);
comak.set_prescribed_coordinates(0,'/jointset/gnd_pelvis/pelvis_tx');
comak.set_prescribed_coordinates(1,'/jointset/gnd_pelvis/pelvis_ty');
comak.set_prescribed_coordinates(2,'/jointset/gnd_pelvis/pelvis_tz');
comak.set_prescribed_coordinates(3,'/jointset/gnd_pelvis/pelvis_tilt');
comak.set_prescribed_coordinates(4,'/jointset/gnd_pelvis/pelvis_list');
comak.set_prescribed_coordinates(5,'/jointset/gnd_pelvis/pelvis_rot');
comak.set_prescribed_coordinates(6,'/jointset/subtalar_r/subt_angle_r');
comak.set_prescribed_coordinates(7,'/jointset/mtp_r/mtp_angle_r');
comak.set_prescribed_coordinates(8,'/jointset/hip_l/hip_flex_l');
comak.set_prescribed_coordinates(9,'/jointset/hip_l/hip_add_l');
comak.set_prescribed_coordinates(10,'/jointset/hip_l/hip_rot_l');
comak.set_prescribed_coordinates(11,'/jointset/pf_l/pf_l_r3');
comak.set_prescribed_coordinates(12,'/jointset/pf_l/pf_l_tx');
comak.set_prescribed_coordinates(13,'/jointset/pf_l/pf_l_ty');
comak.set_prescribed_coordinates(14,'/jointset/knee_l/knee_flex_l');
comak.set_prescribed_coordinates(15,'/jointset/ankle_l/ankle_flex_l');
comak.set_prescribed_coordinates(16,'/jointset/subtalar_l/subt_angle_l');
comak.set_prescribed_coordinates(17,'/jointset/mtp_l/mtp_angle_l');
comak.set_prescribed_coordinates(18,'/jointset/pelvis_torso/lumbar_ext');
comak.set_prescribed_coordinates(19,'/jointset/pelvis_torso/lumbar_latbend');
comak.set_prescribed_coordinates(20,'/jointset/pelvis_torso/lumbar_rot');
comak.set_prescribed_coordinates(21,'/jointset/torso_neckhead/neck_ext');
comak.set_prescribed_coordinates(22,'/jointset/torso_neckhead/neck_latbend');
comak.set_prescribed_coordinates(23,'/jointset/torso_neckhead/neck_rot');
comak.set_prescribed_coordinates(24,'/jointset/acromial_r/arm_add_r');
comak.set_prescribed_coordinates(25,'/jointset/acromial_r/arm_flex_r');
comak.set_prescribed_coordinates(26,'/jointset/acromial_r/arm_rot_r');
comak.set_prescribed_coordinates(27,'/jointset/elbow_r/elbow_flex_r');
comak.set_prescribed_coordinates(28,'/jointset/radioulnar_r/pro_sup_r');
comak.set_prescribed_coordinates(29,'/jointset/radius_hand_r/wrist_flex_r');
comak.set_prescribed_coordinates(30,'/jointset/acromial_l/arm_add_l');
comak.set_prescribed_coordinates(31,'/jointset/acromial_l/arm_flex_l');
comak.set_prescribed_coordinates(32,'/jointset/acromial_l/arm_rot_l');
comak.set_prescribed_coordinates(33,'/jointset/elbow_l/elbow_flex_l');
comak.set_prescribed_coordinates(34,'/jointset/radioulnar_l/pro_sup_l');
comak.set_prescribed_coordinates(35,'/jointset/radius_hand_l/wrist_flex_l');
 
comak.set_primary_coordinates(0,'/jointset/hip_r/hip_flex_r');
comak.set_primary_coordinates(1,'/jointset/hip_r/hip_add_r');
comak.set_primary_coordinates(2,'/jointset/hip_r/hip_rot_r');
comak.set_primary_coordinates(3,'/jointset/knee_r/knee_flex_r');
comak.set_primary_coordinates(4,'/jointset/ankle_r/ankle_flex_r');

secondary_coord_set = COMAKSecondaryCoordinateSet(); 
secondary_coord = COMAKSecondaryCoordinate();

secondary_coord.setName('knee_add_r');
secondary_coord.set_max_change(0.005);
secondary_coord.set_coordinate('/jointset/knee_r/knee_add_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('knee_rot_r');
secondary_coord.set_max_change(0.005);
secondary_coord.set_coordinate('/jointset/knee_r/knee_rot_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('knee_tx_r');
secondary_coord.set_max_change(0.001);
secondary_coord.set_coordinate('/jointset/knee_r/knee_tx_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('knee_ty_r');
secondary_coord.set_max_change(0.001);
secondary_coord.set_coordinate('/jointset/knee_r/knee_ty_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('knee_tz_r');
secondary_coord.set_max_change(0.001);
secondary_coord.set_coordinate('/jointset/knee_r/knee_tz_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('pf_flex_r');
secondary_coord.set_max_change(0.005);
secondary_coord.set_coordinate('/jointset/pf_r/pf_flex_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('pf_rot_r');
secondary_coord.set_max_change(0.005);
secondary_coord.set_coordinate('/jointset/pf_r/pf_rot_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('pf_tilt_r');
secondary_coord.set_max_change(0.005);
secondary_coord.set_coordinate('/jointset/pf_r/pf_tilt_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('pf_tx_r');
secondary_coord.set_max_change(0.001);
secondary_coord.set_coordinate('/jointset/pf_r/pf_tx_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('pf_ty_r');
secondary_coord.set_max_change(0.001);
secondary_coord.set_coordinate('/jointset/pf_r/pf_ty_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

secondary_coord.setName('pf_tz_r');
secondary_coord.set_max_change(0.001);
secondary_coord.set_coordinate('/jointset/pf_r/pf_tz_r');
secondary_coord_set.cloneAndAppend(secondary_coord);

comak.set_COMAKSecondaryCoordinateSet(secondary_coord_set);

comak.set_settle_secondary_coordinates_at_start(true);
comak.set_settle_threshold(1e-4);
comak.set_settle_accuracy(1e-3);
comak.set_settle_internal_step_limit(10000);
comak.set_print_settle_sim_results(true);
comak.set_settle_sim_results_directory(comak_result_dir);
comak.set_settle_sim_results_prefix('walking_tka_settle_sim');
comak.set_max_iterations(15);
comak.set_udot_tolerance(1);
comak.set_udot_worse_case_tolerance(50);
comak.set_unit_udot_epsilon(1e-6);
comak.set_optimization_scale_delta_coord(1);
comak.set_ipopt_diagnostics_level(3);
comak.set_ipopt_max_iterations(5000);
comak.set_ipopt_convergence_tolerance(1e-4);
comak.set_ipopt_constraint_tolerance(1e-4);
comak.set_ipopt_limited_memory_history(200);
comak.set_ipopt_nlp_scaling_max_gradient(10000);
comak.set_ipopt_nlp_scaling_min_value(1e-8);
comak.set_ipopt_obj_scaling_factor(1);
comak.set_activation_exponent(2);
comak.set_contact_energy_weight(0);
comak.set_non_muscle_actuator_weight(1000);
comak.set_model_assembly_accuracy(1e-12);
comak.set_use_visualizer(useVisualizer);

comak.print(['./inputs/' trial '/comak_settings.xml']);

disp('Running COMAK Tool...')
comak.run();

%% Perform Joint Mechanics Analysis
jnt_mech = JointMechanicsTool();
jnt_mech.set_model_file(model_file);
jnt_mech.set_input_states_file([comak_result_dir '/' results_basename '_states.sto']);
jnt_mech.set_use_muscle_physiology(false);
jnt_mech.set_results_file_basename(results_basename);
jnt_mech.set_results_directory(jnt_mech_result_dir);
jnt_mech.set_start_time(start_time);
jnt_mech.set_stop_time(stop_time);
jnt_mech.set_normalize_to_cycle(true);
jnt_mech.set_contacts(0,'all');
jnt_mech.set_contact_outputs(0,'all');
jnt_mech.set_contact_mesh_properties(0,'none');
jnt_mech.set_ligaments(0,'all');
jnt_mech.set_ligament_outputs(0,'all');
jnt_mech.set_muscles(0,'all');
jnt_mech.set_muscle_outputs(0,'all');
jnt_mech.set_attached_geometry_bodies(0,'all');
jnt_mech.set_output_orientation_frame('/bodyset/tibia_r/tibia_implant');
jnt_mech.set_output_position_frame('/bodyset/tibia_r/tibia_implant');
jnt_mech.set_write_vtp_files(true);
jnt_mech.set_vtp_file_format('binary');
jnt_mech.set_write_h5_file(false);
jnt_mech.set_h5_kinematics_data(true);
jnt_mech.set_h5_states_data(true);
jnt_mech.set_use_visualizer(useVisualizer);
jnt_mech.print(['./inputs/' trial '/joint_mechanics_settings.xml']);

disp('Running JointMechanicsTool...');
jnt_mech.run();


%% Plot Simulation Results

% COMAK Inverse Kinematics

settle_file = [ik_result_dir '/' results_basename '_secondary_constraint_settle_states.sto'];
sweep_file = [ik_result_dir '/' results_basename '_secondary_constraint_sweep_states.sto'];



sweep = osimTableToStruct(TimeSeriesTable(sweep_file));
settle = osimTableToStruct(TimeSeriesTable(sweep_file));

coord_names = {'/jointset/pf_r/pf_flex_r','/jointset/pf_r/pf_rot_r','/jointset/pf_r/pf_tilt_r';...
    '/jointset/pf_r/pf_tx_r','/jointset/pf_r/pf_ty_r','/jointset/pf_r/pf_tz_r';...
    '/jointset/knee_r/knee_flex_r','/jointset/knee_r/knee_add_r','/jointset/knee_r/knee_rot_r';...
    '/jointset/knee_r/knee_tx_r','/jointset/knee_r/knee_ty_r','/jointset/knee_r/knee_tz_r'}; 	

groups = {'PF_rot','PF_trans','TF_rot','TF_trans'};

[m,n] = size(coord_names);

% for i = 1:m
%     figure('name',groups{i})
%     for j = 1:n
%         
%         name = coord_names{i,j};
%     
%     
%         ind = find(contains(labels,name));
%         
%         subplot(1,3,j);
% 
%         hold on
% %        scatter(n_data(:,1),n_data(:,ind)*180/pi);
%         plot(data(:,1),data(:,ind));
% 
%         title(name);
%     end
% end


figure('name','Secondary Constrain Functions')
func_set = FunctionSet(secondary_constraint_function_file);

for f = 0:func_set.getSize()-1
    func = GCVSpline.safeDownCast(func_set.get(f));
    
    for n = 0:func.getNumberOfPoints()-1
        X(n+1) = func.getX(n);
        Y(n+1) = func.getY(n);
    end
    
    subplot(4,3,f+1)
    plot(X,Y)
    title(char(func.getName()))
    xlabel('knee flexion [^o]') 
end
% Plot Secondary Constraint Functions