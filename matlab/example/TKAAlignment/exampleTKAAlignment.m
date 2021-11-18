%% example TKA Alignment
%==========================================================================

%% Setup Environment and Folders
% clear;
close all

import org.opensim.modeling.*
Logger.setLevelString('info');

useVisualizer = true;

model_file = '../../../models/knee_tka/grand_challenge/DM/DM.osim';
geometry_path='../../../models/knee_tka/grand_challenge/DM/Geometry';

ModelVisualizer.addDirToGeometrySearchPaths(geometry_path);

if(exist('./inputs','dir')~=7)
    mkdir('./inputs')
end
if(exist('./results','dir')~=7)
    mkdir('./results')
end
if(exist('./results/graphics','dir')~=7)
    mkdir('./results/graphics')
end
if(exist(['./results/ligament-balance' ],'dir')~=7)
    mkdir(['./results/ligament-balance' ])
end
if(exist(['./results/comak-inverse-kinematics' ],'dir')~=7)
    mkdir(['./results/comak-inverse-kinematics' ])
end
if(exist(['./results/comak' ],'dir')~=7)
    mkdir(['./results/comak' ])
end
if(exist('./results/joint-mechanics' ,'dir')~=7)
    mkdir('./results/joint-mechanics' )
end
% sim_names = {'nominal','valgus_3','varus_3','valgus_6','varus_6',};
% femur_x_rot = [0 0 0 0 0] * pi/180;
% tibia_x_rot = [0 3 -3 6 -6] * pi/180;
sim_names = {'varus_6','valgus_6','nominal'};
femur_x_rot = [-3 3 0 0] * pi/180;
tibia_x_rot = [3 -3 0] * pi/180;
nSim = 2;%length(sim_names);

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

%% Setup outputs
results_basename = 'walking_tka';
for i = 1:nSim      
    model_name{i} = ['DM_' sim_names{i}];
    lig_balance_result_dir{i} = ['./results/ligament-balance/' sim_names{i}];
    ik_result_dir{i} = ['./results/comak-inverse-kinematics/' sim_names{i}];
    comak_result_dir{i} = ['./results/comak/' sim_names{i}];
    jnt_mech_result_dir{i} = ['./results/joint-mechanics/' sim_names{i}];
end

%% Create Alignment Models and Perform Ligament Balancing
model = Model(model_file);
state = model.initSystem();
frc_set = model.getForceSet();      


numLigaments = 0;
for i = 0:frc_set.getSize()-1
   frc = frc_set.get(i);
   if(strcmp(frc.getConcreteClassName(),'Blankevoort1991Ligament'))
       lig = Blankevoort1991Ligament.safeDownCast(frc);
       numLigaments = numLigaments + 1;
       default_ref_strain(numLigaments) = lig.getStrain(state);
   end
   
end
    
for i = 1:nSim
    % Change Alignment in Model
    model = Model(model_file);
    state = model.initSystem();
    model.setName(model_name{i});
    
    femur_mesh = Smith2018ContactMesh.safeDownCast(...
        model.getComponent('/contactgeometryset/femur_implant'));

    tibia_mesh = Smith2018ContactMesh.safeDownCast(...
        model.getComponent('/contactgeometryset/tibia_implant'));

    femur_mesh.set_orientation(Vec3(femur_x_rot(i),0,0));
    tibia_mesh.set_orientation(Vec3(tibia_x_rot(i),0,0));
    
    % try to guess final alignment to improve contact at 1st time step
    def_add_value = (tibia_x_rot(i) - femur_x_rot(i))/2;
    model.getCoordinateSet().get('knee_add_r').setDefaultValue(def_add_value);

    % Ligament Balancing
    for j=1:2
        basename = ['ligament_balance_' int2str(j)];

        % Settle Knee into equilibrium position for new implant alignment
        clear lig_forsim;
        
        lig_forsim = ForsimTool();
        lig_forsim.setModel(model);
        lig_forsim.set_results_directory(lig_balance_result_dir{i});
        lig_forsim.set_results_file_basename(basename);
        
        lig_forsim.set_start_time(0);
        lig_forsim.set_stop_time(2);
        lig_forsim.set_integrator_accuracy(1e-3); 
        lig_forsim.set_constant_muscle_control(0.02); %Set all muscles to 2% activation to represent passive state
        lig_forsim.set_use_activation_dynamics(true);
        lig_forsim.set_use_tendon_compliance(false);
        lig_forsim.set_use_muscle_physiology(true);
        lig_forsim.set_equilibrate_muscles(true);

        for s = 1:numSecondaryCoordinates
            lig_forsim.set_unconstrained_coordinates(s-1,['/jointset/' secondary_joints{s} '/' secondary_coordinates{s}]);
        end

        lig_forsim.set_use_visualizer(useVisualizer);
        lig_forsim.print('./inputs/ligament_balance_settings.xml');

        disp('Running ForsimTool to settle knee. iteration: ')
        lig_forsim.run();            

        % Reset Ligament Reference Strains
        lig_balance_states_sto = [lig_balance_result_dir{i} '/' basename '_states.sto'];
        lig_balance_results = osimTableToStruct(TimeSeriesTable(lig_balance_states_sto));           

        for s = 1:numSecondaryCoordinates
            label = ['a_jointset_' secondary_joints{s} '_' secondary_coordinates{s} '_value'];
            value = lig_balance_results.(label);
            value = value(end);
            coord = model.getCoordinateSet().get(secondary_coordinates{s});
            coord.setValue(state,value);
            coord.setDefaultValue(value);
        end

        
        frc_set = model.getForceSet(); 
        nLig = 0;
        for f = 0:frc_set.getSize()-1
           frc = frc_set.get(f);
           if(strcmp(frc.getConcreteClassName(),'Blankevoort1991Ligament'))
               lig = Blankevoort1991Ligament.safeDownCast(frc);
               nLig = nLig + 1;
               lig.setSlackLengthFromReferenceStrain(default_ref_strain(nLig),state);
           end
        end
    end

    model.finalizeConnections();
    new_model_file =['./inputs/DM_' sim_names{i} '.osim'];
    model.print(new_model_file); 
end
    
    %         frc_set = model.getForceSet();
%         for f = 0:frc_set.getSize()-1
%             frc = frc_set.get(f);
%             
%             if(contains(sim_names{i},'MCL'))
%                 lig = Blankevoort1991Ligament.safeDownCast(frc);
%                 if(contains(char(frc.getName()),'MCLs') ||...
%                     contains(char(frc.getName()),'MCLd'))
%                     lig.setSlackLengthFromReferenceStrain(0.0,state);
%                 end
%             end
%         end

%% Perform Simulations
if(true)     
    for i = 1:nSim        
        new_model_file =['./inputs/DM_' sim_names{i} '.osim'];
        
        % COMAK Inverse Kinematics
        comak_ik = COMAKInverseKinematicsTool();
        %comak_ik.setModel(model);

        comak_ik.set_model_file(new_model_file);
        comak_ik.set_results_directory(ik_result_dir{i});
        comak_ik.set_results_prefix(results_basename);
        comak_ik.set_perform_secondary_constraint_sim(true);
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
        comak_ik.set_secondary_constraint_sim_sweep_time(1.0);
        comak_ik.set_secondary_coupled_coordinate_start_value(0);
        comak_ik.set_secondary_coupled_coordinate_stop_value(100);
        comak_ik.set_secondary_constraint_sim_integrator_accuracy(1e-4);
        comak_ik.set_secondary_constraint_sim_internal_step_limit(10000);
        comak_ik.set_secondary_constraint_function_file(...
            [ik_result_dir{i} '/secondary_coordinate_constraint_functions.xml']);
        comak_ik.set_constraint_function_num_interpolation_points(20);
        comak_ik.set_print_secondary_constraint_sim_results(true);
        comak_ik.set_constrained_model_file([ik_result_dir{i} '/ik_constrained_model.osim']);
        comak_ik.set_perform_inverse_kinematics(true);
        comak_ik.set_marker_file('../../../models/knee_tka/grand_challenge/DM/experimental_data/motion_analysis/DM_smooth1_new.trc');

        comak_ik.set_output_motion_file('DM_smooth1_ik.sto');
        comak_ik.set_time_range(0, 2.2);
        comak_ik.set_time_range(1, 3.775);
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
        comak_ik.print(['./inputs/comak_inverse_kinematics_settings.xml']);

        disp('Running COMAKInverseKinematicsTool...')
        comak_ik.run();


        jmt_sweep = JointMechanicsTool();        
        jmt_sweep.set_model_file(new_model_file);
        jmt_sweep.set_use_muscle_physiology(false);
        jmt_sweep.set_results_directory(ik_result_dir{i});
        jmt_sweep.set_normalize_to_cycle(false);
        jmt_sweep.set_input_states_file([ik_result_dir{i} '/' results_basename '_secondary_constraint_sweep_states.sto']);
        jmt_sweep.set_results_file_basename([results_basename '_secondary_constraint_sweep']);
        jmt_sweep.set_start_time(-1);
        jmt_sweep.set_stop_time(-1);
        jmt_sweep.set_contacts(0,'all');
        jmt_sweep.set_contact_outputs(0,'all');
        jmt_sweep.set_ligaments(0,'all');
        jmt_sweep.set_ligament_outputs(0,'total_force');
        jmt_sweep.set_ligament_outputs(1,'strain');
        jmt_sweep.set_muscles(0,'none');
        jmt_sweep.set_muscle_outputs(0,'all');
        jmt_sweep.set_write_vtp_files(false);
        jmt_sweep.set_h5_kinematics_data(true);
        jmt_sweep.set_h5_states_data(false);
        jmt_sweep.set_write_h5_file(true);
        jmt_sweep.set_use_visualizer(useVisualizer);
        jmt_sweep.run();
        
       
        %% Perform COMAK Simulation

        comak = COMAKTool();
        comak.set_model_file(new_model_file);
        comak.set_coordinates_file([ik_result_dir{i} '/DM_smooth1_ik.sto']);
        comak.set_external_loads_file('../../../models/knee_tka/grand_challenge/DM/experimental_data/motion_analysis/DM_smooth1_external_loads.xml'),
        comak.set_results_directory(comak_result_dir{i});
        comak.set_results_prefix(results_basename);
        comak.set_replace_force_set(false);
        comak.set_force_set_file('../../../models/knee_tka/grand_challenge/DM/DM_reserve_actuators.xml');
        comak.set_start_time(2.5);
        comak.set_stop_time(3.775);
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
        comak.set_settle_sim_results_directory(comak_result_dir{i});
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

        comak.print('./inputs/comak_settings.xml');

        disp('Running COMAK Tool...')
         comak.run();

        %% Perform Joint Mechanics Analysis
        jnt_mech = JointMechanicsTool();
        jnt_mech.set_model_file(new_model_file);
        jnt_mech.set_input_states_file([comak_result_dir{i} '/' results_basename '_states.sto']);
        jnt_mech.set_use_muscle_physiology(false);
        jnt_mech.set_input_activations_file([comak_result_dir{i} '/' results_basename '_activations.sto']);
        jnt_mech.set_input_comak_convergence_file([comak_result_dir{i} '/' results_basename '_convergence.sto']);
        jnt_mech.set_results_file_basename(results_basename);
        jnt_mech.set_results_directory(jnt_mech_result_dir{i});
        jnt_mech.set_start_time(2.53);
        jnt_mech.set_stop_time(3.775);
        jnt_mech.set_resample_step_size(-1);
        jnt_mech.set_normalize_to_cycle(true);
        jnt_mech.set_lowpass_filter_frequency(-1);
        jnt_mech.set_print_processed_kinematics(false);
        jnt_mech.set_contacts(0,'all');
        jnt_mech.set_contact_outputs(0,'all');
        jnt_mech.set_contact_mesh_properties(0,'none');
        jnt_mech.set_ligaments(0,'all');
        jnt_mech.set_ligament_outputs(0,'total_force');
        jnt_mech.set_ligament_outputs(1,'strain');
        jnt_mech.set_muscles(0,'none');
        jnt_mech.set_muscle_outputs(0,'all');

        jnt_mech.set_attached_geometry_bodies(0,'none');

        jnt_mech.set_output_orientation_frame('ground');
        jnt_mech.set_output_position_frame('ground');
        jnt_mech.set_write_vtp_files(true);
        jnt_mech.set_write_h5_file(true);
        jnt_mech.set_h5_kinematics_data(true);
        jnt_mech.set_use_visualizer(useVisualizer);
        jnt_mech.print('./inputs/joint_mechanics_settings.xml');

        disp('Running JointMechanicsTool...');
         jnt_mech.run();
    end
end
%% Analyze the results
close all;


for i = 1:nSim
    model_name{i} = ['DM_' sim_names{i}];
%     sweep_files{i} = [ik_result_dir{i} '/' results_basename '_secondary_constraint_sweep.h5'];
    comak_files{i} = [ jnt_mech_result_dir{i} '/' results_basename '.h5'];    
end

% sweep_results = jam_analysis(sweep_files);
comak_results = jam_analysis(comak_files);


line_type = {'-','-.','-.','-*','-*'};
% Plot COMAK IK Settle Simulations
sec_coords = {...
    'knee_flex_r','knee_add_r','knee_rot_r',...
    'knee_tx_r','knee_ty_r','knee_tz_r',...
    'pf_flex_r','pf_rot_r','pf_tilt_r',...
    'pf_tx_r','pf_ty_r','pf_tz_r',...
    };

% Plot COMAK IK Sweep Simulations 
% for i = 1:length(sec_coords)
%     figure('name',sec_coords{i})
%     hold on
%     for n=1:sweep.nFiles
%         plot(sweep.coordinateset.(sec_coords{i}).value(:,n));
%     end
% end
% Plot COMAK IK Joint Angles

% Plot COMAK Secondary Kinematics

for i = 1:length(sec_coords)
    figure('name',sec_coords{i})
    hold on
    for n=1:comak_results.num_files
        plot(comak_results.coordinateset.(sec_coords{i}).value(:,n),line_type{n});
    end

    legend(sim_names)
end
        
% Plot COMAK Ligament Forces 
ligament_names = {...
    'MCLd','MCLs'...
    'LCL',...
    'ACLam','ACLpl',...
    'PCLal','PCLpm',...
    'PT',...
    'ITB',...
    'pCAP'...
    };

for i = 1:length(ligament_names)
    figure('name',ligament_names{i});hold on

    fiber_names = fieldnames(comak_results.forceset.Blankevoort1991Ligament);

    fibers = fiber_names(contains(fiber_names,ligament_names{i}));
    for n=1:comak_results.num_files
        data = 0;
        for k = 1:length(fibers)
            data = data + comak_results.forceset.Blankevoort1991Ligament.(fibers{k}).total_force(:,n);

        end
        plot(data,line_type{n})
        title([ligament_names{i} ' force'])
    end
    legend(sim_names)
end

    
    % Plot COMAK Muscle Forces
if true
    % Plot COMAK Contact Forces
    %PF Contact
    figure('name','PF Contact Force')
    comp = {'X','Y','Z'};
    for i = 1:3
        subplot(1,3,i);hold on
        for n=1:comak_results.num_files
            plot(comak_results.forceset.Smith2018ArticularContactForce.pf_contact.patella_implant.total_contact_force(i,:,n),line_type{n});
        end
        title(['PF Contact ' comp{i}])
    end
    legend(sim_names)

    %TF Contact
    figure('name','TF Contact Force')
    comp = {'X','Y','Z'};
%     for i = 1:3
%         subplot(3,3,i);hold on
%         for n=1:comak_results.nFiles
%             plot(comak_results.forceset.Smith2018ArticularContactForce.knee_contact.tibia_implant.total_contact_force(i,:,n),line_type{n});
%         end
%         title(['Total TF Contact ' comp{i}])
% 
%         subplot(3,3,3+i);hold on
%         for n=1:comak_results.nFiles
%             plot(comak_results.forceset.Smith2018ArticularContactForce.knee_contact.tibia_implant.region(6).regional_contact_force(:,i,n),line_type{n});
%         end
%         title(['Medial TF Contact ' comp{i}])
% 
%         subplot(3,3,6+i);hold on
%         for n=1:comak_results.nFiles
%             plot(comak_results.forceset.Smith2018ArticularContactForce.knee_contact.tibia_implant.region(5).regional_contact_force(:,i,n),line_type{n});
%         end
%         title(['Lateral TF Contact ' comp{i}])
%     end
    for i = 2%1:3
        subplot(1,3,1);hold on
        for n=1:comak_results.num_files
            plot(comak_results.forceset.Smith2018ArticularContactForce.knee_contact.tibia_implant.total_contact_force(:,i,n),line_type{n});
        end
        title(['Total TF Contact ' comp{i}])

        subplot(1,3,2);hold on
        for n=1:comak_results.num_files
            plot(comak_results.forceset.Smith2018ArticularContactForce.knee_contact.tibia_implant.region(6).regional_contact_force(:,i,n),line_type{n});
        end
        title(['Medial TF Contact ' comp{i}])

        subplot(1,3,3);hold on
        for n=1:comak_results.num_files
            plot(comak_results.forceset.Smith2018ArticularContactForce.knee_contact.tibia_implant.region(5).regional_contact_force(:,i,n),line_type{n});
        end
        title(['Lateral TF Contact ' comp{i}])
    end
    legend(sim_names)

    % Plot COMAK Convergence
    figure('name','Comak Convergence')
    subplot(2,1,1);hold on
    for n=1:comak_results.num_files
        plot(comak_results.comak.max_udot_error(:,n),line_type{n});
    end
    title('max udot error')

    subplot(2,1,2);hold on
    for n=1:comak_results.num_files
        plot(comak_results.comak.iterations(:,n),line_type{n});
    end
    title('iterations')
    legend(sim_names)
end