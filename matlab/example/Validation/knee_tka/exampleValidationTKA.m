%% exampleWalkingTKA
%==========================================================================
% Simulate treadmill walking cycles using COMAKTool to predict 
% tibiofemoral kinematics and contact forces to compare against 
% instrumented implant and fluoroscopy measurements.
% 
%
%
%==========================================================================
%% Setup Environment and Folders
clear; close all;
import org.opensim.modeling.*
Logger.setLevelString('info');

useVisualizer = true;

grf_threshold = 15; 
%7 N is used in grand challenge README.rtf
% but 15 N is the lowest that works due to
% grf noise during swing phase

model_file = '../../../../models/knee_tka/grand_challenge/DM/DM.osim';
experimental_data_dir = ...
    '../../../../models/knee_tka/grand_challenge/DM/experimental_data/motion_analysis/';

trial_names = {'DM_ngait_tmf_slow1','DM_ngait_tmf_slow2'};
num_trials = length(trial_names);

if(exist('./inputs','dir')~=7); mkdir('./inputs'); end
if(exist('./results','dir')~=7); mkdir('./results'); end
if(exist(['./results/graphics'],'dir')~=7)
    mkdir(['./results/graphics'])
end

%% Find Start-Stop Times of each cycle
for i = 1:num_trials
    trial(i).name = trial_names{i};
    
    trial(i).trc_file = ...
        [experimental_data_dir '/' trial(i).name '_new.trc']; %'../../../models/knee_tka/grand_challenge/DM/experimental_data/motion_analysis';
    
    trial(i).external_loads_xml_file =  ...
        [experimental_data_dir '/' trial(i).name '_external_loads.xml']; %'../../../models/knee_tka/grand_challenge/DM/experimental_data/motion_analysis';
    
    trial(i).grf_sto_file = ...
        [experimental_data_dir '/' trial(i).name '_grf_filt.sto']; %'../../../models/knee_tka/grand_challenge/DM/experimental_data/motion_analysis';
    
    trial(i).grf.data = osimTableToStruct(TimeSeriesTable(trial(i).grf_sto_file));
    
    trial(i).grf.left_vy = trial(i).grf.data.ground_force_1_vy;
    trial(i).grf.right_vy = trial(i).grf.data.ground_force_2_vy;
    
    % Right Foot GRF
    trial(i).grf.right.thresholded = (grf_threshold >= trial(i).grf.right_vy);
    trial(i).grf.right.crossover = diff(trial(i).grf.right.thresholded);
    
    trial(i).right.heel_strike.frame = find(trial(i).grf.right.crossover==-1);
    trial(i).right.toe_off.frame = find(trial(i).grf.right.crossover==1);
    
    trial(i).right.heel_strike.time = trial(i).grf.data.time(trial(i).right.heel_strike.frame);
    trial(i).right.heel_strike.num = length(trial(i).right.heel_strike.time);
    trial(i).right.toe_off.time = trial(i).grf.data.time(trial(i).right.toe_off.frame);
    trial(i).right.toe_off.num = length(trial(i).right.toe_off.time);
    
    trial(i).right.num_cycles = length(trial(i).right.heel_strike.frame)-1;
    
    for c = 1:trial(i).right.num_cycles
        trial(i).right.cycle(c).full_cycle.start.frame = trial(i).right.heel_strike.frame(c);
        trial(i).right.cycle(c).full_cycle.start.time = trial(i).right.heel_strike.time(c);
        trial(i).right.cycle(c).full_cycle.stop.frame = trial(i).right.heel_strike.frame(c+1);       
        trial(i).right.cycle(c).full_cycle.stop.time = trial(i).right.heel_strike.time(c+1);
        
        trial(i).right.cycle(c).stance.start.frame = trial(i).right.heel_strike.frame(c);
        trial(i).right.cycle(c).stance.start.time = trial(i).right.heel_strike.time(c);
        
        if trial(i).right.toe_off.frame(c) > trial(i).right.heel_strike.frame(c)
            trial(i).right.cycle(c).stance.stop.frame = trial(i).right.toe_off.frame(c);       
            trial(i).right.cycle(c).stance.stop.time = trial(i).right.toe_off.time(c);
        else
            trial(i).right.cycle(c).stance.stop.frame = trial(i).right.toe_off.frame(c+1);       
            trial(i).right.cycle(c).stance.stop.time = trial(i).right.toe_off.time(c+1);
        end
        
        
        
        if trial(i).right.toe_off.frame(c) > trial(i).right.heel_strike.frame(c)
            trial(i).right.cycle(c).swing.start.frame = trial(i).right.toe_off.frame(c);
            trial(i).right.cycle(c).swing.start.time = trial(i).right.toe_off.time(c);            
        else
            trial(i).right.cycle(c).swing.start.frame = trial(i).right.toe_off.frame(c+1);
            trial(i).right.cycle(c).swing.start.time = trial(i).right.toe_off.time(c+1);
        end
        
        trial(i).right.cycle(c).swing.stop.frame = trial(i).right.heel_strike.frame(c+1);       
        trial(i).right.cycle(c).swing.stop.time = trial(i).right.heel_strike.time(c+1);
    end    
    
    % Left Foot GRF
    trial(i).grf.left.thresholded = (grf_threshold >= trial(i).grf.left_vy);
    trial(i).grf.left.crossover = diff(trial(i).grf.left.thresholded);
    
    trial(i).left.heel_strike.frame = find(trial(i).grf.left.crossover==-1);
    trial(i).left.toe_off.frame = find(trial(i).grf.left.crossover==1);
    
    trial(i).left.heel_strike.time = trial(i).grf.data.time(trial(i).left.heel_strike.frame);
    trial(i).left.heel_strike.num = length(trial(i).left.heel_strike.time);
    trial(i).left.toe_off.time = trial(i).grf.data.time(trial(i).left.toe_off.frame);
    trial(i).left.toe_off.num = length(trial(i).left.toe_off.time);
    
    trial(i).left.num_cycles = length(trial(i).left.heel_strike.frame)-1;
          
    
    for c = 1:trial(i).left.num_cycles
        trial(i).left.cycle(c).full_cycle.start.frame = trial(i).left.heel_strike.frame(c);
        trial(i).left.cycle(c).full_cycle.start.time = trial(i).left.heel_strike.time(c);
        trial(i).left.cycle(c).full_cycle.stop.frame = trial(i).left.heel_strike.frame(c+1);       
        trial(i).left.cycle(c).full_cycle.stop.time = trial(i).left.heel_strike.time(c+1);
        
        trial(i).left.cycle(c).stance.start.frame = trial(i).left.heel_strike.frame(c);
        trial(i).left.cycle(c).stance.start.time = trial(i).left.heel_strike.time(c);
        
        if trial(i).left.toe_off.frame(c) > trial(i).left.heel_strike.frame(c)
            trial(i).left.cycle(c).stance.stop.frame = trial(i).left.toe_off.frame(c);       
            trial(i).left.cycle(c).stance.stop.time = trial(i).left.toe_off.time(c);
        else
            trial(i).left.cycle(c).stance.stop.frame = trial(i).left.toe_off.frame(c+1);       
            trial(i).left.cycle(c).stance.stop.time = trial(i).left.toe_off.time(c+1);
        end
        
        
        
        if trial(i).left.toe_off.frame(c) > trial(i).left.heel_strike.frame(c)
            trial(i).left.cycle(c).swing.start.frame = trial(i).left.toe_off.frame(c);
            trial(i).left.cycle(c).swing.start.time = trial(i).left.toe_off.time(c);            
        else
            trial(i).left.cycle(c).swing.start.frame = trial(i).left.toe_off.frame(c+1);
            trial(i).left.cycle(c).swing.start.time = trial(i).left.toe_off.time(c+1);
        end
        
        trial(i).left.cycle(c).swing.stop.frame = trial(i).left.heel_strike.frame(c+1);       
        trial(i).left.cycle(c).swing.stop.time = trial(i).left.heel_strike.time(c+1);
    end
    
    % Plot Thresholded GRF
    
    line_width = 2;
    
    figure('name',[trial(i).name ' Thresholded GRF'])
    % Right
    subplot(2,1,1); hold on
    plot(trial(i).grf.data.time, trial(i).grf.right_vy, 'LineWidth',line_width)
    plot(trial(i).right.heel_strike.time,...
        grf_threshold*ones(trial(i).right.heel_strike.num,1),...
        'o', 'LineWidth',line_width)
    plot(trial(i).right.toe_off.time,...
        grf_threshold*ones(trial(i).right.toe_off.num,1),...
        's', 'LineWidth',line_width)
    ylabel('Force [N]')
    xlabel('Time [s]')
    title('Right GRF')
    
    % Left
    subplot(2,1,2); hold on
    plot(trial(i).grf.data.time, trial(i).grf.left_vy, 'LineWidth',line_width)
    plot(trial(i).left.heel_strike.time, ...
        grf_threshold*ones(trial(i).left.heel_strike.num,1),...
        'o','LineWidth',line_width)
    plot(trial(i).left.toe_off.time, ...
        grf_threshold*ones(trial(i).left.toe_off.num,1),...
        's', 'LineWidth',line_width)
    ylabel('Force [N]')
    xlabel('Time [s]')
    title('Left GRF')    
end

%% Zero the GRF and COP during Swing
for i = 1:num_trials
    
    grf_labels = fieldnames(trial(i).grf.data);
    grf_labels(contains(grf_labels,'time')) =[];
    
    right_grf_labels = grf_labels(contains(grf_labels, '2'));
    left_grf_labels = grf_labels(contains(grf_labels, '1'));
    
    for c = 1:trial(i).right.num_cycles
        swing_start_frame = trial(i).right.cycle(c).swing.start.frame;
        swing_stop_frame = trial(i).right.cycle(c).swing.stop.frame;
        
        for j = 1:length(right_grf_labels)
            trial(i).grf.data.(right_grf_labels{j})...
                (swing_start_frame:swing_stop_frame) = 0;
        end
    end
    
    for c = 1:trial(i).left.num_cycles
        swing_start_frame = trial(i).left.cycle(c).swing.start.frame;
        swing_stop_frame = trial(i).left.cycle(c).swing.stop.frame;
        
        for j = 1:length(left_grf_labels)
            trial(i).grf.data.(left_grf_labels{j})...
                (swing_start_frame:swing_stop_frame) = 0;
        end
    end
    
    % Print new sto file 
     trial(i).zeroed_grf_sto_file = [ './inputs/' trial(i).name '/' trial(i).name '_grf_filt_zeroed.sto'];
%       STOFileAdapter.write(osimTableFromStruct(trial(i).grf.data),trial(i).zeroed_grf_sto_file)
    
     trial(i).zeroed_external_loads_xml_file = [ './inputs/' trial(i).name '/' trial(i).name '_external_loads.xml'];
%        ext_loads = ExternalLoads(trial(i).external_loads_xml_file,true);
%        ext_loads.setDataFileName(trial(i).zeroed_grf_sto_file);
%        ext_loads.print(trial(i).zeroed_external_loads_xml_file);
end


% Plot GRF From All Cycles
figure('name','GRF All Cycles')
for i = 1:num_trials
    for c = 1:trial(i).right.num_cycles
        start_frame = trial(i).right.cycle(c).full_cycle.start.frame;
        stop_frame = trial(i).right.cycle(c).full_cycle.stop.frame;
        
        subplot(4,4,1);hold on
        plot(normcycle(trial(i).grf.data.ground_force_2_vx(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Right Force vx')
        
        subplot(4,4,2);hold on
        plot(normcycle(trial(i).grf.data.ground_force_2_vy(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Right Force vy')
        
        subplot(4,4,3);hold on
        plot(normcycle(trial(i).grf.data.ground_force_2_vz(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Right Force vz')
        
        subplot(4,4,4);hold on
        plot(normcycle(trial(i).grf.data.ground_torque_2_y(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Right Torque y')
        
        subplot(4,4,5);hold on
        plot(normcycle(trial(i).grf.data.ground_force_1_vx(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Left Force vx')
        
        subplot(4,4,6);hold on
        plot(normcycle(trial(i).grf.data.ground_force_1_vy(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Left Force vy')
        
        subplot(4,4,7);hold on
        plot(normcycle(trial(i).grf.data.ground_force_1_vz(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Left Force vz')
        
        subplot(4,4,8);hold on
        plot(normcycle(trial(i).grf.data.ground_torque_1_y(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Left Torque y')
        
        subplot(4,4,9);hold on
        plot(normcycle(trial(i).grf.data.ground_force_2_px(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Right COPx')
        
        subplot(4,4,10);hold on
        plot(normcycle(trial(i).grf.data.ground_force_2_py(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Right COPy')
        
        subplot(4,4,11);hold on
        plot(normcycle(trial(i).grf.data.ground_force_2_pz(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Right COPz')
        
        subplot(4,4,13);hold on
        plot(normcycle(trial(i).grf.data.ground_force_1_px(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('Left COPx')
        
        subplot(4,4,14);hold on
        plot(normcycle(trial(i).grf.data.ground_force_1_py(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('left COPy')
        
        subplot(4,4,15);hold on
        plot(normcycle(trial(i).grf.data.ground_force_1_pz(start_frame:stop_frame)))
        xlabel('Time [s]')
        ylabel('Force [N]')
        title('left COPz')
    end
end

% Plot EMG

% Plot KCF


%% Write OpenSim Input Files for each cycle
% 
% for i = 1:num_trials
%     experimental_data_dir = ...
%         ['../../../../models/knee_tka/grand_challenge/DM/experimental_data/motion_analysis/' trial(i).name '/' cycle_type '/'];
%     
%     trial_trc_table = TimeSeriesTableVec3(trial(i).trc.file);
%     trial_ext_load = ExternalLoads(trial(i).external_loads.file);
%     trial_grf_table = TimeSeriesTable(trial(i).external_loads.file);
%     
%     for c = 1:trial(i).right.num_cycles
%         start_time = trial(i).right.cycle(c).full_cycle.start.time ... 
%                 - cycle_initial_time_pad;
%         stop_time = trial(i).right.cycle(c).full_cycle.stop.time
%         % .trc
%         trc = trial_trc_table.clone();
%                 
%         trc.trim(start_time,stop_time);
%         
%         trial(i).right.cycle(c).full_cycle.trc_file = ...
%             [experimental_data_dir '/' trial(i).name '_' cycle_type '_' int2str(c) '.trc'];
%         
%         TRCFileAdapter.write(trc,trial(i).right.cycle(c).full_cycle.trc_file);
% 
%         % grf.sto
%         grf_table = trial_grf_table.clone();
%         
%         grf_table.trim(start_time,stop_time);
%         
%         trial(i).right.cycle(c).full_cycle.grf_sto_file = ...
%             [experimental_data_dir '/' trial(i).name '_' cycle_type '_' int2str(c) '_grf.sto'];
%         
%         STOFileAdapter.write(trc,trial(i).right.cycle(c).full_cycle.trc_file);
%         
%         % external loads.xml
%     end
% end


%% Perform Inverse Kinematics
% Setup tool and run once for secondary constraint sim 
comak_ik = COMAKInverseKinematicsTool();
comak_ik.set_model_file(model_file);
comak_ik.set_results_directory('./results/comak-inverse-kinematics/');
comak_ik.set_results_prefix('');
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
comak_ik.set_secondary_constraint_sim_sweep_time(3.0);
comak_ik.set_secondary_coupled_coordinate_start_value(0);
comak_ik.set_secondary_coupled_coordinate_stop_value(100);
comak_ik.set_secondary_constraint_sim_integrator_accuracy(1e-3);
comak_ik.set_secondary_constraint_sim_internal_step_limit(10000);
secondary_constraint_function_file = 'secondary_coordinate_constraint_functions.xml';
   
comak_ik.set_secondary_constraint_function_file(['./results/comak-inverse-kinematics/' secondary_constraint_function_file]);
comak_ik.set_constraint_function_num_interpolation_points(20);
comak_ik.set_print_secondary_constraint_sim_results(true);
comak_ik.set_constrained_model_file('ik_constrained_model.osim');

comak_ik.set_perform_inverse_kinematics(false);
% marker_file = trial(i).trc_file;
% comak_ik.set_marker_file(marker_file);
% ik_output_motion_file{i} = [trial(i).name '_ik.sto'];
% comak_ik.set_output_motion_file(ik_output_motion_file{i});
% comak_ik.set_time_range(0, -1);
% comak_ik.set_time_range(1, -1);
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
% tasks(16).Name = 'R_Shank_Superior'; tasks(16).Weight = 1;    
tasks(16).Name = 'R_Shank_Inferior'; tasks(16).Weight = 1;    
tasks(17).Name = 'R_Shank_Lateral'; tasks(17).Weight = 1;    
tasks(18).Name = 'R_Midfoot_Medial'; tasks(18).Weight = 1;    
tasks(19).Name = 'R_Midfoot_Lateral'; tasks(19).Weight = 1;    
tasks(20).Name = 'R_Hindfoot'; tasks(20).Weight = 1;    
tasks(21).Name = 'R_Midfoot_Superior'; tasks(21).Weight = 1;    
tasks(22).Name = 'R_ToeMedial'; tasks(22).Weight = 1;    
tasks(23).Name = 'R_ToeLateral'; tasks(23).Weight = 1;    
tasks(24).Name = 'R_Toe'; tasks(24).Weight = 1;    
tasks(25).Name = 'R_Heel'; tasks(25).Weight = 1;    
% tasks(26).Name = 'L_Patella'; tasks(26).Weight = 1;    
tasks(27-1).Name = 'L_Thigh_Superior'; tasks(27-1).Weight = 1;    
tasks(28-1).Name = 'L_Thigh_Inferior'; tasks(28-1).Weight = 1;    
tasks(29-1).Name = 'L_Thigh_Lateral'; tasks(29-1).Weight = 1;    
tasks(30-1).Name = 'L_Shank_Superior'; tasks(30-1).Weight = 1;    
tasks(31-1).Name = 'L_Shank_Inferior'; tasks(31-1).Weight = 1;    
tasks(32-1).Name = 'L_Shank_Lateral'; tasks(32-1).Weight = 1;    
tasks(33-1).Name = 'L_Midfoot_Medial'; tasks(33-1).Weight = 1;    
tasks(34-1).Name = 'L_Midfoot_Lateral'; tasks(34-1).Weight = 1;    
tasks(35-1).Name = 'L_Hindfoot'; tasks(35-1).Weight = 1;    
tasks(36-1).Name = 'L_Midfoot_Superior'; tasks(36-1).Weight = 1;    
tasks(37-1).Name = 'L_ToeMedial'; tasks(37-1).Weight = 1;    
tasks(38-1).Name = 'L_ToeLateral'; tasks(38-1).Weight = 1;    
tasks(39-1).Name = 'L_Toe'; tasks(39-1).Weight = 1;        
tasks(40-1).Name = 'L_Heel'; tasks(40-1).Weight = 1;

ik_task_set = IKTaskSet();
ik_task = IKMarkerTask();

for j = 1:length(tasks)
    ik_task.setName(tasks(j).Name);
    ik_task.setWeight(tasks(j).Weight);
    ik_task_set.cloneAndAppend(ik_task);
end

comak_ik.set_IKTaskSet(ik_task_set);
comak_ik.print(['./inputs/' trial(i).name '/comak_inverse_kinematics_settings.xml']);

disp(['Running COMAKInverseKinematicsTool: Secondary Constraint Simulation'])
% comak_ik.run();    


for i = 1:num_trials
    if(exist(['./inputs/' trial(i).name],'dir')~=7); mkdir(['./inputs/' trial(i).name]);end
    if(exist(['./results/' trial(i).name],'dir')~=7); mkdir(['./results/' trial(i).name]); end  
    
    results_basename = trial(i).name;
    ik_result_dir{i} = ['./results/' trial(i).name '/comak_inverse_kinematics/']; 
    copyfile(['./results/comak-inverse-kinematics/' secondary_constraint_function_file], ...
        [ik_result_dir{i} secondary_constraint_function_file])
    
%     comak_ik = COMAKInverseKinematicsTool();
    comak_ik.set_model_file(model_file);
    comak_ik.set_results_directory(ik_result_dir(i));
    comak_ik.set_results_prefix(results_basename);
    comak_ik.set_perform_secondary_constraint_sim(false);    
    
%     secondary_constraint_function_file = ...
%         [ik_result_dir{i} '/secondary_coordinate_constraint_functions.xml'];
    comak_ik.set_secondary_constraint_function_file([ik_result_dir{i} secondary_constraint_function_file]);
    comak_ik.set_constrained_model_file('./results/comak-inverse-kinematics/ik_constrained_model.osim');
    
    comak_ik.set_perform_inverse_kinematics(true);
    marker_file = trial(i).trc_file;
    comak_ik.set_marker_file(marker_file);
    ik_output_motion_file{i} = [trial(i).name '_ik.sto'];
    comak_ik.set_output_motion_file(ik_output_motion_file{i});

    comak_ik.print(['./inputs/' trial(i).name '/comak_inverse_kinematics_settings.xml']);

    disp(['Running COMAKInverseKinematicsTool: ' trial(i).name])
%      comak_ik.run();    

    % Inverse Dynamics
    id_result_dir{i} = ['./results/' trial(i).name '/comak_inverse_dynamics/']; 
    id = InverseDynamicsTool();
    id.setResultsDir(id_result_dir{i});
    id.setModelFileName(model_file);
    id.setStartTime(trial(i).grf.data.time(1));
    id.setEndTime(trial(i).grf.data.time(end));
    id.setExcludedForces(ArrayStr('All',1));
    id.setExternalLoadsFileName(trial(i).zeroed_external_loads_xml_file);
    id.setCoordinatesFileName([ik_result_dir{i} '/' ik_output_motion_file{i}]);
    id.setLowpassCutoffFrequency(6);
    id_output_motion_file{i} = [trial(i).name '_id.sto'];
    id.setOutputGenForceFileName([results_basename '/' id_output_motion_file{i}]);
    disp(['Running InverseDynamicsTool: ' trial(i).name])
%       id.run();
end

%% Plot ID Results
id_coords = {'hip_flex_r','knee_flex_r','knee_add_r','ankle_flex_r'};
for i = 1:num_trials
    file = [id_result_dir{i} [trial(i).name '_id.sto']];
    id_results{i} = osimTableToStruct(TimeSeriesTable(file));
    
    figure('name',['ID: ' trial(i).name])
    for j = 1:length(id_coords)
        subplot(4,1,j);hold on;
        plot(id_results{i}.time,id_results{i}.([id_coords{j} '_moment']))
        xlabel('Time [s]')
        ylabel(id_coords{j})
    end
end
        


%% Perform COMAK and Joint Mechanics on each cycle
prescribed_coordinates = {...
    '/jointset/gnd_pelvis/pelvis_tx',...
    '/jointset/gnd_pelvis/pelvis_ty',...
    '/jointset/gnd_pelvis/pelvis_tz',...
    '/jointset/gnd_pelvis/pelvis_tilt',...
    '/jointset/gnd_pelvis/pelvis_list',...
    '/jointset/gnd_pelvis/pelvis_rot',...
    '/jointset/subtalar_r/subt_angle_r',...
    '/jointset/mtp_r/mtp_angle_r',...
    '/jointset/hip_l/hip_flex_l',...
    '/jointset/hip_l/hip_add_l',...
    '/jointset/hip_l/hip_rot_l',...
    '/jointset/pf_l/pf_l_r3',...
    '/jointset/pf_l/pf_l_tx',...
    '/jointset/pf_l/pf_l_ty',...
    '/jointset/knee_l/knee_flex_l',...
    '/jointset/ankle_l/ankle_flex_l',...
    '/jointset/subtalar_l/subt_angle_l',...
    '/jointset/mtp_l/mtp_angle_l',...
    '/jointset/pelvis_torso/lumbar_ext',...
    '/jointset/pelvis_torso/lumbar_latbend',...
    '/jointset/pelvis_torso/lumbar_rot',...
    '/jointset/torso_neckhead/neck_ext',...
    '/jointset/torso_neckhead/neck_latbend',...
    '/jointset/torso_neckhead/neck_rot',...
    '/jointset/acromial_r/arm_add_r',...
    '/jointset/acromial_r/arm_flex_r',...
    '/jointset/acromial_r/arm_rot_r',...
    '/jointset/elbow_r/elbow_flex_r',...
    '/jointset/radioulnar_r/pro_sup_r',...
    '/jointset/radius_hand_r/wrist_flex_r',...
    '/jointset/acromial_l/arm_add_l',...
    '/jointset/acromial_l/arm_flex_l',...
    '/jointset/acromial_l/arm_rot_l',...
    '/jointset/elbow_l/elbow_flex_l',...
    '/jointset/radioulnar_l/pro_sup_l',...
    '/jointset/radius_hand_l/wrist_flex_l'};

primary_coordinates = {...
    '/jointset/hip_r/hip_flex_r',...
    '/jointset/hip_r/hip_add_r',...
    '/jointset/hip_r/hip_rot_r',...
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


for i = 1:num_trials
    for c = 1:trial(i).right.num_cycles
        
        results_basename = [trial(i).name '_fullcycle_' int2str(c)];
        
        initial_time_pad = 0.03;
        
        start_time = trial(i).right.cycle(c).full_cycle.start.time ... 
                - initial_time_pad;
            
        stop_time = trial(i).right.cycle(c).full_cycle.stop.time;

        comak_result_dir = ['./results/' trial(i).name '/comak/']; 
        jnt_mech_result_dir{i} = ['./results/' trial(i).name '/joint-mechanics/']; 
        
        comak = COMAKTool();
        comak.set_model_file(model_file);
        comak.set_coordinates_file([ik_result_dir{i} '/' ik_output_motion_file{i}]);
        comak.set_external_loads_file(trial(i).zeroed_external_loads_xml_file),
        comak.set_results_directory(comak_result_dir);
        comak.set_results_prefix(results_basename);
        comak.set_replace_force_set(false);
        comak.set_force_set_file('../../../../models/knee_tka/grand_challenge/DM/DM_reserve_actuators.xml');
        comak.set_start_time(start_time - 0.1);
        comak.set_stop_time(stop_time);
        comak.set_time_step(0.01);
        comak.set_lowpass_filter_frequency(6);
        comak.set_print_processed_input_kinematics(false);
        
        for j = 1:length(prescribed_coordinates)
            comak.set_prescribed_coordinates(j-1,prescribed_coordinates{j});
        end
        
        for j = 1:length(primary_coordinates)
            comak.set_primary_coordinates(j-1,primary_coordinates{j});
        end

        secondary_coord_set = COMAKSecondaryCoordinateSet();        
        for j = 1:length(secondary_coordinates)
            comak_secondary_coord = COMAKSecondaryCoordinate();
            comak_secondary_coord.setName(secondary_coordinates(j).name);
            comak_secondary_coord.set_max_change(secondary_coordinates(j).max_change);
            comak_secondary_coord.set_coordinate(secondary_coordinates(j).coordinate);
            secondary_coord_set.cloneAndAppend(comak_secondary_coord);
        end
        comak.set_COMAKSecondaryCoordinateSet(secondary_coord_set);
        
        comak.set_settle_secondary_coordinates_at_start(true);
        comak.set_settle_threshold(1e-4);
        comak.set_settle_accuracy(1e-3);
        comak.set_settle_internal_step_limit(10000);
        comak.set_print_settle_sim_results(true);
        
        comak.set_settle_sim_results_directory(comak_result_dir);
        comak.set_settle_sim_results_prefix('validation_tka_settle_sim');
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

        comak.print(['./inputs/' trial(i).name '/fullcycle_' int2str(c) '_comak_settings.xml']);

        disp('Running COMAK Tool...')
%         comak.run();

        % Perform Joint Mechanics Analysis
        jnt_mech = JointMechanicsTool();
        jnt_mech.set_model_file(model_file);
        jnt_mech.set_input_states_file([comak_result_dir '/' results_basename '_states.sto']);
        jnt_mech.set_use_muscle_physiology(false);
        jnt_mech.set_results_file_basename(results_basename);
        jnt_mech.set_results_directory(jnt_mech_result_dir{i});
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
        jnt_mech.set_write_vtp_files(false);
        jnt_mech.set_vtp_file_format('binary');
        jnt_mech.set_write_h5_file(true);
        jnt_mech.set_h5_kinematics_data(true);
        jnt_mech.set_h5_states_data(true);
        jnt_mech.set_use_visualizer(useVisualizer);
        
        % Output kinematics in .h5 file in same reference frame 
        jmft = JointMechanicsFrameTransform();
        jmft.setName('gc_fluoro_frame');
        jmft.set_parent_frame('/bodyset/femur_r/femur_implant');
        jmft.set_child_frame('/bodyset/tibia_r/tibia_implant');
        jmft.set_rotation_type('body');
        jmft.set_rotation_sequence(0,'x');
        jmft.set_rotation_sequence(1,'y');
        jmft.set_rotation_sequence(2,'z');
        jmft.set_output_coordinates(true);
        jmft.set_output_transformation_matrix(false);
        
        jmft_set = JointMechanicsFrameTransformSet();
        jmft_set.cloneAndAppend(jmft);        
        jnt_mech.set_JointMechanicsFrameTransformSet(jmft_set); 
        
        jnt_mech.print(['./inputs/' trial(i).name '/fullcycle_' int2str(c) '_joint_mechanics_settings.xml']);

        disp('Running JointMechanicsTool...');
%          jnt_mech.run();
    end
end




%% Validation Plots
% results\DM_ngait_tmf_slow2\joint-mechanics
% h5_file1 = 'C:/Users/csmith/github/jam-resources/matlab/example/Validation/knee_tka/results/DM_ngait_tmf_slow2/joint-mechanics/DM_ngait_tmf_slow2_fullcycle_1.h5';
% h5_file2 = 'C:/Users/csmith/github/jam-resources/matlab/example/Validation/knee_tka/results/DM_ngait_tmf_slow2/joint-mechanics/DM_ngait_tmf_slow2_fullcycle_2.h5';
fluoro_sto_file = "C:\Users\csmith\github\jam-resources\models\knee_tka\grand_challenge\DM\experimental_data\fluoroscopy\Fluoroscopic Gait Data.sto";
%'../../../../models/knee_tka/grand_challenge/DM/experimental_data/fluoroscopy/Fluoroscopic Gait Data.sto';

numFiles = 0;
for i = 1:1%num_trials
    for c = 1:8%trial(i).right.num_cycles
        numFiles = numFiles + 1;
        h5_file = ['C:/Users/csmith/github/jam-resources/matlab/example/Validation/knee_tka/results/DM_ngait_tmf_slow1/joint-mechanics/' ...
            trial(i).name '_fullcycle_' int2str(c) '.h5' ];
%         h5_file = [jnt_mech_result_dir{i} trial(i).name '_full_cycle_' int2str(c) '.h5' ];
%         h5_file_list{numFiles} = h5_file;
        jam(numFiles) = jam_analysis(h5_file);
    end
end

% jam = jam_analysis(h5_file_list);
fluoro_data = osimTableToStruct(TimeSeriesTable(fluoro_sto_file));


% Knee Kinematics
figure('name','TKA Validation: Tibiofemoral Kinematics')

tf_coords = {'r1','r2','r3','tx','ty','tz'};
fluoro_coords = {'knee_flex_r';'knee_add_r';'knee_rot_r';'knee_ant_r';'knee_sup_r';'knee_lat_r'};

num_coords = length(tf_coords);
for i = 1:num_coords
    subplot(2,3,i);hold on
    if (i<=3)
        for f = 1:numFiles
            plot(jam(f).frametransformsset.coordinates.(tf_coords{i})(:),'LineWidth',line_width)
        end
        plot(fluoro_data.(fluoro_coords{i})* 180 / pi,'k','LineWidth',3)
        
        ylabel('Angle [^o]')
    else
        for f = 1:numFiles
            plot(jam(f).frametransformsset.coordinates.(tf_coords{i})(:),'LineWidth',line_width)
        end
        plot(fluoro_data.(fluoro_coords{i}),'k','LineWidth',3)
        ylabel('Translation [m]')
    end
    
     xlabel('Gait Cycle [%]')
     title(fluoro_coords{i})
end

% legend('cycle 1 COMAK','cycle 2 COMAK','Fluoroscopy Measurement')
% Knee Contact Forces

% EMG

%%




% %% Plot Simulation Results
% 
% % COMAK Inverse Kinematics
% 
% settle_file = [ik_result_dir '/' results_basename '_secondary_constraint_settle_states.sto'];
% sweep_file = [ik_result_dir '/' results_basename '_secondary_constraint_sweep_states.sto'];
% 
% 
% 
% sweep = osimTableToStruct(TimeSeriesTable(sweep_file));
% settle = osimTableToStruct(TimeSeriesTable(sweep_file));
% 
% coord_names = {'/jointset/pf_r/pf_flex_r','/jointset/pf_r/pf_rot_r','/jointset/pf_r/pf_tilt_r';...
%     '/jointset/pf_r/pf_tx_r','/jointset/pf_r/pf_ty_r','/jointset/pf_r/pf_tz_r';...
%     '/jointset/knee_r/knee_flex_r','/jointset/knee_r/knee_add_r','/jointset/knee_r/knee_rot_r';...
%     '/jointset/knee_r/knee_tx_r','/jointset/knee_r/knee_ty_r','/jointset/knee_r/knee_tz_r'}; 	
% 
% groups = {'PF_rot','PF_trans','TF_rot','TF_trans'};
% 
% [m,n] = size(coord_names);
% 
% % for i = 1:m
% %     figure('name',groups{i})
% %     for j = 1:n
% %         
% %         name = coord_names{i,j};
% %     
% %     
% %         ind = find(contains(labels,name));
% %         
% %         subplot(1,3,j);
% % 
% %         hold on
% % %        scatter(n_data(:,1),n_data(:,ind)*180/pi);
% %         plot(data(:,1),data(:,ind));
% % 
% %         title(name);
% %     end
% % end
% 
% 
% figure('name','Secondary Constrain Functions')
% func_set = FunctionSet(secondary_constraint_function_file);
% 
% for f = 0:func_set.getSize()-1
%     func = GCVSpline.safeDownCast(func_set.get(f));
%     
%     for n = 0:func.getNumberOfPoints()-1
%         X(n+1) = func.getX(n);
%         Y(n+1) = func.getY(n);
%     end
%     
%     subplot(4,3,f+1)
%     plot(X,Y)
%     title(char(func.getName()))
%     xlabel('knee flexion [^o]') 
% end
% % Plot Secondary Constraint Functions