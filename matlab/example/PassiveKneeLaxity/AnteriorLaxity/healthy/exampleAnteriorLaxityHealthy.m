%% Example Anterior Laxity
%==========================================================================
% Author: Colin Smith
%
%
%==========================================================================
close all; clear;

import org.opensim.modeling.*
Logger.setLevelString('Info');

useVisualizer = true;

useActivationDynamics = true;
useTendonCompliance = true;
useMusclePhysiology = true;

% Speed vs accuracy (1e-2 to 1e-8)
integratorAccuracy = 1e-2;


if(exist('./inputs','dir')~=7);mkdir('./inputs');end
if(exist('./results','dir')~=7);mkdir('./results');end
if(exist('./results/graphics','dir')~=7);mkdir('./results/graphics');end

model_file = '../../../../../models/knee_healthy/lenhart2015/lenhart2015.osim';

%% Simulation Time
% Simulation consists of four phases:
% settle : allow knee to settle into equilbrium
% flex   : hip and knee flexion
% settle : allow knee to settle into equilbrium 
% force  : ramp up the anterior force
% settle : hold force constant and allow knee to settle into equilbrium 

time_step = 0.01;

settle1_duration = 0.2;
flex_duration = 0.5;
settle2_duration = 0.2;
force_duration = 0.5;
settle3_duration = 0.5;

settle1_time = 0 : time_step : settle1_duration;

flex_time = ...
    settle1_duration + time_step : ...
    time_step : ...
    settle1_duration + flex_duration;

settle2_time = ...
    settle1_duration + flex_duration + time_step : ...
    time_step : ...
    settle1_duration + flex_duration + settle2_duration;

force_time = ...
    settle1_duration + flex_duration + settle2_duration + time_step : ...
    time_step : ...
    settle1_duration + flex_duration + settle2_duration + force_duration;

settle3_time = ...
    settle1_duration + flex_duration + settle2_duration + force_duration + time_step : ...
    time_step : ...
    settle1_duration + flex_duration + settle2_duration + force_duration + settle3_duration;

time = [settle1_time, flex_time, settle2_time, force_time, settle3_time];

time_points = [0,settle1_duration,...
    settle1_duration + flex_duration,...
    settle1_duration + flex_duration + settle2_duration,...
    settle1_duration + flex_duration + settle2_duration + force_duration, ...
    settle1_duration + flex_duration + settle2_duration + force_duration + settle3_duration];

num_settle1_steps = length(settle1_time);
num_flex_steps = length(flex_time);
num_settle2_steps = length(settle2_time);
num_force_steps = length(force_time);
num_settle3_steps = length(settle3_time);
num_steps = length(time);

%% Create Input Files
% Prescribed Coordinates File
%----------------------------
prescribed_coordinates_file = './inputs/prescribed_coordinates.sto';

max_hip_flex = 25;
max_knee_flex = 25;

coord_data.time = time;

hip_flex = [0,0,max_hip_flex,max_hip_flex,max_hip_flex,max_hip_flex];
knee_flex = [0,0,max_knee_flex,max_knee_flex,max_knee_flex,max_knee_flex];

smooth_hip_flex = interp1(time_points, hip_flex, time,'pchip');
smooth_knee_flex = interp1(time_points, knee_flex, time,'pchip');

coord_data.pelvis_tilt = ones(length(time),1)*90;
coord_data.hip_flex_r = smooth_hip_flex';
coord_data.knee_flex_r = smooth_knee_flex';

% Function distributed in OpenSim Resources\Code\Matlab\Utilities
coord_table = osimTableFromStruct(coord_data); 

STOFileAdapter.write(coord_table,prescribed_coordinates_file);

% Plot Prescribed coordinates 
%----------------------------
coord_fig = figure('name','INPUT: Prescribed Coordinates',...
    'Position',  [100, 100, 667, 300]);

subplot(1,3,1);
plot(time,coord_data.pelvis_tilt,'LineWidth',2)
ylim([0.0 100])
xlabel('Time [s]')
ylabel('Angle [^o]')
title('Pelvis Tilt (pelvis\_tilt\_r)')
box off

subplot(1,3,2);
plot(time,coord_data.hip_flex_r,'LineWidth',2)
ylim([0.0 35])
xlabel('Time [s]')
ylabel('Angle [^o]')
title('Hip Flexion (hip\_flex\_r)')
box off

subplot(1,3,3);
plot(time,coord_data.knee_flex_r,'LineWidth',2)
ylim([0.0 35])
xlabel('Time [s]')
ylabel('Angle [^o]')
title('Knee Flexion (knee\_flex\_r)')
box off

saveas(coord_fig,'./results/graphics/prescribed_coordinates.png')

% External Loads Files
%---------------------
force_magnitude = 100; % 100 N anterior force, similar to KT-1000 arthrometer
force_point_height = -0.1; %Apply at the tibial tuberosity height

% write .sto file
external_loads_sto_file = 'external_loads.sto';

force_vx = [0,0,0,0,force_magnitude,force_magnitude];
smooth_force_vx = interp1(time_points, force_vx, time,'pchip');

force_data.time = time;
force_data.tibia_proximal_r_force_vx = smooth_force_vx';
force_data.tibia_proximal_r_force_vy = zeros(num_steps,1);
force_data.tibia_proximal_r_force_vz = zeros(num_steps,1);
force_data.tibia_proximal_r_force_px = zeros(num_steps,1);
force_data.tibia_proximal_r_force_py = ...
                                    ones(num_steps,1) * force_point_height;
force_data.tibia_proximal_r_force_pz = zeros(num_steps,1);
force_data.tibia_proximal_r_torque_x = zeros(num_steps,1);
force_data.tibia_proximal_r_torque_y = zeros(num_steps,1);
force_data.tibia_proximal_r_torque_z = zeros(num_steps,1);

% Function distributed in OpenSim Resources\Code\Matlab\Utilities
force_table = osimTableFromStruct(force_data); 

force_table.addTableMetaDataString(...
    'header','Anterior Tibial External Force')

STOFileAdapter.write(force_table,['./inputs/' external_loads_sto_file]);

% write .xml file
external_loads_xml_file = './inputs/external_loads.xml';
ext_force = ExternalForce();
ext_force.setName('AnteriorForce');
ext_force.set_applied_to_body('tibia_proximal_r');
ext_force.set_force_expressed_in_body('tibia_proximal_r');
ext_force.set_point_expressed_in_body('tibia_proximal_r');
ext_force.set_force_identifier('tibia_proximal_r_force_v');
ext_force.set_point_identifier('tibia_proximal_r_force_p');
ext_force.set_torque_identifier('tibia_proximal_r_torque_');

ext_loads = ExternalLoads();
ext_loads.setDataFileName(external_loads_sto_file);
ext_loads.cloneAndAppend(ext_force);
ext_loads.print(external_loads_xml_file );

% Plot External Loads
%--------------------
ext_loads_fig = figure('name','INPUT: External Loads', ...
                'Position',  [100, 100, 333, 300]);
plot(time,force_data.tibia_proximal_r_force_vx,'LineWidth',2)
ylim([0.0 100])
xlabel('Time [s]')
ylabel('Anterior Force [N]')
title('External Loads on the Tibia')
box off

saveas(ext_loads_fig,'./results/graphics/external_loads.png')

% Perform Simulation with ForsimTool
Healthy
forsim_result_dir = './results/forsim';
basename = 'healthy';

forsim = ForsimTool();
forsim.set_model_file(model_file);
forsim.set_results_directory(forsim_result_dir);
forsim.set_results_file_basename(basename);
forsim.set_start_time(-1);
forsim.set_stop_time(-1);
forsim.set_integrator_accuracy(integratorAccuracy);
%Set all muscles to 2% activation to represent passive state
forsim.set_constant_muscle_control(0.02); 
forsim.set_use_activation_dynamics(useActivationDynamics);
forsim.set_use_tendon_compliance(useTendonCompliance);
forsim.set_use_muscle_physiology(useMusclePhysiology)
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
forsim.set_prescribed_coordinates_file(prescribed_coordinates_file);
forsim.set_external_loads_file(external_loads_xml_file);
forsim.set_use_visualizer(useVisualizer);
forsim.print('./inputs/forsim_settings.xml');

disp('Running Forsim Tool...')
forsim.run();

%% Perform Analysis with JointMechanicsTool
%Healthy
jnt_mech_result_dir = './results/joint_mechanics';


jnt_mech = JointMechanicsTool();
jnt_mech.set_model_file(model_file);
jnt_mech.set_input_states_file(...
    [forsim_result_dir '/' basename '_states.sto']);
jnt_mech.set_results_file_basename(basename);
jnt_mech.set_results_directory(jnt_mech_result_dir);
jnt_mech.set_use_activation_dynamics(useActivationDynamics);
jnt_mech.set_use_tendon_compliance(useTendonCompliance);
jnt_mech.set_use_muscle_physiology(useMusclePhysiology);
jnt_mech.set_start_time(-1);
jnt_mech.set_stop_time(-1);
jnt_mech.set_normalize_to_cycle(false);
jnt_mech.set_contacts(0,'all');
jnt_mech.set_ligaments(0,'all');
jnt_mech.set_muscles(0,'all');
jnt_mech.set_muscle_outputs(0,'all');
jnt_mech.set_attached_geometry_bodies(0,'/bodyset/femur_distal_r');
jnt_mech.set_attached_geometry_bodies(1,'/bodyset/tibia_proximal_r');
jnt_mech.set_attached_geometry_bodies(2,'/bodyset/patella_r');
jnt_mech.set_output_orientation_frame('ground');
jnt_mech.set_output_position_frame('ground');
jnt_mech.set_write_vtp_files(false);
jnt_mech.set_write_h5_file(true);
jnt_mech.set_h5_kinematics_data(true);
jnt_mech.set_h5_states_data(true);
jnt_mech.set_use_visualizer(useVisualizer);
jnt_mech.print('./inputs/joint_mechanics_settings.xml');

disp('Running JointMechanicsTool...');
jnt_mech.run();

%% Plot Results

jam = jam_analysis({[jnt_mech_result_dir '/' basename '.h5']});


   
% TF Kinematics
coords = {'knee_flex_r','knee_add_r','knee_rot_r','knee_tx_r','knee_ty_r','knee_tz_r'};
num_coords = length(coords);
tf_kinematics_fig = figure('name','TF Kinematics');

for i = 1:num_coords
    subplot(2,3,i); hold on
    plot(jam.time, jam.coordinateset.(coords{i}).value,'LineWidth',2);

    xlabel('Time [s]')
    ylabel(strrep(coords{i},'_',' '))
end
saveas(tf_kinematics_fig,'./results/graphics/tibiofemoral_kinematics.png')

% PF Kinematics
coords = {'pf_flex_r','pf_rot_r','pf_tilt_r','pf_tx_r','pf_ty_r','pf_tz_r'};
num_coords = length(coords);
pf_kinematics_fig = figure('name','PF Kinematics');

for i = 1:num_coords
    subplot(2,3,i); hold on
    plot(jam.time, jam.coordinateset.(coords{i}).value,'LineWidth',2);

    xlabel('Time [s]')
    ylabel(strrep(coords{i},'_',' '))

end
saveas(pf_kinematics_fig,'./results/graphics/patellofemoral_kinematics.png')

% Plot Ligament Forces 
ligament_names = {'ACLam','ACLpl','PT','MCL','ITB1','LCL','PCL','pCAP'};
lig_sym = {'-sk','-^k','g','r','b','m','c','y'};

lig_force_fig = figure('name','Ligament Loading');

for i = 1:length(ligament_names)
    subplot(3,3,i);hold on;
    fiber_names = fieldnames(jam.forceset.Blankevoort1991Ligament);
    fibers = fiber_names(contains(fiber_names,ligament_names{i}));       

    data = 0;
    for k = 1:length(fibers)
        data = data + jam.forceset.Blankevoort1991Ligament.(fibers{k}).('total_force');
    end
    plot(jam.time, data, 'LineWidth',2)
    title(ligament_names{i})
    xlabel('Time [s]')
    ylabel('Force [N]')
end
saveas(lig_force_fig,'./results/graphics/ligament_forces.png')


% Plot Muscle Forces
% for m = 0:model.getMuscles().getSize()-1
%     muscles{m+1} = char(msl_set.get(m).getName());
% end
% msl_frc_fig = figure('name', 'Muscle Forces');
% hold on;
% h=gca;
% for i = 1:length(muscles)
%     jam.plot_muscle_output(gca,muscles{i},'actuation');
% end
muscles = {'vasmed_r','vaslat_r','vasint_r','recfem_r',...
    'bfsh_r','bflh_r','semimem_r','semiten_r',...
    'gasmed_r','gaslat_r','soleus_r','tibpost_r'};

msl_force_fig = figure('name', 'Muscle Forces');

hold on;
h=gca;
for i = 1:length(muscles)
    subplot(3,4,i);
    jam.plot_muscle_output(gca,muscles{i},'actuation');
    title(strrep(muscles{i},'_',' '))
    xlabel('Time [s]')
    ylabel('Force [N]')
end
saveas(msl_force_fig,'./results/graphics/muscle_forces.png')


%% Plot Contact Force
cnt_names = {'tf_contact','pf_contact'};
mesh_names = {'tibia_cartilage','patella_cartilage'};
num_cnt = length(cnt_names);

cmp_names = {'Force X [N]','Force Y [N]','Force Z [N]'};

cnt_force_fig = figure('name','Contact Forces');
s=1;

for j = 1:3 %xyz force components
    for i = 1:num_cnt 
        subplot(3,num_cnt,s);hold on;
        plot(jam.time, jam.forceset.Smith2018ArticularContactForce.(cnt_names{i}).(mesh_names{i}).total_contact_force(:,j))

        ylabel(cmp_names{j})
        xlabel('Time [s]')
        
        if(j==1)
            title(strrep(cnt_names{i},'_',' '))
        end
        s=s+1;
    end    
end
saveas(cnt_force_fig,'./results/graphics/contact_forces.png')