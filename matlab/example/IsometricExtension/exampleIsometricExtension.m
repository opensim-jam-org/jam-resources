%% Example Isometric Extension
%==========================================================================
% Author: Colin Smith
%
%Simulation consists of four phases:
%settle1: allow knee to settle into equilbrium
%flex: hip and knee flexion
%settle2: allow knee to settle into equilbrium 
%force: ramp up the muscle forces
%==========================================================================

%% Setup 
import org.opensim.modeling.*
close all; clear

% Settings
Logger.setLevelString('Debug');

useVisualizer = true;
model_file = '../../../models/knee_healthy/smith2019/smith2019.osim';

time.step = 0.01;
time.settle1_duration = 0.01;
time.flex_duration = 0.1;
time.settle2_duration = 0.01;
time.force_duration = 0.01;

inputs.max_hip_flex = 0;
inputs.max_knee_flex = 0;

inputs.max_control = 0.25;
inputs.max_activation = 0.25;
inputs.max_force = 25;

% Create Directories
if(exist('./inputs','dir')~=7); mkdir('./inputs'); end
if(exist('./results','dir')~=7); mkdir('./results'); end

graphics_dir = './results/graphics';
if(exist(graphics_dir,'dir')~=7); mkdir(graphics_dir); end

%% Simulation Time
time.settle1_time = 0 : time.step : time.settle1_duration;
time.flex_time = time.settle1_duration + time.step : time.step : time.settle1_duration + time.flex_duration;
time.settle2_time = time.settle1_duration + time.flex_duration + time.step : time.step : time.settle1_duration + time.flex_duration + time.settle2_duration;
time.force_time = time.settle1_duration + time.flex_duration + time.settle2_duration + time.step : time.step : time.settle1_duration + time.flex_duration + time.settle2_duration + time.force_duration;
time.values = [time.settle1_time, time.flex_time, time.settle2_time, time.force_time];

time_points = [0,time.settle1_duration,...
    time.settle1_duration + time.flex_duration,...
    time.settle1_duration + time.flex_duration + time.settle2_duration,...
    time.settle1_duration + time.flex_duration + time.settle2_duration + time.force_duration];

time.num_settle1_steps = length(time.settle1_time);
time.num_flex_steps = length(time.flex_time);
time.num_settle2_steps = length(time.settle2_time);
time.num_force_steps = length(time.force_time);
time.num_steps = length(time.values);


%% Prescribed Coordinates File
prescribed_coordinates_file = './inputs/prescribed_coordinates.sto';



inputs.hip_flex = [0,0,inputs.max_hip_flex ,inputs.max_hip_flex ,inputs.max_hip_flex ];
inputs.knee_flex = [0,0,inputs.max_hip_flex ,inputs.max_hip_flex ,inputs.max_hip_flex ];

inputs.smooth_hip_flex = interp1(time_points, inputs.hip_flex, time.values,'pchip');
inputs.smooth_knee_flex = interp1(time_points, inputs.knee_flex, time.values,'pchip');

coord_data.hip_flex_r = inputs.smooth_hip_flex';
coord_data.knee_flex_r = inputs.smooth_knee_flex';
coord_data.time = time.values;

coord_table = osimTableFromStruct(coord_data); %% Function distributed in OpenSim 4.0 resources
STOFileAdapter.write(coord_table,prescribed_coordinates_file);

% Prescribed coordinates plot
coord_fig = figure('name','prescribed_coordinates','Position',  [100, 100, 667, 300]);

subplot(1,2,1);
plot(time.values,coord_data.hip_flex_r,'LineWidth',2)
ylim([0.0 40])
xlabel('Time [s]')
ylabel('Angle [^o]')
title('Hip Flexion (hip\_flex\_r)')
box off

subplot(1,2,2);
plot(time.values,coord_data.knee_flex_r,'LineWidth',2)
ylim([0.0 40])
xlabel('Time [s]')
ylabel('Angle [^o]')
title('Knee Flexion (knee\_flex\_r)')
box off

saveas(coord_fig,[graphics_dir '/prescribed_coordinates.png'])

%% Muscle File
muscle_input_sto_file = './inputs/muscle_inputs.sto';



msl_data.time = time.values;

inputs.control_points = [0,0,0,0,inputs.max_control];
inputs.activation_points = [0,0,0,0,inputs.max_activation];
inputs.force_points = [0,0,0,0,inputs.max_force];

inputs.smooth_control = interp1(time_points,inputs.control_points, time.values,'pchip');
inputs.smooth_activation = interp1(time_points,inputs.activation_points, time.values,'pchip');
inputs.smooth_force = interp1(time_points,inputs.force_points, time.values,'pchip');

msl_data.vasmed_r_control = inputs.smooth_control';

msl_data.vaslat_r_activation = inputs.smooth_activation';

msl_data.vasint_r_force = inputs.smooth_force';

msl_table = osimTableFromStruct(msl_data); %% Function distributed in OpenSim 4.0 resources

STOFileAdapter.write(msl_table,muscle_input_sto_file);

% Muscle input plots
msl_fig = figure('name','muscle_inputs','Position',  [100, 100, 1000, 300]);

subplot(1,3,1);
plot(time.values,msl_data.vasmed_r_control,'LineWidth',2)
ylim([0.0 1.0])
xlabel('Time [s]')
ylabel('Control')
title('Vastus Medialis (vasmed\_r)')
box off

subplot(1,3,2);
plot(time.values,msl_data.vaslat_r_activation,'LineWidth',2)
ylim([0.0 1.0])
xlabel('Time [s]')
ylabel('Activation')
title('Vastus Lateralis (vaslat\_r)')
box off

subplot(1,3,3);
plot(time.values,msl_data.vasint_r_force,'LineWidth',2)
ylim([0.0 100.0])
xlabel('Time [s]')
ylabel('Force [N]')
title('Vastus Intermedius (vasint\_r)')
box off

saveas(msl_fig,[graphics_dir '/muscle_inputs.png'])

%% Run Simulation
forsim_result_dir = './results/forsim';
jnt_mech_result_dir = './results/joint-mechanics';
basename = 'healthy';

model = Model(model_file);

forsim = ForsimTool();
forsim.setModel(model);
forsim.set_results_directory(forsim_result_dir);
forsim.set_results_file_basename(basename);
forsim.set_stop_time(-1);
forsim.set_integrator_accuracy(1e-2);%Note this should be 1e-6 for research
%Set all muscles to 2% activation to represent passive state
forsim.set_constant_muscle_control(0.02); 
forsim.set_use_activation_dynamics(true);
forsim.set_use_tendon_compliance(true);
forsim.set_use_muscle_physiology(true);
forsim.set_equilibrate_muscles(true);
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
forsim.set_unconstrained_coordinates(11,'/jointset/meniscus_medial_r/meniscus_medial_flex_r');
forsim.set_unconstrained_coordinates(12,'/jointset/meniscus_medial_r/meniscus_medial_add_r');
forsim.set_unconstrained_coordinates(13,'/jointset/meniscus_medial_r/meniscus_medial_rot_r');
forsim.set_unconstrained_coordinates(14,'/jointset/meniscus_medial_r/meniscus_medial_tx_r');
forsim.set_unconstrained_coordinates(15,'/jointset/meniscus_medial_r/meniscus_medial_ty_r');
forsim.set_unconstrained_coordinates(16,'/jointset/meniscus_medial_r/meniscus_medial_tz_r');
forsim.set_unconstrained_coordinates(17,'/jointset/meniscus_lateral_r/meniscus_lateral_flex_r');
forsim.set_unconstrained_coordinates(18,'/jointset/meniscus_lateral_r/meniscus_lateral_add_r');
forsim.set_unconstrained_coordinates(19,'/jointset/meniscus_lateral_r/meniscus_lateral_rot_r');
forsim.set_unconstrained_coordinates(20,'/jointset/meniscus_lateral_r/meniscus_lateral_tx_r');
forsim.set_unconstrained_coordinates(21,'/jointset/meniscus_lateral_r/meniscus_lateral_ty_r');
forsim.set_unconstrained_coordinates(22,'/jointset/meniscus_lateral_r/meniscus_lateral_tz_r');
forsim.set_prescribed_coordinates_file(prescribed_coordinates_file);

forsim.set_actuator_input_file(muscle_input_sto_file);
forsim.set_use_visualizer(useVisualizer);

analysis_set = AnalysisSet();

frc_reporter = ForceReporter();
frc_reporter.setName('ForceReporter');
analysis_set.cloneAndAppend(frc_reporter);


msl_analysis = MuscleAnalysis();
msl_analysis.setName('MuscleAnalysis');

analysis_set.cloneAndAppend(msl_analysis);
forsim.set_AnalysisSet(analysis_set);

forsim.print('./inputs/healthy_forsim_settings.xml');
disp('Running Forsim Tool...')
forsim.run();

jnt_mech = JointMechanicsTool();
jnt_mech.setModel(model);
jnt_mech.set_use_activation_dynamics(true);
jnt_mech.set_use_tendon_compliance(true);
jnt_mech.set_use_muscle_physiology(true);
jnt_mech.set_input_states_file(...
    [forsim_result_dir '/' basename '_states.sto']);
jnt_mech.set_input_activations_file(...
    [forsim_result_dir '/' basename '_activations.sto']);
jnt_mech.set_input_forces_file(...
    [forsim_result_dir '/' basename '_forces.sto']);
jnt_mech.set_results_file_basename(basename);
jnt_mech.set_results_directory(jnt_mech_result_dir);
jnt_mech.set_stop_time(-1);
jnt_mech.set_normalize_to_cycle(false);
jnt_mech.set_contacts(0,'all');
jnt_mech.set_ligaments(0,'all');
% jnt_mech.set_muscles(0,'all');
jnt_mech.set_muscles(0,'/forceset/vasmed_r');
jnt_mech.set_muscles(1,'/forceset/vaslat_r');
jnt_mech.set_muscles(2,'/forceset/vasint_r');
jnt_mech.set_muscle_outputs(0,'all');
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
jnt_mech.print('./inputs/healthy_joint_mechanics_settings.xml');

disp('Running JointMechanicsTool...');
jnt_mech.run();

%% Analyze Results
states_sto = [forsim_result_dir '/' basename '_states.sto'];
results.states = osimTableToStruct(TimeSeriesTable(states_sto));

frc_report_sto = [forsim_result_dir '/' basename '_ForceReporter_forces.sto'];
results.frc = osimTableToStruct(TimeSeriesTable(frc_report_sto));

results.msl_TendonForce_sto = ...
    [forsim_result_dir '/' basename '_MuscleAnalysis_TendonForce.sto'];
results.msl_ActiveFiberForce_sto = ...
    [forsim_result_dir '/' basename '_MuscleAnalysis_ActiveFiberForce.sto'];
results.msl_PassiveFiberForce_sto = ...
    [forsim_result_dir '/' basename '_MuscleAnalysis_PassiveFiberForce.sto'];
results.msl_TendonLength_sto = ...
    [forsim_result_dir '/' basename '_MuscleAnalysis_TendonLength.sto'];
results.msl_FiberLength_sto = ...
    [forsim_result_dir '/' basename '_MuscleAnalysis_FiberLength.sto'];
results.msl_NormalizedFiberLength_sto = ...
    [forsim_result_dir '/' basename '_MuscleAnalysis_NormalizedFiberLength.sto'];
results.msl_Length_sto = ...
    [forsim_result_dir '/' basename '_MuscleAnalysis_Length.sto'];

results.msl.TendonForce = osimTableToStruct(TimeSeriesTable(results.msl_TendonForce_sto));
results.msl.ActiveFiberForce = osimTableToStruct(TimeSeriesTable(results.msl_ActiveFiberForce_sto));
results.msl.PassiveFiberForce = osimTableToStruct(TimeSeriesTable(results.msl_PassiveFiberForce_sto));
results.msl.TendonLength = osimTableToStruct(TimeSeriesTable(results.msl_TendonLength_sto));
results.msl.FiberLength = osimTableToStruct(TimeSeriesTable(results.msl_FiberLength_sto));
results.msl.NormalizedFiberLength = osimTableToStruct(TimeSeriesTable(results.msl_NormalizedFiberLength_sto));
results.msl.Length = osimTableToStruct(TimeSeriesTable(results.msl_Length_sto));

jam = jam_analysis('smith2019',{[jnt_mech_result_dir '/' basename '.h5']});
%%
% Muscle Results
figure('name','Results: Muscles')

% Control
subplot(5,4,1);hold on;
plot(jam.forceset.Muscle.vasmed_r.excitation)
title('Vastus Medialis')
ylabel('Control')

subplot(5,4,2);hold on;
plot(jam.forceset.Muscle.vaslat_r.excitation)

title('Vastus Lateralis')

subplot(5,4,3);hold on;
plot(jam.forceset.Muscle.vasint_r.excitation)
title('Vastus Intermedius')


% Activation
subplot(5,4,5);hold on;
plot(results.states.a_forceset_vasmed_r_activation,'*')
plot(jam.forceset.Muscle.vasmed_r.activation)
ylabel('Activation')

subplot(5,4,6);hold on;
%plot(results.states.a_forceset_vaslat_r_activation,'*')
plot(jam.forceset.Muscle.vaslat_r.activation)

subplot(5,4,7);hold on;
plot(results.states.a_forceset_vasint_r_activation,'*')
plot(jam.forceset.Muscle.vasint_r.activation)

subplot(5,4,8);hold on;
p1 = plot(0,0,'*');
p2 = plot(0,0);
set(gca,'Visible','off')
leg=legend('states','h5','Location','north');
set(leg, 'Position', get(gca,'Position') );

subplot(5,4,9);hold on;
plot(results.frc.vasmed_r,'*')
plot(results.msl.TendonForce.vasmed_r,'o')
plot(results.msl.ActiveFiberForce.vasmed_r,'o')
plot(results.msl.PassiveFiberForce.vasmed_r,'o')
plot(jam.forceset.Muscle.vasmed_r.passive_fiber_force)
plot(jam.forceset.Muscle.vasmed_r.active_fiber_force)
plot(jam.forceset.Muscle.vasmed_r.tendon_force)
plot(jam.forceset.Muscle.vasmed_r.tension,'^')
plot(jam.forceset.Muscle.vasmed_r.actuation,'s')
ylabel('Force [N]')

subplot(5,4,10);hold on;
plot(results.frc.vaslat_r,'*')
plot(results.msl.TendonForce.vaslat_r,'o')
plot(results.msl.ActiveFiberForce.vaslat_r,'o')
plot(results.msl.PassiveFiberForce.vaslat_r,'o')
plot(jam.forceset.Muscle.vaslat_r.passive_fiber_force)
plot(jam.forceset.Muscle.vaslat_r.active_fiber_force)
plot(jam.forceset.Muscle.vaslat_r.tendon_force)
plot(jam.forceset.Muscle.vaslat_r.tension,'^')
plot(jam.forceset.Muscle.vaslat_r.actuation,'s')

subplot(5,4,11);hold on;
plot(results.frc.vasint_r,'*')
plot(results.msl.TendonForce.vasint_r,'o')
plot(results.msl.ActiveFiberForce.vasint_r,'o')
plot(results.msl.PassiveFiberForce.vasint_r,'o')
plot(jam.forceset.Muscle.vasint_r.passive_fiber_force)
plot(jam.forceset.Muscle.vasint_r.active_fiber_force)
plot(jam.forceset.Muscle.vasint_r.tendon_force)
plot(jam.forceset.Muscle.vasint_r.tension,'^')
plot(jam.forceset.Muscle.vasint_r.actuation,'s')

subplot(5,4,12);hold on;
plot(0,0,'*')
plot(0,0,'o')
plot(0,0,'o')
plot(0,0,'o')
plot(0,0)
plot(0,0)
plot(0,0)
plot(0,0,'^')
plot(0,0,'s')
set(gca,'Visible','off')
leg=legend('frc report','msl-analysis Tendon','msl-analysis Active Fiber','msl-analysis Passive Fiber','h5 Passive Fiber','h5 Active Fiber','h5 Tendon','h5 Tension','h5 Actuation','Location','north');
set(leg, 'Position', get(gca,'Position') );

% Length
subplot(5,4,13);hold on;
plot(results.states.a_forceset_vasmed_r_fiber_length,'*')
plot(results.msl.TendonLength.vasmed_r,'o')
plot(results.msl.FiberLength.vasmed_r,'o')
plot(results.msl.Length.vasmed_r,'o')
plot(jam.forceset.Muscle.vasmed_r.tendon_length)
plot(jam.forceset.Muscle.vasmed_r.fiber_length)
plot(jam.forceset.Muscle.vasmed_r.length)
ylabel('Length')


subplot(5,4,14);hold on;
plot(results.states.a_forceset_vaslat_r_fiber_length,'*')
plot(results.msl.TendonLength.vaslat_r,'o')
plot(results.msl.FiberLength.vaslat_r,'o')
plot(results.msl.Length.vaslat_r,'o')
plot(jam.forceset.Muscle.vaslat_r.tendon_length)
plot(jam.forceset.Muscle.vaslat_r.fiber_length)
plot(jam.forceset.Muscle.vaslat_r.length)

subplot(5,4,15);hold on;
plot(results.states.a_forceset_vasint_r_fiber_length,'*')
plot(results.msl.TendonLength.vasint_r,'o')
plot(results.msl.FiberLength.vasint_r,'o')
plot(results.msl.Length.vasint_r,'o')
plot(jam.forceset.Muscle.vasint_r.tendon_length)
plot(jam.forceset.Muscle.vasint_r.fiber_length)
plot(jam.forceset.Muscle.vasint_r.length)


subplot(5,4,16);hold on;
plot(0,0,'*')
plot(0,0,'o')
plot(0,0,'o')
plot(0,0,'o')
plot(0,0)
plot(0,0)
plot(0,0)
set(gca,'Visible','off')
leg=legend('states','msl-analysis Tendon Length','msl-analysis Fiber Length','msl-analysis Length','h5 Tendon Length','h5 Fiber Length','h5 Length','Location','north');
set(leg, 'Position', get(gca,'Position') );

% leg=legend('states','msl-analysis Tendon Length','msl-analysis Fiber Length','msl-analysis Length','h5 Tendon Length','h5 Fiber Length','h5 Length','Location','eastoutside');
% set(leg, 'Units','normal');
% set(leg, 'Position', get(leg,'Position') .* [1.2 1 1 1]);
% Normalized Length
subplot(5,4,17);hold on;
plot(results.msl.NormalizedFiberLength.vasmed_r,'*')
plot(jam.forceset.Muscle.vasmed_r.normalized_fiber_length)
ylabel('Normalized Fiber Length')
xlabel('Time [s]')

subplot(5,4,18);hold on;
plot(results.msl.NormalizedFiberLength.vaslat_r,'*')
plot(jam.forceset.Muscle.vaslat_r.normalized_fiber_length)
xlabel('Time [s]')

subplot(5,4,19);hold on;
plot(results.msl.NormalizedFiberLength.vasint_r,'*')
plot(jam.forceset.Muscle.vasint_r.normalized_fiber_length)
xlabel('Time [s]')

subplot(5,4,20);hold on;
p1 = plot(0,0,'*');
p2 = plot(0,0);
set(gca,'Visible','off')
leg = legend('msl analysis','h5 file','Location','north');
set(leg, 'Position', get(gca,'Position') );


% set(leg, 'Units','normal');
% set(leg, 'Position', get(leg,'Position') .* [1.2 1 1 1]);
% Secondary Kinematics