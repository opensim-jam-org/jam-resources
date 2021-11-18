%% Simple Muscle Demo
%==========================================================================
% Author: Colin Smith
%
%==========================================================================

%% Settings
clear; close all;
import org.opensim.modeling.*

useVisualizer = true;

time.step = 0.01;
time.hold1_duration = 0.1;
time.slide_duration = 0.1;
time.hold2_duration = 0.05;

inputs.start_displacement = 0.95;
inputs.max_displacement = 1.05;


inputs.max_control = 0.3;
inputs.max_activation = 0.3;
inputs.max_force = 25;

integratorAccuracy = 1e-10;

if(exist('./inputs','dir')~=7); mkdir('./inputs'); end
if(exist('./results','dir')~=7); mkdir('./results'); end
if(exist('./results/default','dir')~=7); mkdir('./results/default'); end

graphics_dir = './results/graphics';
if(exist(graphics_dir,'dir')~=7); mkdir(graphics_dir); end



%% Build Model
model = Model();
model.setName('SimpleMuscleDemo')

% Define Bodies and Joints in the Model
% Get a reference to the model's ground body
ground = model.getGround();

% Attach Anchor geometry to the Ground
% Add offset frames so we can position the geometry
% anchor1Offset = PhysicalOffsetFrame('anchor1', ground, Transform(Vec3(0,0.0,0.0)));
ground.attachGeometry(Brick(Vec3(0.05,0.05,0.6)));

% Add the frames and Geometry
% model.addComponent(anchor1Offset)

% Instantiate a Body with mass, inertia, and a display geometry
block = Body();
block.setName('Block');
block.setMass(20);
block.setMassCenter(Vec3(0));
block.setInertia(Inertia(0.133,0.133,0.133,0,0,0));
% Add display geometry for the block
block.attachGeometry(Brick(Vec3(0.05,0.05,0.6)));


% Instantiate a Slider Joint that connects the block and ground.
blockSideLength      = 0.1;
locationInParentVec3 = Vec3(0, blockSideLength/2, 0);
blockToGround        = SliderJoint('blockToGround', ...
                            ground, Vec3(0), Vec3(0), ...
                            block, Vec3(0), Vec3(0));

% Set bounds on the 6 coordinates of the Free Joint.
positionRange = [0, 2];
blockToGround.upd_coordinates(0).setRange(positionRange);
blockToGround.upd_coordinates(0).setName('displacement');
blockToGround.upd_coordinates(0).setDefaultValue(1);
% Add the block body and joint to the model
model.addBody(block);
model.addJoint(blockToGround);

% Define Muscles in the Model
% Define parameters for a Muscle
maxIsometricForce  = 100.0;
optimalFiberLength = 0.6;
tendonSlackLength  = 0.4;
pennationAngle 	   = 15.0 * pi/180;

% Instantiate a Muscle
muscle1 = Millard2012EquilibriumMuscle();
muscle1.setName('muscle1')
muscle1.setMaxIsometricForce(maxIsometricForce)
muscle1.setOptimalFiberLength(optimalFiberLength)
muscle1.setTendonSlackLength(tendonSlackLength);
muscle1.setPennationAngleAtOptimalFiberLength(pennationAngle)
muscle1.addNewPathPoint('muscle1-point1', ground, Vec3(0.0,0.0,-0.3))
muscle1.addNewPathPoint('muscle1-point2', block, Vec3(0.0,0.0,-0.3))

muscle2 = Millard2012EquilibriumMuscle();
muscle2.setName('muscle2');
muscle2.setMaxIsometricForce(maxIsometricForce)
muscle2.setOptimalFiberLength(optimalFiberLength)
muscle2.setTendonSlackLength(tendonSlackLength)
muscle2.setPennationAngleAtOptimalFiberLength(pennationAngle)
muscle2.addNewPathPoint('muscle2-point1', ground, Vec3(0.0,0.0,-0.1))
muscle2.addNewPathPoint('muscle2-point2', block, Vec3(0.0,0.0,-0.1))

muscle3 = Millard2012EquilibriumMuscle();
muscle3.setName('muscle3')
muscle3.setMaxIsometricForce(maxIsometricForce)
muscle3.setOptimalFiberLength(optimalFiberLength)
muscle3.setTendonSlackLength(tendonSlackLength);
muscle3.setPennationAngleAtOptimalFiberLength(pennationAngle)
muscle3.addNewPathPoint('muscle3-point1', ground, Vec3(0.0,0.0,0.1))
muscle3.addNewPathPoint('muscle3-point2', block, Vec3(0.0,0.0,0.1))

muscle4 = Millard2012EquilibriumMuscle();
muscle4.setName('muscle4')
muscle4.setMaxIsometricForce(maxIsometricForce)
muscle4.setOptimalFiberLength(optimalFiberLength)
muscle4.setTendonSlackLength(tendonSlackLength);
muscle4.setPennationAngleAtOptimalFiberLength(pennationAngle)
muscle4.addNewPathPoint('muscle4-point1', ground, Vec3(0.0,0.0,0.3))
muscle4.addNewPathPoint('muscle4-point2', block, Vec3(0.0,0.0,0.3))

% Add the muscles (as forces) to the model
model.addForce(muscle1)
model.addForce(muscle2);
model.addForce(muscle3);
model.addForce(muscle4);

% Finalize connections so that sockets connectees are correctly saved
model.finalizeConnections();

% Print the model to a XML file (.osim)
model_file = './inputs/simple_muscle_demo.osim';
model.print(model_file);

%% Simulation Time
time.hold1_time = 0 : time.step : time.hold1_duration;
time.slide_time = time.hold1_duration + time.step : time.step : time.hold1_duration + time.slide_duration;
time.hold2_time = time.hold1_duration + time.slide_duration + time.step : time.step : time.hold1_duration + time.slide_duration + time.hold2_duration;
time.values = [time.hold1_time, time.slide_time, time.hold2_time];

time_points = [0,time.hold1_duration,...
    time.hold1_duration + time.slide_duration,...
    time.hold1_duration + time.slide_duration + time.hold2_duration];

time.num_hold1_steps = length(time.hold1_time);
time.num_slide_steps = length(time.slide_time);
time.num_hold2_steps = length(time.hold2_time);
time.num_steps = length(time.values);


%% Prescribed Coordinates File
prescribed_coordinates_file = './inputs/prescribed_coordinates.sto';



inputs.displacement = [inputs.start_displacement,inputs.start_displacement,inputs.max_displacement ,inputs.max_displacement];

inputs.smooth_displacement = interp1(time_points, inputs.displacement, time.values,'pchip');


coord_data.displacement = inputs.smooth_displacement';
coord_data.time = time.values;

coord_table = osimTableFromStruct(coord_data); %% Function distributed in OpenSim 4.0 resources
STOFileAdapter.write(coord_table,prescribed_coordinates_file);

% Prescribed coordinates plot
coord_fig = figure('name','prescribed_coordinates','Position',  [100, 100, 330, 300]);

%subplot(1,2,1);
plot(time.values,coord_data.displacement,'LineWidth',2)
ylim([0.0 2.0])
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Slider Coordinate (displacement)')
box off

saveas(coord_fig,[graphics_dir '/prescribed_coordinates.png'])

%% Write actuator input file
muscle_input_sto_file = './inputs/muscle_inputs.sto';

msl_data.time = time.values;

inputs.control_points = [inputs.max_control,inputs.max_control,inputs.max_control,inputs.max_control];
inputs.activation_points = [inputs.max_activation,inputs.max_activation,inputs.max_activation,inputs.max_activation];
inputs.force_points = [inputs.max_force,inputs.max_force,inputs.max_force,inputs.max_force];

inputs.smooth_control = interp1(time_points,inputs.control_points, time.values,'pchip');
inputs.smooth_activation = interp1(time_points,inputs.activation_points, time.values,'pchip');
inputs.smooth_force = interp1(time_points,inputs.force_points, time.values,'pchip');

msl_data.muscle2_control = inputs.smooth_control';

msl_data.muscle3_activation = inputs.smooth_activation';

msl_data.muscle4_force = inputs.smooth_force';

msl_table = osimTableFromStruct(msl_data); % Function distributed in OpenSim 4.0 resources

STOFileAdapter.write(msl_table,muscle_input_sto_file);

% Muscle input plots
msl_fig = figure('name','muscle_inputs','Position',  [100, 100, 1000, 300]);

subplot(1,3,1);
plot(time.values,msl_data.muscle2_control,'LineWidth',2)
ylim([0.0 1.0])
xlabel('Time [s]')
ylabel('Control')
title('Muscle 2')
box off

subplot(1,3,2);
plot(time.values,msl_data.muscle3_activation,'LineWidth',2)
ylim([0.0 1.0])
xlabel('Time [s]')
ylabel('Activation')
title('Muscle 3')
box off

subplot(1,3,3);
plot(time.values,msl_data.muscle4_force,'LineWidth',2)
ylim([0.0 100.0])
xlabel('Time [s]')
ylabel('Force [N]')
title('Muscle 4')
box off

saveas(msl_fig,[graphics_dir '/muscle_inputs.png'])

%% Perfom Simulation - Default
useMusclePhysiology = true;
useActivationDynamics = true;
useTendonCompliance = true;

default_basename = 'default';
default_forsim_result_dir = './results/default/forsim';
default_jnt_mech_result_dir = './results/default/joint_mechanics';

forsim = ForsimTool();
forsim.set_model_file(model_file);
forsim.set_results_directory(default_forsim_result_dir);
forsim.set_results_file_basename(default_basename);
forsim.set_stop_time(-1);
forsim.set_integrator_accuracy(integratorAccuracy); %Note this should be 1e-6 for research
%forsim.set_constant_muscle_control(0.02); %Set all muscles to 2% activation to represent passive state

forsim.set_use_activation_dynamics(true);
forsim.set_use_tendon_compliance(true);
forsim.set_use_muscle_physiology(true);

forsim.set_equilibrate_muscles(true);
forsim.set_prescribed_coordinates_file(prescribed_coordinates_file);
forsim.set_actuator_input_file(muscle_input_sto_file);
forsim.set_use_visualizer(useVisualizer);
forsim.print('./inputs/forsim_settings_default.xml');

disp('Running Forsim Tool...')
forsim.run();

jnt_mech = JointMechanicsTool();
jnt_mech.set_model_file(model_file);
jnt_mech.set_input_states_file([default_forsim_result_dir '/' default_basename '_states.sto']);
jnt_mech.set_input_activations_file([default_forsim_result_dir '/' default_basename '_activations.sto']);
jnt_mech.set_input_forces_file([default_forsim_result_dir '/' default_basename '_forces.sto']);
jnt_mech.set_results_file_basename(default_basename);
jnt_mech.set_results_directory(default_jnt_mech_result_dir);

jnt_mech.set_use_activation_dynamics(useActivationDynamics);
jnt_mech.set_use_tendon_compliance(useTendonCompliance);
jnt_mech.set_use_muscle_physiology(useMusclePhysiology);

jnt_mech.set_muscles(0,'all');
jnt_mech.set_muscle_outputs(0,'all');
jnt_mech.set_attached_geometry_bodies(0,'/ground');
jnt_mech.set_attached_geometry_bodies(1,'/bodyset/Block');
jnt_mech.set_write_vtp_files(true);
jnt_mech.set_write_h5_file(true);
jnt_mech.set_h5_kinematics_data(true);
jnt_mech.set_h5_states_data(true);
jnt_mech.set_use_visualizer(useVisualizer);
jnt_mech.print('./inputs/joint_mechanics_settings_default.xml');

disp('Running JointMechanicsTool...');
jnt_mech.run();

%% Perform Simulation - Muscle Physiology
% useMusclePhysiology = false;
% useActivationDynamics = true;
% useTendonCompliance = true;
% 
% mp_basename = 'muscle_physiology';
% mp_forsim_result_dir = './results/muscle_physiology/forsim';
% mp_jnt_mech_result_dir = './results/muscle_physiology/joint_mechanics';
% 
% forsim.set_results_directory(mp_basename);
% forsim.set_results_file_basename(mp_forsim_result_dir)
% forsim.set_use_activation_dynamics(useActivationDynamics);
% forsim.set_use_tendon_compliance(useTendonCompliance);
% forsim.set_use_muscle_physiology(useMusclePhysiology);
% forsim.print('./inputs/forsim_settings_muscle_physiology.xml');
% 
% disp('Running Forsim Tool...')
% forsim.run();
% 
% jnt_mech.set_results_directory(mp_basename);
% jnt_mech.set_results_file_basename(mp_jnt_mech_result_dir);
% jnt_mech.set_use_activation_dynamics(useActivationDynamics);
% jnt_mech.set_use_tendon_compliance(useTendonCompliance);
% jnt_mech.set_use_muscle_physiology(useMusclePhysiology);
% jnt_mech.set_write_vtp_files(false);
% jnt_mech.print('./inputs/joint_mechanics_settings_muscle_physiology.xml');
% 
% disp('Running JointMechanicsTool...');
% jnt_mech.run();
% 
% %% Perform Simulation - Activation Dynamics
% useMusclePhysiology = true;
% useActivationDynamics = false;
% useTendonCompliance = true;
% 
% ad_basename = 'activation_dynamics';
% ad_forsim_result_dir = './results/activation_dynamics/forsim';
% ad_jnt_mech_result_dir = './results/activation_dynamics/joint_mechanics';
% 
% forsim.set_results_directory(ad_basename);
% forsim.set_results_file_basename(ad_forsim_result_dir)
% forsim.set_use_activation_dynamics(useActivationDynamics);
% forsim.set_use_tendon_compliance(useTendonCompliance);
% forsim.set_use_muscle_physiology(useMusclePhysiology);
% forsim.print('./inputs/forsim_settings_activation_dynamics.xml');
% 
% disp('Running Forsim Tool...')
% forsim.run();
% 
% jnt_mech.set_results_directory(ad_basename);
% jnt_mech.set_results_file_basename(ad_jnt_mech_result_dir);
% jnt_mech.set_use_activation_dynamics(useActivationDynamics);
% jnt_mech.set_use_tendon_compliance(useTendonCompliance);
% jnt_mech.set_use_muscle_physiology(useMusclePhysiology);
% jnt_mech.set_write_vtp_files(false);
% jnt_mech.print('./inputs/joint_mechanics_settings_activation_dynamics.xml');
% 
% disp('Running JointMechanicsTool...');
% jnt_mech.run();
% 
% %% Perform Simulation - Tendon Compliance
% useMusclePhysiology = true;
% useActivationDynamics = true;
% useTendonCompliance = false;
% 
% tc_basename = 'tendon_compliance';
% tc_forsim_result_dir = './results/tendon_compliance/forsim';
% tc_jnt_mech_result_dir = './results/tendon_compliance/joint_mechanics';
% 
% forsim.set_results_directory(tc_basename);
% forsim.set_results_file_basename(tc_forsim_result_dir)
% forsim.set_use_activation_dynamics(useActivationDynamics);
% forsim.set_use_tendon_compliance(useTendonCompliance);
% forsim.set_use_muscle_physiology(useMusclePhysiology);
% forsim.print('./inputs/forsim_settings_tendon_compliance.xml');
% 
% disp('Running Forsim Tool...')
% forsim.run();
% 
% jnt_mech.set_results_directory(tc_basename);
% jnt_mech.set_results_file_basename(tc_jnt_mech_result_dir);
% jnt_mech.set_use_activation_dynamics(useActivationDynamics);
% jnt_mech.set_use_tendon_compliance(useTendonCompliance);
% jnt_mech.set_use_muscle_physiology(useMusclePhysiology);
% jnt_mech.set_write_vtp_files(false);
% jnt_mech.print('./inputs/joint_mechanics_settings_tendon_compliance.xml');
% 
% disp('Running JointMechanicsTool...');
% jnt_mech.run();
% 
%% Analyze Results
h5_list = {...
    [default_jnt_mech_result_dir '/' default_basename '.h5']};%, ...
%     [mp_jnt_mech_result_dir '/' mp_basename '.h5'], ...
%     [ad_jnt_mech_result_dir '/' ad_basename '.h5'], ...
%     [tc_jnt_mech_result_dir '/' tc_basename '.h5']};

jam = jam_analysis(h5_list);

muscle1_line = 'k';
muscle2_line = 'b';
muscle3_line = 'r';
muscle4_line = 'c';

line_width = 2;

% Default Results
figure('name','Results: Default')

% Activation
subplot(4,3,1);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.activation,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.activation,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.activation,muscle3_line, 'LineWidth', 2)
% plot(jam.time, jam.forceset.Muscle.muscle4.activation,muscle4_line, 'LineWidth', 2)
title('Activation')
ylabel('Activation')

% Normalized Fiber Length
subplot(4,3,2);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.normalized_fiber_length,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.normalized_fiber_length,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.normalized_fiber_length,muscle3_line, 'LineWidth', 2)
% plot(jam.time, jam.forceset.Muscle.muscle4.normalized_fiber_length,muscle4_line, 'LineWidth', 2)
title('Normalized Fiber Length')
ylabel('Norm Length')

% Normalized Fiber Velocity
subplot(4,3,3);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.normalized_fiber_velocity,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.normalized_fiber_velocity,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.normalized_fiber_velocity,muscle3_line, 'LineWidth', 2)
plot(0,0,muscle4_line, 'LineWidth', 2)
% plot(jam.time, jam.forceset.Muscle.muscle4.normalized_fiber_length,muscle4_line, 'LineWidth', 2)
title('Normalized Fiber Velocity')
ylabel('Norm Velocity')
legend('Muscle1','Muscle2 (control)','Muscle3 (activation)','Muscle4 (force)')

% MTU Length
subplot(4,3,4);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.length,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.length,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.length,muscle3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle4.length,muscle4_line, 'LineWidth', 2)
ylabel('Length [m]')
title('Total MTU')

% Fiber Length
subplot(4,3,5);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.fiber_length,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.fiber_length,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.fiber_length,muscle3_line, 'LineWidth', 2)
% plot(jam.time, jam.forceset.Muscle.muscle4.fiber_length,muscle4_line, 'LineWidth', 2)
title('Muscle Fiber')
ylabel('Length [m]')

% Tendon Length
subplot(4,3,6);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.tendon_length,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.tendon_length,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.tendon_length,muscle3_line, 'LineWidth', 2)
% plot(jam.time, jam.forceset.Muscle.muscle4.tendon_length,muscle4_line, 'LineWidth', 2)
title('Tendon')
ylabel('Length [m]')

% MTU velocity
subplot(4,3,7);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.speed,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.speed,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.speed,muscle3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle4.speed,muscle4_line, 'LineWidth', 2)
ylabel('Velocity [m/s]')
title('Total MTU')

% Fiber Velocity
subplot(4,3,8);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.fiber_velocity,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.fiber_velocity,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.fiber_velocity,muscle3_line, 'LineWidth', 2)
% plot(jam.time, jam.forceset.Muscle.muscle4.fiber_velocity,muscle4_line, 'LineWidth', 2)
title('Muscle Fiber')
ylabel('Velocity [m/s]')

% Tendon Velocity
subplot(4,3,9);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.tendon_velocity,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.tendon_velocity,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.tendon_velocity,muscle3_line, 'LineWidth', 2)
% plot(jam.time, jam.forceset.Muscle.muscle4.tendon_velocity,muscle4_line, 'LineWidth', 2)
title('Tendon')
ylabel('Velocity [m/s]')

% Active Force
subplot(4,3,10);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.active_fiber_force,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.active_fiber_force,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.active_fiber_force,muscle3_line, 'LineWidth', 2)
% plot(jam.time, jam.forceset.Muscle.muscle4.active_fiber_force,muscle4_line, 'LineWidth', 2)
title('Active Fiber')
ylabel('Force [N]')
xlabel('Time [s]')

% Passive Force
subplot(4,3,11);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.passive_fiber_force,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.passive_fiber_force,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.passive_fiber_force,muscle3_line, 'LineWidth', 2)
% plot(jam.time, jam.forceset.Muscle.muscle4.passive_fiber_force,muscle4_line, 'LineWidth', 2)
title('Passive Fiber')
ylabel('Force [N]')
xlabel('Time [s]')

% Tendon Force
subplot(4,3,12);hold on;
plot(jam.time, jam.forceset.Muscle.muscle1.tendon_force,muscle1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle2.tendon_force,muscle2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle3.tendon_force,muscle3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Muscle.muscle4.actuation,muscle4_line, 'LineWidth', 2)
title('Tendon (Total MTU)')
ylabel('Force [N]')
xlabel('Time [s]')




