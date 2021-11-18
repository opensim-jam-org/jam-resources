%% Simple Ligament Demo
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
time.hold2_duration = 0.2;

inputs.start_displacement = 0.1;
inputs.max_displacement = 0.12;


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

ground = model.getGround();
ground.attachGeometry(Brick(Vec3(0.01,0.01,0.06)));

% Instantiate a Body with mass, inertia, and a display geometry
block = Body();
block.setName('Block');
block.setMass(20);
block.setMassCenter(Vec3(0));
block.setInertia(Inertia(0.133,0.133,0.133,0,0,0));
block.attachGeometry(Brick(Vec3(0.01,0.01,0.06)));

blockToGround = SliderJoint('blockToGround', ...
                    ground, Vec3(0), Vec3(0), ...
                    block, Vec3(0), Vec3(0));

% Set bounds on the 6 coordinates of the Free Joint.
positionRange = [0, 0.2];
blockToGround.upd_coordinates(0).setRange(positionRange);
blockToGround.upd_coordinates(0).setName('displacement');
blockToGround.upd_coordinates(0).setDefaultValue(0.1);
% Add the block body and joint to the model
model.addBody(block);
model.addJoint(blockToGround);

% Add Ligaments
lig1_pt1 = Vec3(0.0, 0.0, -0.03);
lig1_pt2 = Vec3(0.0, 0.0, -0.03);
lig1 = Blankevoort1991Ligament('ligament1',ground,lig1_pt1,block,lig1_pt2);
model.addForce(lig1);

lig2_pt1 = Vec3(0.0, 0.0, -0.01);
lig2_pt2 = Vec3(0.0, 0.0, -0.01);
lig2 = Blankevoort1991Ligament('ligament2',ground,lig2_pt1,block,lig2_pt2);
model.addForce(lig2);

lig3_pt1 = Vec3(0.0, 0.0, 0.01);
lig3_pt2 = Vec3(0.0, 0.0, 0.01);
lig3 = Blankevoort1991Ligament('ligament3',ground,lig3_pt1,block,lig3_pt2);
model.addForce(lig3);

lig4_pt1 = Vec3(0.0, 0.0, 0.03);
lig4_pt2 = Vec3(0.0, 0.0, 0.03);
lig4 = Blankevoort1991Ligament('ligament4',ground,lig4_pt1,block,lig4_pt2);
model.addForce(lig4);

state = model.initSystem();

lig1.set_slack_length(0.1);
lig1.set_linear_stiffness(100);

lig2.set_slack_length(0.1);
lig2.set_linear_stiffness(100);
lig2.set_damping_coefficient(0.006);

% Need to set linear stiffness before using setSlackLengthFromXX functions
lig3.set_linear_stiffness(100);
lig3.setSlackLengthFromReferenceStrain(0.03,state);

lig4.set_linear_stiffness(100);
lig4.setSlackLengthFromReferenceForce(5,state);


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
% ylim([0.0 2.0])
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Slider Coordinate (displacement)')
box off

saveas(coord_fig,[graphics_dir '/prescribed_coordinates.png'])

%% Perfom Simulation - Default
basename = 'ligament_demo';
forsim_result_dir = './results/forsim';
jnt_mech_result_dir = './results/joint_mechanics';

forsim = ForsimTool();
forsim.set_model_file(model_file);
forsim.set_results_directory(forsim_result_dir);
forsim.set_results_file_basename(basename);
forsim.set_stop_time(-1);
forsim.set_integrator_accuracy(integratorAccuracy);
forsim.set_maximum_time_step(1e-4);
forsim.set_prescribed_coordinates_file(prescribed_coordinates_file);
forsim.set_use_visualizer(useVisualizer);
forsim.print('./inputs/forsim_settings_default.xml');

disp('Running Forsim Tool...')
forsim.run();

jnt_mech = JointMechanicsTool();
jnt_mech.set_model_file(model_file);
jnt_mech.set_input_states_file([forsim_result_dir '/' basename '_states.sto']);
jnt_mech.set_results_file_basename(basename);
jnt_mech.set_results_directory(jnt_mech_result_dir);
jnt_mech.set_ligaments(0,'all');
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


%% Analyze Results
h5_list = {...
    [jnt_mech_result_dir '/' basename '.h5']};%, ...
%     [mp_jnt_mech_result_dir '/' mp_basename '.h5'], ...
%     [ad_jnt_mech_result_dir '/' ad_basename '.h5'], ...
%     [tc_jnt_mech_result_dir '/' tc_basename '.h5']};

jam = jam_analysis(h5_list);

lig1_line = 'k';
lig2_line = 'b';
lig3_line = 'r';
lig4_line = 'c';

line_width = 2;

% Default Results
figure('name','Results')

% Length
subplot(3,3,1);hold on;
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament1.length,lig1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament2.length,lig2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament3.length,lig3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament4.length, lig4_line, 'LineWidth', 2)
title('Length')
ylabel('Length [m]')

% Strain
subplot(3,3,2);hold on;
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament1.strain,lig1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament2.strain,lig2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament3.strain,lig3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament4.strain, lig4_line, 'LineWidth', 2)
title('Strain')
ylabel('Strain')

% Spring Force
subplot(3,3,3);hold on;
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament1.spring_force,lig1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament2.spring_force,lig2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament3.spring_force,lig3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament4.spring_force, lig4_line, 'LineWidth', 2)
title('Spring Force')
ylabel('Force [N]')
legend('ligament1','ligament2 (damping)','ligament3 (ref strain)','ligament4 (ref force)','Location','southeast')

% Lengthening Speed
subplot(3,3,4);hold on;
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament1.lengthening_speed,lig1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament2.lengthening_speed,lig2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament3.lengthening_speed,lig3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament4.lengthening_speed, lig4_line, 'LineWidth', 2)
title('Lengthening Speed')
ylabel('Speed [m/s]')

% Strain Rate
subplot(3,3,5);hold on;
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament1.strain,lig1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament2.strain,lig2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament3.strain,lig3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament4.strain, lig4_line, 'LineWidth', 2)
title('Strain Rate')
ylabel('Strain Rate [1/s]')

% Damping Force
subplot(3,3,6);hold on;
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament1.damping_force,lig1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament2.damping_force,lig2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament3.damping_force,lig3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament4.damping_force, lig4_line, 'LineWidth', 2)
title('Damping Force')
ylabel('Force [N]')

% Potential Energy 
subplot(3,3,7);hold on;
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament1.potential_energy,lig1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament2.potential_energy,lig2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament3.potential_energy,lig3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament4.potential_energy, lig4_line, 'LineWidth', 2)
title('Potential Energy')
ylabel('Energy [J]')
xlabel('Time [s]')

% Displacement 
subplot(3,3,8);hold on;
plot(jam.time, jam.coordinateset.displacement.value,lig1_line, 'LineWidth', 2)
plot(jam.time, jam.coordinateset.displacement.value,lig2_line, 'LineWidth', 2)
plot(jam.time, jam.coordinateset.displacement.value,lig3_line, 'LineWidth', 2)
plot(jam.time, jam.coordinateset.displacement.value, lig4_line, 'LineWidth', 2)
title('Displacement')
ylabel('Distance [m]')
xlabel('Time [s]')

% Total Force
subplot(3,3,9);hold on;
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament1.total_force,lig1_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament2.total_force,lig2_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament3.total_force,lig3_line, 'LineWidth', 2)
plot(jam.time, jam.forceset.Blankevoort1991Ligament.ligament4.total_force, lig4_line, 'LineWidth', 2)
title('Total Force')
ylabel('Force [N]')
xlabel('Time [s]')




