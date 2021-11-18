%% Example Visualize Kinematics
%==========================================================================
%
%==========================================================================
clear; close all;
%Setup Environment and Folders
import org.opensim.modeling.*
Logger.setLevelString('Info');

useVisualizer = true;

if(exist('./inputs','dir')~=7)
    mkdir('./inputs')
end
if(exist('./inputs/Geometry','dir')~=7)
    mkdir('./inputs/Geometry')
end
if(exist('./results','dir')~=7)
    mkdir('./results')
end
if(exist('./results/passive','dir')~=7)
    mkdir('./results/passive')
end
if(exist('./results/inertia','dir')~=7)
    mkdir('./results/inertia')
end
if(exist('./results/elastic','dir')~=7)
    mkdir('./results/elastic')
end
if(exist('./results/graphics','dir')~=7)
    mkdir('./results/graphics')
end

%% Build Model

% Copy Geometries to example dir
geometry_dir = '../../../../models/knee_healthy/lenhart2015/Geometry/';
local_geometry_dir = './inputs/Geometry/';

femur_bone_mesh_file = 'lenhart2015-R-femur-bone.stl';
tibia_bone_mesh_file = 'lenhart2015-R-tibia-bone.stl';
patella_bone_mesh_file = 'lenhart2015-R-patella-bone.stl';
femur_cartilage_mesh_file = 'lenhart2015-R-femur-cartilage.stl';
tibia_cartilage_mesh_file = 'lenhart2015-R-tibia-cartilage.stl';
patella_cartilage_mesh_file = 'lenhart2015-R-patella-cartilage.stl';

% femur_cartilage_mesh_file = 'femur-cart_tune_scaled_align.stl';
% tibia_cartilage_mesh_file = 'tibia-cart_tune_scaled_align.stl';

% femur_bone_mesh_file = 'ACLC01_R_Femur.stl';
% tibia_bone_mesh_file = 'ACLC01_R_Tibia.stl';
% patella_bone_mesh_file = 'ACLC01_R_Patella.stl';
% femur_cartilage_mesh_file = 'ACLC01_R_FemurCartilage.stl';
% tibia_cartilage_mesh_file = 'ACLC01_R_TibiaCartilage.stl';
% patella_cartilage_mesh_file = 'ACLC01_R_PatellaCartilage.stl';

% femur_bone_mesh_file = 'femur-cart_tune.stl';
% tibia_bone_mesh_file = 'tibia-cart_tune.stl';
% patella_bone_mesh_file = 'ACLC01_R_Patella.stl';
% femur_cartilage_mesh_file = 'femur-cart_tune.stl';
% tibia_cartilage_mesh_file = 'tibia-cart_tune.stl';
% patella_cartilage_mesh_file = 'ACLC01_R_PatellaCartilage.stl';


copyfile([geometry_dir femur_bone_mesh_file],[local_geometry_dir femur_bone_mesh_file])
copyfile([geometry_dir tibia_bone_mesh_file],[local_geometry_dir tibia_bone_mesh_file])
copyfile([geometry_dir patella_bone_mesh_file],[local_geometry_dir patella_bone_mesh_file])

copyfile([geometry_dir femur_cartilage_mesh_file],[local_geometry_dir femur_cartilage_mesh_file])
copyfile([geometry_dir tibia_cartilage_mesh_file],[local_geometry_dir tibia_cartilage_mesh_file])
copyfile([geometry_dir patella_cartilage_mesh_file],[local_geometry_dir patella_cartilage_mesh_file])

ModelVisualizer.addDirToGeometrySearchPaths(local_geometry_dir);

model = Model();
modelName = 'knee';
model.setName(modelName);

% Add Bodies
% Mass properties are arbitrary as the model will be kinematically driven,
% so no dynamic simulation is performed
massCenter = Vec3(0.0);                       
inertia = Inertia(0.1, 0.1, 0.1, 0, 0, 0);

femur  = Body();
femur.setName('femur_distal_r');
femur.setMass(1.0);                                      
femur.setMassCenter(massCenter);                        
femur.setInertia(inertia); 
femur.attachGeometry(Mesh(femur_bone_mesh_file));
model.addBody(femur);

tibia  = Body();
tibia.setName('tibia_proximal_r');
tibia.setMass(1.0);                                      
tibia.setMassCenter(massCenter);                        
tibia.setInertia(inertia);
tibia.attachGeometry(Mesh(tibia_bone_mesh_file));
model.addBody(tibia);

patella  = Body();
patella.setName('patella_r');
patella.setMass(1.0);                                      
patella.setMassCenter(massCenter);                        
patella.setInertia(inertia); 
patella.attachGeometry(Mesh(patella_bone_mesh_file));
model.addBody(patella);


%Add Joints

ground_transform = SpatialTransform();

ground_transform.get_rotation1().set_coordinates(0,'ground_femur_rx_r');
ground_transform.get_rotation1().set_axis(Vec3(1,0,0));
ground_transform.get_rotation1().set_function(LinearFunction());

ground_transform.get_rotation2().set_coordinates(0,'ground_femur_ry_r');
ground_transform.get_rotation2().set_axis(Vec3(0,1,0));
ground_transform.get_rotation2().set_function(LinearFunction());

ground_transform.get_rotation3().set_coordinates(0,'ground_femur_rz_r');
ground_transform.get_rotation3().set_axis(Vec3(0,0,1));
ground_transform.get_rotation3().set_function(LinearFunction());

ground_transform.get_translation1().set_coordinates(0,'ground_femur_tx_r');
ground_transform.get_translation1().set_axis(Vec3(1,0,0));
ground_transform.get_translation1().set_function(LinearFunction());

ground_transform.get_translation2().set_coordinates(0,'ground_femur_ty_r');
ground_transform.get_translation2().set_axis(Vec3(0,1,0));
ground_transform.get_translation2().set_function(LinearFunction());

ground_transform.get_translation3().set_coordinates(0,'ground_femur_tz_r');
ground_transform.get_translation3().set_axis(Vec3(0,0,1));
ground_transform.get_translation3().set_function(LinearFunction());

ground_femur_r = CustomJoint('ground_femur_r',...
    model.getGround(),Vec3(0.0),Vec3(0.0),femur,Vec3(0.0),Vec3(0.0),...
    ground_transform);

model.addJoint(ground_femur_r);

knee_transform = SpatialTransform();

knee_transform.get_rotation1().set_coordinates(0,'knee_flex_r');
knee_transform.get_rotation1().set_axis(Vec3(0,0,1));
knee_transform.get_rotation1().set_function(LinearFunction(-1,0));

knee_transform.get_rotation2().set_coordinates(0,'knee_add_r');
knee_transform.get_rotation2().set_axis(Vec3(1,0,0));
knee_transform.get_rotation2().set_function(LinearFunction());

knee_transform.get_rotation3().set_coordinates(0,'knee_rot_r');
knee_transform.get_rotation3().set_axis(Vec3(0,1,0));
knee_transform.get_rotation3().set_function(LinearFunction());

knee_transform.get_translation1().set_coordinates(0,'knee_ant_r');
knee_transform.get_translation1().set_axis(Vec3(1,0,0));
knee_transform.get_translation1().set_function(LinearFunction());

knee_transform.get_translation2().set_coordinates(0,'knee_sup_r');
knee_transform.get_translation2().set_axis(Vec3(0,1,0));
knee_transform.get_translation2().set_function(LinearFunction());

knee_transform.get_translation3().set_coordinates(0,'knee_lat_r');
knee_transform.get_translation3().set_axis(Vec3(0,0,1));
knee_transform.get_translation3().set_function(LinearFunction());

knee_r = CustomJoint('knee_r',...
    femur,Vec3(0.0),Vec3(0.0),tibia,Vec3(0.0),Vec3(0.0),...
    knee_transform);

model.addJoint(knee_r);

knee_flex = model.getCoordinateSet().get('knee_flex_r');
knee_flex.set_range(0,-pi/4);
knee_flex.set_range(1,120*pi/180);
knee_flex.set_clamped(true);

knee_add = model.getCoordinateSet().get('knee_add_r');
knee_add.set_range(0,-10*pi/180);
knee_add.set_range(1,10*pi/180);
knee_add.set_clamped(true);

knee_rot = model.getCoordinateSet().get('knee_rot_r');
knee_rot.set_range(0,-30*pi/180);
knee_rot.set_range(1,30*pi/180);
knee_rot.set_clamped(true)

knee_tx = model.getCoordinateSet().get('knee_ant_r');
knee_tx.set_range(0,-0.2);
knee_tx.set_range(1,0.2);
knee_tx.set_clamped(true);

knee_ty = model.getCoordinateSet().get('knee_sup_r');
knee_ty.set_range(0,-0.2);
knee_ty.set_range(1,0.2);
knee_ty.set_clamped(true);

knee_tz = model.getCoordinateSet().get('knee_lat_r');
knee_tz.set_range(0,-0.2);
knee_tz.set_range(1,0.2);
knee_tz.set_clamped(true);

pf_transform = SpatialTransform();

pf_transform.get_rotation1().set_coordinates(0,'pf_flex_r');
pf_transform.get_rotation1().set_axis(Vec3(0,0,1));
pf_transform.get_rotation1().set_function(LinearFunction());

pf_transform.get_rotation2().set_coordinates(0,'pf_rot_r');
pf_transform.get_rotation2().set_axis(Vec3(1,0,0));
pf_transform.get_rotation2().set_function(LinearFunction());

pf_transform.get_rotation3().set_coordinates(0,'pf_tilt_r');
pf_transform.get_rotation3().set_axis(Vec3(0,1,0));
pf_transform.get_rotation3().set_function(LinearFunction());

pf_transform.get_translation1().set_coordinates(0,'pf_ant_r');
pf_transform.get_translation1().set_axis(Vec3(1,0,0));
pf_transform.get_translation1().set_function(LinearFunction());

pf_transform.get_translation2().set_coordinates(0,'pf_sup_r');
pf_transform.get_translation2().set_axis(Vec3(0,1,0));
pf_transform.get_translation2().set_function(LinearFunction());

pf_transform.get_translation3().set_coordinates(0,'pf_lat_r');
pf_transform.get_translation3().set_axis(Vec3(0,0,1));
pf_transform.get_translation3().set_function(LinearFunction());

pf_r = CustomJoint('pf_r',...
    femur,Vec3(0.0),Vec3(0.0),patella,Vec3(0.0),Vec3(0.0),...
    pf_transform);

model.addJoint(pf_r);

pf_flex = model.getCoordinateSet().get('pf_flex_r');
pf_flex.set_range(0,-pi/4);
pf_flex.set_range(1,120*pi/180);
pf_flex.set_clamped(true);

pf_add = model.getCoordinateSet().get('pf_rot_r');
pf_add.set_range(0,-10*pi/180);
pf_add.set_range(1,10*pi/180);
pf_add.set_clamped(true);

pf_rot = model.getCoordinateSet().get('pf_tilt_r');
pf_rot.set_range(0,-30*pi/180);
pf_rot.set_range(1,30*pi/180);
pf_rot.set_clamped(true)

pf_tx = model.getCoordinateSet().get('pf_ant_r');
pf_tx.set_range(0,-0.2);
pf_tx.set_range(1,0.2);
pf_tx.set_clamped(true);

pf_ty = model.getCoordinateSet().get('pf_sup_r');
pf_ty.set_range(0,-0.2);
pf_ty.set_range(1,0.2);
pf_ty.set_clamped(true);

pf_tz = model.getCoordinateSet().get('pf_lat_r');
pf_tz.set_range(0,-0.2);
pf_tz.set_range(1,0.2);
pf_tz.set_clamped(true);

femur_cnt_mesh = Smith2018ContactMesh('femur_cartilage',...
                    femur_cartilage_mesh_file,femur);
model.addContactGeometry(femur_cnt_mesh);
                

tibia_cnt_mesh = Smith2018ContactMesh('tibia_cartilage',...
                    tibia_cartilage_mesh_file,tibia);
model.addContactGeometry(tibia_cnt_mesh);


patella_cnt_mesh = Smith2018ContactMesh('patella_cartilage',...
                    patella_cartilage_mesh_file,patella);
model.addContactGeometry(patella_cnt_mesh);

tf_contact = Smith2018ArticularContactForce('tf_contact',...
                femur_cnt_mesh,tibia_cnt_mesh);
model.addForce(tf_contact);

pf_contact = Smith2018ArticularContactForce('pf_contact',...
                femur_cnt_mesh,patella_cnt_mesh);
model.addForce(pf_contact);

%Add MCL Ligament
mcl_wrap = WrapEllipsoid();
mcl_wrap.set_active(true);
mcl_wrap.set_dimensions(Vec3(0.03, 0.018, 0.012));
mcl_wrap.set_translation(Vec3(0.001, -0.0346, -0.0187));
mcl_wrap.set_xyz_body_rotation(Vec3(-0.105, -0.035, 0.182));
mcl_wrap.set_quadrant('all');
tibia.addWrapObject(mcl_wrap);

MCL1_pt1 = Vec3(0.0115, 0.0075, -0.0405);
MCL1_pt2 = Vec3(0.0149, -0.0618, -0.0153);
MCL1 = Blankevoort1991Ligament('MCL1',femur,MCL1_pt1,tibia,MCL1_pt2);
model.addForce(MCL1);
MCL1.get_GeometryPath().addPathWrap(mcl_wrap);

MCL2_pt1 = Vec3(0.0052, 0.0083, -0.0406);
MCL2_pt2 = Vec3(0.0076, -0.0611, -0.0191);
MCL2 = Blankevoort1991Ligament('MCL2',femur,MCL2_pt1,tibia,MCL2_pt2);
model.addForce(MCL2);
MCL2.get_GeometryPath().addPathWrap(mcl_wrap);

MCL3_pt1 = Vec3(-0.0004, 0.0072, -0.0388);
MCL3_pt2 = Vec3(0.0012, -0.0619, -0.0188);
MCL3 = Blankevoort1991Ligament('MCL3',femur,MCL3_pt1,tibia,MCL3_pt2);
model.addForce(MCL3);
MCL3.get_GeometryPath().addPathWrap(mcl_wrap);

state = model.initSystem();

MCL1.setSlackLengthFromReferenceStrain(-0.03,state);
MCL1.set_linear_stiffness(1000);
MCL2.setSlackLengthFromReferenceStrain(0.0,state);
MCL2.set_linear_stiffness(1000);
MCL3.setSlackLengthFromReferenceStrain(0.03,state);
MCL3.set_linear_stiffness(1000);

%Print model file
model.finalizeConnections();
model_file = './inputs/knee.osim';
model.print(model_file);

%% Perform Analysis with JointMechanicsTool
basename = 'passive';
passive_result_dir = './results/passive/joint-mechanics';

jnt_mech = JointMechanicsTool();
jnt_mech.set_model_file(model_file);
jnt_mech.set_input_transforms_file(...
   '../../../../models/knee_healthy/experimental_data/dynamic_mri/passive_body_transformations.sto');

% jnt_mech.set_input_transforms_file('C:\Users\csmith\github\jam-resources\matlab\example\VisualizeKinematics\results\passive\joint-mechanics\passive_frame_transforms_in_ground.sto');
% jnt_mech.set_input_states_file(...
%     '../../../models/knee_healthy/experimental_data/dynamic_mri/passive_body_transformations.sto');

jnt_mech.set_results_file_basename('passive');
jnt_mech.set_results_directory(passive_result_dir);
jnt_mech.set_start_time(-1);
jnt_mech.set_stop_time(-1);
jnt_mech.set_normalize_to_cycle(false);
jnt_mech.set_contacts(0,'all');
jnt_mech.set_ligaments(0,'all');
jnt_mech.set_muscles(0,'none');
% jnt_mech.set_attached_geometry_bodies(0,'/bodyset/femur_distal_r');
% jnt_mech.set_attached_geometry_bodies(1,'/bodyset/tibia_proximal_r');
% jnt_mech.set_attached_geometry_bodies(2,'/bodyset/patella_r');
jnt_mech.set_output_orientation_frame('ground');
jnt_mech.set_output_position_frame('ground');
jnt_mech.set_write_vtp_files(true);
jnt_mech.set_write_h5_file(true);
jnt_mech.set_h5_kinematics_data(true);
jnt_mech.set_h5_states_data(true);
jnt_mech.set_use_visualizer(true);
% jnt_mech.set_write_transforms_file(true)
% jnt_mech.set_output_transforms_file_type('sto');
jnt_mech.print(['./inputs/' basename '_joint_mechanics_settings.xml']);

disp('Running JointMechanicsTool...');
% jnt_mech.run();

% Inertia
basename = 'inertia';
inertia_result_dir = './results/inertia/joint-mechanics';

jnt_mech.set_input_transforms_file(...
   '../../../../models/knee_healthy/experimental_data/dynamic_mri/inertia_body_transformations.sto');
jnt_mech.set_results_file_basename(basename);
jnt_mech.set_results_directory(inertia_result_dir);
jnt_mech.set_write_vtp_files(true);
jnt_mech.set_write_h5_file(true);
jnt_mech.set_h5_kinematics_data(true);
jnt_mech.set_h5_states_data(true);
jnt_mech.set_use_visualizer(useVisualizer);
jnt_mech.print(['./inputs/' basename '_joint_mechanics_settings.xml']);

disp('Running JointMechanicsTool...');
% jnt_mech.run();

% Elastic
basename = 'elastic';
elastic_result_dir = './results/elastic/joint-mechanics';

jnt_mech.set_input_transforms_file(...
   '../../../../models/knee_healthy/experimental_data/dynamic_mri/elastic_body_transformations.sto');
jnt_mech.set_results_file_basename(basename);
jnt_mech.set_results_directory(elastic_result_dir);
jnt_mech.set_write_vtp_files(true);
jnt_mech.set_write_h5_file(true);
jnt_mech.set_h5_kinematics_data(true);
jnt_mech.set_h5_states_data(true);
jnt_mech.set_use_visualizer(useVisualizer);
jnt_mech.print(['./inputs/' basename '_joint_mechanics_settings.xml']);

disp('Running JointMechanicsTool...');
% jnt_mech.run();

%% Plot Results
close all;

modelName='model';
passive_h5 = jam_analysis(modelName,{[passive_result_dir '/passive.h5']});
modelName='knee';
inertia_h5 = jam_analysis(modelName,{[inertia_result_dir '/inertia.h5']});
elastic_h5 = jam_analysis(modelName,{[elastic_result_dir '/elastic.h5']});

line_width = 2;

% TF Kinematics
figure('name','TF Kinematics');
coords = {'knee_flex_r','knee_add_r','knee_rot_r','knee_ant_r','knee_sup_r','knee_lat_r'};
for i = 1:length(coords)
    subplot(2,3,i);hold on;
    plot(passive_h5.coordinateset.(coords{i}).value,'LineWidth',line_width);
    plot(inertia_h5.coordinateset.(coords{i}).value,'LineWidth',line_width);
    plot(elastic_h5.coordinateset.(coords{i}).value,'LineWidth',line_width);
    xlabel('Frame')
    if(i<4)
        ylabel('Angle [^o]')
    else
        ylabel('Translation [m]')
    end
    title(strrep(coords{i},'_',' '))
end    
legend('passive','elastic','inertia')

% PF Kinematics
figure('name','PF Kinematics');
coords = {'pf_flex_r','pf_rot_r','pf_tilt_r','pf_ant_r','pf_sup_r','pf_lat_r'};
for i = 1:length(coords)
    subplot(2,3,i);hold on;
    plot(passive_h5.coordinateset.(coords{i}).value,'LineWidth',line_width);
    plot(inertia_h5.coordinateset.(coords{i}).value,'LineWidth',line_width);
    plot(elastic_h5.coordinateset.(coords{i}).value,'LineWidth',line_width);
    xlabel('Frame')
    if(i<4)
        ylabel('Angle [^o]')
    else
        ylabel('Translation [m]')
    end
    title(strrep(coords{i},'_',' '))
end    
legend('passive','elastic','inertia')

% Center of Proximity
figure('name','Center of Proximity time');
coords = {'knee_flex_r','knee_add_r','knee_rot_r','knee_ant_r','knee_sup_r','knee_lat_r'};

    
    
passive_COP_med = ...
    passive_h5.forceset.Smith2018ArticularContactForce.tf_contact.tibia_cartilage.region(5).regional_center_of_proximity;
passive_COP_lat = ...
    passive_h5.forceset.Smith2018ArticularContactForce.tf_contact.tibia_cartilage.region(6).regional_center_of_proximity;

inertia_COP_med = ...
    inertia_h5.forceset.Smith2018ArticularContactForce.tf_contact.tibia_cartilage.region(5).regional_center_of_proximity;
inertia_COP_lat = ...
    inertia_h5.forceset.Smith2018ArticularContactForce.tf_contact.tibia_cartilage.region(6).regional_center_of_proximity;

elastic_COP_med = ...
    elastic_h5.forceset.Smith2018ArticularContactForce.tf_contact.tibia_cartilage.region(5).regional_center_of_proximity;
elastic_COP_lat = ...
    elastic_h5.forceset.Smith2018ArticularContactForce.tf_contact.tibia_cartilage.region(6).regional_center_of_proximity;

passive_COP_med(passive_COP_med==-1)=nan;
passive_COP_lat(passive_COP_lat==-1)=nan;
inertia_COP_med(inertia_COP_med==-1)=nan;
inertia_COP_lat(inertia_COP_lat==-1)=nan;
elastic_COP_med(elastic_COP_med==-1)=nan;
elastic_COP_lat(elastic_COP_lat==-1)=nan;

index = reshape(1:6, 2,3).';
comp ={'anterior','superior','lateral'};

for i = 1:3
    subplot(3,2,index(i));hold on;
    plot(passive_COP_med(:,i),'LineWidth',line_width);
    plot(inertia_COP_med(:,i),'LineWidth',line_width);
    plot(elastic_COP_med(:,i),'LineWidth',line_width);
end 
title('Medial Tibia Plateau COP')
xlabel('Frames')
ylabel([comp{i} ' [m]'])

for i = 1:3
    subplot(3,2,index(i+3));hold on;
    plot(passive_COP_lat(:,i),'LineWidth',line_width);
    plot(inertia_COP_lat(:,i),'LineWidth',line_width);
    plot(elastic_COP_lat(:,i),'LineWidth',line_width);
end 
legend('passive','elastic','inertia')


figure('name','Center of Proximity');
subplot(1,2,1);hold on;
% plot(passive_COP_med(:,2),passive_COP_med(:,1),'LineWidth',line_width);
% plot(inertia_COP_med(:,2),inertia_COP_med(:,1),'LineWidth',line_width);
plot(elastic_COP_med(:,3),elastic_COP_med(:,1),'LineWidth',line_width);
title('Medial Tibia COP')
xlabel('Lateral [m]')
ylabel('Anterior [m]')

subplot(1,2,2);hold on;
% plot(passive_COP_lat(:,2),passive_COP_lat(:,1),'LineWidth',line_width);
% plot(inertia_COP_lat(:,2),inertia_COP_lat(:,1),'LineWidth',line_width);
plot(elastic_COP_lat(:,3),elastic_COP_lat(:,1),'LineWidth',line_width);
title('Lateral Tibia COP')
xlabel('Lateral [m]')
ylabel('Anterior [m]')
legend('passive','elastic','inertia')

% Ligament Strains
% for 1:3 
%     subplot(3,1,i);hold on;
%     plot(,'LineWidth',line_width)
% end