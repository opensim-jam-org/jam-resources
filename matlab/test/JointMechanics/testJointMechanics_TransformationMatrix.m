%% test JointMechanics transformation matrix input and output
%==========================================================================
import org.opensim.modeling.*

useVisualizer = true;

%model_file = '../../../models/knee_healthy/lenhart2015/lenhart2015.osim';
model_file = "C:\Users\csmith\github\jam-resources\matlab\example\VisualizeMeasuredKinematics\VisualizeDynamicMRI\inputs\knee.osim";

% Copy Geometries to example dir
geometry_dir = '../../../models/knee_healthy/lenhart2015/Geometry/';
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

femur_cnt_mesh = Smith2018ContactMesh('femur_cartilage',...
                    femur_cartilage_mesh_file,femur);
model.addContactGeometry(femur_cnt_mesh);

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

%% TF Joint
tibia  = Body();
tibia.setName('tibia_proximal_r');
tibia.setMass(1.0);                                      
tibia.setMassCenter(massCenter);                        
tibia.setInertia(inertia);
tibia.attachGeometry(Mesh(tibia_bone_mesh_file));
model.addBody(tibia);

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

knee_tx = model.getCoordinateSet().get('knee_tx_r');
knee_tx.set_range(0,-0.2);
knee_tx.set_range(1,0.2);
knee_tx.set_clamped(true);

knee_ty = model.getCoordinateSet().get('knee_ty_r');
knee_ty.set_range(0,-0.2);
knee_ty.set_range(1,0.2);
knee_ty.set_clamped(true);

knee_tz = model.getCoordinateSet().get('knee_tz_r');
knee_tz.set_range(0,-0.2);
knee_tz.set_range(1,0.2);
knee_tz.set_clamped(true);




                

tibia_cnt_mesh = Smith2018ContactMesh('tibia_cartilage',...
                    tibia_cartilage_mesh_file,tibia);
model.addContactGeometry(tibia_cnt_mesh);

tf_contact = Smith2018ArticularContactForce('tf_contact',...
                femur_cnt_mesh,tibia_cnt_mesh);
model.addForce(tf_contact);

%% PF Joint

patella  = Body();
patella.setName('patella_r');
patella.setMass(1.0);                                      
patella.setMassCenter(massCenter);                        
patella.setInertia(inertia); 
patella.attachGeometry(Mesh(patella_bone_mesh_file));
model.addBody(patella);

pf_transform = SpatialTransform();

pf_transform.get_rotation1().set_coordinates(0,'pf_flex_r');
pf_transform.get_rotation1().set_axis(Vec3(0,0,1));
pf_transform.get_rotation1().set_function(LinearFunction(-1,0));

pf_transform.get_rotation2().set_coordinates(0,'pf_rot_r');
pf_transform.get_rotation2().set_axis(Vec3(1,0,0));
pf_transform.get_rotation2().set_function(LinearFunction());

pf_transform.get_rotation3().set_coordinates(0,'pf_tilt_r');
pf_transform.get_rotation3().set_axis(Vec3(0,1,0));
pf_transform.get_rotation3().set_function(LinearFunction());

pf_transform.get_translation1().set_coordinates(0,'pf_tx_r');
pf_transform.get_translation1().set_axis(Vec3(1,0,0));
pf_transform.get_translation1().set_function(LinearFunction());

pf_transform.get_translation2().set_coordinates(0,'pf_ty_r');
pf_transform.get_translation2().set_axis(Vec3(0,1,0));
pf_transform.get_translation2().set_function(LinearFunction());

pf_transform.get_translation3().set_coordinates(0,'pf_tz_r');
pf_transform.get_translation3().set_axis(Vec3(0,0,1));
pf_transform.get_translation3().set_function(LinearFunction());

pf_r = CustomJoint('pf_r',...
    femur,Vec3(0.0),Vec3(0.0),patella,Vec3(0.0),Vec3(0.0),...
    pf_transform);

model.addJoint(pf_r);

pf_flex = model.getCoordinateSet().get('pf_flex_r');
pf_flex.set_range(1,-pi/4);
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

patella_cnt_mesh = Smith2018ContactMesh('patella_cartilage',...
                    patella_cartilage_mesh_file,patella);
model.addContactGeometry(patella_cnt_mesh);
pf_contact = Smith2018ArticularContactForce('pf_contact',...
                femur_cnt_mesh,patella_cnt_mesh);

model.addForce(pf_contact);

%Print model file
model.finalizeConnections();
model_file = './inputs/knee.osim';
model.print(model_file);

%%


jm1 = JointMechanicsTool();
jm1.set_model_file(model_file);
input_states_file = './inputs/passive_flexion_states.sto';
jm1.set_input_states_file(input_states_file);
jm1.set_results_directory('./results1');
jm1.set_results_file_basename('result1');
jm1.set_write_transforms_file(true);
jm1.set_output_transforms_file_type('sto')
jm1.set_write_h5_file(false);
jm1.set_write_vtp_files(false);
% jm1.run()

jm2 = JointMechanicsTool();
jm2.set_model_file(model_file);
jm2.set_input_transforms_file('./results1/result1_frame_transforms_in_ground.sto');
jm2.set_results_directory('./results2');
jm2.set_results_file_basename('result2');
jm2.set_write_transforms_file(false);
jm2.set_write_h5_file(true);
jm2.set_write_vtp_files(false);
jm2.set_use_visualizer(useVisualizer);
jm2.print("jm2_settings.xml");
% jm2.run()

%% 
input_states = osimTableToStruct(TimeSeriesTable(input_states_file));

jam = jam_analysis('./results2/result2.h5');
%%
% Plot TF Coordinates
tf_coords = {'knee_flex_r','knee_add_r','knee_rot_r','knee_tx_r','knee_ty_r','knee_tz_r'};
figure('name','TF Kinematics');

for i = 1:length(tf_coords)
    subplot(2,3,i);hold on;
    plot(input_states.(['a_jointset_knee_r_' tf_coords{i} '_value'])*180/pi);
    plot(jam.coordinateset.(tf_coords{i}).value)
    title(tf_coords{i})
end