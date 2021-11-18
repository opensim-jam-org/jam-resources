%% Example Build Knee Model
%===============================================================================

%Setup Environment and Folders
import org.opensim.modeling.*
Logger.setLevelString('Info');


useVisualizer = true;


if(exist('./results/Geometry','dir')~=7)
    mkdir('./results/Geometry')
end
if(exist('./results','dir')~=7)
    mkdir('./results')
end


%% Build Model

% Copy Geometries to example dir
geometry_dir = '../../../models/knee_healthy/lenhart2015/Geometry/';
local_geometry_dir = './results/Geometry/';

femur_bone_mesh_file = 'lenhart2015-R-femur-bone.stl';
tibia_bone_mesh_file = 'lenhart2015-R-tibia-bone.stl';

femur_cartilage_mesh_file = 'lenhart2015-R-femur-cartilage.stl';
tibia_cartilage_mesh_file = 'lenhart2015-R-tibia-cartilage.stl';

copyfile([geometry_dir femur_bone_mesh_file],[local_geometry_dir femur_bone_mesh_file])
copyfile([geometry_dir tibia_bone_mesh_file],[local_geometry_dir tibia_bone_mesh_file])
copyfile([geometry_dir femur_cartilage_mesh_file],[local_geometry_dir femur_cartilage_mesh_file])
copyfile([geometry_dir tibia_cartilage_mesh_file],[local_geometry_dir tibia_cartilage_mesh_file])


ModelVisualizer.addDirToGeometrySearchPaths(local_geometry_dir);

% Create Model
model = Model();
modelName = 'knee';
model.setName(modelName);

% Add Bodies
% Mass properties are arbitrary as simple example
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

%Add Joints

ground_transform = SpatialTransform();

ground_transform.get_rotation1().set_coordinates(0,'ground_femur_rx_r');
ground_transform.get_rotation1().set_axis(Vec3(1,0,1));
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

femur_cnt_mesh = Smith2018ContactMesh('femur_cartilage',...
                    femur_cartilage_mesh_file,femur);
femur_cnt_mesh.set_elastic_modulus(5e6);
femur_cnt_mesh.set_poissons_ratio(0.45);
femur_cnt_mesh.set_thickness(0.03);                
model.addContactGeometry(femur_cnt_mesh);
                

tibia_cnt_mesh = Smith2018ContactMesh('tibia_cartilage',...
                    tibia_cartilage_mesh_file,tibia);
tibia_cnt_mesh.set_elastic_modulus(5e6);
tibia_cnt_mesh.set_poissons_ratio(0.45);
tibia_cnt_mesh.set_thickness(0.03);
model.addContactGeometry(tibia_cnt_mesh);


tf_contact = Smith2018ArticularContactForce('tf_contact',...
                femur_cnt_mesh,tibia_cnt_mesh);
model.addForce(tf_contact);



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

MCL1.set_slack_length(0.03);
MCL1.set_linear_stiffness(1000);
MCL2.setSlackLengthFromReferenceStrain(0.03,state);
MCL2.set_linear_stiffness(1000);
MCL3.setSlackLengthFromReferenceForce(20,state);
MCL3.set_linear_stiffness(1000);

%Print model file
model_file = './results/knee.osim';
model.print(model_file);
