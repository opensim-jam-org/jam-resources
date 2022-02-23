%% Test Joint Mechanics Output 
import org.opensim.modeling.*

useVisualizer = true;

model_file = '../../../models/knee_tka\grand_challenge\DM\DM.osim';
states_file = './tka_passive_flexion_states.sto';
results_basename = 'test';
jnt_mech_result_dir = './results';

jnt_mech = JointMechanicsTool();
jnt_mech.set_model_file(model_file);
jnt_mech.set_input_states_file(states_file);
jnt_mech.set_use_muscle_physiology(false);
jnt_mech.set_results_file_basename(results_basename);
jnt_mech.set_results_directory(jnt_mech_result_dir);
jnt_mech.set_start_time(0);
jnt_mech.set_stop_time(-1);
jnt_mech.set_normalize_to_cycle(false);
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
jmft.set_parent_frame('/contactgeometryset/femur_implant/mesh_frame');
jmft.set_child_frame('/contactgeometryset/tibia_implant/mesh_frame');
jmft.set_rotation_type('body');
jmft.set_rotation_sequence(0,'x');
jmft.set_rotation_sequence(1,'y');
jmft.set_rotation_sequence(2,'z');
jmft.set_output_coordinates(true);
jmft.set_output_transformation_matrix(false);

jmft_set = JointMechanicsFrameTransformSet();
jmft_set.cloneAndAppend(jmft);        
jnt_mech.set_JointMechanicsFrameTransformSet(jmft_set); 

jnt_mech.print(['joint_mechanics_settings.xml']);

disp('Running JointMechanicsTool...');
jnt_mech.run();