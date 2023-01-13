import os
import opensim as osim

models_path = os.path.abspath('../../../../models')
model_file = os.path.join(models_path, 'knee_healthy', 'lenhart2015', 'lenhart2015.osim')

comak_result_dir = "/root/jam-resources/python/examples/Walking/Walking/results/comak";
results_basename = "walking";
save_xml_path = "/root/jam-resources/python/examples/Walking/Walking/inputs/comak_settings.xml";

results_basename = 'walking'
comak_result_dir = './results/comak'
jnt_mech_result_dir = './results/joint-mechanics'

## Perform Joint Mechanics Analysis
jnt_mech = osim.JointMechanicsTool();
jnt_mech.set_model_file(model_file);
jnt_mech.set_input_states_file(os.path.join(comak_result_dir, results_basename + '_states.sto'));
jnt_mech.set_use_muscle_physiology(False);
jnt_mech.set_results_file_basename(results_basename);
jnt_mech.set_results_directory(jnt_mech_result_dir);
jnt_mech.set_start_time(1.16);
jnt_mech.set_stop_time(-1);
jnt_mech.set_resample_step_size(-1);
jnt_mech.set_normalize_to_cycle(True);
jnt_mech.set_lowpass_filter_frequency(-1);
jnt_mech.set_print_processed_kinematics(False);
jnt_mech.set_contacts(0,'all');
jnt_mech.set_contact_outputs(0,'all');
jnt_mech.set_contact_mesh_properties(0,'none');
jnt_mech.set_ligaments(0,'all');
jnt_mech.set_ligament_outputs(0,'all');
jnt_mech.set_muscles(0,'all');
jnt_mech.set_muscle_outputs(0,'all');

jnt_mech.set_attached_geometry_bodies(0,'all');

jnt_mech.set_output_orientation_frame('ground');
jnt_mech.set_output_position_frame('ground');
jnt_mech.set_write_vtp_files(True);
jnt_mech.set_vtp_file_format('binary');
jnt_mech.set_write_h5_file(False);
jnt_mech.set_h5_kinematics_data(True);
jnt_mech.set_h5_states_data(True);
jnt_mech.set_write_transforms_file(False);
jnt_mech.set_output_transforms_file_type('sto');
jnt_mech.set_use_visualizer(False);
jnt_mech.setDebugLevel(0);

analysis_set = osim.AnalysisSet();

frc_reporter = osim.ForceReporter();
frc_reporter.setName('ForceReporter');

analysis_set.cloneAndAppend(frc_reporter);
jnt_mech.set_AnalysisSet(analysis_set);
jnt_mech.printToXML('./inputs/joint_mechanics_settings.xml');

print('Running JointMechanicsTool...');
jnt_mech.run();
print('Finished JointMechanicsTool!')