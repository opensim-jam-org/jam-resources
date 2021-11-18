%% Build DM Model
import org.opensim.modeling.*

model = Model("C:\Users\csmith\github\jam-resources\models\knee_tka\grand_challenge\DM\DM_old.osim");

model.setName('DM');

secondary_coord = {...
    'knee_add_r',...
    'knee_rot_r',...
    'knee_tx_r',...
    'knee_ty_r',...
    'knee_tz_r',...
    'pf_flex_r',...
    'pf_rot_r',...
    'pf_tilt_r',...
    'pf_tx_r',...
    'pf_ty_r',...
    'pf_tz_r',...
    };
secondary_viscosity = [5,5,100,100,100,2,2,2,10,10,10];

%% Add Secondary Damping
for i = 1:length(secondary_coord)
    damp_frc = SpringGeneralizedForce(secondary_coord{i});
    damp_frc.setName(secondary_coord{i});
    damp_frc.setViscosity(secondary_viscosity(i));
    
    model.addForce(damp_frc);
end
state = model.initSystem();
%% Set Ligament Slack Length
for i=0:model.getForceSet().getSize()-1
    frc = model.getForceSet().get(i);
    if(strcmp(frc.getConcreteClassName(),'Blankevoort1991Ligament'))
        lig = Blankevoort1991Ligament.safeDownCast(frc);
        
        if(contains(char(lig.getName()),'pCAP'))
            lig.setSlackLengthFromReferenceStrain(0.08,state);
        elseif(contains(char(lig.getName()),'ITB'))
            lig.setSlackLengthFromReferenceStrain(0.02,state);
        elseif(contains(char(lig.getName()),'MCLd'))
            lig.setSlackLengthFromReferenceStrain(0.03,state);
        elseif(contains(char(lig.getName()),'MCLs'))
            lig.setSlackLengthFromReferenceStrain(0.03,state);
        elseif(contains(char(lig.getName()),'PCLal'))
            lig.setSlackLengthFromReferenceStrain(0.01,state);
        elseif(contains(char(lig.getName()),'PCLpm'))
            lig.setSlackLengthFromReferenceStrain(-0.06,state);
        elseif(contains(char(lig.getName()),'LCL'))
            lig.setSlackLengthFromReferenceStrain(0.06,state);
        elseif(contains(char(lig.getName()),'PT'))
            lig.setSlackLengthFromReferenceStrain(0.02,state);
        elseif(contains(char(lig.getName()),'lPFL'))
            lig.setSlackLengthFromReferenceStrain(0.01,state);
        elseif(contains(char(lig.getName()),'mPFL'))
            lig.setSlackLengthFromReferenceStrain(0.-0.05,state);
        elseif(contains(char(lig.getName()),'pMC'))
            lig.setSlackLengthFromReferenceStrain(0.05,state);
        elseif(contains(char(lig.getName()),'PFL'))
            lig.setSlackLengthFromReferenceStrain(-0.01,state);    
        end
    end
end
%% Add Contact Geometry
femur_implant_mesh_file = 'femur_component_surface_gc_m1.stl';
patella_implant_mesh_file = 'patella_button_surface_v2.stl';
tibia_implant_mesh_file = 'tibia_insert_surface_gc_p1.stl';

femur_cnt_mesh = Smith2018ContactMesh('femur_implant',...
                    femur_implant_mesh_file,model.getBodySet.get('femur_r'));
femur_cnt_mesh.set_thickness(0.003);
femur_cnt_mesh.set_elastic_modulus(463e6);
femur_cnt_mesh.set_poissons_ratio(0.46);
femur_cnt_mesh.set_scale_factors(Vec3(1.055));
model.addContactGeometry(femur_cnt_mesh);
                

tibia_cnt_mesh = Smith2018ContactMesh('tibia_implant',...
                    tibia_implant_mesh_file,model.getBodySet.get('tibia_proximal_r'));

tibia_cnt_mesh.set_thickness(0.003);
tibia_cnt_mesh.set_elastic_modulus(463e6);
tibia_cnt_mesh.set_poissons_ratio(0.46);   
tibia_cnt_mesh.set_scale_factors(Vec3(1.055));
model.addContactGeometry(tibia_cnt_mesh);


patella_cnt_mesh = Smith2018ContactMesh('patella_implant',...
                    patella_implant_mesh_file,model.getBodySet.get('patella_r'));
patella_cnt_mesh.set_thickness(0.003);
patella_cnt_mesh.set_elastic_modulus(463e6);
patella_cnt_mesh.set_poissons_ratio(0.46);  
patella_cnt_mesh.set_scale_factors(Vec3(1.055));
                model.addContactGeometry(patella_cnt_mesh);

tf_contact = Smith2018ArticularContactForce('knee_contact',...
                femur_cnt_mesh,tibia_cnt_mesh);
tf_contact.set_elastic_foundation_formulation('linear');
tf_contact.set_use_lumped_contact_model(true);
model.addForce(tf_contact);

pf_contact = Smith2018ArticularContactForce('pf_contact',...
                femur_cnt_mesh,patella_cnt_mesh);
model.addForce(pf_contact);
%% Change Muscles to Millard
model.initSystem();

new_frc_set = ForceSet();%model.getForceSet());

for i=0:model.getForceSet().getSize()-1
    frc = model.getForceSet().get(i);
    
    if(strcmp(frc.getConcreteClassName(),'Thelen2003Muscle'))
        msl = Thelen2003Muscle.safeDownCast(frc);
        
        new_msl = Millard2012EquilibriumMuscle(...
            msl.getName(),msl.getMaxIsometricForce(),...
            msl.getOptimalFiberLength(), msl.getTendonSlackLength(),...
            msl.getPennationAngleAtOptimalFiberLength());
        
        new_msl.set_GeometryPath(msl.get_GeometryPath());
        
        new_frc_set.cloneAndAppend(new_msl);

    else
           new_frc_set.cloneAndAppend(frc);
    end
end
model.getForceSet().clearAndDestroy();
model.initSystem();
for i=0:new_frc_set.getSize()-1
    model.getForceSet().cloneAndAppend(new_frc_set.get(i));
end

model.print("C:\Users\csmith\github\jam-resources\models\knee_tka\grand_challenge\DM\DM.osim");

%% Create Reserve Actuator Files
state = model.initSystem();

%% Get the number of coordinates  and a handle to the coordainte set
coordSet = model.getCoordinateSet();
nCoord = coordSet.getSize();

%% Instantiate some empty vec3's for later.
massCenter = Vec3();
axisValues = Vec3();

%% Instantiate an empty Force set
forceSet = ForceSet();

%% Set the optimal force
optimalForce = 1;

%% Start going through the coordinates, creating an actuator for each
for iCoord = 0 : nCoord - 1

    % get a reference to the current coordinate
    coordinate = coordSet.get(iCoord);
    % If the coodinate is constrained (locked or prescribed), don't
    % add an actuator
    if coordinate.isConstrained(state)
        continue
    end

    % get the joint, parent and child names for the coordiante
    joint = coordinate.getJoint();
    parentName = joint.getParentFrame().getName();
    childName = joint.getChildFrame().getName();

    % If the coordinates parent body is connected to ground, we need to
    % add residual actuators (torque or point).
    if strcmp(parentName, model.getGround.getName() )

        % Custom and Free Joints have three translational and three
        % rotational coordinates.
        if strcmp(joint.getConcreteClassName(), 'CustomJoint') || strcmp(joint.getConcreteClassName(), 'FreeJoint')
               % get the coordainte motion type
               motion = char(coordinate.getMotionType());
               % to get the axis value for the coordinate, we need to drill
               % down into the coordinate transform axis
               eval(['concreteJoint = ' char(joint.getConcreteClassName()) '.safeDownCast(joint);'])
               sptr = concreteJoint.getSpatialTransform();
               for ip = 0 : 5
                  if strcmp(char(sptr.getCoordinateNames().get(ip)), char(coordinate.getName))
                        sptr.getTransformAxis(ip).getAxis(axisValues);
                        break
                  end
               end


               % make a torque actuator if a rotational coordinate
               if strcmp(motion, 'Rotational')
                   newActuator = TorqueActuator(joint.getParentFrame(),...
                                         joint.getParentFrame(),...
                                         axisValues);

               % make a point actuator if a translational coordinate.
             elseif strcmp(motion, 'translational')
                    % make a new Point actuator
                    newActuator = PointActuator();
                    % set the body
                    newActuator.set_body(char(joint.getChildFrame().getName()))
                    % set point that forces acts at
                    newActuator.set_point(massCenter)
                    % the point is expressed in the local
                    newActuator.set_point_is_global(false)
                    % set the direction that actuator will act in
                    newActuator.set_direction(axisValues)
                    % the force is expressed in the global
                    newActuator.set_force_is_global(true)
              else % something else that we don't support right now
                    newActuator = CoordinateActuator();
              end
        else % if the joint type is not free or custom, just add coordinate actuators
                % make a new coordinate actuator for that coordinate
                newActuator = CoordinateActuator();
        end

    else % the coordinate is not connected to ground, and can just be a
         % coordinate actuator.
         newActuator = CoordinateActuator( char(coordinate.getName) );
    end

    % set the optimal force for that coordinate
    newActuator.setOptimalForce(optimalForce);
    % set the actuator name
    newActuator.setName( coordinate.getName() );
    % set min and max controls
    newActuator.setMaxControl(Inf)
    newActuator.setMinControl(-Inf)

    % append the new actuator onto the empty force set
    forceSet.cloneAndAppend(newActuator);
end

forceSet.print('DM_reserve_actuators.xml');