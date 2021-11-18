%% Adjust Passive Fiber Force-Length
%==========================================================================
% Author: Colin Smith
%
% This script changes the passive fiber force-length curve to reduce
% passive muscle forces when the knee is in deep flexion. 
%
% These changes match those in the orignal SIMM model made by Darryl 
% Thelen.
%
%
%
%==========================================================================
import org.opensim.modeling.*

model_file = "./knee_healthy/lenhart2015/lenhart2015.osim";
% model_file = "./knee_healthy/smith2019/smith2019.osim";
% model_file = "./knee_healthy/smith2019/smith2019.osim";
msl_to_adjust = {'gaslat_r','gasmed_r','recfem_r','soleus_r','vasint_r','vaslat_r','vasmed_r'};

model = Model(model_file);

msl_set = model.getMuscles();

for m = 0:msl_set.getSize()-1
    msl = msl_set.get(m);
    
    if (~contains(msl_to_adjust,char(msl.getName())))
        continue;
    end
    
    if(strcmp(char(msl.getConcreteClassName()),"Millard2012EquilibriumMuscle"))
        msl = Millard2012EquilibriumMuscle.safeDownCast(msl);
        
        FL = msl.get_FiberForceLengthCurve();
        FL.set_strain_at_one_norm_force(2.0);
    end
end

model.print(model_file);
fprintf('Wrote updated model file: \n%s\n',model_file)