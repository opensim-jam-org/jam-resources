close all
results_dir = './data';
model = 'smith2019';

dir_info = dir(results_dir);
dir_info=dir_info(~ismember({dir_info.name},{'.','..'}));

nResultFiles = size(dir_info,1);

file_list = cell(nResultFiles,1);
for i = 1:nResultFiles
    subresult_dir = [results_dir '/' dir_info(i).name];
    file = [subresult_dir '/walking.h5'];        
    file_list{i}= file;                
end

jam = jam_analysis(model,file_list);


%% Plot Secondary Coordinates
sec_coords = {...
    'knee_flex_r','knee_add_r','knee_rot_r',...
    'knee_tx_r','knee_ty_r','knee_tz_r',...
    };
figure('name','Knee Kinematics')
for i = 1:length(sec_coords)    
    subplot(2,3,i);hold on
    for n=1:jam.nFiles
        plot(0:100,jam.coordinateset.(sec_coords{i}).value(:,n));
    end
    xlim([0,100])
end

%% Plot Ligament Forces
ligament_names = {...
    'MCLd','MCLs'...
    'LCL',...
    'ACLam','ACLpl',...
    'PCL','PCL',...
    'PT',...
    'ITB',...
    'pCAP'...
    };

for i = 1:length(ligament_names)
    figure('name',ligament_names{i})
    
    fiber_names = fieldnames(jam.forceset.Blankevoort1991Ligament);
    
    fibers = fiber_names(contains(fiber_names,ligament_names{i}));
    for n=1:jam.nFiles
        data = 0;
        for k = 1:length(fibers)
            data = data + jam.forceset.Blankevoort1991Ligament.(fibers{k}).total_force;
            
        end
        plot(data)
        title([ligament_names{i} ' force'])
    end
end

%% Plot Muscle Forces
muscles = {'vasmed_r','recfem_r'};

for i = 1:length(muscles)
    figure('name', muscles{i});
    hold on;
    jam.plot_muscle_output(gca,muscles{i},'activation');

end

%% Plot Contact Forces
figure('name','TF Cartilage-Cartilage Contact Force')
comp = {'X','Y','Z'};
for i = 1:3
    subplot(3,3,i);hold on
    for n=1:jam.nFiles
        plot(jam.forceset.Smith2018ArticularContactForce.tf_contact.tibia_cartilage.total_contact_force(:,i,n));
    end
    title(['Total TF Cart-Cart Contact ' comp{i}])
    
    subplot(3,3,3+i);hold on
    for n=1:jam.nFiles
        plot(jam.forceset.Smith2018ArticularContactForce.tf_contact.tibia_cartilage.region(5).regional_contact_force(:,i,n));
    end
    title(['Medial TF Cart-Cart Contact ' comp{i}])
    
    subplot(3,3,6+i);hold on
    for n=1:jam.nFiles
        plot(jam.forceset.Smith2018ArticularContactForce.tf_contact.tibia_cartilage.region(6).regional_contact_force(:,i,n));
    end
    title(['Lateral TF Cart-Cart Contact ' comp{i}])
end

%% Plot COMAK convergence
figure('name','Comak Convergence')
subplot(2,1,1);hold on
for n=1:jam.nFiles
    plot(jam.comak.max_udot_error(:,n));
end
title('max udot error')

subplot(2,1,2);hold on
for n=1:jam.nFiles
    plot(jam.comak.iterations(:,n));
end
title('iterations')