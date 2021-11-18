%% Generate OpenSim Input Files
%==========================================================================
%
% Treadmill trials DM_ngait_tmf_slow1_new.trc and DM_ngait_tmf_slow2_new.trc
% are a different .trc format. Manually delete extra tabs at end of header
% lines, then open in mokka and save as .trc to get them into a .trc format
% that opensim can read in this script. 
%
%==========================================================================
import org.opensim.modeling.*
clear;
data_dir = 'P:\Projects\LMB_grand_challenge_database\read_only\DataforSixthCompetition\'; 
overground_gait_video_dir = [data_dir 'Synchronized Motion Data\Overground Gait Trials\Video Motion Data'];
treadmill_gait_video_dir = [data_dir 'Synchronized Motion Data\Treadmill Gait Trials\Video Motion Data'];
out_folder ='./motion_analysis';


for d = 2:2%nDir
    if (d==1) 
        dir_info = dir(overground_gait_video_dir);        
    elseif (d==2)
        dir_info = dir(treadmill_gait_video_dir);
    end
    
    dir_info = dir_info(~ismember({dir_info.name},{'.','..'}));
    
    nFiles = size(dir_info,1);
    
    for i = 1:nFiles
        clear new_data;
        split_name = strsplit(dir_info(i).name,'_');
        file = [dir_info(i).folder '\' dir_info(i).name];

        
        if (contains(split_name{end},'grf') || ...
            (contains(split_name{end-1},'grf') && contains(split_name{end},'filt')))
            % Create GRF sto file
            fid = fopen(file);
            labels = strsplit(fgetl(fid), ',');
            fclose(fid);
            data = csvread(file,2);

            labels{ismember(labels,'time(sec)')}='time';

            nFP=0;
            FP_id=cell(0);
            nFrame = size(data,1);

            for j=1:length(labels)
                if(contains(labels{j},'time'))
                    new_data(:,j) = data(:,j);
                    continue;
                end
                id_num = labels{j}(end);
                if(contains(labels{j},'Fx'))
                    labels{j} = ['ground_force_' id_num '_vx'];
                    new_data(:,j) = -data(:,j);
                elseif(contains(labels{j},'Fy'))
                    labels{j} = ['ground_force_' id_num '_vy'];
                    new_data(:,j) = data(:,j+1);
                elseif(contains(labels{j},'Fz'))
                    labels{j} = ['ground_force_' id_num '_vz'];
                    new_data(:,j) = data(:,j-1);
                elseif(contains(labels{j},'COPx'))
                    labels{j} = ['ground_force_' id_num '_px'];
                    new_data(:,j) = -data(:,j)/1000;
                elseif(contains(labels{j},'COPy'))
                    labels{j} = ['ground_force_' id_num '_py'];
                    new_data(:,j) = data(:,j+1)/1000;
                elseif(contains(labels{j},'COPz'))
                    labels{j} = ['ground_force_' id_num '_pz'];
                    new_data(:,j) = data(:,j-1)/1000;
                    
                    % CAREFUL CHECK THIS WORKS FOR OVERGROUND
%                 elseif(contains(labels{j},'Tx'))
%                     labels{j} = ['ground_torque_' id_num '_x'];
%                     new_data(:,j) = -data(:,j);
%                 elseif(contains(labels{j},'Ty'))
%                     labels{j} = ['ground_torque_' id_num '_y'];
%                     new_data(:,j) = data(:,j+1);
%                 elseif(contains(labels{j},'Tz'))
%                     labels{j} = ['ground_torque_' id_num '_z'];
%                     new_data(:,j) = data(:,j-1);
                elseif(contains(labels{j},'Tz'))
                    labels{j} = ['ground_torque_' id_num '_y'];
                    new_data(:,j) = data(:,j)/1000;
                end

                if(~any(contains(FP_id,id_num)))
                    nFP = nFP+1;
                    FP_id{nFP} = id_num;
                end
            end



            for j=1:length(labels)
                data_struct.(labels{j}) = new_data(:,j);
            end

            for j=1:nFP
                data_struct.(['ground_torque_' FP_id{j} '_x']) = zeros(nFrame,1);
                data_struct.(['ground_torque_' FP_id{j} '_z']) = zeros(nFrame,1);
            end

            data_table = osimTableFromStruct(data_struct);
            data_table.addTableMetaDataString('nRows',num2str(data_table.getNumRows()));
            data_table.addTableMetaDataString('nColumns',num2str(data_table.getNumColumns()+1));
            STOFileAdapter.write(data_table,[out_folder '/' dir_info(i).name(1:end-3) 'sto']);

            %Create External Loads Files
            %ExternalLoads = 
        end

        if contains(split_name{end},'trc')
            if contains(split_name{end},'org'); continue; end
            trc = TRCFileAdapter();
            trc_table = TimeSeriesTableVec3(file);
            
            trc_struct = osimTableToStruct(trc_table);
            
            fields = fieldnames(trc_struct);
            for j = 1:length(fields)
                if contains(fields{j},'time')
                    continue;
                end
               p = trc_struct.(fields{j});
               
               p = fillmissing(p,'pchip','EndValues','none');
               trc_struct.(fields{j}) = [-p(:,1) p(:,3) p(:,2)];
            end
            
            % Account for miss sync in .trc file
            if contains(dir_info(i).name,'tmf_slow1')
                trc_struct.time = trc_struct.time + 0.439997;
            elseif contains(dir_info(i).name,'tmf_slow2')
                trc_struct.time = trc_struct.time + 0.289997;
            end
            new_table = osimTableFromStruct(trc_struct);
            new_table.addTableMetaDataString('DataRate', num2str(120));
            new_table.addTableMetaDataString('Units', 'mm');
            trc.write(new_table,char([out_folder '/' dir_info(i).name]));
            
        end
    end
end