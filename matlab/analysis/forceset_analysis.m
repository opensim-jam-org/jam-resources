classdef forceset_analysis
    properties
        data
    end
    methods
        function obj = forceset_analysis(h5_file_list)
        %==================================================================
        %CLASS CONSTRUCTER
        %------------------------------------------------------------------
        %file : full path to .h5 or .sto file to read and analyze 
        %(include extension)
        %==================================================================
            
            for n=1:length(h5_file_list)
                h5_file = h5_file_list{n};

                forceset_path = ['/' model_name '/forceset'];
                info = h5info(h5_file,forceset_path);

                nComponents = length(info.Groups);

                for i=1:nComponents
                    label = strsplit(info.Groups(i).Name,'/');
                    comp = label{4};                

                    nForce = length(info.Groups(i).Groups);

                    for j=1:nForce
                        label = strsplit(info.Groups(i).Groups(j).Name,'/');
                        force = label{5}; 


                        nDataSet = length(info.Groups(i).Groups(j).Datasets);

                        for k=1:nDataSet
                            param = info.Groups(i).Groups(j).Datasets(k).Name;
                            data_set_name = [ info.Groups(i).Groups(j).Name '/' param];
                            obj.data.(comp).(force).(param)(:,n) = h5read(h5_file,data_set_name);
                        end
                    end
                end
            end
        end
    
    
        function h = plot_muscle_param(obj,muscle_names,param_name)
        %==================================================================
        % PLOT_MUSCLE_PARAM
        %------------------------------------------------------------------
        %
        %
        %==================================================================      

%             if strcmp(muscle_names{1},'all')
%                 muscle_names = obj.names;
%             end

            h = figure;

            nMsl = length(muscle_names);

            hold on;
            for i = 1:nMsl                   
                plot(obj.data.Muscle.(muscle_names{i}).(param_name),'LineWidth',2);                
            end
            %xlabel('time [s]')
            ylabel(param_name)
            legend(muscle_names)
        end
        
        function plot_muscle_param_groups(obj,model,param_name)
        %==================================================================
        % PLOT_MUSCLE_PARAM
        %------------------------------------------------------------------
        %
        %
        %==================================================================            
            for g = 1:model.getForceSet().getNumGroups
                group = model.getForceSet().getGroup(g-1);
                
                for m = 1:group.getMembers.getSize()
                    muscle_names{m} = char(group.getMembers.get(m-1));            
                end
                %figure('name',[char(group.getName()) ' : ' param_name])
                
                h = obj.plot_muscle_param(muscle_names,param_name);
                set(h,'name',[char(group.getName()) ' : ' param_name])
                title([char(group.getName()) ' : ' param_name])
            end
        end
    end
end
        
%         
%         
%         
% 
%             info = h5info(h5_file,'/Muscles');
% 
%             obj.nMuscles = length(info.Groups);
% 
%             obj.names = cell(obj.nMuscles,1);
% 
%             for i = 1:obj.nMuscles
%                 label = strsplit(info.Groups(i).Name,'/');
%                 obj.names{i} = label{3};
% 
%                 nDataSet = length(info.Groups(i).Datasets);
% 
%                 for j = 1:nDataSet 
%                     param = info.Groups(i).Datasets(j).Name;
%                     data_set_name = ['/Muscles/' obj.names{i} '/' param];
%                     obj.data.(obj.names{i}).(param) = h5read(h5_file,data_set_name);
%                 end
%             end           
%         end
%     end
% end