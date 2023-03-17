function [expt, Nexpt, dataTable] = GetVascExpt(dataTablePath ) % dataCol,  dataTable, row
    mainDir = 'D:\2photon\';    
    dataCol = struct('mouse',1, 'date',2, 'FOV',3, 'runs',4, 'CSD',5, 'superficial',6, 'penetrating',7, 'Done',8); %'color',6
    dataTable = readcell( dataTablePath );  
    dataTable(1,:) = [];
    dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);
    missingData = cellfun(@all, cellfun(@ismissing, dataTable, 'UniformOutput',false), 'UniformOutput',true );
    Nexpt = size(dataTable, 1);
    expt = repmat( struct('mouse','', 'date',[], 'fov',NaN, 'name','', 'dir','', 'runs',[], 'Nruns',NaN, 'csd',NaN, 'superficial',NaN, 'penetrating',NaN, 'done',NaN), Nexpt, 1); % 

    %xDone = find([dataTable{:,dataCol.done}] > 0);
    %x3D = intersect( xDone, find([dataTable{:,dataCol.volume}] == 1 ));
    %x2D = intersect( xDone, find( [dataTable{:,dataCol.volume}] == 0 ));
    %xPresent = 2; %x3D; %xAll; % x3Dcsd; 
    for row = 1:Nexpt
        expt(row).mouse = dataTable{row,dataCol.mouse};
        expt(row).date = dataTable{row,dataCol.date};  if isnumeric(expt(row).date), expt(row).date = num2str(expt(row).date); end
        if missingData(row,dataCol.FOV)
            expt(row).fov = 1;
        else
            expt(row).fov = dataTable{row,dataCol.FOV};
            if ischar(expt(row).fov), expt(row).fov = str2double(expt(row).fov); end
        end
        expt(row).dir = sprintf('%s%s\\%s_%s\\', mainDir, expt(row).mouse, expt(row).date, expt(row).mouse); %sprintf('%s%s\\%s_FOV%i_%s\\', mainDir, expt(row).mouse, expt(row).date, expt(row).fov, expt(row).mouse);% 
        expt(row).name = sprintf('%s_%s_FOV%i', expt(row).mouse, expt(row).date, expt(row).fov); %sprintf('%s_%s', expt(row).mouse, expt(row).date);
        %{
        if expt(row).fov == 1
            %D:\2photon\CGRP02\201110_CGRP02
            expt(row).dir = sprintf('%s%s\\%s_%s\\', mainDir, expt(row).mouse, expt(row).date, expt(row).mouse); %sprintf('%s%s\\%s_%s\\', mainDir, expt(x).mouse, expt(x).date, expt(x).mouse);
            expt(row).name = sprintf('%s_%s', expt(row).mouse, expt(row).date);
        else
            expt(row).dir = sprintf('%s%s\\%s_FOV%i_%s\\', mainDir, expt(row).mouse, expt(row).date, expt(row).fov, expt(row).mouse);
            expt(row).name = sprintf('%s_%s_FOV%i', expt(row).mouse, expt(row).date, expt(row).fov);
        end
        %}

        fprintf('\n\n%s', expt(row).name )
        if isnumeric(dataTable{row,dataCol.runs}) 
            expt(row).runs = dataTable{row,dataCol.runs}; 
        elseif ischar(dataTable{row,dataCol.runs})
            expt(row).runs = str2num( dataTable{row,dataCol.runs} );  %#ok<*ST2NM>
        else
            error('Runs format needs to be numeric or character');
        end
        expt(row).Nruns = numel(expt(row).runs); 
        
        % Get central planes for superficial and penetrating vessels
        expt(row).superficial = dataTable{row, dataCol.superficial};
        expt(row).penetrating = dataTable{row, dataCol.penetrating};
    end
end