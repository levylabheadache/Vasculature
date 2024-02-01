clear; clc; close all;
dataDir = 'D:\2photon\Simone\Simone_Macrophages\'; % 'D:\2photon\Simone\'; %'C:\2photon';
dataSet = 'Vasculature'; %'Afferents'; %  'Neutrophil_Simone'; %'Pollen'; %'Astrocyte'; %    'NGC'; % 'Neutrophil'; % 'Vasculature %'Macrophage'

% Parse data table
dataTablePath = 'R:\Levy Lab\2photon\ImagingDatasets_Simone_240124.xlsx'; % 'R:\Levy Lab\2photon\ImagingDatasetsSimone2.xlsx'; %'D:\MATLAB\ImagingDatasets.xlsx'; % 'D:\MATLAB\NGCdata.xlsx';  Simone
dataTable = readcell(dataTablePath, 'sheet',dataSet);  % 'NGC', ''
colNames = dataTable(1,:); dataTable(1,:) = [];
dataCol = struct('mouse',find(contains(colNames, 'Mouse')), 'date',find(contains(colNames, 'Date')), 'FOV',find(contains(colNames, 'FOV')), 'vascChan',find(contains(colNames, 'VascChan')),...
    'volume',find(contains(colNames, 'Volume')), 'run',find(contains(colNames, 'Runs')), 'Ztop',find(contains(colNames, 'Zbot')), 'Zbot',find(contains(colNames, 'Ztop')), 'csd',find(contains(colNames, 'CSD')), 'ref',find(contains(colNames, 'Ref')), 'done',find(contains(colNames, 'Done')));
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);

% Initialize variables
expt = cell(1,Nexpt); 
runInfo = cell(1,Nexpt); 
Tscan = cell(1,Nexpt); 
loco = cell(1,Nexpt); %vessels(x) = cell(1,Nexpt);
vesselROI = cell(1,Nexpt); 
NvesselROI = cell(1,Nexpt); 
tifStackMax = cell(1,Nexpt);

% Choose which subset to  process
%xGLM = [18,22,24,30:32];
xPresent = 268; % xGLM; flip(100:102); %45:47; % [66:69]; %6;  62,64,
Npresent = numel(xPresent);

for x = xPresent  %30 %x2D % x2Dcsd % x3D %% 51
    % Experiment/run level metadata
    [expt{x}, runInfo{x}, regParams, projParam] = ParseDataTable(dataTable, x, dataCol, dataDir); % , regParams, projParam
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x});

    % Get locomotion data
    for runs = flip(expt{x}.runs) % expt{x}.runs %
        loco{x}(runs) = GetLocoData( runInfo{x}(runs), 'show',true );
    end
   
    %Get locomotion state
    loco{x} = GetLocoState(expt{x}, loco{x}, 'dir',strcat(dataDir, expt{x}.mouse,'\'), 'name',expt{x}.mouse, 'var','velocity', 'show',true); %
    
    %[regParam.refRun, regParam.refScan] = DetermineReference(expt{x}, Tscan{x}, loco{x});
    %expt{x}.refRun = regParam.refRun; % vestigal at this point?

    % load _projParam.mat if 'Unrecognized field name 'color' 
    % load('D:\2photon\Simone\Simone_Vasculature\WT0119\201221_FOV2\Projections\WT0119_201221_FOV2_projParam.mat', 'projParam')
    
    % Analyze vasculature
    catInfo(x) = ConcatenateRunInfo(expt{x}, runInfo{x}, 'overwrite',false);  % 'suffix','sbxcat',   
    [~,deform{x}] = GetDeformCat3D( expt{x}, catInfo(x), 'show',true, 'overwrite',false, 'window',find(Tscan{x}{1}<=32,1,'last') ); % true false
    projParam = GenerateExptProjections(expt{x}, catInfo(x), Tscan{x}); %  projParam
    [vesselROI{x}, NvesselROI{x}, tifStackMax{x}] = SegmentVasculature(expt{x}, projParam, 'overwrite',false, 'review',false );
    vesselROI{x} = GetVesselDiameter(expt{x}, projParam, vesselROI{x}, 'smooth',4, 'overwrite',false, 'show',true);
    
       % figure of vessel diameter changes over time
    for Z = find(NvesselROI{x})
        colorMat = distinguishable_colors(NvesselROI{x}(Z));
        close all;
        DiameterTime = figure;
        subplot(1,2,1);
        imshow(tifStackMax{x}{Z}, []);
        impixelinfo;
        for roi = 1:NvesselROI{x}(Z) 
            drawpolygon('Position', vesselROI{x}{Z}(roi).boxPosition.xy, 'color',colorMat(roi,:));  
        end
        title( sprintf('%s Z projection %i', expt{x}.name, Z), 'Interpreter','none' )

        tempDiam = [vesselROI{x}{Z}.diameter];
        subplot(1,2,2)
        for roi = 1:NvesselROI{x}(Z)
            plot(vesselROI{x}{Z}(roi).diameter.um_lp, 'color',colorMat(roi,:)); hold on;
            line( [1,projParam.totBin], vesselROI{x}{Z}(roi).diameter.um_max*[1,1], 'color',colorMat(roi,:), 'linestyle','--');
        end
        %plot( vertcat(tempDiam.um_gauss)')  % vertcat(Tscan{x}{:}), %plot(vesselROI{Z}.diameter.um, 'k')
        
        xlabel('Scan') %xlabel('Time (sec)'); 
        ylabel('Diameter (\mum)');
        axis tight; axis square;

        % save figure
        figPath = sprintf('%s%s_DiameterTime', expt{x}.dir, expt{x}.name);
    if ~exist(figPath, 'file') || overwrite
        fprintf('\nSaving %s', figPath);
        saveas(DiameterTime, figPath)
    end
        pause;
    end
    
     

%Combine vessels of different Z planes
    vesselROI{x} = [vesselROI{x}{~cellfun(@isempty, vesselROI{x})}];

    %{
    
    % Sort vessels by subtype
    [~,~,roiTypeInd] = unique({vesselROIpool.vesselType}); % A,D,V
    roiTypeSort = []; typeLims = nan(1,3); k = 0;
    for s = [2,1,3] %  D,A,V
        k = k+1;
        roiTypeSort = [roiTypeSort, find(roiTypeInd == s)'];  
        typeLims(k) = roiTypeSort(end);
    end 
    Ntype = [typeLims(1), diff(typeLims)];
    vesselROIpool = vesselROIpool(roiTypeSort);
    vessTypes = {'Dural','Artery','Vein'};
    vessTypes = vessTypes(Ntype > 0);

    diamPool = [vesselROIpool.diameter];
    allDiam = cat(1, diamPool.um_gauss)';
    allDiamZ = zscore(allDiam, [], 1);
    
    leftOpt =  {[0.1,0.02], [0.1,0.05], [0.1,0.1]};  % {[vert, horz], [bottom, top], [left, right]}
    rightOpt = {[0.05,0.02], [0.1,0.05], [0.1,0.1]};  % {[vert, horz], [bottom, top], [left, right]}
    close all; clearvars sp
    figure('WindowState','maximized');
    sp(1) = subtightplot(4,2,1:2:5,leftOpt{:});
    imagesc(allDiamZ'); caxis([-5,5]); colormap(bluewhitered); CB = colorbar; CB.Location = 'NorthOutside';
    set(gca,'Ytick', typeLims+0.5, 'YtickLabel',vessTypes(Ntype > 0), 'TickDir','out'); 
    xlabel('Scan'); ylabel('Vessel ROI')
    impixelinfo; 

    sp(2) = subtightplot(4,2,7,leftOpt{:});
    plot(vertcat(Tscan{x}{:}), vertcat(loco{x}(runs).Vdown)') %plot(vesselROI{Z}.diameter.um, 'k')
    xlabel('Time (sec)'); ylabel('Velocity (cm/s)');
    axis tight; % axis square;

    for Z = 1:projParam.Nz
        colorMat = distinguishable_colors(NvesselROI{x}(Z));
        subtightplot(projParam.Nz,2,2*Z,rightOpt{:});
        imshow(tifStackMax{x}{Z}, []);
        for roi = 1:NvesselROI{x}(Z), drawpolygon('Position', vesselROI{x}{Z}(roi).boxPosition.xy, 'color',colorMat(roi,:));  end
        title( sprintf('%s Z projection %i', expt{x}.name, Z), 'Interpreter','none' )
    end

    plot(vesselROIpool(1).diameter )

    %{
    subplot(2,2,2);
    imagesc(corr(allDiamZ)); axis square; caxis([-1,1]); colormap(bluewhitered); impixelinfo; 
    set(gca,'Xtick', typeLims+0.5, 'XtickLabel',vessTypes(Ntype > 0),'Ytick', typeLims+0.5, 'YtickLabel',vessTypes(Ntype > 0), 'TickDir','out'); 
    %}
    %}

end

%%
allVessROI = [vesselROI{xPresent}];
allVessROI.missing_frac
