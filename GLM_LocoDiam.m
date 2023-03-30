%% Use GLM to assess contribution of different variables

% Clear any previous variables in the Workspace and Command Window to start fresh
clear; clc; close all;

% TODO -- Set the directory of where animal folders are located
dataDir =  'D:\2photon\Simone\Simone_Vasculature\'; %  'D:\2photon\Simone\Simone_Macrophages\'; %  'Z:\2photon\Simone\Simone_Macrophages\'; %  'D:\2photon\'; %  'D:\2photon\Simone\'; % 

% PARSE DATA TABLE 

% TODO --- Set excel sheet
dataSet = 'Vasculature'; %'Macrophage'; 'AffCSD'; 'Pollen'; 'Vasculature'; %  'Astrocyte'; %  'Anatomy'; %  'Neutrophil_Simone'; %  'NGC'; % 'Neutrophil'; % 'Afferents'
[regParam, projParam] = DefaultProcessingParams(dataSet); % get default parameters for processing various types of data

regParam.method = 'affine'; %rigid 
regParam.name = 'affine'; %rigid  

% TODO --- Set data spreadsheet directory
dataTablePath = 'R:\Levy Lab\2photon\ImagingDatasets.xlsx'; % 'R:\Levy Lab\2photon\ImagingDatasetsSimone2.xlsx';
dataTable = readcell(dataTablePath, 'sheet',dataSet);  % 'NGC', ''
colNames = dataTable(1,:); dataTable(1,:) = [];
dataCol = struct('mouse',find(contains(colNames, 'Mouse')), 'date',find(contains(colNames, 'Date')), 'FOV',find(contains(colNames, 'FOV')), 'vascChan',find(contains(colNames, 'VascChan')),...
    'volume',find(contains(colNames, 'Volume')), 'run',find(contains(colNames, 'Runs')), 'Ztop',find(contains(colNames, 'Zbot')), 'Zbot',find(contains(colNames, 'Ztop')), 'csd',find(contains(colNames, 'CSD')), ...
    'ref',find(contains(colNames, 'Ref')), 'edges',find(contains(colNames, 'Edge')), 'Zproj',find(contains(colNames, 'Zproj')), 'done',find(contains(colNames, 'Done')));
Nexpt = size(dataTable, 1);
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false);

% Initialize variables
expt = cell(1,Nexpt); 
runInfo = cell(1,Nexpt); 
Tscan = cell(1,Nexpt); 
loco = cell(1,Nexpt); % Tcat = cell(1,Nexpt);
vesselROI = cell(1,Nexpt); 
NvesselROI = cell(1,Nexpt); 
tifStackMax = cell(1,Nexpt);

% Use GLM to assess contribution of different variables
locoDiam_pred = cell(1,Nexpt); %predictors
locoDiam_resp = cell(1,Nexpt); %response
locoDiam_opts = cell(1,Nexpt); %options
locoDiam_result = cell(1,Nexpt); 
locoDiam_summary = cell(1,Nexpt);

% TODO --- Specify xPresent - row number(X) within excel sheet
xPresent = 74; %[18,22,24,30:32]; % flip(100:102); %45:47; % [66:69];
Npresent = numel(xPresent);
overwrite = false;

GLMname = 'locoDiam';
figDir = 'D:\MATLAB\Figures\GLM_Vasculature\';  % CSD figures\
mkdir( figDir );

% Set GLM rate
GLMrate = 15.49/1; %15.49/16 %number of planes, it matches the slowest sampling rate

for x = xPresent % x3D % 

    % Parse data table
    [expt{x}, runInfo{x}, regParam, projParam] = ParseDataTable(dataTable, x, dataCol, dataDir, regParam, projParam);
    [Tscan{x}, runInfo{x}] = GetTime(runInfo{x});  % , Tcat{x}

    % Load vascular data
    [vesselROI{x}, NvesselROI{x}, tifStackMax{x}] = SegmentVasculature(expt{x}, projParam, 'overwrite',false, 'review',false );

    % Set GLMparallel options
    % housekeeping
    locoDiam_opts{x}.name = sprintf('%s_%s', expt{x}.name, GLMname); %strcat(expt{x}.name, , '_preCSDglm');
    locoDiam_opts{x}.rShow = NaN;
    locoDiam_opts{x}.figDir = ''; % figDir;

    % Signal processing parameters
    locoDiam_opts{x}.trainFrac = 0.75; % 1; %
    locoDiam_opts{x}.Ncycle = 20;
    locoDiam_opts{x}.minDev = 0.05; 
    locoDiam_opts{x}.minDevFrac = 0.1;
    locoDiam_opts{x}.maxP = 0.05;
    locoDiam_opts{x}.Nshuff = 0;  
    locoDiam_opts{x}.minShuff = 15; 
    locoDiam_opts{x}.window = [-4,4]; % [0,0]; % [-0.5, 0.5]; %consider temporal shifts this many seconds after/before response
    locoDiam_opts{x}.lopo = true; %false; %LOPO = Leave One Predictor Out

    % Downsample data to GLMrate target 
    locoDiam_opts{x}.frameRate = expt{x}.scanRate;  % GLMrate; %
    locoDiam_opts{x}.binSize = max([1,round(expt{x}.scanRate/GLMrate)]); 
    locoDiam_opts{x}.minShuffFrame = round( locoDiam_opts{x}.frameRate*locoDiam_opts{x}.minShuff );
    windowFrame = round(locoDiam_opts{x}.window*locoDiam_opts{x}.frameRate); %window(sec)*frameRate
    locoDiam_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    locoDiam_opts{x}.maxShift = max( abs(windowFrame) );
    locoDiam_opts{x}.Nshift = numel( locoDiam_opts{x}.shiftFrame );  %Nshift = preCSDOpts(x).Nshift;
    locoDiam_opts{x}.lags = locoDiam_opts{x}.shiftFrame/locoDiam_opts{x}.frameRate; %[-sec,+sec]
    locoDiam_opts{x}.xVar = 'Time';

    % GLMnet parameters - don't change without a good reason
    locoDiam_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    locoDiam_opts{x}.CVfold = 10;
    locoDiam_opts{x}.nlamda = 1000;
    locoDiam_opts{x}.maxit = 5*10^5;
    locoDiam_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    locoDiam_opts{x}.standardize = true; 

    % Get locomotion data
    for runs = expt{x}.runs % flip(expt{x}.runs) % 
        loco{x}(runs) = GetLocoData( runInfo{x}(runs), 'show',false ); 
        %plot(loco{x}(regParam.refRun).Vdown)
    end

    %Get locomotion state
    if any(cellfun(@isempty, {loco{x}.stateDown})) %isempty(loco{x}.stateDown)
            try
                loco{x} = GetLocoState(expt{x}, loco{x}, 'dir',strcat(dataDir, expt{x}.mouse,'\'), 'name',expt{x}.mouse, 'var','velocity', 'show',true); %
            catch
                fprintf('\nGetLocoState failed for %s', expt{x}.name)
            end
    end

    % Concatenate input variables pre-CSD
    % Define locomotion predictors
    tempVelocityCat = BinDownMean( vertcat(loco{x}(expt{x}.preRuns).Vdown), locoDiam_opts{x}.binSize );
    tempAccelCat = BinDownMean( abs(vertcat(loco{x}(expt{x}.preRuns).Adown)), locoDiam_opts{x}.binSize ); 
    tempStateCat = BinDownMean( vertcat(loco{x}(expt{x}.preRuns).stateDown), locoDiam_opts{x}.binSize );
    
    locoDiam_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); 
    locoDiam_pred{x}.data = [tempVelocityCat, tempAccelCat, tempStateCat];  
    locoDiam_pred{x}.name = {'Velocity', '|Accel|', 'State'}; % ,'Speed',  'Str-Exp', 'Str-Comp',
    locoDiam_pred{x}.N = size(locoDiam_pred{x}.data,2);
    for p = flip(1:locoDiam_pred{x}.N) 
        locoDiam_pred{x}.lopo.name{p} = ['No ',locoDiam_pred{x}.name{p}]; 
    end
    
    % Set up leave-one-family-out
    locoDiam_pred{x}.fam.col = {}; %{1:4, 5:7}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    locoDiam_pred{x}.fam.N = numel(locoDiam_pred{x}.fam.col); 
    locoDiam_pred{x}.fam.name = {}; %{'All'};%  'Onset Time',
    
%     Define Response
%     Vascular diamater data
%     vesselROIpool = [vesselROI{x}{~cellfun(@isempty, vesselROI{x})}];
%     bigVesselInd = find(strcmpi({vesselROI{x}{:}.vesselType}, {'A'}) | strcmpi({vesselROI{x}{:}.vesselType}, {'D'})); %Find index of cells containing string A(rtery) and D(ura)
%     diamPool = [vesselROI{x}(bigVesselInd).diameter]; 

    % Define Response
    vesselROIpool = [vesselROI{x}{:}];
    diamPool = [vesselROIpool.diameter];
    allDiam = cat(1, diamPool.um_lp)';
    allDiamZ = zscore(allDiam, [], 1);
    diamResp = BinDownMean( allDiam, locoDiam_opts{x}.binSize ); % allDiamZ
    locoDiam_resp{x}.data = diamResp; 
    locoDiam_resp{x}.N = size(locoDiam_resp{x}.data, 2); 
    locoDiam_resp{x}.name = sprintfc('Diameter %i', 1:locoDiam_resp{x}.N);
 
    % Remove scans with missing data 
    nanFrame = find(any(isnan([locoDiam_pred{x}.data, locoDiam_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    locoDiam_pred{x}.data(nanFrame,:) = []; 
    locoDiam_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    locoDiam_opts{x}.load = true; % false; % 
    locoDiam_opts{x}.saveRoot = expt{x}.dir; %''; %
    [locoDiam_result{x}, locoDiam_summary{x}, ~, locoDiam_pred{x}, locoDiam_resp{x}] = GLMparallel(locoDiam_pred{x}, locoDiam_resp{x}, locoDiam_opts{x}); 
    %locoDiam_summary{x} = SummarizeGLM(locoDiam_result{x}, locoDiam_pred{x}, locoDiam_resp{x}, locoDiam_opts{x});
end


%%
for x = xPresent
    locoDiam_opts{x}.rShow = 1:sum(NvesselROI{x}); %1:locoDiam_resp{x}.N; % 1:LocoDeform_resp{x}.N; %NaN;
    locoDiam_opts{x}.xVar = 'Time';
    ViewGLM(locoDiam_pred{x}, locoDiam_resp{x}, locoDiam_opts{x}, locoDiam_result{x}, locoDiam_summary{x}); %GLMresultFig = 
end


%% Compare GLM to data for each experiment
close all; clearvars sp SP;
PreGLMresults = figure('WindowState','maximized', 'color','w');
opt = {[0.02,0.07], [0.06,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
rightOpt = {[0.1,0.07], [0.1,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
jitterWidth = 0.45;
xAngle = 30;
Nrow = locoDiam_pred{xPresent(1)}(1).N+1; Ncol = 3;
spGrid = reshape( 1:Nrow*Ncol, Ncol, Nrow )';
for x = xPresent
    sp(locoDiam_pred{x}.N+1) = subtightplot(locoDiam_pred{x}.N+1, 3, 1:2, opt{:});
    imagesc( locoDiam_resp{x}.data' );
    ylabel('Diameter', 'Interpreter','none');
    title( sprintf('%s', expt{x}.name), 'Interpreter','none');
    set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'Xtick',[]); % , 'Ytick',onStruct(x).fluor.responder
    text( repmat(size(locoDiam_resp{x}.data,1)+1, locoDiam_summary{x}.Ngood, 1), locoDiam_summary{x}.rGood+0.5, '*', 'VerticalAlignment','middle', 'FontSize',8);
    impixelinfo;
    
    for v = 1:locoDiam_pred{x}.N
        sp(v) = subtightplot(Nrow, Ncol, spGrid(v+1,1:2), opt{:});
        plot( locoDiam_pred{x}.data(:,v) ); hold on;
        ylabel(locoDiam_pred{x}.name{v}, 'Interpreter','none');
        xlim([-Inf,Inf]);
        if v < locoDiam_pred{x}.N
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'XtickLabel',[]);
        else
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off');
        end
    end
    xlabel('Scan');
    
    subtightplot(3,3,3, rightOpt{:});
    bar([locoDiam_summary{x}.Ngood]/expt{x}.Nroi ); % numel(onStruct(x).fluor.responder),   , numel(rLocoPreFit{x})
    set(gca,'Xtick',1, 'XtickLabel',{'Fit'}, 'box','off'); % 'Loco','Fit','Both'  :3
    ylabel('Fraction of ROI');
    ylim([0,1]);
    
    subtightplot(3,3,6, rightOpt{:});
    JitterPlot( locoDiam_summary{x}.lopo.devFrac(:,locoDiam_summary{x}.rGood)', jitterWidth ); hold on;
    line([0,locoDiam_pred{x}.N+1], [1,1], 'color','k', 'lineStyle','--');
    xlim([0,locoDiam_pred{x}.N+1]); ylim([0,Inf]); 
    ylabel('Fraction of total deviance'); title('Leave One Predictor Out (well-fit units only)');
    set(gca, 'Xtick',1:locoDiam_pred{x}.N,  'XtickLabel', locoDiam_summary{x}.lopo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    
    
    subtightplot(3,3,9, rightOpt{:});
    JitterPlot( locoDiam_summary{x}.lofo.devFrac(:,locoDiam_summary{x}.rGood)', jitterWidth ); hold on;
    line([0,locoDiam_pred{x}.fam.N]+0.5, [1,1], 'color','k', 'lineStyle','--');
    xlim([0,locoDiam_pred{x}.fam.N]+0.5);
    ylabel('Fraction of total deviance'); title('Leave One Family Out (well-fit units only)');
    set(gca, 'Xtick',1:locoDiam_pred{x}.fam.N,  'XtickLabel', locoDiam_summary{x}.lofo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    ylim([0,Inf]);
    
    linkaxes(sp,'x');
    % {
    figPath = sprintf('%s%s_Deviance.tif', figDir, GLMname);
    if exist(figPath,'file'), delete(figPath); end
    fprintf('\nSaving %s', figPath);
    %export_fig( figPath, '-pdf', '-painters','-q101', '-append', LocoSensitivePrePost ); pause(1);
    %print(PreGLMresults, figPath, '-dtiff' ); 
    pause%(1);   
    clf;
    %}
    %pause; clf;
end
%%
for x = find(~cellfun(@isempty, locoDiam_result)) %  xPresent
    locoDiam_opts{x}.rShow = locoDiam_summary{x}.rGood; % 2; %1:7; %1:LocoDeform_resp{x}.N; %NaN; % 1:LocoDeform_resp{x}.N; %NaN;
    locoDiam_opts{x}.xVar = 'Time';
    ViewGLM(locoDiam_pred{x}, locoDiam_resp{x}, locoDiam_opts{x}, locoDiam_result{x}, locoDiam_summary{x}); %GLMresultFig = 
end

%% Divide units into Insensitiveensitive, Mixed, loco-only or deformation-only units
locoDiam_Nsubtype = []; k = 0;
for x = intersect( find(~cellfun(@isempty, locoDiam_result)), xPresent ) %xPresent
   
    k = k+1;
    locoDiam_Nsubtype(k,:) = [locoDiam_summary{x}.nInsensitive, locoDiam_summary{x}.nMixed, locoDiam_summary{x}.nDeform, locoDiam_summary{x}.nLoco]; % /expt{x}.Nroi
end
locoDiam_Nsubtype(k+1,:) = sum(locoDiam_Nsubtype, 1);
locoDiam_subtypeFrac = locoDiam_Nsubtype./repmat( sum(locoDiam_Nsubtype,2), 1, 4);

bar(locoDiam_subtypeFrac,'stacked')

%% Show single examples of each subtype
for x = 30 %xPresent
    locoDiam_opts{x}.xVar = 'Time';
    %{
    if locoDiam_summary{x}.nMixed > 0
        locoDiam_opts{x}.rShow = locoDiam_summary{x}.rMixed;
        ViewGLM(locoDiam_pred{x}, locoDiam_resp{x}, locoDiam_opts{x}, locoDiam_result{x}, locoDiam_summary{x});
    end
    %}
    if locoDiam_summary{x}.nDeform > 0
        locoDiam_opts{x}.rShow = locoDiam_summary{x}.rDeform;
        ViewGLM(locoDiam_pred{x}, locoDiam_resp{x}, locoDiam_opts{x}, locoDiam_result{x}, locoDiam_summary{x});
    end
    %{
    if locoDiam_summary{x}.nLoco > 0
        locoDiam_opts{x}.rShow = locoDiam_summary{x}.rLoco;
        ViewGLM(locoDiam_pred{x}, locoDiam_resp{x}, locoDiam_opts{x}, locoDiam_result{x}, locoDiam_summary{x});
    end
    %}
end

%% Pool results across experiments
preCSDdevPool = []; goodDevPool = [];
preCSDdevFracPool = []; %lofoDevFracPool = [];
for x = xPresent%find(~cellfun(@isempty, rLocoPreFit))
    if ~isempty( locoDiam_summary{x}.rGood )
        preCSDdevPool = [preCSDdevPool, locoDiam_summary{x}.dev]; 
        goodDevPool = [goodDevPool, locoDiam_summary{x}.dev( locoDiam_summary{x}.rGood )];
        preCSDdevFracPool = [preCSDdevFracPool, vertcat(locoDiam_summary{x}.lopo.devFrac(:,locoDiam_summary{x}.rGood), locoDiam_summary{x}.lofo.devFrac(:,locoDiam_summary{x}.rGood) )]; % rLocoPreFit{x}, :
    end
end
%% Summarize deviance explained
opt = {[0.02,0.07], [0.1,0.07], [0.09,0.09]};  % {[vert, horz], [bottom, top], [left, right] }
DevianceFig = figure('WindowState','maximized', 'color','w');
k = 1; clearvars h;
subtightplot(1,3,1,opt{:});
for x = xPresent
    [Ftemp, Xtemp] = ecdf( locoDiam_summary{x}.dev ); hold on;
    h(k) = plot(Xtemp, Ftemp, 'color',0.7*[1,1,1] );
    k = k + 1;
end
[Fdev, Xdev] = ecdf( preCSDdevPool );
h(k) = plot( Xdev, Fdev, 'color','k', 'LineWidth',2 ); 
axis square;
legend(h, {expt(xPresent).name, 'Pooled', 'Threshold'}, 'Location','southEast', 'Interpreter','none', 'AutoUpdate',false );
xlim([0, 0.6]);
line(locoDiam_opts{xPresent(1)}.minDev*[1,1], [0,1], 'Color','r', 'LineStyle','--'); % h(k+1) = 
xlabel('Deviance Explained'); ylabel('Fraction of Units');
title( sprintf('%s Fit Results', GLMname), 'Interpreter','none' );

subtightplot(1,3,2,opt{:});
JitterPlot( 1 - preCSDdevFracPool', 0.5, 'ErrorCap',10, 'monochrome',0.6); hold on;
line([0,size(preCSDdevFracPool,1)+1], [0,0], 'color','k');
axis square;
ylabel('Relative Explanatory Value'); %ylabel('Cumulative Distribution');
ylim([-1,1]);
set(gca,'Xtick', 1:size(preCSDdevFracPool,1), 'XtickLabel',[locoDiam_pred{xPresent(1)}.name, locoDiam_pred{xPresent(1)}.fam.name]);
xtickangle(30);
title('Well-Fit Units');

subtightplot(1,3,3,opt{:});
bar(locoDiam_subtypeFrac,'stacked');
ylabel('Fraction of Units');
barPos = get(gca, 'Position');
xlim([0, size(locoDiam_subtypeFrac,1)+1]);
axis square;
title('Subtype Breakdown');
legend('Insensitive','Mixed','Deform-only','Loco-only', 'Location','EastOutside');
set(gca, 'Xtick',1:size(locoDiam_subtypeFrac,1), 'XtickLabel',{expt(xPresent).name, 'Pooled'}, 'TickLabelInterpreter','none', 'FontSize',10, 'Position',barPos );
xtickangle(30);

% Save the figure
figPath = sprintf('%s%s_DevianceResults.tif', figDir, GLMname);
if exist(figPath,'file'), delete(figPath); end
fprintf('\nSaving %s', figPath);
print(DevianceFig, figPath, '-dtiff', '-r300'); %pause(1); clf;

%% Show mean coeff by type
Ntype = 4;
subtypeColor = distinguishable_colors(Ntype);
close all;
SubtypeCoeffFig = figure('WindowState','maximized', 'color','w');
opt = {[0.02,0.07], [0.1,0.1], [0.1,0.06]};  % {[vert, horz], [bottom, top], [left, right] }
LW = 1.5;
colororder( subtypeColor ) 
for x = xPresent
    tempInsensitiveCoeff = cat(3, locoDiam_result{x}( locoDiam_summary{x}.rIns ).coeff );
    meanInsensitiveCoeff = mean(tempInsensitiveCoeff, 3, 'omitnan' );
    
    tempMixedCoeff = cat(3, locoDiam_result{x}( locoDiam_summary{x}.rMixed ).coeff );
    meanMixedCoeff = mean(tempMixedCoeff, 3, 'omitnan' );
    
    tempDeformCoeff = cat(3, locoDiam_result{x}( locoDiam_summary{x}.rDeform ).coeff );
    meanDeformCoeff = mean(tempDeformCoeff, 3, 'omitnan' );
    
    tempLocoCoeff = cat(3, locoDiam_result{x}( locoDiam_summary{x}.rLoco ).coeff );
    meanLocoCoeff = mean(tempLocoCoeff, 3, 'omitnan' );
    
    sp(1) = subtightplot(1,4,1,opt{:});
    if locoDiam_summary{x}.nIns > 0
        plot(locoDiam_opts{x}.lags,  meanInsCoeff, 'LineWidth',LW );
    end
    axis square;
    tempPos = get(gca,'Position');
    xlabel('Lag (s)'); ylabel('Coefficient'); 
    title( sprintf('Insensitive (n = %i)', locoDiam_summary{x}.nIns) );
    legend(locoDiam_pred{x}.name, 'Location','NorthWest', 'AutoUpdate',false)
    set(gca,'Position',tempPos);
    
    sp(2) = subtightplot(1,4,2,opt{:});
    if locoDiam_summary{x}.nMixed > 0
        plot(locoDiam_opts{x}.lags,  meanMixedCoeff, 'LineWidth',LW ); 
    end
    axis square;
    title( sprintf('Mixed (n = %i)', locoDiam_summary{x}.nMixed) );
    xlabel('Lag (s)'); 
    
    sp(3) = subtightplot(1,4,3,opt{:});
    if locoDiam_summary{x}.nDeform > 0
        plot(locoDiam_opts{x}.lags,  meanDeformCoeff, 'LineWidth',LW ); 
    end
    axis square;
    title( sprintf('Deformation-dependent (n = %i)', locoDiam_summary{x}.nDeform) );
    xlabel('Lag (s)'); %title('Deformation-dependent'); % ylabel('Coefficient');
    
    sp(4) = subtightplot(1,4,4,opt{:});
    if locoDiam_summary{x}.nLoco > 0
        plot(locoDiam_opts{x}.lags,  meanLocoCoeff, 'LineWidth',LW );
    end
    axis square;
    xlabel('Lag (s)'); 
    title(sprintf('Locomotion-dependent (n = %i)', locoDiam_summary{x}.nLoco)); % ylabel('Coefficient'); 
    linkaxes(sp,'xy');
    
    pause;
    
    % Save the figure
    figPath = sprintf('%s%s_%s_SubtypeCoeff.tif', figDir, GLMname, expt{x}.name);
    if exist(figPath,'file'), delete(figPath); end
    %print(SubtypeCoeffFig, figPath, '-dtiff', '-r300' );  fprintf('\nSaved %s\n', figPath);

    clf;
end

%% Show  coeff by subtype
subtype = {'Insensitive', 'Mixed', 'Deform', 'Loco'};
Nsubtype = 4;
pooledCoeff = struct('Insensitive',[], 'Mixed',[], 'Deform',[], 'Loco',[]);
for x = xPresent
    pooledCoeff.Insensitive = cat(3, pooledCoeff.Insensitive, locoDiam_result{x}( locoDiam_summary{x}.rIns ).coeff );
    pooledCoeff.Mixed = cat(3, pooledCoeff.Mixed, locoDiam_result{x}( locoDiam_summary{x}.rMixed).coeff );
    pooledCoeff.Deform = cat(3, pooledCoeff.Deform, locoDiam_result{x}( locoDiam_summary{x}.rDeform ).coeff );
    pooledCoeff.Loco = cat(3, pooledCoeff.Loco, locoDiam_result{x}( locoDiam_summary{x}.rLoco ).coeff );
    
end
pooledLags = locoDiam_opts{xPresent(1)}.lags;

Ncol = locoDiam_pred{xPresent(1)}.N;
k = 0;
close all;
opt = {[0.03,0.03], [0.07,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] }
SubtypeCoeffFig = figure('WindowState','maximized', 'color','w');
for row = 1:Ntype
    for col = 1:locoDiam_pred{xPresent(1)}.N
        k = k + 1;
        subtightplot(Ntype, Ncol, k, opt{:});
        plot( pooledLags, squeeze(pooledCoeff.(subtype{row})(:,col,:) ), 'color',[0,0,0,0.1] );
        axis square;
        if row == 1, title( locoDiam_pred{xPresent(1)}.name{col} ); end
        if col == 1, ylabel( sprintf('%s coeff', subtype{row} )); end
        if row == Nsubtype, xlabel('Lag (s)'); end
    end
    %pause;
end
% Save the figure
figPath = sprintf('%s%s_SubtypeCoeff.tif', figDir, GLMname);
if exist(figPath,'file'), delete(figPath); end
print(SubtypeCoeffFig, figPath, '-dtiff', '-r300' );  fprintf('\nSaved %s\n', figPath);

%%
devSummMat = cell(1,Nexpt); devSummPool = [];
opt = {[0.02,0.05], [0.02,0.02], [0.06,0.03]};  % {[vert, horz], [bottom, top], [left, right] }
GLMresultFig = figure('WindowState','maximized', 'color','w');

for x = xPresent
    cla;
    devSummMat{x} = [locoDiam_summary{x}.dev; locoDiam_summary{x}.lopo.dev; locoDiam_summary{x}.lofo.dev];
    %devSummPool = [devSummPool, devSummMat{x}];
    % {
    subtightplot(1,1,1,opt{:});
    imagesc( devSummMat{x} ); %imagesc( devSummPool ); %
    axis image;
    set(gca, 'Ytick', 1:size(devSummMat{x}, 1), 'YtickLabel', [{'All'}; locoDiam_summary{x}.lopo.name(:); locoDiam_summary{x}.lofo.name(:)], ...
        'Xtick',1:10:locoDiam_resp{x}.N, 'TickDir','out', 'TickLength',[0.003,0], 'FontSize',8); % 'XtickLabel',locoDiam_resp{x}.name
    title(sprintf('%s: GLM Fit Summary', locoDiam_opts{x}.name), 'Interpreter','none');
    %xtickangle(30);
    
    CB = colorbar; CB.Label.String = 'Deviance Explained';
    impixelinfo;
    pause;
    %}
end

devSummCat = cat(2, devSummMat{:});
devSummMed = median(devSummCat, 2, 'omitnan' );

GLMdevFig = figure('WindowState','maximized', 'color','w');
imagesc( devSummMed' );
set(gca,'XtickLabel',  [{'All'}; locoDiam_summary{x}.lopo.name(:); locoDiam_summary{x}.lofo.name(:)], ...
    'YtickLabel',locoDiam_resp{x}.name, 'TickDir','out', 'TickLength',[0.003,0]);
title('GLM Fit Summary (All 3D Data)', 'Interpreter','none');
xtickangle(30);
axis image;
CB = colorbar; CB.Label.String = 'Median Deviance Explained';
figPath = sprintf('%s%s_DevSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
print( GLMdevFig, figPath, '-dtiff');
impixelinfo;

%% Plot coefficient values for each predictor, unit and experiment
close all; clearvars sp SP;
PreGLMcoeff = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.04], [0.09,0.05], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
for v = 1:locoDiam_pred{x}.N
    c = 0;
    for x = xPresent %find(~cellfun(@isempty, rLocoPreFit)) %
        c = c+1;
        if locoDiam_summary{x}.Ngood > 0
            subtightplot(1, locoDiam_pred{x}.N, v, opt{:});
            plot(c, locoDiam_summary{x}.peakCoeff(locoDiam_summary{x}.rGood,v), 'k.' ); hold on; %rLocoPreFit{x}
        end
    end
    line([0,c+1], [0,0], 'color','k');
    title( locoDiam_pred{x}.name{v}, 'Interpreter','none' );
    if v == 1, ylabel('Peak Coefficient'); end
    xlim([0,c+1]);
    set(gca, 'Xtick', 1:c, 'XtickLabel', {expt(xPresent).name}, 'TickLabelInterpreter','none' );
    xtickangle(45);
    axis square;
end

%% Plot peak coefficient vs latency values for each predictor, good unit, and experiment
close all; clearvars sp SP;
PreGLMcoeff = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.04], [0.09,0.05], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
exptColor = distinguishable_colors(numel(xPresent));
for v = 1:locoDiam_pred{x}.N
    c = 0;
    for x = xPresent %find(~cellfun(@isempty, rLocoPreFit)) %
        c = c+1;
        if locoDiam_summary{x}.Ngood > 0
            subtightplot(1, locoDiam_pred{x}.N, v, opt{:});
            plot(locoDiam_summary{x}.peakLag(locoDiam_summary{x}.rGood,v), locoDiam_summary{x}.peakCoeff(locoDiam_summary{x}.rGood,v), '.', 'color',exptColor(c,:) ); hold on; %rLocoPreFit{x}
        end
    end
    title( locoDiam_pred{xPresent(1)}.name{v}, 'Interpreter','none' );
    if v == 1, ylabel('Coefficient'); end
    xlabel('Lag (s)');
    xlim([-6,6]);
    axis square;
end

%% Identify the best-fit units from each experiment
for x = xPresent
    [devSort, rDevSort] = sort( [locoDiam_result{x}.dev], 'descend' );
    %[maxDevExp, rMaxDevExp] = max( [locoDiam_result{x}.dev] );
    WriteROIproj(expt{x}, ROI{x}, 'edges',segParams{x}.edges, 'overwrite',true, 'rSet',rDevSort(1:10), 'buffer',20*[1,1,1,1]); % ROI{x} =   
    %ViewResults3D( expt{x}, Tscan{x}, deform{x}, loco{x}, fluor{x}, allVars, ROI{x}, 'cat',true, 'limits', viewLims ); 
end

%% Check mechanosensitive units for sigmoidal stimulus response curve to various forms of deformation PRE-CSD
responderMat_pre = cell(1,Nexpt);
sigmResp_pre = repmat( struct('speed',[], 'Nspeed',NaN, 'speedFrac',NaN, 'trans',[], 'Ntrans',NaN, 'transFrac',NaN, 'scale',[], 'Nscale',NaN, 'scaleFrac',NaN, 'stretch',[], 'Nstretch',NaN, 'stretchFrac',NaN, ...
    'shear',[], 'Nshear',NaN, 'shearFrac',NaN, 'shearRate',[], 'NshearRate',NaN, 'shearRateFrac',NaN, 'poly',[], 'Npoly',NaN, 'polyFrac',NaN), 1, Nexpt);
tic
for x = x2D
    rCheck = 1:expt{x}.Nroi; % [locoDiam_summary{x}.rDeform, locoDiam_summary{x}.rMixed];
    Ncheck = numel(rCheck);
    responderMat_pre{x} = nan(5, expt{x}.Nroi);
    responderMat_pre{x}(:,rCheck) = 0;
    
    fluorPre = [fluor{x}(expt{x}.preRuns).dFF];
    fluorPre = vertcat(fluorPre.ROI);
    deformPre = deform{x}(expt{x}.preRuns);
    scaleMagPre = vertcat(deformPre.scaleMag);
    scaleMagPre = mean(scaleMagPre(:,segParams{x}.zProj), 2, 'omitnan');
    %{
    speedPre = loco{x}(expt{x}.preRuns);
    speedPre = vertcat(speedPre.speedDown);
    transMagPre = vertcat(deformPre.transMag);
    transMagPre = mean(transMagPre(:,segParams{x}.zProj), 2, 'omitnan');
    stretchMagPre = vertcat(deformPre.stretchMag);
    stretchMagPre = mean(stretchMagPre(:,segParams{x}.zProj), 2, 'omitnan');
    shearMagPre = vertcat(deformPre.shearMag);
    shearMagPre = mean(shearMagPre(:,segParams{x}.zProj), 2, 'omitnan');
    shearRateMagPre = vertcat(deformPre.DshearMag);
    shearRateMagPre = mean(shearRateMagPre(:,segParams{x}.zProj), 2, 'omitnan');
    %zshiftMag 
    %}
    % Scaling
    [scaleStim_pre{x}, scaleResp_pre{x}, scaleResult_pre{x}, scaleSumm_pre{x}, scaleFit_pre{x}] = StimResponse(scaleMagPre, fluorPre, 0, 10, 0, 'fit',true, 'show',false); % , 'show',false  dffResp{x}
    sigmResp_pre(x).scale = intersect(rCheck, scaleSumm_pre{x}.rGood);
    sigmResp_pre(x).Nscale = numel(sigmResp_pre(x).scale);
    sigmResp_pre(x).scaleFrac = sigmResp_pre(x).Nscale/Ncheck;
    responderMat_pre{x}(2,sigmResp_pre(x).scale) = 1;
    sigmResp_pre(x).scaleThresh = [scaleResult{x}(sigmResp_pre(x).scale).thresh];
    sigmResp_pre(x).scaleSlope = [scaleResult{x}(sigmResp_pre(x).scale).rate];
    %{
    % Speed
    [speedStim{x}, speedResp{x}, speedResult{x}, speedSumm{x}, speedFit{x}] = StimResponse(speedPre, fluorPre, 0, 10, 0, 'fit',true); %   dffResp{x} , 'show',true
    sigmResp_pre(x).speed = intersect([locoDiam_summary{x}.rLoco, locoDiam_summary{x}.rMixed], speedSumm{x}.rGood);
    sigmResp_pre(x).Nspeed = numel(sigmResp_pre(x).speed);
    sigmResp_pre(x).speedFrac = sigmResp_pre(x).Nspeed/numel([locoDiam_summary{x}.rLoco, locoDiam_summary{x}.rMixed]);
    sigmResp_pre(x).speedThresh = [speedResult{x}(sigmResp_pre(x).speed).thresh];
    sigmResp_pre(x).speedSlope = [speedResult{x}(sigmResp_pre(x).speed).rate];
    % Translation
    [transStim{x}, transResp{x}, transResult{x}, transSumm{x}, transFit{x}] = StimResponse(transMagPre, fluorPre, 0, 10, 0, 'fit',true); % , 'show',false  dffResp{x}
    sigmResp_pre(x).trans = intersect(rCheck, transSumm{x}.rGood);
    sigmResp_pre(x).Ntrans = numel(sigmResp_pre(x).trans);
    sigmResp_pre(x).transFrac = sigmResp_pre(x).Ntrans/Ncheck;
    responderMat_pre{x}(1,sigmResp_pre(x).trans) = 1;
    sigmResp_pre(x).transThresh = [transResult{x}(sigmResp_pre(x).trans).thresh];
    sigmResp_pre(x).transSlope = [transResult{x}(sigmResp_pre(x).trans).rate];

    % Stretch
    [stretchStim{x}, stretchResp{x}, stretchResult{x}, stretchSumm{x}, stretchFit{x}] = StimResponse(stretchMagPre, fluorPre, 0, 10, 0, 'fit',true); % , 'show',false
    sigmResp_pre(x).stretch = intersect(rCheck, stretchSumm{x}.rGood);
    sigmResp_pre(x).Nstretch = numel(sigmResp_pre(x).stretch);
    sigmResp_pre(x).stretchFrac = sigmResp_pre(x).Nstretch/Ncheck;
    responderMat_pre{x}(3,sigmResp_pre(x).stretch) = 1;
    sigmResp_pre(x).stretchThresh = [stretchResult{x}(sigmResp_pre(x).stretch).thresh];
    sigmResp_pre(x).stretchSlope = [stretchResult{x}(sigmResp_pre(x).stretch).rate];
    % Shear
    [shearStim{x}, shearResp{x}, shearResult{x}, shearSumm{x}, shearFit{x}] = StimResponse(shearMagPre, fluorPre, 0, 10, 0, 'fit',true); % , 'show',false
    sigmResp_pre(x).shear = intersect(rCheck, shearSumm{x}.rGood);
    sigmResp_pre(x).Nshear = numel(sigmResp_pre(x).shear);
    sigmResp_pre(x).shearFrac = sigmResp_pre(x).Nshear/Ncheck;
    responderMat_pre{x}(4,sigmResp_pre(x).shear) = 1;
    sigmResp_pre(x).shearThresh = [shearResult{x}(sigmResp_pre(x).shear).thresh];
    sigmResp_pre(x).shearSlope = [shearResult{x}(sigmResp_pre(x).shear).rate];
    % Shear rate
    [shearRateStim{x}, shearRateResp{x}, shearRateResult{x}, shearRateSumm{x}, shearRateFit{x}] = StimResponse(shearRateMagPre, fluorPre, 0, 10, 0, 'fit',true); % , 'show',false
    sigmResp_pre(x).shearRate = intersect(rCheck, shearRateSumm{x}.rGood);
    sigmResp_pre(x).NshearRate = numel(sigmResp_pre(x).shearRate);
    sigmResp_pre(x).shearRateFrac = sigmResp_pre(x).NshearRate/Ncheck;
    responderMat_pre{x}(5,sigmResp_pre(x).shearRate) = 1;
    sigmResp_pre(x).shearRateThresh = [shearRateResult{x}(sigmResp_pre(x).shearRate).thresh];
    sigmResp_pre(x).shearRateSlope = [shearRateResult{x}(sigmResp_pre(x).shearRate).rate];  
    
    % Poly-responders
    sigmResp_pre(x).poly = find( sum(responderMat_pre{x}, 1) > 1);
    sigmResp_pre(x).Npoly = numel(sigmResp_pre(x).poly);
    sigmResp_pre(x).polyFrac = sigmResp_pre(x).Npoly/Ncheck;
    %}
    toc
    %bar( [deformResponders(x).transFrac, deformResponders(x).scaleFrac, deformResponders(x).stretchFrac, deformResponders(x).shearFrac] )
end

%% Check mechanosensitive units for sigmoidal stimulus response curve to various forms of deformation POST-CSD
postCSD_cutoff = 15; % how many minutes after CSD to start grabbing data?
responderMat_post = cell(1,Nexpt);
sigmResp_post = repmat( struct('speed',[], 'Nspeed',NaN, 'speedFrac',NaN, 'trans',[], 'Ntrans',NaN, 'transFrac',NaN, 'scale',[], 'Nscale',NaN, 'scaleFrac',NaN, 'stretch',[], 'Nstretch',NaN, 'stretchFrac',NaN, ...
    'shear',[], 'Nshear',NaN, 'shearFrac',NaN, 'shearRate',[], 'NshearRate',NaN, 'shearRateFrac',NaN, 'poly',[], 'Npoly',NaN, 'polyFrac',NaN), 1, Nexpt);
tic
for x = x3Dcsd %xPresent
    postRuns = 3:4;
    Tpost = vertcat(Tscan{x}{postRuns});
    Tpost = Tpost - Tpost(1);
    postScans = find( Tpost > 60*postCSD_cutoff );
    fluorPost = [fluor{x}(postRuns).dFF];
    fluorPost = vertcat(fluorPost.ROI);
    fluorPost = fluorPost(postScans,:);
    deformPost = deform{x}(postRuns);
    scaleMagPost = vertcat(deformPost.scaleMag);
    scaleMagPost = scaleMagPost(postScans,:);
    scaleMagPost = mean(scaleMagPost(:,segParams{x}.zProj), 2, 'omitnan');
    %{
    speedPost = loco{x}(postRuns);
    speedPost = vertcat(speedPost.speedDown);
    speedPost = speedPost(postScans,:);
    
    transMagPost = vertcat(deformPost.transMag);
    transMagPost = transMagPost(postScans,:);

    stretchMagPost = vertcat(deformPost.stretchMag);
    stretchMagPost = stretchMagPost(postScans,:);
    shearMagPost = vertcat(deformPost.shearMag);
    shearMagPost = shearMagPost(postScans,:);
    shearRateMagPost = vertcat(deformPost.DshearMag);
    shearRateMagPost = shearRateMagPost(postScans,:);
    %}
    rCheck = 1:expt{x}.Nroi; %  [locoDiam_summary{x}.rDeform, locoDiam_summary{x}.rMixed];
    Ncheck = numel(rCheck);
    responderMat_post{x} = nan(5, expt{x}.Nroi);
    responderMat_post{x}(:,rCheck) = 0;
    %{
    % Speed
    [speedStim{x}, speedResp{x}, speedResult{x}, speedSumm{x}, speedFit{x}] = StimResponse(speedPost, fluorPost, 0, 10, 0, 'fit',true); %   dffResp{x} , 'show',true
    sigmResp_post(x).speed = intersect([locoDiam_summary{x}.rLoco, locoDiam_summary{x}.rMixed], speedSumm{x}.rGood);
    sigmResp_post(x).Nspeed = numel(sigmResp_post(x).speed);
    sigmResp_post(x).speedFrac = sigmResp_post(x).Nspeed/numel([locoDiam_summary{x}.rLoco, locoDiam_summary{x}.rMixed]);
    sigmResp_post(x).speedThresh = [speedResult{x}(sigmResp_post(x).speed).thresh];
    sigmResp_post(x).speedSlope = [speedResult{x}(sigmResp_post(x).speed).rate];
    % Translation
    [transStim{x}, transResp{x}, transResult{x}, transSumm{x}, transFit{x}] = StimResponse(transMagPost, fluorPost, 0, 10, 0, 'fit',true); % , 'show',false  dffResp{x}
    sigmResp_post(x).trans = intersect(rCheck, transSumm{x}.rGood);
    sigmResp_post(x).Ntrans = numel(sigmResp_post(x).trans);
    sigmResp_post(x).transFrac = sigmResp_post(x).Ntrans/Ncheck;
    responderMat_post{x}(1,sigmResp_post(x).trans) = 1;
    sigmResp_post(x).transThresh = [transResult{x}(sigmResp_post(x).trans).thresh];
    sigmResp_post(x).transSlope = [transResult{x}(sigmResp_post(x).trans).rate];
    %}
    % Scaling
    [scaleStim_post{x}, scaleResp_post{x}, scaleResult_post{x}, scaleSumm_post{x}, scaleFit_post{x}] = ...
        StimResponse(scaleMagPost, fluorPost, 0, 10, 0, 'fit',true, 'show',false); % , 'show',false  dffResp{x}
    sigmResp_post(x).scale = intersect(rCheck, scaleSumm{x}.rGood);
    sigmResp_post(x).Nscale = numel(sigmResp_post(x).scale);
    sigmResp_post(x).scaleFrac = sigmResp_post(x).Nscale/Ncheck;
    responderMat_post{x}(2,sigmResp_post(x).scale) = 1;
    sigmResp_post(x).scaleThresh = [scaleResult{x}(sigmResp_post(x).scale).thresh];
    sigmResp_post(x).scaleSlope = [scaleResult{x}(sigmResp_post(x).scale).rate];
    %{
    % Stretch
    [stretchStim{x}, stretchResp{x}, stretchResult{x}, stretchSumm{x}, stretchFit{x}] = StimResponse(stretchMagPost, fluorPost, 0, 10, 0, 'fit',true); % , 'show',false
    sigmResp_post(x).stretch = intersect(rCheck, stretchSumm{x}.rGood);
    sigmResp_post(x).Nstretch = numel(sigmResp_post(x).stretch);
    sigmResp_post(x).stretchFrac = sigmResp_post(x).Nstretch/Ncheck;
    responderMat_post{x}(3,sigmResp_post(x).stretch) = 1;
    sigmResp_post(x).stretchThresh = [stretchResult{x}(sigmResp_post(x).stretch).thresh];
    sigmResp_post(x).stretchSlope = [stretchResult{x}(sigmResp_post(x).stretch).rate];
    % Shear
    [shearStim{x}, shearResp{x}, shearResult{x}, shearSumm{x}, shearFit{x}] = StimResponse(shearMagPost, fluorPost, 0, 10, 0, 'fit',true); % , 'show',false
    sigmResp_post(x).shear = intersect(rCheck, shearSumm{x}.rGood);
    sigmResp_post(x).Nshear = numel(sigmResp_post(x).shear);
    sigmResp_post(x).shearFrac = sigmResp_post(x).Nshear/Ncheck;
    responderMat_post{x}(4,sigmResp_post(x).shear) = 1;
    sigmResp_post(x).shearThresh = [shearResult{x}(sigmResp_post(x).shear).thresh];
    sigmResp_post(x).shearSlope = [shearResult{x}(sigmResp_post(x).shear).rate];
    % Shear rate
    [shearRateStim{x}, shearRateResp{x}, shearRateResult{x}, shearRateSumm{x}, shearRateFit{x}] = StimResponse(shearRateMagPost, fluorPost, 0, 10, 0, 'fit',true); % , 'show',false
    sigmResp_post(x).shearRate = intersect(rCheck, shearRateSumm{x}.rGood);
    sigmResp_post(x).NshearRate = numel(sigmResp_post(x).shearRate);
    sigmResp_post(x).shearRateFrac = sigmResp_post(x).NshearRate/Ncheck;
    responderMat_post{x}(5,sigmResp_post(x).shearRate) = 1;
    sigmResp_post(x).shearRateThresh = [shearRateResult{x}(sigmResp_post(x).shearRate).thresh];
    sigmResp_post(x).shearRateSlope = [shearRateResult{x}(sigmResp_post(x).shearRate).rate];  
    
    % Poly-responders
    sigmResp_post(x).poly = find( sum(responderMat_post{x}, 1) > 1);
    sigmResp_post(x).Npoly = numel(sigmResp_post(x).poly);
    sigmResp_post(x).polyFrac = sigmResp_post(x).Npoly/Ncheck;
    %}
    
    toc
    %bar( [deformResponders(x).transFrac, deformResponders(x).scaleFrac, deformResponders(x).stretchFrac, deformResponders(x).shearFrac] )
end

%%
for x = x3Dcsd
    for r = sigmResp_pre(x).scale %1
        % Plot pre-CSD sigmoid
        errorbar( scaleStim_pre{x}.mean, scaleResp_pre{x}.mean(:,r), scaleResp_pre{x}.sem(:,r), 'LineStyle','none', 'Color','k' ); hold on;
        xRange = linspace(0,scaleStim_pre{x}.mean(end))';  %(0:0.01:scaleStim_pre{x}.mean(end))';
        plot( xRange, predict(scaleFit_pre{x}{r}, xRange), 'b' );
        % Plot post-CSD sigmoids
        errorbar( scaleStim_post{x}.mean, scaleResp_post{x}.mean(:,r), scaleResp_post{x}.sem(:,r), 'LineStyle','none', 'Color','k' ); hold on;
        xRange = linspace(0,scaleStim_post{x}.mean(end))';  %(0:0.01:scaleStim_post{x}.mean(end))';
        plot( xRange, predict(scaleFit_post{x}{r}, xRange), 'r' );

        axis square;
        pause;
        cla;
    end
end

%% Are mechanosensitive units segregated by layer?
zIns = cell(1,Nexpt);  zMech = cell(1,Nexpt);
for x = x3D
    %{

    plot(locoDiam_summary{x}.lofo.devFrac(1,:), zROI,  '.');
    ylabel('Z Position'); xlabel('Fraction of deviance explained'); title('Deformation');
    axis square;
    pause;
    cla;
    %}
    zROI = vertcat(ROI{x}.cent);
    zROI = zROI(:,3);
    zIns{x} = zROI(locoDiam_summary{x}.rIns);
    zMixed = zROI(locoDiam_summary{x}.rMixed);
    zDeform = zROI(locoDiam_summary{x}.rDeform);
    zMech{x} = [zMixed;zDeform];
    zLoco = zROI(locoDiam_summary{x}.rLoco);

    %{
    zType = cell2padmat({zROI, zIns{x}, zMech{x}, zDeform, zMixed, zLoco} );
    violin(zType); hold on;
    %boxplot(zType); hold on;
    
    pZ = ranksum(zIns{x}, zMech{x});
    line([2,3], expt{x}.Nplane*[1,1], 'color','k'); 
    text(2.5, expt{x}.Nplane+1, num2str(pZ), 'HorizontalAlignment','center' )
    title( sprintf('x = %i: %s', x, expt{x}.name), 'Interpreter','none' );
    set(gca, 'Xtick',1:5, 'XtickLabel', ...
        {sprintf('All (n=%i)',expt{x}.Nroi), ...
        sprintf('Ins (n=%i)',locoDiam_summary{x}.nIns), ...
        sprintf('Mech (n=%i)',numel(zMech{x})), ...
        sprintf('Deform (n=%i)',locoDiam_summary{x}.nDeform), ...
        sprintf('Mixed (n=%i)',locoDiam_summary{x}.nMixed), ...
        sprintf('Loco (n=%i)',locoDiam_summary{x}.nLoco)}   )
    xtickangle(30);
    ylim([1, expt{x}.Nplane+2]);
    pause; cla;
    %}
end
zMed = [cellfun(@median, zIns(x3D))', cellfun(@median, zMech(x3D))'];
plot( [1,2], zMed, 'color',0.7*[1,1,1]); hold on;
plot( [1,2], zMed, '.', 'markersize',10);
xlim([0,3]);
signrank(zMed(:,1), zMed(:,2) )
