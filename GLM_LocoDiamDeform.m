%% Use GLM to assess contribution of different variables
locoDiamDeform_pred = cell(1,Nexpt); 
locoDiamDeform_resp = cell(1,Nexpt); 
locoDiamDeform_opts = cell(1,Nexpt); 
locoDiamDeform_result = cell(1,Nexpt); 
locoDiamDeform_summary = cell(1,Nexpt);
GLMname = 'locoDiamDeform_V';
GLMrate = 1; %15.49/30;
for x = xPresent % x3D % 
    projParam = GenerateExptProjections(expt{x}, catInfo(x), Tscan{x}); % get projection parameters
    % GLMparallel options
    locoDiamDeform_opts{x}.name = sprintf('%s_%s', expt{x}.name, GLMname); %strcat(expt{x}.name, , '_preCSDglm');
    locoDiamDeform_opts{x}.rShow = NaN;
    locoDiamDeform_opts{x}.figDir = ''; % figDir;
    locoDiamDeform_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    locoDiamDeform_opts{x}.standardize = false; 
    locoDiamDeform_opts{x}.trainFrac = 0.75; % 1; %
    locoDiamDeform_opts{x}.Ncycle = 20;
    locoDiamDeform_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    locoDiamDeform_opts{x}.CVfold = 10;
    locoDiamDeform_opts{x}.nlamda = 1000;
    locoDiamDeform_opts{x}.maxit = 5*10^5;
    locoDiamDeform_opts{x}.minDev = 0.05; 
    locoDiamDeform_opts{x}.minDevFrac = 0.1;
    locoDiamDeform_opts{x}.maxP = 0.05;
    locoDiamDeform_opts{x}.Nshuff = 0;  
    locoDiamDeform_opts{x}.minShuff = 15; 
    locoDiamDeform_opts{x}.window = [-60,60]; % [0,0]; % [-0.5, 0.5]; % 
    locoDiamDeform_opts{x}.lopo = true; %false; %
    % downsample data to GLMrate target 
    locoDiamDeform_opts{x}.binSize = max([1,round(expt{x}.scanRate/GLMrate)]);  %projParam.bin; %;
    locoDiamDeform_opts{x}.frameRate = expt{x}.scanRate/locoDiamDeform_opts{x}.binSize; % GLMrate; %projParam.rate_bin; %;  %  %
    locoDiamDeform_opts{x}.minShuffFrame = round( locoDiamDeform_opts{x}.frameRate*locoDiamDeform_opts{x}.minShuff );
    windowFrame = round(locoDiamDeform_opts{x}.window*locoDiamDeform_opts{x}.frameRate);
    locoDiamDeform_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    locoDiamDeform_opts{x}.maxShift = max( abs(windowFrame) );
    locoDiamDeform_opts{x}.Nshift = numel( locoDiamDeform_opts{x}.shiftFrame );  %Nshift = preCSDOpts(x).Nshift;
    locoDiamDeform_opts{x}.lags = locoDiamDeform_opts{x}.shiftFrame/locoDiamDeform_opts{x}.frameRate;
    locoDiamDeform_opts{x}.xVar = 'Time';
    % Massage the data
    % vasc diamater data
    %vesselROIpool = [vesselROI{x}{~cellfun(@isempty, vesselROI{x})}];
    bigVesselInd = find(strcmpi({vesselROI{x}.vesselType}, {'A'}) | strcmpi({vesselROI{x}.vesselType}, {'D'}));
    diamPool = [vesselROI{x}(bigVesselInd).diameter];
    allDiam = cat(1, diamPool.um_lp)';
    allDiamZ = zscore(allDiam, [], 1);
    diamData = BinDownMean( allDiamZ, locoDiamDeform_opts{x}.binSize ); % allDiamZ; % allDiamZ
    % locomotion data
    velocityData = BinDownMean( vertcat(loco{x}(expt{x}.preRuns).Vdown), locoDiamDeform_opts{x}.binSize );
    accelData = BinDownMean( abs(vertcat(loco{x}(expt{x}.preRuns).Adown)), locoDiamDeform_opts{x}.binSize ); 
    stateData = BinDownMean( vertcat(loco{x}(expt{x}.preRuns).stateDown), locoDiamDeform_opts{x}.binSize );
    % deformation data
    zUse = [projParam.z{:}];
    transMag = vertcat(deform{x}(expt{x}.preRuns).transMag);
    transData = prctile(transMag(:,zUse),80,2);
    transData(normalize(transData) > 5) = NaN;
    transData = BinDownMean( transData, locoDiamDeform_opts{x}.binSize );
    scaleMag = vertcat(deform{x}(expt{x}.preRuns).scaleMag);
    scaleData = prctile(scaleMag(:,zUse),80,2);
    scaleData(normalize(scaleData) > 5) = NaN;
    scaleData = BinDownMean( scaleData, locoDiamDeform_opts{x}.binSize );
    shearMag = vertcat(deform{x}(expt{x}.preRuns).shearMag);
    shearData = prctile(shearMag(:,zUse),80,2);
    shearData(normalize(shearData) > 5) = NaN;
    shearData = BinDownMean( shearData, locoDiamDeform_opts{x}.binSize );
    %{
    transSpdMag = vertcat(deform{x}(expt{x}.preRuns).DtransMag);
    transSpdData = prctile(transSpdMag(:,zUse),80,2);
    transSpdData(normalize(transSpdData) > 5) = NaN;
    transSpdData = BinDownMean( transSpdData, locoDiamDeform_opts{x}.binSize );
    stretchMag = vertcat(deform{x}(expt{x}.preRuns).stretchMag);
    stretchData = prctile(stretchMag(:,zUse),80,2);
    stretchData(normalize(stretchData) > 5) = NaN;
    stretchData = BinDownMean( stretchData, locoDiamDeform_opts{x}.binSize );
    shearRateMag = vertcat(deform{x}(expt{x}.preRuns).DshearMag);
    shearRateData = prctile(shearRateMag(:,zUse),80,2);
    shearRateData(normalize(shearRateData) > 5) = NaN;
    shearRateData = BinDownMean( shearRateData, locoDiamDeform_opts{x}.binSize );
    %}
    tempShift = vertcat(deform{x}(expt{x}.preRuns).shiftZ);
    shiftData = prctile(tempShift(:,zUse),80,2);
    shiftData(abs(normalize(shiftData)) > 5) = NaN;
    shiftData = BinDownMean(shiftData, locoDiamDeform_opts{x}.binSize);

    % DEFINE PREDICTORS
    locoDiamDeform_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); 
    locoDiamDeform_pred{x}.data = [velocityData, diamData]; % accelData, stateData, 
    locoDiamDeform_pred{x}.name = [{'Velocity'}, sprintfc('Diam %i', 1:size(diamData,2))]; %  '|Accel|', 'State'
    locoDiamDeform_pred{x}.N = size(locoDiamDeform_pred{x}.data,2);
    for p = flip(1:locoDiamDeform_pred{x}.N), locoDiamDeform_pred{x}.lopo.name{p} = ['No ',locoDiamDeform_pred{x}.name{p}]; end
    % Set up leave-one-family-out
    firstDiamInd = find(contains(locoDiamDeform_pred{x}.name, 'Diam'), 1);
    locoDiamDeform_pred{x}.fam.col = {1:firstDiamInd-1, firstDiamInd:locoDiamDeform_pred{x}.N}; %{1:4, 5:7}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    locoDiamDeform_pred{x}.fam.N = numel(locoDiamDeform_pred{x}.fam.col); 
    locoDiamDeform_pred{x}.fam.name = {'Loco','Diam'}; 

    % DEFINE RESPONSE
    locoDiamDeform_resp{x}.data = [transData, scaleData, shearData, shiftData]; %diamData;  transSpdData, stretchData, shearRateData, 
    locoDiamDeform_resp{x}.N = size(locoDiamDeform_resp{x}.data, 2); 
    locoDiamDeform_resp{x}.name = {'|Translation|', '|Scale|', '|Shear|', 'Z Shift'}; % , 'TransSpd', 'Stretch', 'Shear Rate'
    % Remove scans with missing data 
    nanFrame = find(any(isnan([locoDiamDeform_pred{x}.data, locoDiamDeform_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    locoDiamDeform_pred{x}.data(nanFrame,:) = []; locoDiamDeform_resp{x}.data(nanFrame,:) = [];

    % RUN THE GLM
    locoDiamDeform_opts{x}.load = false; % true; %   
    locoDiamDeform_opts{x}.saveRoot = expt{x}.dir; %''; %
    [locoDiamDeform_result{x}, locoDiamDeform_summary{x}, ~, locoDiamDeform_pred{x}, locoDiamDeform_resp{x}] = GLMparallel(locoDiamDeform_pred{x}, locoDiamDeform_resp{x}, locoDiamDeform_opts{x}); 
    %locoDiamDeform_summary{x} = SummarizeGLM(locoDiamDeform_result{x}, locoDiamDeform_pred{x}, locoDiamDeform_resp{x}, locoDiamDeform_opts{x});
end
%%
for x = xPresent
    locoDiamDeform_opts{x}.rShow = 2; % 1:locoDiamDeform_resp{x}.N; % 3; %1:locoDiamDeform_resp{x}.N; % 1:LocoDeform_resp{x}.N; %NaN;
    locoDiamDeform_opts{x}.xVar = 'Time';
    ViewGLM(locoDiamDeform_pred{x}, locoDiamDeform_resp{x}, locoDiamDeform_opts{x}, locoDiamDeform_result{x}, locoDiamDeform_summary{x}); %GLMresultFig = 
end

%%
minDev = locoDiamDeform_opts{xPresent(1)}.minDev;
tempSumm = [locoDiamDeform_summary{xPresent}];
lagVec = locoDiamDeform_opts{xPresent(1)}.lags;
Nsumm = numel(tempSumm);
Nresp = locoDiamDeform_resp{xPresent(1)}.N;
% Relative value of different individual predictors (LOPO) and families of predictors (LOFO)
fullDev = nan(Nsumm,Nresp); goodFits = false(Nsumm,Nresp); lopoVal = nan(firstDiamInd-1,Nresp,Nsumm); lofoVal = nan(2,Nresp,Nsumm); goodCoeff = cell(1,Nsumm);
for s = 1:numel(tempSumm)
    fullDev(s,:) = tempSumm(s).dev;
    goodFits(s,:) = tempSumm(s).dev > minDev; % which forms of deformation were well-fit?
    tempVal = 100*(1-tempSumm(s).lopo.devFrac);
    lopoVal(:,goodFits(s,:),s) = tempVal(1:firstDiamInd-1,goodFits(s,:)); % locomotion data only
    
    tempVal = 100*(1-tempSumm(s).lofo.devFrac);
    lofoVal(:,goodFits(s,:),s) = tempVal(:,goodFits(s,:));

    tempCoeff = cat(3, locoDiamDeform_result{xPresent(s)}.coeff);
    goodCoeff{s} = nan(size(tempCoeff));
    goodCoeff{s}(:,:,goodFits(s,:)) = tempCoeff(:,:,goodFits(s,:));
end

figDir = 'D:\MATLAB\Figures\VascDeformGrant\';

% Which locomotion variables contribute most to good fits for differnt forms of deformation?
locoDiam_locoVarsLOPO = figure('Units','inches', 'OuterPosition',[2,2,6,4]);
for resp = 1:Nresp
    subplot(1,Nresp,resp)
    if firstDiamInd > 2
        JitterPlot(squeeze(lopoVal(:,resp,:))', 'new',false, 'bar',true) % squeeze(lopoVal(:,resp,:))'
    else
        JitterPlot(squeeze(lopoVal(:,resp,:)), 'new',false, 'bar',true) % squeeze(lopoVal(:,resp,:))'
    end
    title(locoDiamDeform_resp{xPresent(1)}.name(resp));
    if resp == 1, ylabel('Exp Value'); end
    xlim([0,firstDiamInd])
    set(gca,'Xtick',1:firstDiamInd-1, 'XtickLabel',locoDiamDeform_pred{xPresent(1)}.name(1:firstDiamInd-1))
    xtickangle(30)
    axis square;
    %pause;
end
figPath = sprintf('%slocoDiam_locoVarsLOPO_%s.pdf', figDir, GLMname);
exportgraphics(locoDiam_locoVarsLOPO, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

locoDiam_DevianceAndLOFO = figure('Units','inches', 'OuterPosition',[2,2,6,4]);
sp(1) = subplot(1,2,1);
JitterPlot(fullDev, 'new',false, 'bar',true)
%bar(mean(fullDev, 'omitnan')); hold on;
plot([0,8], minDev*[1,1], 'color','r')
ylabel('Deviance explained');
set(gca,'Xtick',1:Nresp, 'XtickLabel',locoDiamDeform_resp{xPresent(1)}.name)
axis square;

sp(2) = subplot(1,2,2);
bar(mean(lofoVal, 3, 'omitnan')'); 
ylabel('Explanatory value (%)')
set(gca,'Xtick',1:Nresp, 'XtickLabel',locoDiamDeform_resp{xPresent(1)}.name)
linkaxes(sp,'x')
xlim([0,Nresp+1])
axis square;
legend('Loco', 'Diam', 'Location','northeastoutside');
figPath = sprintf('%slocoDiam_DevianceAndLOFO_%s.pdf', figDir, GLMname);
exportgraphics(locoDiam_DevianceAndLOFO, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

locoDiam_velocityCoeffByResp = figure('Units','inches', 'OuterPosition',[2,2,9,6]);
for resp = 1:Nresp
    sp(resp) = subplot(1,Nresp,resp);
    for s = 1:Nsumm
        plot(locoDiamDeform_opts{xPresent(s)}.lags, goodCoeff{s}(:,1,resp)); hold on; % velocity coefficiends
        title(locoDiamDeform_resp{xPresent(1)}.name(resp)) % sprintf()
        if resp == 1, ylabel('Velocity Coefficients'); end
    end
    xlabel('Delay (s)')
    set(gca, 'Xtick',-60:30:60)
    xtickangle(0);
    axis square;
end
figPath = sprintf('%slocoDiam_velocityCoeffByResp_%s.pdf', figDir, GLMname);
exportgraphics(locoDiam_velocityCoeffByResp, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);


locoDiam_diamCoeffByResp = figure('Units','inches', 'OuterPosition',[2,2,9,6]);
for resp = 1:Nresp
    sp(resp) = subplot(1,Nresp,resp);
    for s = 1:Nsumm
        plot(locoDiamDeform_opts{xPresent(s)}.lags, mean(goodCoeff{s}(:,firstDiamInd:end,resp),2) ); hold on; % velocity
        title(locoDiamDeform_resp{xPresent(1)}.name(resp)) % sprintf()
        if resp == 1, ylabel('Mean diameter coefficients'); end
    end
    xlabel('Delay (s)')
    set(gca, 'Xtick',-60:30:60)
    xtickangle(0);
    axis square;
end
figPath = sprintf('%slocoDiam_diamCoeffByResp_%s.pdf', figDir, GLMname);
exportgraphics(locoDiam_diamCoeffByResp, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

% How much does each predictor influence deformation
figure
for s = 1:Nsumm
    cla;
    bar( [1-tempSumm(s).lopo.devFrac(:,2)] )
    set(gca,'Xtick',1:locoDiamDeform_pred{xPresent(s)}.N, 'XtickLabel',locoDiamDeform_pred{xPresent(s)}.name)
    pause; 
end

%% Example of well fit GLM
x = 30;
Tglm = ((0:length(locoDiamDeform_pred{x}.data(:,1))-1)/locoDiamDeform_opts{x}.frameRate)/60;

locoDiam_scaleGLMexample = figure('WindowState','maximized'); % 'Units','inches', 'OuterPosition',[2,2,9,6]
sp(1) = subplot(4,1,1);
plot(Tglm, locoDiamDeform_resp{x}.data(:,2), 'color','c' ); hold on;
plot(Tglm, locoDiamDeform_result{x}(2).prediction, 'color','k' );
plot(Tglm, locoDiamDeform_result{x}(2).lofo.prediction(:,1), 'color','b')
plot(Tglm, locoDiamDeform_result{x}(2).lofo.prediction(:,2), 'color','r')
legend('Data','Full model','Diameter only','Locomotion only')
ylabel('Scaling (um)');
axis tight; 

sp(2) = subplot(4,1,2);
plot(Tglm, locoDiamDeform_pred{x}.data(:,1) ); axis tight;  % ylim([-Inf,Inf])
ylabel('Velocity (cm/s)');

sp(3) = subplot(4,1,3);
plot(Tglm, locoDiamDeform_pred{x}.data(:,2) ); ylim([1,2.1])
set(gca, 'Ytick',[1,2])
ylabel('Loco. state');

[~,maxInd] = max(1-locoDiamDeform_result{x}(2).lopo.devFrac(firstDiamInd:end));
maxInd = maxInd + firstDiamInd - 1;
sp(4) = subplot(4,1,4);
plot(Tglm, locoDiamDeform_pred{x}.data(:,maxInd) );  % , 'Color',[0,0,0,0.2]
linkaxes(sp,'x')
axis tight;
xlabel('Time (min)')
ylabel('Diameter (z-score)');
xlim([0,20]);
figPath = sprintf('%slocoDiam_scaleGLMexample_%s_%s.pdf', figDir, expt{x}.name, GLMname);
exportgraphics(locoDiam_scaleGLMexample, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

%% Vascular coefficients for scaling and shearing
figure;
resp = 2;
for s = find(fullDev(:,resp) > minDev)'
    %plot(locoDiamDeform_opts{xPresent(s)}.lags, goodCoeff{s}(:,4:end,resp)); hold on; 
    plot(locoDiamDeform_opts{xPresent(s)}.lags, mean(goodCoeff{s}(:,4:end,resp),2) ); hold on;
    line(locoDiamDeform_opts{xPresent(s)}.lags([1,end]), [0,0], 'color','k'); hold off;
    axis square;
    pause;
end



