function vesselROI = EstimateVesselDiameter(vesselROI, umPerPixel, scanRate, varargin) % MScanData, theFrames
IP = inputParser;
addRequired( IP, 'vesselROI', @isstruct )
addRequired( IP, 'umPerPixel', @isnumeric )
addRequired( IP, 'scanRate', @isnumeric )
addParameter( IP, 'motion', false, @islogical )
addParameter( IP, 'smooth', 5, @isnumeric )
addParameter( IP, 'minR2', 0.7, @isnumeric ) %0.8 %SCN
addParameter( IP, 'minDiam', 3, @isnumeric )
addParameter( IP, 'lp', NaN, @isnumeric ) % low-pass cutoff frequency (for Gaussian filtering)
addParameter( IP, 'show', true, @islogical )
parse( IP, vesselROI, umPerPixel, scanRate, varargin{:} ); % mouse, exptDate,
%motionCorrect = IP.Results.motion;
smoothingWidth = IP.Results.smooth;
lpCutoff = IP.Results.lp;
minR2 = IP.Results.minR2;
minDiam = IP.Results.minDiam;
minDiamPix = minDiam/umPerPixel;
show = IP.Results.show;
lpFilt = MakeGaussFilt( scanRate, 10, 0, 1, 'cutoff',lpCutoff, 'show',true );  % scanRate/2

% Calculate diameter using FWHM and get the baseline diameter (Add in a 5 pixel median filter)
Nscan = size(vesselROI(1).projection, 1);
Nroi = numel(vesselROI);
for roi = 1:Nroi
    % Estimate diameter of the max projection, and use those results to restrict fitting for individual scans
    %[vesselROI(roi).diameter.pix_max, vesselROI(roi).diameter.lims] = CalcFWHM( vesselROI(roi).projection_max(1,:), 'smooth',smoothingWidth, 'show',true); % medfilt1(vesselROI(roi).projection_max(1, :), 5)

    [maxFitResult, maxInit, stepLower, stepUpper] = CalculateDiameter( vesselROI(roi).projection_max(1,:), 'range',vesselROI(roi).radonRange, 'smooth',smoothingWidth, 'show',true ); % 'boxRadon',vesselROI(roi).boxRadon, 
    vesselROI(roi).diameter.pix_max = maxFitResult.width;
    stepInit(1) = maxFitResult.onset_amp; % onset amplitude
    stepInit(2) = maxFitResult.onset_rate; % onset rate
    stepInit(3) = maxFitResult.onset_x; % Onset position
    stepInit(4) = maxFitResult.offset_amp; % offset amplitude
    stepInit(5) = maxFitResult.offset_rate; % offset rate
    stepInit(6) = maxFitResult.offset_x; % Offset position
    stepInit(7) = maxFitResult.constant; % Constant term
    %stepInit(8) = maxFitResult.gauss_amp; %
    %stepInit(9) = maxFitResult.gauss_x; %
    %stepInit(10) = maxFitResult.gauss_std; %

    % Estimate the diameter for each scan
    tempFitResult = cell(1,Nscan);
    tic
    parfor scan = 1:Nscan
        tempFitResult{scan} = CalculateDiameter( vesselROI(roi).projection(scan,:), 'range',vesselROI(roi).radonRange, 'smooth',smoothingWidth, 'stepInit',stepInit, 'stepLower',stepLower, 'stepUpper',stepUpper, 'show',true); % pause;
    end
    tempFitResult = [tempFitResult{:}]; toc   
    vesselROI(roi).diameter.pix_raw = [tempFitResult.width];
    vesselROI(roi).diameter.R2 = [tempFitResult.R2];
   
%{
    figure; 
    sp(1) = subplot(2,2,1); plot(vesselROI(roi).diameter.pix_raw); hold on;
    line(Nscan*[0,1], minDiamPix*[1,1], 'color','r', 'lineStyle','--')
    sp(2) = subplot(2,2,3); plot(vesselROI(roi).diameter.R2); hold on;
    line(Nscan*[0,1], minR2*[1,1], 'color','r', 'lineStyle','--')
    linkaxes(sp, 'x')
    xlim([1,Nscan]);
    subplot(2,2,[2,4]); plot(vesselROI(roi).diameter.R2, vesselROI(roi).diameter.pix_raw, '.')
%}
    
    % Suppress badly fit data points
    vesselROI(roi).missing_scan = find(vesselROI(roi).diameter.R2 < minR2 | isnan(vesselROI(roi).diameter.pix_raw) | vesselROI(roi).diameter.pix_raw < minDiamPix);  % find(isnan(vesselROI(roi).diameter.pix_raw));
    vesselROI(roi).Nmissing = numel(vesselROI(roi).missing_scan);
    vesselROI(roi).missing_frac = vesselROI(roi).Nmissing/Nscan;
    vesselROI(roi).diameter.pix_raw(vesselROI(roi).missing_scan) = NaN;
%{    
    figure; 
    subplot(2,2,1); plot(vesselROI(roi).diameter.pix_raw); hold on; %sp(1) = 
    plot(vesselROI(roi).diameter.pix_corr);
    line(Nscan*[0,1], minDiamPix*[1,1], 'color','r', 'lineStyle','--')
%}

    % interpolate bad/missing scans
    vesselROI(roi).diameter.pix_corr = vesselROI(roi).diameter.pix_raw;
    if vesselROI(roi).Nmissing > 0
        goodScan = setdiff(1:Nscan, vesselROI(roi).missing_scan);
        vesselROI(roi).diameter.pix_corr(vesselROI(roi).missing_scan) = interp1( goodScan, vesselROI(roi).diameter.pix_raw(goodScan), vesselROI(roi).missing_scan, 'spline' );
    end

  %{
    figure; 
    sp(1) = subplot(3,1,1); plot(vesselROI(roi).diameter.R2); hold on;
    line([1,Nscan], minR2*[1,1], 'color','r')
    ylabel('R^2');

    sp(2) = subplot(3,1,2); plot(vesselROI(roi).diameter.pix_raw); hold on;
    line([1,Nscan], vesselROI(roi).diameter.pix_max*[1,1], 'color','r')
    line([1,Nscan], minDiamPix*[1,1], 'color','r')
    plot(vesselROI(roi).missing_scan, vesselROI(roi).diameter.pix_raw(vesselROI(roi).missing_scan), 'x')
    ylabel('Diameter (pix)'); title('Raw');
    sp(3) = subplot(3,1,3); plot(vesselROI(roi).diameter.pix_corr); ylabel('Diameter (pix)'); title('Corrected');
    linkaxes(sp, 'x')
    xlim([1,Nscan])
%}

    % Convert pix -> um
    vesselROI(roi).diameter.um_max = umPerPixel*vesselROI(roi).diameter.pix_max;
    vesselROI(roi).diameter.um = umPerPixel*vesselROI(roi).diameter.pix_corr; 
    vesselROI(roi).diameter.um_modal = prctile(vesselROI(roi).diameter.um, 50); %mode(vesselROI(roi).diameter.um);  % (setdiff(1:Nscan, vesselROI(roi).missing_scan))

    % Lowpass filter vessel diameter 
    if ~isnan(lpCutoff)
        vesselROI(roi).diameter.um_lp = filtfilt(lpFilt, 1, vesselROI(roi).diameter.um);
    else
        vesselROI(roi).diameter.um_lp = nan(1, Nscan);
    end

    if show
        figure;
        plot(1:Nscan, vesselROI(roi).diameter.um, 'c', 'linewidth',2); hold on;
        plot(1:Nscan,vesselROI(roi).diameter.um_lp, 'k', 'linewidth',1); hold on;
        line([1,Nscan], vesselROI(roi).diameter.um_max*[1,1], 'color','r')
        line([1,Nscan], vesselROI(roi).diameter.um_modal*[1,1], 'color','b')
        ylabel('Diameter (um)'); xlabel('Scan');
    end
end
end



    %{
    %vesselROI(roi).diameter.lims = vesselROI(roi).diameter.lims + smoothingWidth/2*[-1,1];
    tic
    for scan = flip(1:Nscan)
        %CalculateDiameter( vesselROI(roi).projection(scan,:), 'smooth',smoothingWidth, 'show',true ); % 'edgeLims',vesselROI(roi).diameter.lims, 
        %vesselROI(roi).diameter.pix_raw(scan) = CalcFWHM( vesselROI(roi).projection(scan, :), 'smooth',smoothingWidth, 'show',true ); % 'edgeLims',vesselROI(roi).diameter.lims, 
        tempFitResult(scan) = CalculateDiameter( vesselROI(roi).projection(scan,:), 'smooth',smoothingWidth, 'stepInit',stepInit, 'stepLower',stepLower, 'stepUpper',stepUpper, 'show',false); 
    end
    toc
    %}

    % Remove motion artifacts
    %{
    if motionCorrect
        diffThresh = 2;
        if strcmpi(vesselROI(roi).vesselType, 'D') || strcmpi(vesselROI(roi).vesselType, 'V')
            rawThresh = 0.3;
        else
            rawThresh = 0.5;
        end
        vesselROI(roi).diameter.um = RemoveVesselMotion(vesselROI(roi).diameter.um, vesselROI(roi).diameter.um_modal, diffThresh, rawThresh, false);
    end
    %}
%{
    figure;
    subplot(1,2,1); plot(vesselROI(roi).projection_max(1,:)); 
    subplot(1,2,2); plot(vesselROI(roi).boxRadon )
%}
    %{
    fitMat = [vertcat(tempFitResult.onset_amp), vertcat(tempFitResult.onset_x), vertcat(tempFitResult.onset_rate), vertcat(tempFitResult.offset_amp), ...
        vertcat(tempFitResult.offset_x), vertcat(tempFitResult.offset_rate), vertcat(tempFitResult.constant), vertcat(tempFitResult.R2), vertcat(tempFitResult.width)];
    fitNames = ["X_L", "Amp_L", "Rate_L", "X_R", "Amp_R", "Rate_R","Const","R2", "Width"]
    
    imagesc(fitMat'); impixelinfo;
    imagesc(zscore(fitMat)'); caxis([-4,4]); colormap bluewhitered

    tiledlayout('flow')
    for col = 1:size(fitMat,2)
        sp(col) = nexttile;
        plot(fitMat(:,col)); hold on;
        plot(vesselROI(roi).missing_scan, fitMat(vesselROI(roi).missing_scan,col), '.'); 

        ylabel(fitNames(col), 'Interpreter','none')
    end
    linkaxes(sp, 'x')
    %}
