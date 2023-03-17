function [fitResult, stepInit, stepLower, stepUpper] = CalculateDiameter(vesselProfile_raw, varargin)
IP = inputParser;
addRequired( IP, 'data', @isnumeric )
%addParameter(IP, 'boxRadon',[], @isnumeric)
addParameter(IP, 'range', [], @isnumeric)
addParameter(IP, 'stepInit', [], @isnumeric)
addParameter(IP, 'stepUpper', [], @isnumeric)
addParameter(IP, 'stepLower', [], @isnumeric)
addParameter( IP, 'threshold', NaN, @isnumeric )
addParameter( IP, 'smooth', 1, @isnumeric )
addParameter( IP, 'maxTau', 1, @isnumeric )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'pause', 0, @isnumeric )
parse( IP, vesselProfile_raw, varargin{:} ); 
fitRange = IP.Results.range;
%boxRadon = IP.Results.boxRadon;
smoothing = IP.Results.smooth;
threshold = IP.Results.threshold;
maxTau = IP.Results.maxTau;
stepInit = IP.Results.stepInit;
stepLower = IP.Results.stepLower;
stepUpper = IP.Results.stepUpper;
show = IP.Results.show;
pauseTime = IP.Results.pause;

fitResult = struct('R2',NaN, 'onset_amp',NaN, 'onset_x',NaN, 'onset_rate',NaN, 'offset_amp',NaN, 'offset_x',NaN, 'offset_rate',NaN, 'constant',NaN, 'prediction',[], 'width',NaN); 
vesselProfile_raw = double(vesselProfile_raw(:));     % make sure this is column, and cast to double
if isempty(fitRange) 
    fitRange = find(vesselProfile_raw >0); 
else
    fitRange = fitRange(1):fitRange(end); 
end
if ~isempty(fitRange)
    % smooth data, if appropriate
    if smoothing > 1
        vesselProfile_smooth = medfilt1(vesselProfile_raw, smoothing, 'omitnan', 'truncate'); % conv2(data, rectwin(smoothing)./ smoothing, 'valid');
        %plot(vesselProfile_raw); hold on; plot(vesselProfile_smooth);
    else
        vesselProfile_smooth = vesselProfile_raw;
    end

    % Data to be modeled
    xData = fitRange(:); % 
    yData = vesselProfile_smooth(fitRange); % restrict to useful part of the profile.    ./boxRadon(fitRange)
    if numel(unique(yData)) == 1
        fprintf('\nSkipped fitting - constant data'); 
        return; 
    end
    % figure; plot(xData, yData)
    modelStep = @(b,x)( b(1)./(1+exp(-(x-b(3))/b(2))) + b(4)./(1+exp((x-b(6))/b(5))) + b(7)); %
    %modelStepGauss = @(b,x)( b(1)./(1+exp(-(x-b(3))/b(2))) + b(4)./(1+exp((x-b(6))/b(5))) + b(7) +  b(8)*exp( -(x-b(9)).^2/(2*b(10)^2)  ) );
    %modelGauss = @(b,x)( b(1)*exp( -(x-b(2)).^2/(2*b(3)^2)  ) );   
    % b = [1,1,-5,1,1,5,0, -2, 0, 0.01]
    % Determine initial parameters
    if isempty(stepInit)
        % use thresholding to roughly separate background from foreground
        if isnan(threshold)
            otsuLevel = multithresh(yData, 1);
            offset = median(yData(yData < otsuLevel)); % find the baseline
            baseStd = std(yData(yData < otsuLevel));
            threshold = max(yData - offset)/2 + offset;  % threshold is half max, taking offset into account  max(vesselProfile)/2;
        end

        % Set up sigmoidal step function fit
        medDiff = median(yData(yData > otsuLevel)) - offset;
        stepInit(1) = medDiff; % onset amplitude
        stepInit(2) = 2; % onset time constant
        stepInit(3) = xData(find(yData > threshold, 1)); % Onset time
        stepInit(4) = medDiff; % offset amplitude
        stepInit(5) = 2; % offset time constant
        stepInit(6) = xData(find(yData > threshold, 1, 'last')); % Offset time
        stepInit(7) = offset; % Constant term
        %stepInit(8) = 0; % amplitude
        %stepInit(9) = median(xData); %0; % center
        %stepInit(10) = 1; % std dev
    end

    % Set upper and lower bounds on the parameters: important to allow for possibility of long delay in onset time
    quarterRange = (max(xData)-min(xData))/4;
    if isempty(stepLower)
        stepLower(1) = 0; % Onset amplitude
        stepLower(2) = 0; % onset rate constant
        stepLower(3) = xData(2); % Onset x position
        stepLower(4) = 0; % offset amplitude
        stepLower(5) = 0;
        stepLower(6) = median(xData);
        stepLower(7) = 0; % constant
        %stepLower(8) = -medDiff; % amplitude
        %stepLower(9) = min(xData)+quarterRange; % center
        %stepLower(10) = 0; % std dev
    end
    if isempty(stepUpper)
        %maxTau = xData(end)-xData(1);
        stepUpper(1) = 2*medDiff; % Onset amplitude
        stepUpper(2) = maxTau; % onset rate constant
        stepUpper(3) = median(xData); % Onset x position
        stepUpper(4) = 2*medDiff; % offset amplitude
        stepUpper(5) = maxTau; % offset rate constant
        stepUpper(6) = xData(end-1); % offset x position
        stepUpper(7) = threshold; % constant term % offset-2*baseStd
        %stepUpper(8) = 0; % amplitude
        %stepUpper(9) = max(xData)-quarterRange; % center
        %stepUpper(10) = quarterRange; % std dev
    end
    
    if ~all(isnan(yData))
        lsqfitOpts = struct('Display','off', 'Algorithm','trust-region-reflective'); % fitOpts = statset('Display','off', 'MaxIter',1000); 
        [stepFit, resnorm] = lsqcurvefit(modelStep, stepInit, xData, yData, stepLower, stepUpper, lsqfitOpts); %     modelStepGauss
        fitResult.R2 = 1-resnorm/sum((yData - mean(yData)).^2); %R2; %fiber(f).(waveType).traceFit{k,roi}.Rsquared.Adjusted; % goodness of fit
        fitResult.onset_amp = stepFit(1);
        fitResult.onset_x = stepFit(3); 
        fitResult.onset_rate = stepFit(2);
        fitResult.offset_amp = stepFit(4);
        fitResult.offset_x = stepFit(6);
        fitResult.offset_rate = stepFit(5); 
        fitResult.constant = stepFit(7);
        fitResult.prediction = modelStep(stepFit, xData); % modelStepGauss(stepFit, xData); %
        fitResult.width = fitResult.offset_x - fitResult.onset_x;
        %fitResult.gauss_x = stepFit(9);
        %fitResult.gauss_std = stepFit(10);
        %fitResult.gauss_amp = stepFit(8);


    else
        fprintf('\nSkipped fitting - no data available');
    end

    % plot illustrating the process
    if show
        dataMax = max(vesselProfile_raw(xData));
        dataMin = min(yData);
        cla
        h(1) = plot(xData, vesselProfile_raw(xData)  ); hold on; % ./boxRadon(fitRange)
        h(2) = plot(xData, yData);
        h(3) = plot(xData, fitResult.prediction);
        % Isolate the step function and gaussian components
        %stepIso = [stepFit(1:7), zeros(1,3)];
        %gaussIso = [zeros(1,7), stepFit(8:10)];
        %plot(xData, modelStepGauss(stepIso, xData), '-.'); hold on;
        %plot(xData, modelStepGauss(gaussIso, xData), '--');
        xlim([xData(1), xData(end)]);
        ylim([dataMin, dataMax])
        legend(h, {'Raw','Smoothed','Fit'}, 'AutoUpdate','off')
        line(fitResult.onset_x*[1,1], dataMax*[0,1], 'color','k','linestyle','--')
        line(fitResult.offset_x*[1,1], dataMax*[0,1], 'color','k','linestyle','--')
        title(sprintf('width = %2.1f pix (R^2 = %2.3f. Rates = [%2.1f, %2.1f] pix)', fitResult.width, fitResult.R2, fitResult.onset_rate, fitResult.offset_rate));
        pause(pauseTime); % 0.2pause; %
    end
else
    fprintf('Blank input data - returning NaNs')
end
%{
onsetModel =  @(b,x)(b(1)./(1+exp(-(x-b(3))/b(2))));
offsetModel =  @(b,x)(b(1)./(1+exp((x-b(3))/b(2))));
testX = 1:100;
%plot(testX, onsetModel([1000,2,testX(10)], testX) + offsetModel([1000,50,530], testX))
plot(testX, onsetModel([100,1,testX(50)], testX))
onsetModel([1000,2,500], xData)
%}
end