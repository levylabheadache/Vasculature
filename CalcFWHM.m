function [width, point, threshold] = CalcFWHM(vesselProfile_raw, varargin)
% Calc vasc full-width at half-max
IP = inputParser;
addRequired( IP, 'data', @isnumeric )
addParameter( IP, 'threshold', NaN, @isnumeric )
addParameter( IP, 'edgeLims', [-Inf, Inf], @isnumeric )
addParameter( IP, 'smooth', 1, @isnumeric )
addParameter( IP, 'show', false, @islogical )
parse( IP, vesselProfile_raw, varargin{:} ); % mouse, exptDate,
smoothing = IP.Results.smooth;
threshold = IP.Results.threshold;
edgeLims = IP.Results.edgeLims;
show = IP.Results.show;
width = NaN; point = nan(1,2);

if all(~vesselProfile_raw)
    warning('vessel profile is blank, returning NaNs');
else
    vesselProfile_raw = double(vesselProfile_raw(:));     % make sure this is column, and cast to double
    nonzeroInd = find(vesselProfile_raw >0);
    nonzeroLim = nonzeroInd([1,end])';
    
    % smooth data, if appropriate
    if smoothing > 1
        vesselProfile_smooth = medfilt1(vesselProfile_raw, smoothing); % conv2(data, rectwin(smoothing)./ smoothing, 'valid');
        %plot(vesselProfile_raw); hold on; plot(vesselProfile_smooth); 
        nonzeroLim = nonzeroLim + (smoothing/2)*[1,-1];
    else
        vesselProfile_smooth = vesselProfile_raw;
        nonzeroLim = nonzeroLim + [1,-1];
    end
    if isnan(threshold)
        smoothWindowVals = vesselProfile_smooth(nonzeroLim(1):nonzeroLim(2)); % restrict to useful part of the profile
        otsuLevel = multithresh(smoothWindowVals, 1);
        offset = median(smoothWindowVals(smoothWindowVals < otsuLevel)); % find the baseline
        threshold = max(smoothWindowVals - offset)/2 + offset;  % threshold is half max, taking offset into account  max(vesselProfile)/2;
    end

    tempCC = bwconncomp(vesselProfile_smooth > threshold); % all the indices where the data is above half max
    if tempCC.NumObjects > 0
        CCsize = cellfun(@numel, tempCC.PixelIdxList); % require continuity
        %[maxCCsize, maxCCind] = max(CCsize);
        %goodCC = find(CCsize > smoothing);
        if min(CCsize) >= smoothing
            aboveI = vertcat(tempCC.PixelIdxList{CCsize >= smoothing}); %tempCC.PixelIdxList{maxCCind};
        else
            return;
        end
    else
        return; % nothing was above threshold!
    end
    firstI = aboveI(1);                 % index of the first point above threshold
    lastI = aboveI(end);                % index of the last point above threshold

    if (firstI-1 < 1) || (lastI+1) > length(vesselProfile_smooth)
        return % interpolation would result in error
    end

    % use linear interpolation to get a more accurate picture of where the max was
    % find value difference between the point and the threshold value,
    % and scale this by the difference between integer points ...
    pointOffset(1) = (threshold-vesselProfile_smooth(firstI-1)) / (vesselProfile_smooth(firstI)-vesselProfile_smooth(firstI-1));
    pointOffset(2) = (threshold-vesselProfile_smooth(lastI)) / (vesselProfile_smooth(lastI+1)-vesselProfile_smooth(lastI));
    point(1) = firstI- 1 + pointOffset(1);
    point(2) = lastI + pointOffset(2);
    % Prevent going beyond edges defined by max projection
    if point(1) < edgeLims(1), point(1) = edgeLims(1);  end
    if point(2) > edgeLims(2), point(2) = edgeLims(2);  end
    width = point(2)-point(1);

    % plot illustrating the process
    if show
        dataMax = max(vesselProfile_smooth);
        cla
        plot(vesselProfile_raw); hold on;
        plot(vesselProfile_smooth); 
        xlim(nonzeroLim);
        line([1,size(vesselProfile_smooth,1)], offset*[1,1], 'color','k','linestyle','--');
        line([1,size(vesselProfile_smooth,1)], threshold*[1,1], 'color','r','linestyle','--');
        line(point(1)*[1,1], dataMax*[0,1], 'color','k','linestyle','--')
        line(point(2)*[1,1], dataMax*[0,1], 'color','k','linestyle','--')
        if ~isinf(edgeLims(1)), line(edgeLims(1)*[1,1], dataMax*[0,1], 'color','k'); end
        if ~isinf(edgeLims(2)), line(edgeLims(2)*[1,1], dataMax*[0,1], 'color','k'); end
        plot([firstI-1, lastI], vesselProfile_smooth([firstI-1, lastI]), 'o')
        title(sprintf('width = %2.1f pix', width));
        pause(); % 0.2
    end
end
end