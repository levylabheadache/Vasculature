function [MScanData] = FWHM_MovieProjection(vesselROI, umPerPixel) % MScanData, theFrames

% Calculate diameter using FWHM and get the baseline diameter (Add in a 5 pixel median filter)
for roi = 1:numel(vesselROI)
    for scan = 1:size(vesselROI(roi).projection, 1)
        vesselROI(roi).diameter.pix(scan) = CalcFWHM( medfilt1(vesselROI(roi).projection(scan, :), 5) );
    end
    vesselROI(roi).diameter.um = umPerPixel*vesselROI(roi).rawVesselDiameter; % xFactor
    [holdHist, d] = hist(vesselROI(roi).tempVesselDiameter, 0:.25:100);
    [~, maxD] = max(holdHist);
    vesselROI(roi).diameter.um_modal = d(maxD);
end

%{
for f = min(theFrames):max(theFrames)
    MScanData.data.rawVesselDiameter(f) = CalcFWHM( medfilt1(MScanData.notes.vesselROI.projection(f, :), 5) );
end
MScanData.data.tempVesselDiameter = MScanData.data.rawVesselDiameter*MScanData.notes.micronsPerPixel; % MScanData.notes.xFactor
[holdHist, d] = hist(MScanData.data.tempVesselDiameter, 0:.25:100);
[~, maxD] = max(holdHist);
MScanData.notes.vesselROI.modalFixedDiameter = d(maxD);
%}

end