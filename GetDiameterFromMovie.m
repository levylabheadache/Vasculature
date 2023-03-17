function [MScanData] = GetDiameterFromMovie(MScanData, fileID)
% Opens the tiff file and gets the  vessel projections from the defined polygons
MScanData.notes.firstFrame = imread(fileID, 'TIFF', 'Index', 1);
%fftFirstFrame = fft2(double(MScanData.notes.firstFrame));
X = repmat(1:MScanData.notes.xSize, MScanData.notes.ySize, 1);
Y = repmat((1:MScanData.notes.ySize)', 1, MScanData.notes.xSize);
MScanData.notes.vesselROI.projectionAngle = atand(diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 1))/diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 2)));
MScanData.notes.pixelShift = zeros(4, MScanData.notes.endframe-MScanData.notes.startframe+1);
for frame = MScanData.notes.startframe:MScanData.notes.endframe
    rawFrame = imread(fileID, 'TIFF', 'Index', frame);
    %fftRawFrame = fft2(double(rawFrame));
    %[MScanData.notes.pixelShift(:, theFrame), ~] = dftregistration(fftFirstFrame, fftRawFrame, 1); %DftRegistration(fftFirstFrame, fftRawFrame, 1);
    inpolyFrame = inpolygon(X + MScanData.notes.pixelShift(3, frame), Y + MScanData.notes.pixelShift(4, frame), MScanData.notes.vesselROI.boxPosition.xy(:, 1), MScanData.notes.vesselROI.boxPosition.xy(:, 2));
    boundedrawFrame = rawFrame.*uint16(inpolyFrame);
    MScanData.notes.vesselROI.projection(frame, :) = radon(boundedrawFrame, MScanData.notes.vesselROI.projectionAngle);
end

end