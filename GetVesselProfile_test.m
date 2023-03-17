function vesselROI = GetVesselProfile(vesselROI, tifPath)
% Opens the tiff file and gets the  vessel projections from the defined polygons

Nroi = numel(vesselROI);

runTif = loadtiff(tifPath{1}); 
[ySize, xSize, Nscan] = size(runTif);
X = repmat(1:xSize, ySize, 1);
Y = repmat((1:ySize)', 1, xSize);
Nruns = size(tifPath, 1);

inpolyScan = cell(1,Nroi);
for roi = 1:Nroi
    inpolyScan{roi} = inpolygon(X, Y, vesselROI(roi).boxPosition.xy(:,1), vesselROI(roi).boxPosition.xy(:,2));
end

for runs = 1:Nruns
    if run > 1, runTif = loadtiff(tifPath{runs}); end
    for roi = 1:Nroi
        for scan = 1:Nscan
            boundedScan = tifPath(:,:,scan).*uint16(inpolyScan{roi}); % imshow(boundedScan, [])
            vesselROI(roi).projection(scan, :) = radon(boundedScan, vesselROI(roi).projectionAngle);
        end
    end
end
end