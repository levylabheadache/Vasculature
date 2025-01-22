function vesselROI = GetVesselProfile(vesselROI, tifStack, projName)
% Opens the tiff file and gets the  vessel projections from the defined polygons
if nargin > 2 
    projName = sprintf('projection_%s', projName);
else
    projName = 'projection';
end
Nroi = numel(vesselROI);
[ySize, xSize, Nscan] = size(tifStack);
X = repmat(1:xSize, ySize, 1);
Y = repmat((1:ySize)', 1, xSize);
%Nruns = size(tifStack, 1)
radonFactor = 0.67;
for roi = 1:Nroi
    inpolyScan = inpolygon(X, Y, vesselROI(roi).boxPosition.xy(:,1), vesselROI(roi).boxPosition.xy(:,2));
    if Nscan == 1 % max projection only
        vesselROI(roi).boxIm = tifStack(:,:,1).*uint16(inpolyScan); %*uint16(inpolyScan)% figure; imshow(vesselROI(roi).boxIm, []); impixelinfo;
        vesselROI(roi).boxRadon = radon(inpolyScan, vesselROI(roi).projectionAngle); % use this to set the limits for subsequent curve fitting     figure; plot( vesselROI(roi).boxRadon )
        % Restrict analysis to range (in radon projection space) where there are a useful number of pixels (ie  >= radonFactor of max)    
        vesselROI(roi).radonRange = find(vesselROI(roi).boxRadon >= radonFactor*max(vesselROI(roi).boxRadon), 1):find(vesselROI(roi).boxRadon >= radonFactor*max(vesselROI(roi).boxRadon), 1, 'last');
    end
    for scan = 1:Nscan
        boxedScan = tifStack(:,:,scan).*uint16(inpolyScan); %*uint16(inpolyScan) % imshow(boxedScan, [])
        boxedRadon = radon(boxedScan, vesselROI(roi).projectionAngle); % note: radon transform integrates (sums) over all pixels, not an average
        boxedRadon([1:vesselROI(roi).radonRange(1)-1, vesselROI(roi).radonRange(end)+1:end]) = NaN; % suppress values outside the useful range
        vesselROI(roi).(projName)(scan,:) = boxedRadon./vesselROI(roi).boxRadon; % convert from sum to average
    end
end
end

%{
% Find the edges that are aligned/orthogonal to the vessel 
vesselVec = diff(vesselROI(roi).vesselLine.position.xy)

edgeVec = nan(4,2);
for v = 1:3
    edgeVec(v,:) = vesselROI(roi).boxPosition.xy(v+1,:) - vesselROI(roi).boxPosition.xy(v,:);
end
edgeVec(4,:) = vesselROI(roi).boxPosition.xy(4,:) - vesselROI(roi).boxPosition.xy(1,:);

edgeAngleDiff = nan(4,1);
for v = 1:4
    edgeAngleDiff(v) = rad2deg(acos(dot(vesselVec, edgeVec(v,:))/(norm(vesselVec)*norm(vesselVec))));
end
[~,vSort] = sort(abs(edgeAngleDiff - 90));
%}

%{
I = double(imread('peppers.png')); % 
X = reshape(I,size(I,1)*size(I,2),3);
coeff = pca(X);
Itransformed = X*coeff;
Ipc1 = reshape(Itransformed(:,1),size(I,1),size(I,2));
Ipc2 = reshape(Itransformed(:,2),size(I,1),size(I,2));
Ipc3 = reshape(Itransformed(:,3),size(I,1),size(I,2));
figure, imshow(Ipc1,[]);
figure, imshow(Ipc2,[]);
figure, imshow(Ipc3,[]);
%}
%{
[y,x] = ind2sub(size(vesselROI(roi).boxIm),  find(vesselROI(roi).boxIm > 0))
imXY = [x,y];
weightVec = double(vesselROI(roi).boxIm(vesselROI(roi).boxIm > 0));
imPC = pca(imXY, 'Weights',weightVec)

imPC(:,1)
pc1angle = rad2deg( cart2pol(  imPC(1,1), imPC(2,1) ) )
pc2angle = rad2deg( cart2pol(  imPC(1,2), imPC(2,2) ) )
%}
