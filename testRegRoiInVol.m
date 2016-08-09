subject = 'J115';
theDate = '2015-11-03';
location = 'L01';
nS = getNameStruct(subject, theDate, location);
refLocSumm = load(nS.referenceLocalizationFileName);
stackIs = imageSeries([jlgDataDir '/' refLocSumm.stackPath]);
stack = stackIs.images;

roiWThresh = .005;
nRoiColors = 10;

whichClusters = 1:2;
whichBlocks = 2:6:36; %;1:36;


%%
[roiS.day1103, xyzrcCP] = registerRois(subject, '2015-11-03', location, ...
    'whichBlocks',whichBlocks, 'whichClusters', whichClusters);

%%
[roiS.day0925, xyzrcCP] = registerRois(subject, '2015-09-25', location, ...
    'whichBlocks',whichBlocks, 'whichClusters', whichClusters);
roiSday0925 = roiS.day0925
%%
roiDays = fields(roiS);
nDays = length(roiDays);
whichClustersToPlot = 1;

%
figure;

for pp = 1:2
    subplot(1,2,pp);
    hold on
    
    
    zSliceToPlot = 22;
    zSlice = normalizeToZeroOne(double(stack(:,:,zSliceToPlot)));
    [xx, yy] = meshgrid(1:size(stack,2),1:size(stack,1));
    surf(xx,yy,ones(size(xx)).*zSliceToPlot, zSlice);
    colormap bone
    freezeColors;
    
    colormap(hsv(nRoiColors));
    roiColorVals = linspace(0, 1, nRoiColors);
    for dd = 1:nDays
        dayFieldName = roiDays{dd};
        thisRoiS = roiS.(dayFieldName);
        for cc=whichClustersToPlot,
            for bb=whichBlocks,
                nRois = size(thisRoiS(cc,bb).w,3);
                for ww=1:nRois
                    thisRoi = thisRoiS(cc,bb).w(:,:,ww);
                    thisRoi(thisRoi>roiWThresh) = 1;
                    thisRoi(thisRoi~=1) = nan;
                    thisRoiColorNum = thisRoi .* roiColorVals(mod(ww,nRoiColors) + 1);
                    thisRoiColorNum = thisRoi .* roiColorVals(dd);
                    surf(thisRoiS(cc,bb).x,thisRoiS(cc,bb).y,round(thisRoiS(cc,bb).z), thisRoiColorNum);
                end
            end
        end
    end
    
    shading flat
    
    xlabel('x')
    axis vis3d
    alpha(.7)
    set(gca,'ydir','rev','zdir','rev')
end

