function h = blocksByTimePlots(subj,movieDate,location)
if isempty(location)
    location = 'L01';
end
theRefLocDir = referenceLocalizationDir(subj, movieDate, location);
s = matfile(fullfile(theRefLocDir, 'summary.mat'));
par = s.params;
mBlocks = par.mByNBlocks(1);
nBlocks = par.mByNBlocks(2);
xyzrcoPeak = s.xyzrcoPeak;
h = figure; 
[hh] = makeSubplots(h,2,3,0.15,0.05,[0.005 0.05 .99 .95]);

movieLength = size(xyzrcoPeak,3);
titles = 'xyzrco';
for dd = 1:size(xyzrcoPeak,1)
    ax = hh(dd);
    imagesc(permute(xyzrcoPeak(dd,:,:),[2 3 1]),'parent',ax);
    ylabel(ax, titles(dd)); 
    colormap(ax,colormapRedBlue);
    colorbar(ax);
    if dd == 3,
        caxis(ax,[par.whichSlices(1) par.whichSlices(end)]);
    end
    if dd < size(xyzrcoPeak,1)-1
        set(ax,'XTick',[])
    end
end
xlabel(ax,'frame #')
set(h, 'position', [560    23   745   925]);
linkaxes(hh)