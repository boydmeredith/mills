function h = blocksByTimeHeatPlots(subj,movieDate,location,varargin)
p = inputParser;
addOptional(p, 'doShow', 'on');
addOptional(p, 'doSave', false);
parse(p, varargin{:});

doShow = p.Results.doShow;
doSave = p.Results.doSave;

if doShow, doShow = 'on'; else doShow = 'off'; end
theRefLocDir = referenceLocalizationDir(subj, movieDate, location);
s = load(fullfile(theRefLocDir, 'summary.mat'),'xyzrcoPeak','params');
par = s.params;
mBlocks = par.mByNBlocks(1);
nBlocks = par.mByNBlocks(2);
xyzrcoPeak = s.xyzrcoPeak;
h = figure('visible',doShow); 
[hh] = makeSubplots(h,2,3,0.15,0.1,[0.005 0.05 .99 .95]);

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
linkaxes(hh)

set(h, 'position',[53 5 1220 700], 'paperpositionmode','manual',...
        'paperunits','inches','paperposition',[0 0 11 8.5],'papersize',[ 11 8.5]);

if doSave,
    print(h, fullfile(theRefLocDir, 'blocksByTimePlots.pdf'),'-dpdf','-opengl');
end
