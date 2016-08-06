function h = blocksByTimeLinePlots(subj,movieDate,location,varargin)

p = inputParser;
addOptional(p,'meanCenter',false);
addOptional(p,'plotMean',true);
addOptional(p,'nBlocksPerPlot',6);
addOptional(p,'variablesToPlot',1:6);
addOptional(p,'doShow','on');
addOptional(p,'doSave',false);
parse(p,varargin{:});


if isempty(location)
    location = 'L01';
end
theRefLocDir = referenceLocalizationDir(subj, movieDate, location);
s = load(fullfile(theRefLocDir, 'summary.mat'),'xyzrcoPeak','params');
par = s.params;
mBlocks = par.mByNBlocks(1);
nBlocks = par.mByNBlocks(2);
xyzrcoPeak = s.xyzrcoPeak;
h = figure('visible',p.Results.doShow);
%[hh] = makeSubplots(h,2,3,0.15,0.1,[0.005 0.05 .99 .95]);
nPlots = mBlocks*nBlocks/p.Results.nBlocksPerPlot;

[hh, hM] = makeSubplots(h, nPlots, length(p.Results.variablesToPlot), 0.2, .3, [.01 .02 .92 .98]);
hC = makeSubplots(h,1,1,0,0,[.94 .3 .02 .3]);
movieLength = size(xyzrcoPeak,3);
cmap = colormapRedBlue;
colorSet = cmap(round(linspace(1,length(colormapRedBlue),p.Results.nBlocksPerPlot)),:);
titles = 'xyzrco';


for dd = p.Results.variablesToPlot
    
    
    
    
    %     hold(ax,'all')
    
    if p.Results.meanCenter && dd~=6
        valsToPlot =  bsxfun(@minus, squeeze(xyzrcoPeak(dd,:,:))', mean(xyzrcoPeak(dd,:,:),3));
    else
        valsToPlot = squeeze(xyzrcoPeak(dd,:,:))';
    end
    %plot(valsToPlot,'parent',ax);
    %ylabel(ax, titles(dd));
    
    for pp=1:nPlots
        ax = hM(find(p.Results.variablesToPlot==dd),pp);        set(ax,'ColorOrder',colorSet,'color',[1 1 1]*.7);
        grid(ax,'on')
        hold(ax,'all');
        whichBlocksToPlot = (pp-1)*p.Results.nBlocksPerPlot+[1:p.Results.nBlocksPerPlot];
        plot(ax,valsToPlot(:,whichBlocksToPlot),'linewidth',1)
        ylim(ax,[min(valsToPlot(:)) max(valsToPlot(:))]);
        
        if pp==1
            ylabel(ax,titles(dd))
        end
        if dd==p.Results.variablesToPlot(1)
            title(ax,sprintf('blocks %02d-%02d',whichBlocksToPlot(1), whichBlocksToPlot(end)),'fontsize',15);
        elseif dd==p.Results.variablesToPlot(end)
            xlabel(ax,'frame #')
        end
        
        
    end
    if p.Results.plotMean
        plot(ax,squeeze(mean(xyzrcoPeak(dd,:,:),2)),'k','linewidth',2);
    end
    if dd==3
        set(ax,'ydir','rev')
    end
    
    %
    %     if dd == 3,
    %         caxis(ax,[par.whichSlices(1) par.whichSlices(end)]);
    %     end
    %     if dd < size(xyzrcoPeak,1)-1
    %         set(ax,'XTick',[])
    %     end
end

set(h, 'position',[1 1 1220 900], 'paperpositionmode','manual',...
    'paperunits','inches','paperposition',[0 0 11 8.5],'papersize',[ 11 8.5]);

if p.Results.doSave,
    print(h, fullfile(theRefLocDir, 'blocksByTimePlots.pdf'),'-dpdf','-opengl');
end

imagesc([1:p.Results.nBlocksPerPlot]','parent',hC(1)); colormap(hC(1),colorSet)


