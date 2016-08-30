function plotDriftAllSessions(subject,varargin)

p=inputParser;
addParameter(p,'meanCenter',false);
addParameter(p,'plotVar',3);
addParameter(p,'ignoreEdges',false);
addParameter(p,'whichBlocks',[]);
parse(p,varargin{:});

if isempty(p.Results.whichBlocks)
    whichBlocks = 1:36;
else
    whichBlocks = p.Results.whichBlocks;
end

plotVarNames = 'xyzrco';

sumFiles=ls(fullfile(jlgDataDir, subject, ...
    '*/L01_referenceLocalization/summary.mat'));
sumFiles = strsplit(sumFiles,'\n');
sumFiles(strcmp(sumFiles,''))=[];
nSummaries = length(sumFiles);
fig = figure;
nCols = 7;
nRows = ceil(nSummaries/nCols);
hM = makeSubplots(fig, nCols, nRows, .2 , .4, [0.03 .03 .95 .95]);
set(fig,'position',[1 1 1500 1000]);
yLabel = [plotVarNames(p.Results.plotVar)];
if p.Results.meanCenter, yLabel = [yLabel ' (mean-centered)']; end;

interiorIx = setdiff(1:36,[1:6 7 12 13 18 19  24 25 30 31:36]);

nBlocks=36;
if p.Results.ignoreEdges
    nBlocks=length(interiorIx);
end
cmap = colormapRedBlue;

colorSet = cmap(round(linspace(1,length(colormapRedBlue),nBlocks)),:);

plotVal = [];
dateVal = [];

for ss = 1:nSummaries
    
    

    thisSumFile = sumFiles{ss};
    summ = load(thisSumFile,'xyzrcoPeak','params');
    ax = hM(ss);
    
    valToPlot = squeeze(summ.xyzrcoPeak(p.Results.plotVar,whichBlocks,:));
    if p.Results.ignoreEdges
        valToPlot = valToPlot(interiorIx,:);
    end
    
    if p.Results.meanCenter
        valToPlot =  bsxfun(@minus, valToPlot, mean(valToPlot,2));
    end
    set(ax,'ColorOrder',colorSet,'color',[1 1 1]*.7,'xlim',[1 size(valToPlot,2)]);
    hold(ax,'all');
    plot(ax, valToPlot','linewidth',1);
    
    title(ax,[subject ' ' summ.params.movieDate],'fontsize',10);
    if mod(ss,nCols)==1, ylabel(ax,yLabel,'fontsize',12); end
    if ss>nCols*(nRows-1) , xlabel(ax,'frame #'); end;
    if ~p.Results.meanCenter
        switch p.Results.plotVar
            case {1 , 2},
                set(ax,'ylim',[1 512])
            case 3,
                set(ax,'ydir','rev','ylim',[1 51])
            case 5,
                set(ax,'ylim',[0 1]);
        end
    end
    
    dateVal = [dateVal; datenum(strrep(strrep(sumFiles{ss},[fullfile(jlgDataDir, subject) '/'],''),'/L01_referenceLocalization/summary.mat',''))];
    plotVal = [plotVal; valToPlot(1,1)];
end

linkaxes(hM(1:nSummaries));


figure;plot(dateVal,plotVal,'.-')
