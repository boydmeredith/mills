function inspectAllDatasets(varargin)
f=figure;
p=inputParser;
addParamValue(p,'subjects', {'J114' 'J115' 'J116' 'J117' 'J118' 'J119' 'J122' 'J123'});
parse(p,varargin{:});
nSubjects = length(p.Results.subjects);
subjects = p.Results.subjects;

for ss = 1:nSubjects
    ax = subplot(nSubjects,1,ss,'parent',f);
    try 
        sumFiles=ls(fullfile(jlgDataDir, subjects{ss}, ...
        '*/L01_referenceLocalization/summary.mat'));
    catch
        continue;
    end
    
    sumFiles = strsplit(sumFiles,'\n');
    
    counter = 0;
    for ff = 1:length(sumFiles)
        if isempty(sumFiles{ff}), continue, end
        thisSum = load(sumFiles{ff},'xyzrcoPeak', 'params');
        
        mnBlocks = thisSum.params.mByNBlocks;
        
        x = squeeze(thisSum.xyzrcoPeak(1,:,:));
        y = squeeze(thisSum.xyzrcoPeak(2,:,:));
        z = squeeze(thisSum.xyzrcoPeak(3,:,:));
        r = squeeze(thisSum.xyzrcoPeak(4,:,:));
        c = squeeze(thisSum.xyzrcoPeak(5,:,:));
        o = squeeze(thisSum.xyzrcoPeak(6,:,:));
        
        hold(ax, 'on')

        % Z value over time each day
%         scatter(ax, counter+(1:size(z,2)), mean(z),'.'); set(ax,'ydir','rev','ylim',[1 51])
%         counter = counter+length(z)+1;
%         plot(ax, counter*[1 1], [1 51],':k')
        % C value over time each day
        scatter(ax, counter+(1:size(c,2)), mean(c),'.'); set(ax, 'ylim',[0 1])
        counter = counter+length(z)+1;
        plot(ax, counter*[1 1], [0 1],':k')
        % Y value over day
        %scatter(ax, repmat(ff,length(y(:,1)),1), y(:,1), 10,[1 1 1]*.2, 'filled' );
        
    end
    title(ax,subjects{ss})
    
end
xlabel(ax,'session #')
%ylabel(ax,'y')

