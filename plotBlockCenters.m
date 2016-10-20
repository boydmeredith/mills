function h = plotBlockCenters(subj,theDate,location,varargin)

p = inputParser;
addParamValue(p, 'ax', []);
parse(p, varargin{:});

h = p.Results.ax;
if isempty(h)
    f = figure;
    h = axes('parent',f);
end
cla(h);


summaryFname = [referenceLocalizationDir(subj, theDate, location) '/summary.mat'];

summary = load(summaryFname, 'xyzrcoPeak','params');

edgeBlocks = getEdgeBlocks(summary.params.mByNBlocks);

plot(h,reshape(summary.xyzrcoPeak(1,edgeBlocks,:),[],1),...
    reshape(summary.xyzrcoPeak(2,edgeBlocks,:),[],1),'r.');
hold(h,'on')
plot(h,reshape(summary.xyzrcoPeak(1,~edgeBlocks,:),[],1),...
    reshape(summary.xyzrcoPeak(2,~edgeBlocks,:),[],1),'b.');


xlim(h,[-50 550])
ylim(h,[-50 550])
axis(h, 'image');