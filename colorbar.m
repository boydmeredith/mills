function cb=colorbar(varargin)
%COLORBAR Display color bar (color scale)
%   COLORBAR appends a colorbar to the current axes in the default (right)
%   location
%
%   COLORBAR('peer',AX) creates a colorbar associated with axes AX
%   instead of the current axes.
%
%   COLORBAR(...,LOCATION) appends a colorbar in the specified location
%   relative to the axes.  LOCATION may be any one of the following strings:
%       'North'              inside plot box near top
%       'South'              inside bottom
%       'East'               inside right
%       'West'               inside left
%       'NorthOutside'       outside plot box near top
%       'SouthOutside'       outside bottom
%       'EastOutside'        outside right
%       'WestOutside'        outside left
%
%   COLORBAR(...,P/V Pairs) specifies additional property name/value pairs
%   for colorbar
%
%   COLORBAR('off'), COLORBAR('hide'), and COLORBAR('delete') delete all
%   colorbars associated with the current axes.
%
%   COLORBAR(H,'off'), COLORBAR(H,'hide'), and COLORBAR(H,'delete') delete
%   the colorbar specified by H.
%
%   H = COLORBAR(...) returns a handle to the colorbar axes, which is a
%   child of the current figure. If a colorbar exists, a new one is still \
%   created.
%
%   See also COLORMAP.

%   Unsupported APIs for internal use:
%   LOCATION strings can be abbreviated N, SO, etc or lower case.
%   COLORBAR(...,'horiz') is the same as COLORBAR(...,'SouthOutside')
%   COLORBAR(...,'vert') is the same as COLORBAR(...,'EastOutside')
%   COLORBAR(...,'off') and COLORBAR(...,'delete') deletes the colorbar, if any
%   COLORBAR(AX) puts an R13 colorbar into axes AX.

%   Copyright 1984-2015 The MathWorks, Inc.
%   $Revision: 1.1.4.17

args = varargin;

% Continue warning that the v6 form will go away in the future.
if (nargin >0 && ischar(args{1}) && strcmp(args{1},'v6'))
    warning(message('MATLAB:colorbar:DeprecatedV6Argument'));
    args(1) = [];
end

% check for colorbar(cbar, ...)
cbarin = [];
if ~isempty(args) && ~ischar(args{1}) && ishghandle(args{1}) && isa(args{1}, 'matlab.graphics.illustration.ColorBar')
    cbarin = args{1};    
    args(1) = [];
    if isempty(args) || ~ismember(args{1},{'off','hide','delete'})
        error(message('MATLAB:colorbar:InvalidHandleInput'));
    end
end

do_explicit_delete = false;
peeraxes = [];
location = 'eastoutside';
locationInArgList = false;
pvpairs = {};

% hard code the enumeration values until we can query the datatype directly
locations = {'North','South','East','West','NorthOutside','SouthOutside','EastOutside','WestOutside','manual'};
locationAbbrevs = cell(1,length(locations));
for k=1:length(locations)
    str = locations{k};
    locationAbbrevs{k} = str(str>='A' & str<='Z');
end

% Allow peer axes to be frist argument without 'peer' key passed in
% g818842
if ~isempty(args) && all(ishghandle(args{1},'axes'))
    args = ['peer' args];
end

% n <= narg && isempty(pvpairs)
% done when varargin is empty or reach PV pairs
done = isempty(args);
while ~done
    arg = args{1};
    if ischar(arg)
        if strcmpi(arg,'peer')
            if length(args) > 1
                peeraxes = handle(args{2});
                args(1:2) = [];
            end
            if ~ishghandle(peeraxes,'axes')
                error(message('MATLAB:colorbar:InvalidPeer'));
            end
        elseif ismember(lower(arg), {'off','hide','delete'})
            % found 'off', 'hide', or 'delete'
            do_explicit_delete = true;
            if isempty(peeraxes)
                fig = get(0,'CurrentFigure');
                if ~isempty(fig)
                    peeraxes = get(fig,'CurrentAxes');
                end
            end
            args(1) = [];
        elseif strcmp(arg(1),'h')
            % found 'horizontal'
            location = 'southoutside';
            args(1) = [];
        elseif strcmp(arg(1),'v')
            % found 'vertical'
            location = 'eastoutside';
            args(1) = [];
        elseif (any(strcmpi(arg, locations)) || any(strcmpi(arg, locationAbbrevs)))
            % found a Location in long or short form,
            % colorbar(...,LOCATION)
            location = lower(arg);
            locationInArgList = true;
            % look up the long form location string if needed
            abbrev = find(strcmpi(location, locationAbbrevs));
            if ~isempty(abbrev)
                location = lower(locations{abbrev});
            end
            args(1) = [];
        else
            % assume prop-value pairs have started
            pvpairs = args;
            % find last instance of Location in PV pairs and remove all
            % instances
            [loc, pvpairs] = findLocationInPVPairs(pvpairs);
            if ~isempty(loc)
                location = loc;
                locationInArgList = true;
            end
            break
        end
    else
        error(message('MATLAB:colorbar:InvalidInput'))
    end
done = isempty(args);
end

if do_explicit_delete
    if ~isempty(cbarin)
        delete(cbarin);
    else
        if ~isempty(peeraxes)
            deleteMatchingColorbars(peeraxes);
        end
    end
    if nargout > 0
        error(message('MATLAB:colorbar:TooManyOutputs'));
    end
else
    
    % check pv args
    if ~isempty(pvpairs)
        msg = check_pv_args(pvpairs);
        if ~isempty(msg)
            error(msg)
        end
    end
    
    % create an axes and figure if necessary
    if isempty(peeraxes)
        peeraxes = gca;
    end

    % delete any existing colorbars with the same peer and location
    deleteMatchingColorbars(peeraxes, location);
    
    % construct colorbar
    cbar = matlab.graphics.illustration.ColorBar;
    
    % set the peer axes
    cbar.Axes = peeraxes;
    
    % match Axes color for box and ruler
    cbar.Color_I = get(peeraxes.Parent,'DefaultAxesXColor');
    
    % set location, preserve mode if Location not passed in by user.
    if locationInArgList
        cbar.Location = location;
    else
        cbar.Location_I = location;
    end
    
    % set PV pairs
    if ~isempty(pvpairs)
        set(cbar,pvpairs{:});
    end

    % return colorbar handle
    if nargout > 0
        cb = cbar;
    end
end

%----------------------------------------------------%
function deleteMatchingColorbars(peeraxes, varargin)

doMatchLocation = false;
if ~isempty(varargin)
    location = varargin{1};
    doMatchLocation = true;
end

fig = ancestor(peeraxes,'figure');
cbar = findobj(fig,'-class','matlab.graphics.illustration.ColorBar');
for k=1:length(cbar)
    if isequal(cbar(k).Axes, peeraxes) % matching peer
        if doMatchLocation
            if strcmpi(cbar(k).Location, location)
                delete(cbar(k));
            end
        else
            delete(cbar(k));
        end
    end
end

%---------------------------------------------
function [loc, pvpairs] = findLocationInPVPairs(pvpairs)

loc = [];
removeThese = [];
for i = 1:(length(pvpairs)-1)
    if strcmpi(pvpairs{i},'Location')
        loc = pvpairs{i+1};
        removeThese = [removeThese i i+1]; %#ok<AGROW>
    end
end
pvpairs(removeThese) = [];

%---------------------------------------------
function msg = check_pv_args(args)

msg = [];
n=length(args);
if mod(n,2)==0 % must be even
    for i=1:2:n
        metaClass = ?matlab.graphics.illustration.ColorBar;
        propNames = cellfun(@(x) (x.Name), metaClass.Properties, 'UniformOutput', false);
        if ~ischar(args{i})
            msg = message('MATLAB:colorbar:InvalidPropertyName');
        elseif ~any(strcmpi(propNames,args{i}))
            msg = message('MATLAB:colorbar:UnknownProperty', args{i});
        elseif strcmpi(args{i},'Parent')
            if ~isscalar(args{i+1}) || ~ishghandle(args{i+1})
                msg = message('MATLAB:colorbar:ParentNotHandle');
            elseif ~isa(handle(args{i+1}),'matlab.ui.container.CanvasContainer')
                msg = message('MATLAB:colorbar:InvalidParent',get(args{i+1},'Type'));
            end
        end
    end
else
    msg = message('MATLAB:colorbar:InvalidInput');
end
