function quiverC3D(x, y, z, u, v, w, varargin)
%quiverC3D creates a 3D quiver plot and adds a color coding. The color coding is
%given by the magnitudes of the component vectors. Large values result in colors 
%from the upper end of the used colormap. 
% 
%   INPUT:
%       x - array, x components of initial points
%       y - array, y components of initial points
%       z - array, z components of initial points
%       u - array, x components of arrows
%       v - array, y components of arrows
%       w - array, z components of arrows
%       scale - number, if > 0, automatic scaling is used to avoid overlap
%           of adjacent arrows. If scale = 0, the scaling is disabled.
%       MaxNumArrows - a positive integer (non-integer should work as well)
%           number limiting the maximum number of plotted arrows. Since vectors
%           of length zero aren't plotted and the matrices are sampled
%           uniformly, it's possible that fewer arrows will be visible in the
%           end. If MaxNumArrows is not given its value is set to inf.
% 
%   OUTPUT:
%       none
% 
% WARNING!: This function might not create arrows with the same length as
%   would be obtained using quiver3(x, y, z, u, v, w). I do not know
%   whether the scaling method is the same is in quiver3.
% 
% 
% --------------------------------------------------------------------------------------
% 
%   EXAMPLE:
%       load wind
%       quiverC3D(x,y,z,u,v,w,1,1000);
%       view(-45,20)
%       set(gca, 'Color', 'black');
%       axis tight
%       axis equal
%   
% --------------------------------------------------------------------------------------
% 
%   See also: QUIVER3, LINESPEC, COLORMAP.
% 

%% ***** Parse arguments ***** 

p = inputParser;
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addRequired(p, 'z', @isnumeric);
addRequired(p, 'u', @isnumeric);
addRequired(p, 'v', @isnumeric);
addRequired(p, 'w', @isnumeric);
addOptional(p, 'scale', 1, @isscalar);
addOptional(p, 'MaxNumArrows', inf, @isscalar);
addParameter(p, 'Normalize', false);
addParameter(p, 'LineWidth', 0.5, @isscalar);
addParameter(p, 'MaxHeadSize', 1, @isscalar);
addParameter(p, 'Colorbar', false, @islogical);
addParameter(p, 'Parent', gca(), @isgraphics);

parse(p, x, y, z, u, v, w, varargin{:});

scale = p.Results.scale;
maxnumarrows = p.Results.MaxNumArrows;
normalize = p.Results.Normalize;
lw = p.Results.LineWidth;
hs = p.Results.MaxHeadSize;
ax = p.Results.Parent;

assert(scale >= 0);

if ~isequal(size(x), size(y), size(z), size(u), size(v), size(w))
    error('X, Y, Z, U, V, W have to be arrays of the same size.');
end

%% ***** Initialization ***** 
if numel(u) > maxnumarrows
    if isvector(u)
        N = ceil(numel(u)/maxnumarrows);
    else
        N = ceil(nthroot(numel(u)/maxnumarrows, ndims(u)));
    end
    
    if isvector(x)
        idx = {1:N:length(x)};
    else
        idx = arrayfun(@(x) 1:N:x, size(x), 'UniformOutput', false);
    end
    
    x = x(idx{:});
    y = y(idx{:});
    z = z(idx{:});
    u = u(idx{:});
    v = v(idx{:});
    w = w(idx{:});
end


%% ***** Colormap definition ***** 
C = colormap(ax);
ncolors = size(C, 1);
I = sqrt(u.^2 + v.^2 + w.^2);

% Commented by AnV July 2, 2022 
%if verLessThan('matlab','9.6')
%    Imax = max(I(:), [], 'omitnan');
%    Imin = min(I(:), [], 'omitnan');
%else
    Imax = nanmax(I(:), []);
    Imin = nanmin(I(:), []);
%end



% Assume that the colormap matrix contains the colors in its rows
Ic = round(I/Imax*ncolors);
Ic(Ic == 0) = 1;


if normalize
    u = u./I;
    v = v./I;
    w = w./I;
else
    % Scale to avoid overlaps
    if scale > 0
        [~, dxi] = uniquetol(x);
        dx = diff(x(sort(dxi)));
        [~, dyi] = uniquetol(y);
        dy = diff(y(sort(dyi)));
        [~, dzi] = uniquetol(z);
        dz = diff(z(sort(dzi)));
        
        % Average difference 
        dr_mean = mean(abs([dx; dy; dz]));
        
        lenscale = scale*dr_mean/Imax;
    else
        lenscale = 1;
    end
end

u = u*lenscale;
v = v*lenscale;
w = w*lenscale;

%% ***** Plotting ***** 
ishold_flag = ishold(ax);
hold(ax, 'on');

if numel(u) > ncolors
    % Group all values which belong to the same color value
    for k = 1:ncolors
        mask = (Ic == k);    
        quiver3(ax, x(mask), y(mask), z(mask), u(mask), v(mask), w(mask), 0, ...
            'Color', C(k, :), 'LineWidth', lw, 'MaxHeadSize', hs);
    end
else
    % If there are so few values, it each arrow is plotted individually 
    for k = 1:numel(u)
        quiver3(ax, x(k), y(k), z(k), u(k), v(k), w(k), 0, 'Color', C(Ic(k), :), ...
            'LineWidth', lw, 'MaxHeadSize', hs);
    end
end

if ~ishold_flag
    hold(ax, 'off');
end

if p.Results.Colorbar
    cb = colorbar(ax);
    tlabels = cb.TickLabels;
    tlabels = sprintfc('%1.3g ',linspace(Imin, Imax, length(tlabels)));
    cb.TickLabels = tlabels;
end

MyColorBar = colorbar;
%TickL = get(MyColorBar,'Ticklabels');
%TickL = sprintfc('%1.3g ',((1:length(TickL)) - 1)*Imax/length(TickL));
%set(MyColorBar, 'TickLabels', TickL); 

end 


function [z,ii,jj] = uniquetol(x,tol,varargin)
%UNIQUETOL Unique element within a tolerance.
%   [Y,I,J] = UNIQUETOL(X,TOL) is very similar to UNIQUE, but allows an
%   additional tolerance input, TOL. TOL can be taken as the total absolute
%   difference between similar elements. TOL must be a none negative
%   scalar. If not provided, TOL is assumed to be 0, which makes UNIQUETOL
%   identical to UNIQUE.
%
%   UNIQUETOL(...,'ROWS')
%   UNIQUETOL(...,'FIRST')
%   UNIQUETOL(...,'LAST')
%   These expressions are identical to the UNIQUE counterparts.
%
%   See also UNIQUE.

% Siyi Deng; 03-19-2010; 05-15-2010; 10-29-2010;

if size(x,1) == 1, x = x(:); end
if nargin < 2 || isempty(tol) || tol == 0
    [z,ii,jj] = unique(x,varargin{:});
    return;
end
[y,ii,jj] = unique(x,varargin{:}); 
if size(x,2) > 1
    [~,ord] = sort(sum(x.^2,1),2,'descend');
    [y,io] = sortrows(y,ord);
    [~,jo] = sort(io);
    ii = ii(io);
    jj = jo(jj);
end
d = sum(abs(diff(y,1,1)),2);
isTol = [true;d > tol];
z = y(isTol,:);
bin = cumsum(isTol); % [n,bin] = histc(y,z);
jj = bin(jj);
ii = ii(isTol);

end % UNIQUETOL;

