% -------------------------------------------------------------------------
%  quiverC3DNew.m  --  Variant of quiverC3D with explicit nanmax/nanmin and
%                       an embedded uniquetol fallback.
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  NOTE: the original file in this repository defined a function named
%        quiverC3D (matching the *other* helper) which caused a name clash
%        on the MATLAB path. The function is renamed to match the file
%        name (quiverC3DNew) so both helpers can coexist. Behaviour and
%        parameters are otherwise identical to the upstream version.
% -------------------------------------------------------------------------
function quiverC3DNew(x, y, z, u, v, w, varargin)
%QUIVERC3DNEW 3-D quiver plot with magnitude-coloured arrows.
%   See quiverC3D for full documentation; this variant differs only in:
%     * Direct nanmax / nanmin (instead of verLessThan dispatch).
%     * Local fallback uniquetol implementation appended below.

    %% --- Parse arguments -----------------------------------------------
    p = inputParser;
    addRequired(p, 'x', @isnumeric);
    addRequired(p, 'y', @isnumeric);
    addRequired(p, 'z', @isnumeric);
    addRequired(p, 'u', @isnumeric);
    addRequired(p, 'v', @isnumeric);
    addRequired(p, 'w', @isnumeric);
    addOptional(p, 'scale',        1,   @isscalar);
    addOptional(p, 'MaxNumArrows', inf, @isscalar);
    addParameter(p, 'Normalize',   false);
    addParameter(p, 'LineWidth',   0.5, @isscalar);
    addParameter(p, 'MaxHeadSize', 1,   @isscalar);
    addParameter(p, 'Colorbar',    false, @islogical);
    addParameter(p, 'Parent',      gca(), @isgraphics);
    parse(p, x, y, z, u, v, w, varargin{:});

    scale        = p.Results.scale;
    maxnumarrows = p.Results.MaxNumArrows;
    normalize    = p.Results.Normalize;
    lw           = p.Results.LineWidth;
    hs           = p.Results.MaxHeadSize;
    ax           = p.Results.Parent;

    assert(scale >= 0);
    if ~isequal(size(x), size(y), size(z), size(u), size(v), size(w))
        error('X, Y, Z, U, V, W have to be arrays of the same size.');
    end

    %% --- Subsample if there are too many arrows ------------------------
    if numel(u) > maxnumarrows
        if isvector(u)
            N = ceil(numel(u) / maxnumarrows);
        else
            N = ceil(nthroot(numel(u) / maxnumarrows, ndims(u)));
        end
        if isvector(x)
            idx = {1:N:length(x)};
        else
            idx = arrayfun(@(s) 1:N:s, size(x), 'UniformOutput', false);
        end
        x = x(idx{:}); y = y(idx{:}); z = z(idx{:});
        u = u(idx{:}); v = v(idx{:}); w = w(idx{:});
    end

    %% --- Magnitude-based colour mapping --------------------------------
    C       = colormap(ax);
    ncolors = size(C, 1);
    I       = sqrt(u.^2 + v.^2 + w.^2);

    Imax = nanmax(I(:), []);
    Imin = nanmin(I(:), []);   %#ok<NASGU>  % kept for optional colourbar

    Ic = round(I/Imax*ncolors);
    Ic(Ic == 0) = 1;

    if normalize
        u = u./I;  v = v./I;  w = w./I;
    else
        if scale > 0
            [~, dxi] = uniquetol(x);  dx = diff(x(sort(dxi)));
            [~, dyi] = uniquetol(y);  dy = diff(y(sort(dyi)));
            [~, dzi] = uniquetol(z);  dz = diff(z(sort(dzi)));
            dr_mean  = mean(abs([dx; dy; dz]));
            lenscale = scale * dr_mean / Imax;
        else
            lenscale = 1;
        end
    end
    u = u*lenscale;  v = v*lenscale;  w = w*lenscale;

    %% --- Plot ----------------------------------------------------------
    ishold_flag = ishold(ax);
    hold(ax, 'on');

    if numel(u) > ncolors
        for k = 1:ncolors
            mask = (Ic == k);
            quiver3(ax, x(mask), y(mask), z(mask), ...
                u(mask), v(mask), w(mask), 0, ...
                'Color', C(k, :), 'LineWidth', lw, 'MaxHeadSize', hs);
        end
    else
        for k = 1:numel(u)
            quiver3(ax, x(k), y(k), z(k), u(k), v(k), w(k), 0, ...
                'Color', C(Ic(k), :), 'LineWidth', lw, 'MaxHeadSize', hs);
        end
    end

    if ~ishold_flag
        hold(ax, 'off');
    end

    if p.Results.Colorbar
        colorbar(ax);
    end
end


function [z, ii, jj] = uniquetol(x, tol, varargin)
%UNIQUETOL Local fallback (Siyi Deng, 2010) used when the built-in is not
%   available. Tolerance-based unique. With tol omitted or 0 it is
%   identical to UNIQUE.
    if size(x, 1) == 1
        x = x(:);
    end
    if nargin < 2 || isempty(tol) || tol == 0
        [z, ii, jj] = unique(x, varargin{:});
        return;
    end
    [y, ii, jj] = unique(x, varargin{:});
    if size(x, 2) > 1
        [~, ord]  = sort(sum(x.^2, 1), 2, 'descend');
        [y, io]   = sortrows(y, ord);
        [~, jo]   = sort(io);
        ii        = ii(io);
        jj        = jo(jj);
    end
    d     = sum(abs(diff(y, 1, 1)), 2);
    isTol = [true; d > tol];
    z     = y(isTol, :);
    bin   = cumsum(isTol);
    jj    = bin(jj);
    ii    = ii(isTol);
end
