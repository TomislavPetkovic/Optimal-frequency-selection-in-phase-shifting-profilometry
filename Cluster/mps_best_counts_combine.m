function [out_min, out_max] = mps_best_counts_combine(in_min, in_max, varargin)
% MPS_BEST_COUNTS_COMBINE Combines results of multiple searches.
%   [out_min, out_max] = MPS_BEST_COUNTS_COMBINE(in_min, in_max, ...)
%   combines result of multiple exhaustive searches into one result.
%
%   See also MPS_BEST_COUNTS_MIN_MAX.

% $Revision: 1.0 $  $Date: 2022/04/13 $
% $Author(s): Tomislav Petkovic $

% Get parameters.
N = nargin;
M = numel(varargin);

assert( N == 2 + M );
assert( mod(M, 2) == 0 );

% Start with first two results.
out_min = in_min;
out_max = in_max;
assert( out_min{1}(1) <= out_max{1}(1) );

% Process all remaining inputs.
for i = 1 : M
    distance_i = varargin{i}{1}(1);
    
    if out_min{1}(1) > distance_i
        out_min = varargin{i};
    end
    
    if out_max{1}(1) < distance_i
        out_max = varargin{i};
    end
end

assert( out_min{1}(1) <= out_max{1}(1) );