function [out_min, out_max] = mps_best_counts_min_max(counts_in, q)
% MPS_BEST_COUNTS_MIN_MAX Searches for constellation with extremal half-distances.
%   [out_min, out_max] = MPS_BEST_COUNTS_MIN_MAX(counts_in) preforms a
%   sequential search of all constellations stored in counts_in and returns
%   two extremal constellations, one having a minimal and another having a
%   maxmal half-distance between constellation points.
%
%   See also MPS_BEST_COUNTS_COMBINE.

% $Revision: 1.1 $  $Date: 2022/04/20 $
% $Author(s): Tomislav Petkovic $

assert( ismatrix(counts_in) );

% Get dimensions.
N = size(counts_in, 1);
K = size(counts_in, 2);
assert( 2 <= N );
assert( 1 <= K );

% Preallocate intermediate values.
lambda_max = [];
width_max = [];
counts_max = [];
distance_max = -Inf;

lambda_min = [];
width_min = [];
counts_min = [];
distance_min = Inf;

% Fetch data queue and define pooling step.
have_queue = 1 < nargin;
if have_queue
    step = floor(K/10);
    step_min = 50000;
    if step > step_min; step = step_min; end
    next_i = step;
    last_i = K - next_i;
end

% Process all inputs.
for i = 1 : K
    
    counts1 = double( counts_in(:, i) );
    
    [lambda, W1] = mps_periods_from_fringe_counts(counts1);
    [counts2, W2] = mps_fringe_counts_from_periods(lambda);
    assert( all(counts1 == counts2) && (W1 == W2) );
    
    [O, Xk, Xv] = mps_get_projection_matrix_and_centers(lambda, W1);
    d_min = mps_get_minimal_distance(Xk, Xv);
    
    if d_min < distance_min
        lambda_min = lambda;
        width_min = W1;
        counts_min = counts1;
        distance_min = d_min;
    end
    
    if d_min > distance_max
        lambda_max = lambda;
        width_max = W1;
        counts_max = counts1;
        distance_max = d_min;
    end
    
    if have_queue && (i == next_i)
        send(q, step);
        last_i = i;
        next_i = next_i + step;
    end
end

if have_queue && (0 < K - last_i)
    send(q, K - last_i);
end

% Assign outputs. Always assign the half-distance distance as the first
% element in the cell array as it is used later to combine multiple outputs
% of this function.
out_min = {distance_min, lambda_min, width_min, counts_min};
out_max = {distance_max, lambda_max, width_max, counts_max};