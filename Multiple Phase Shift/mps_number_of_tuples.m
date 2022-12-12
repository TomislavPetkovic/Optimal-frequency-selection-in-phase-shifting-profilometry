function [N, width] = mps_number_of_tuples(lambda, width)
% MPS_NUMBER_OF_TUPLES Returns the total number of period tuples.
%   N = MPS_NUMBER_OF_TUPLES(lambda, width) returns the total number of
%   period tuples M for input periods lambda. All elements of lambda must
%   be whole numbers. If width is omitted or empty then it is set to the
%   least common multiple of all elements of lambda.
%
%   [N, width] = MPS_NUMBER_OF_TUPLES(lambda) also returns maximal width.
%
%   The number of tuples N is computed using inclusion-exclusion principle.
%
%   See also MPS_LCM, MPS_GET_PERIOD_TUPLES.

% $Revision: 1.0 $  $Date: 2016/03/30 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(0, 2);

% Sanity check.
assert( isnumeric(lambda) && (1 < numel(lambda)) && all(lambda(:) == round(lambda(:))) );

% Set default width.
width_max = mps_lcm(lambda);
if (2 > nargin) || isempty(width); width = width_max; end
assert( isnumeric(width) && (1 == numel(width)) && (round(width) == width) );

if width > width_max
    error(['Invalid input! Supplied width must be a whole number in [0,' num2str(width_max) '] interval']);
end
assert( (0 <= width) && (width <= width_max) );

% Compute the number of tuples using inclusion-exclusion formula.
n = numel(lambda);
N = 0;
for i = 1 : n
    p = lcm_of_tuples(lambda, i);
    N = N - (-1)^i * sum( floor(width ./ p) );
end

assert(0 < N);

end



function t = lcm_of_tuples(x, N)
% LCM_OF_TUPLES Generate the least common multiple of all tuples.
%   t = LCM_OF_TUPLES(x, N) generates combinations of elements of x of the
%   order N and computes the least common multiple for all tuples, e.g. if
%   N = 2 then t contains least common multiples of all pairs [x(i) x(j)]
%   where indices i and j differ, if N = 3 then least common multiples are
%   computed for all triples [x(i) x(j) x(k)] where indices differ etc.

% First generate all tuples as combinations over x.
C = nchoosek(reshape(x, 1, numel(x)), N);
n = size(C, 1);
m = size(C, 2);

% Then compute LCM of tuples. Note that if all elements of x are relatively
% prime then LCM may be replaced by a product, i.e. t is simply prod(C, 2).
if 1 < m
    t = zeros(n, 1);
    for i = 1 : n; t(i) = mps_lcm(C(i, :)); end
else
    t = C;
end

end