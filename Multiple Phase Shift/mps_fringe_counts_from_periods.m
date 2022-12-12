function [count, width] = mps_fringe_counts_from_periods(lambda, width)
% MPS_FRINGE_COUNTS_FROM_PERIODS Returns fringe counts from periods.
%   count = MPS_FRINGE_COUNTS_FROM_PERIODS(lambda) returns an array of
%   whole numbers count which holds fringe counts for fringe patterns with
%   periods stored in array lamgda. All numbers stored in lambda must be
%   whole.
%
%   [count, width] = MPS_PERIODS_FROM_FRINGE_COUNTS(lambda) also returns
%   maximal (virtual) width of the projector screen.
%
%   count = MPS_PERIODS_FROM_FRINGE_COUNTS(lambda, width) works for
%   non-integer lambas also.
%
%   See also MPS_PERIOD_FROM_FRINGE_COUNT, MPS_LCM.

% $Revision: 1.2 $  $Date: 2017/06/19 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(0, 2);

% Check inputs.
all_integers = all( round(lambda) == lambda );
have_width = (1 < nargin) && ~isempty(width);

if all_integers
    
    % Check if inputs are relatively prime integers and reduce them if they
    % are not.
    gcm_lambda = mps_gcd(lambda);
    if 1 < gcm_lambda
        warning('Periods are not relatively prime integers.');
        lambda = lambda ./ gcm_lambda;
    end
    
    % For integer lambda the maximal allowed width is the least common
    % multiple of numbers stored in lambda.
    lcm_lambda = mps_lcm(lambda);
    if ~have_width; width = lcm_lambda; end
    
    assert( isnumeric(width) && (1 == numel(width)) && (0 < width) && (width <= lcm_lambda) );
    
else
    
    % For real periods the maximal allowed width is not well defined
    % therefore we warn the user that width should be selected. If the
    % width is not supplied then compute one as least common multiple of
    % integer periods only.
    if (1 == nargin) || ~have_width
        warning('Width should be defined for non-integer periods!');
        
        lcm_lambda = mps_lcm( lambda(round(lambda) == lambda) );
        width = lcm_lambda;
    end
    
    assert( isnumeric(width) && (1 == numel(width)) && (0 < width) );
    
end

% Compute fringe counts.
count = width ./ lambda;

% Reduce fringe counts.
if all_integers
    
    gcd_count = mps_gcd(count);
    if 1 < gcd_count
        warning('Fringe counts are not relatively prime integers!');
        count = count ./ gcm_count;
    end
    
    lcm_count = mps_lcm(count);
    if width > lcm_count
        warning('Reducing width to LCM of fringe counts!');
        width = lcm_count;
    end
    
    assert(lcm_count == lcm_lambda);
    assert( all( round(count) == count ) );
end