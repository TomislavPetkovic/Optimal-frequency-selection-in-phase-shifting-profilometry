function [k, width, idx] = mps_get_spurious_tuples(lambda, width, limit_deg)
% MPS_GET_SPURIOUS_TUPLES Returns spurious period tuples.
%   [k, width, idx] = MPS_GET_SPURIOUS_TUPLES(lambda, width, limit_deg)
%   returns all spurious tuples k for input period lambdas.
%
%   The input lambda defines wavelengths. All elements of lambda must be
%   whole numbers. If there are M elements in lambda then k is N x M matrix
%   whose each row contains one of M-tuples of spurious periods.
%
%   The input width uses the width as the maximal allowed unwrapped value.
%   If left empty then a theoretical maximum is used.
%
%   The input limit_deg defines expected error in wrapped phase in degrees.
%   If omitted or left empty then the value of 10 degrees is used.
%
%   See also MPS_GET_PERIOD_TUPLES and MPS_GET_WRAPPED_TUPLES.

% $Revision: 1.0 $  $Date: 2022/07/30 $
% $Author(s): Tomislav Petkovic $

narginchk(2, 3);
nargoutchk(0, 3);

% Set default limit to 5 degrees.
if (3 > nargin) || isempty(limit_deg); limit_deg = 90; end
assert( isnumeric(limit_deg) && (1 == numel(limit_deg)) && (0 < limit_deg) && (limit_deg < 360) );

% Set default width.
lambda = reshape(lambda, 1, numel(lambda));
if all( round(lambda) == lambda )
    width_max = mps_lcm(lambda);
    if (2 > nargin) || isempty(width); width = width_max; end
    assert( (0 <= width) && (width <= width_max) );
end
assert( isnumeric(width) && (1 == numel(width)) && (0 < width) );

% Create outputs.
k = zeros(0, numel(lambda), 'int32');
j = 1;
idx = zeros(0, 1, 'int32');

% Get period tuples.
if 2 > nargin; width = []; end
[k_in, width, w_min, w_max] = mps_get_period_tuples(lambda, width);
N = size(k_in, 1);

% First compute lengths of all intervals over which period-order numbers
% are constant. Then we scale the selected limit on wrapped phase error
% from degrees to the common variable x using the largest wavelength as the
% error in the largest wavelength is the worst case. Finally, we find
% indices of period order numbers after which there is an interval shorter
% than the computed limit.
delta = diff([w_min; w_max(end)]);
delta_max = max(lambda) * limit_deg / 360;
s1 = find(delta < delta_max);

% Add all spurious period-order numbers to the output.
if ~isempty(s1)
    
    % Add first-order spurious tuples.
    for i = 1 : numel(s1)
        n = s1(i);
        if (2 < n) && (n < N - 1)
            % The period order tuple after one at n is too close to one
            % at n which means noise may perturb the measurement
            % producing a spurious perior tuple with theoretically
            % impossible period order numbers.            
            k(j,:) = k_in(n-1,:) + k_in(n+1,:) - k_in(n,:);
            idx(j,1) = n;
            j = j + 1;
        end
    end
    
    s2 = find(1 == diff(s1));
    if ~isempty(s2)
        
        % Add second-order spurious tuples.        
        for i = 1 : numel(s2)
            n = s1(s2(i));
            if (2 < n) && (n < N - 2)
                % We have to consecutive close perior-order tuples which
                % means we have to two add additional perturbations to
                % spurious period tuples.
                k(j,:) = k_in(n-1,:) + k_in(n+2,:) - k_in(n+1,:);
                idx(j,1) = n;
                j = j + 1;
                
                k(j,:) = k_in(n-1,:) + k_in(n+2,:) - k_in(n,:);
                idx(j,1) = n;
                j = j + 1;                
            end
        end
        
        s3 = find(1 == diff(s2));
        if any(s3)
            warning('The code to generate spurious period order tuples of order higher than 2 is not implemented!');
        end
    end
end