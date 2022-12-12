function lcm = mps_lcm(a, b)
% MPS_LCM Computes the least common multiple of two whole numbers.
%   lcm = MPS_LCM(a, b) returns the least common multiple of two whole
%   numbers a and b.
%
%   lcm = MPS_LCM(x) where x is a vector of whole numbers returns the least
%   common multiple of all numbers in x. Note that all elements in x must
%   be whole numbers.
%
%   See also MPS_GCD, MPS_EXTENDED_GCD.

% $Revision: 1.2 $  $Date: 2022/04/08 $
% $Author(s): Tomislav Petkovic $

narginchk(1, 2);
nargoutchk(0, 1);

% If there are two inputs then compute the least common multiple directly.
% Otherwise, if there is only one input argument with more than two
% elements call MPS_LCM iteratively as the least common multiple may be
% computed recursively.
if 2 == nargin
    
    assert( isnumeric(a) && (1 == numel(a)) );
    assert( isnumeric(b) && (1 == numel(b)) );
    
    % Test if inputs are whole numbers.
    a_is_whole = (round(a) == a);
    b_is_whole = (round(b) == b);
    if (false == a_is_whole); a = round(a); end
    if (false == b_is_whole); b = round(b); end
    if (false == a_is_whole) || (false == b_is_whole); warning('Rounding input(s) to the closest whole number!'); end
    
    % Test for a possibility of integer overflow.
    assert( isa(a, class(b)) );
    if isa(a, 'double')
        limit = pow2(53); % 53 binary digits
        assert( (a < limit) && (b < limit) );
    elseif isa(a, 'single')
        limit = pow2(24); % 24 binary digits
        assert( (a < limit) && (b < limit) );
    elseif isa(a, 'int64')
        limit = int64(9223372036854775807);
    elseif isa(a, 'int32')
        limit = int32(2147483647);
        assert( (a < limit) && (b < limit) );
    elseif isa(a, 'int16')
        limit = int16(32767);
        assert( (a < limit) && (b < limit) );
    elseif isa(a, 'int8')
        limit = int8(127);
        assert( (a < limit) && (b < limit) );
    elseif isa(a, 'uint64')
        limit = uint64(18446744073709551615);
    elseif isa(a, 'uint32')
        limit = uint32(4294967295);
        assert( (a < limit) && (b < limit) );
    elseif isa(a, 'uint16')
        limit = uint16(65535);
        assert( (a < limit) && (b < limit) );
    elseif isa(a, 'uint8')
        limit = uint8(255);
        assert( (a < limit) && (b < limit) )        
    else
        warning('Unsupported data type! Skipping integer overflow tests.');
    end
    
    gcd = mps_gcd(a, b);
    
    if double(a) * double(b) < double(limit)
        lcm = abs(a * b) / gcd;
    else
        gcd = mps_gcd(a, b);
        fac_a = abs(a) / gcd;
        fac_b = abs(b) / gcd;
        assert( (round(fac_a) == fac_a) && (round(fac_b) == fac_b) );
        lcm = fac_a * fac_b * gcd;
        assert(lcm <= limit);
    end
    
else
    
    assert( isnumeric(a) && all(a(:) == round(a(:))) );
    
    if (1 < numel(a))
        lcm = mps_lcm(a(1), a(2));
        for i = 3 : numel(a); lcm = mps_lcm(lcm, a(i)); end
    elseif 1 == numel(a)
        lcm = a(1);
    else
        lcm = [];
    end
    
end