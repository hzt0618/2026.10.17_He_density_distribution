function out = cross(varargin)
% CROSS Compatibility wrapper.
% - No input args: return legacy metastable-He cross-section profile.
% - With input args: perform vector cross-product behavior compatible with
%   MATLAB cross(A,B[,dim]) without relying on builtin(''cross'',...).

if nargin == 0
    Energy = [19.82,19.9,20.1,20.3,20.5,20.7,20.9,21.1,21.3,21.5, ...
        21.7,21.9,22.1,22.3,22.5,22.7,22.9,30,50,100];    % eV
    ctab = [0.00,2.15,5.09,6.10,5.40,3.85,4.88,5.01,4.73,4.80, ...
        4.92,5.09,5.27,5.23,6.32,5.46,5.53,8.80,4.20,0.51] * 1e-18; % cm^2
    energy = linspace(19.82, 100, 101);
    out = interp1(Energy, ctab, energy, 'pchip', 'extrap');
    return;
end

if nargin < 2 || nargin > 3
    error('cross:InvalidNumInputs', 'cross requires 0, 2, or 3 inputs.');
end

A = varargin{1};
B = varargin{2};
if ~isnumeric(A) || ~isnumeric(B)
    error('cross:InvalidType', 'Inputs must be numeric arrays.');
end
if ~isequal(size(A), size(B))
    error('cross:SizeMismatch', 'A and B must have the same size.');
end

if nargin == 3
    dim = varargin{3};
    if ~isscalar(dim) || dim < 1 || dim ~= floor(dim)
        error('cross:InvalidDim', 'DIM must be a positive integer scalar.');
    end
else
    dim = find(size(A) == 3, 1, 'first');
    if isempty(dim)
        error('cross:NoDimensionOfThree', ...
            'A and B must have a dimension of length 3.');
    end
end

if dim > ndims(A) || size(A, dim) ~= 3
    error('cross:NoDimensionOfThree', ...
        'Selected DIM must index a dimension of length 3.');
end

idx = repmat({':'}, 1, ndims(A));
idx1 = idx; idx2 = idx; idx3 = idx;
idx1{dim} = 1;
idx2{dim} = 2;
idx3{dim} = 3;

A1 = A(idx1{:}); A2 = A(idx2{:}); A3 = A(idx3{:});
B1 = B(idx1{:}); B2 = B(idx2{:}); B3 = B(idx3{:});

C1 = A2 .* B3 - A3 .* B2;
C2 = A3 .* B1 - A1 .* B3;
C3 = A1 .* B2 - A2 .* B1;

out = zeros(size(A), 'like', A + B);
out(idx1{:}) = C1;
out(idx2{:}) = C2;
out(idx3{:}) = C3;
end
