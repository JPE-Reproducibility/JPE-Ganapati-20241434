% compute cosine basis of dimension d at x

function [XX, DX] = cosine_basis(x, dim)

N = length(x);

%   initialize for recursion
XX  = zeros(N, dim);
XX(:, 1) = 1.0;
for ii = 2:dim
    XX(:, ii) = sqrt(2) * cos((ii-1) * pi * x);
end

if nargout >= 2
    DX = zeros(N, dim);
    for ii = 2:dim
        DX(:, ii) = -sqrt(2) * (ii - 1) * pi * sin((ii-1) * pi * x);
    end
end

