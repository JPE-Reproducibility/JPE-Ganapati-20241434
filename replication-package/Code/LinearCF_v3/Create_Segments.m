function [vector] = Create_Segments(tau_1,tau_2,n_segments)
%CREATE_SEGMENTS Create a vector of segments
%   tau_1 - nxn
%   tau_2 - nxn
%   n_segments, integer >= 1

% Code for testing
% n_segments= 3
% tau_1 = zeros(4,4)
% tau_2 = ones(4,4)
% [tau_1 tau_2]

    steps   = permute(((1:(n_segments+1))-1)/n_segments,[3,1,2]);
    A       = repmat(steps,[size(tau_1) size(steps,2)]);
    B       = repmat((tau_2-tau_1),1,1,n_segments+1);
    vector  = tau_1+A.*B;

end

