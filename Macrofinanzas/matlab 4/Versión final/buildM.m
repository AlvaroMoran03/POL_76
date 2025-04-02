 function [M] = buildM(x, mu, sig)
% Construct the M matrix using three vectors (dM, dD, dU), corresponding to
% the main diagonal, the diagonal below the main one, the diagonal above
% the main one:
% M = [dM(1) dU(1)    0     0  ....      0       0      0
%      dD(2) dM(2) dU(2)    0  ....      0       0      0  
%         0  dD(3) dM(3) dU(3) ....      0       0      0 
%      ....  ....  ....  ....  ....   ....    ....    ....
%         0     0     0     0  .... dD(N-1) dM(N-1) dU(N-1) 
%         0     0     0     0  ....      0    dD(N)   dM(N)]
% =========================================================================
% Input:
% x - grid (N-by-1), equally spaced
% mu - drift term (N-by-1)
% sig - volatility term (N-by-1)

N = length(x); % Grid size
dx = x(2)-x(1); % Grid step
dx2 = dx^2; % Grid step squared

% Constrict the diagonals
dD = -min(mu, 0)/dx + sig.^2/(2*dx2); 
dM = -max(mu, 0)/dx + min(mu, 0)/dx - sig.^2/dx2;
dU = max(mu, 0)/dx + sig.^2/(2*dx2);

% Construct the M matrix
M = spdiags([dD dM dU],[1 0 -1],N,N)';