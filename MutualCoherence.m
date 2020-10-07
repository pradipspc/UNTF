function mu = MutualCoherence(F)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mutual coherence of a frame F %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m, N] = size(F);
F = F./sqrt(sum(abs(F).^2));
% F=normc(F);
% Gram matrix
G = F'*F;
% Absolute values of the Gram matrix
absgram = abs(G);

% Increase diagonal elements to exclude them
% from calculations
H = 10*eye(N);
absgram = absgram+H;
% off-diagonal abs entries
muvals = absgram(find(absgram<10));

% maximum abs entry
[mu mp]= max(muvals);

% disp(['MUTUAL COHERENCE OF ',num2str(m), 'X', num2str(N), ' MATRIX IS = ', num2str(mu), '(', num2str(mp), ')']);

