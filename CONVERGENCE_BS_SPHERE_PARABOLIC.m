% DESCRIPTION - tests the convergence for a parabolic BS problem on the 3D
% sphere using VEM on a polyhedral mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars

L2errors = zeros(1,4);
meshsizes = zeros(1,4);
L2rates = zeros(1,3);

T = 2^(-2); % final time

for it=1:3
    tau = 2^(-2-2*it);
    switch it
        case 1
            meshfile = 'sphere5.mat';
        case 2
            meshfile = 'sphere9.mat';
        case 3
            meshfile = 'sphere17.mat';
        case 4
            meshfile = 'sphere33.mat';
    end
    [L2errors(it), meshsizes(it)] = parabolic3d_bs_sphere(meshfile, T, tau, false);
end

for it=1:3
   L2rates(it) = log(L2errors(it+1)/L2errors(it)) / log(meshsizes(it+1)/meshsizes(it));
end