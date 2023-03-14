% DESCRIPTION - Tests the generation and local matrix computation for the
% 2D VEM element considered in [Hitchhiker's guide to the virtual element
% method]. The computed matrices match the result known in closed form in
% the mentioned paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = [0 0 0;
     3 0 0;
     3 2 0;
     3/2 4 0;
     0 4 0];

P0 = [1 1 0];

E = element2d(P, P0);

figure
set(gcf,'Color', 'white')
fill(P(:,1), P(:,2), 1)
set(gca,'FontSize',18)
axis equal