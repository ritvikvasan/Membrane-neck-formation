% Inputs:
%   alpha - dimensionless patch area
%   mesh - meshing for the domain, runs from 0 to 1, i.e. 0:0.01:1
%   lambda - membrane tension at the boundary, in units of pN/nm
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   R0 - nondimensionalization length

% Output:
%   initSol - initialized solution array


function initSol = endoInit(alpha, mesh, lambda, k0, R0)

ds=0.0001;    % Small deviation from x=0 (necessary to avoid division by 0)

t = alpha*mesh;     % area mesh points

lam = lambda*R0^2/k0;   % dimensionless membrane tension

initSol = [ds + sqrt(2*t)           % x
           zeros(1, length(t))      % y
           zeros(1, length(t))      % psi
           zeros(1, length(t))      % h
           zeros(1, length(t))      % l  
           lam*ones(1, length(t))
           sqrt(2*(4*t + 40))%sqrt(2*(3*t + 50))           % x
           zeros(1, length(t))      % y
           zeros(1, length(t))      % psi
           zeros(1, length(t))      % h
           zeros(1, length(t))      % l  
           lam*ones(1, length(t))]; % lambda