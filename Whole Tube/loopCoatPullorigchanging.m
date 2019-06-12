% Input(s):
%   alpha - dimensionless patch area
%   mesh - meshing for the domain, runs from 0 to 1, i.e. 0:0.01:1
%   lambda - membrane tension at the boundary, in units of pN/nm
%   acoat0 - dimensionless coat area
%   rF - radius of the area of applied force (on a flat membrane), in units of nm
%   zpRng - range of pole z-positions to loop over, in units of nm
%   f0 - Initial guess for the applied force, in units of pN
%   k0 - bending rigidity of bare membrane, in units of pN*nm
%   dk - ratio between rigidity of coated membrane and bare membrane, dk = k_coated/k_bare
%   P - pressure difference from outside to inside, in units of pN/nm^2
%   gamma - sharpness of transition from coated to bare membrane, i.e. tanh(gamma*x)
%   C0 - preferred curvature of the coat
%   R0 - nondimensionalization length
%   initSol - initial guess for the solution

% Output(s):
%   FvsZp - array of computed force vs. position
%   coatPullSol - solution array
%   compZpRng - computed range of pole z-positions, in units of nm

function [FvsZp, coatPullSol, compZpRng] = loopCoatPullorigchanging(alpha, mesh, lambda, acoat, rF, zpRng, f0, k0, dk, P, gamma, C0, R0, initSol)

t=alpha*mesh;   % area mesh points

if isempty(initSol)
    initSol = endoInitchanging(alpha, mesh, lambda, k0, R0);        % initial guess
end

FvsZp = zeros(2, length(zpRng));    % initialize FvsZp matrix
coatPullSol = zeros(12, length(mesh), length(zpRng));   % initialize solution matrix
compZpRng = zpRng;

% display a status bar for the calculation
h = waitbar(0,sprintf('Calculating... z_p = %0.1f nm/%0.1f nm', zpRng(1), zpRng(end)));

figure; % open a figure for the intermediate solutions

% loop over the zpRng vector
for ii = 1:length(zpRng)
    
    lastwarn('');
   
    % update the status bar
    waitbar(ii/length(zpRng), h, sprintf('Calculating... z_p = %0.1f nm/%0.1f nm', zpRng(ii), zpRng(end)))
    

    try
    
    % solve for the iith value of zpRng
    [t,Sol,f] = coatPullorigchanging(alpha, mesh, lambda, acoat, rF, zpRng(ii), f0, k0, dk, P, gamma, C0, R0, initSol);
    
    [warnStr, warnID] = lastwarn;
    
    % catches errors from loopCoatPull
    if strcmp(warnID, 'MATLAB:bvp4c:RelTolNotMet')
        
        error(warnStr)
        
    end
    
    catch ME
        
        display(ME.message);
        
        coatPullSol = coatPullSol(:,:,1:ii-1);
        
        FvsZp = FvsZp(:,1:ii-1);
        
        compZpRng = compZpRng(1:ii-1);
        
        break;  % breaks out of the loop
        
    end
    
    % assign iith solution
    coatPullSol(:,:,ii) = Sol;
    
    % assign iith value of FvsZp
    FvsZp(:,ii) = [Sol(2,1)*R0, f];
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
end

close(h)    % close status bar

display(sprintf('Final solution: z_p = %0.0f nm', compZpRng(end)));

% plot the resultant profile of the membrane
coatArea = [0 acoat];
actArea = [0 ((rF/R0)^2)/2];
plotTitle = sprintf('Membrane profile, \\lambda = %0.3f pN/nm,  P = %0.3f pN/nm^2, A_{coat} = %0.3f nm^2', lambda, P, acoat*2*pi*R0^2);
%xLim = [-20 20];

plotMemProfileArea(Sol, t, R0, coatArea, actArea, [], [], [], plotTitle, 0)
