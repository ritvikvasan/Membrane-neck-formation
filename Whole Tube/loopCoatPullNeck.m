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

%load coatPullorig_breakatt=50_morepoints,initSol of last 

function [FvsZp, coatPullSol, compZpRng] = loopCoatPullNeck(alpha, mesh, lambda, acoat, rF, xprng, f0, k0, dk, P, gamma, C0, R0, initSol, Fz, zp)

t=alpha*mesh;   % area mesh points

if isempty(initSol)
    initSol = endoInit(alpha, mesh, lambda, k0, R0);        % initial guess
end

FvsZp = zeros(4, length(xprng));    % initialize FvsZp matrix
coatPullSol = zeros(12, length(mesh), length(xprng));   % initialize solution matrix
compZpRng = xprng;
hinit = initSol(2,end);

% display a status bar for the calculation
h = waitbar(0,sprintf('Calculating... z_p = %0.1f nm/%0.1f nm', xprng(1), xprng(end)));

figure; % open a figure for the intermediate solutions

% loop over the zpRng vector
for ii = 1:length(xprng)
    
    lastwarn('');
   
    % update the status bar
    waitbar(ii/length(xprng), h, sprintf('Calculating... z_p = %0.1f nm/%0.1f nm', xprng(ii), xprng(end)))
    

    try
    
    % solve for the iith value of zpRng
    [t,Sol,faxial, fradial, Z_c] = coatPullneck(alpha, mesh, lambda, 0, rF, xprng(ii), f0, k0, dk, P, gamma, C0, R0, initSol, Fz, zp);
    %[t,Sol,faxial, fradial] = coatPullnecktest(alpha, mesh, lambda, 0, rF, xprng(ii), f0, k0, dk, P, gamma, C0, R0, initSol, Fz, zp);

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
    FvsZp(:,ii) = [Sol(1,2001)*R0, faxial, fradial, Z_c];
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
    t = alpha*mesh;
    ttt = [t t*(200-alpha)/alpha + alpha];
    %ttt = [t t*(800-alpha)/alpha + alpha];
    ff = 1:1:4002;
    temp = [Sol(2,:) Sol(8,:)];
    %alpha = ttt(find(temp<=hinit,1));
    al = ff(find(temp<=hinit,1));
    alpha = ttt(al);
    
    
    %temp = t(find(Sol(2,:)<=hinit,1));
   % if isempty(temp)
      %  temp = t(find(Sol(8,:)<=hinit,1));
       % alpha = alpha + temp;
   % else 
     %   alpha = temp;
    %end
 
end

close(h)    % close status bar

display(sprintf('Final solution: z_p = %0.0f nm', compZpRng(end)));

% plot the resultant profile of the membrane
coatArea = [alpha-5 alpha];
actArea = [0 ((rF/R0)^2)/2];
plotTitle = sprintf('Membrane profile, \\lambda = %0.3f pN/nm,  P = %0.3f pN/nm^2, A_{coat} = %0.3f nm^2', lambda, P, acoat*2*pi*R0^2);
%xLim = [-20 20];

plotMemProfileArea(Sol, t, R0, coatArea, actArea, [], [], [], plotTitle, 0)
