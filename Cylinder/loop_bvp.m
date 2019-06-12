function [loopsol_2,s] = loop_bvp(R0,Z0,mesh, a0rng, C0, F0, ZF) 
x = mesh*Z0/R0;
aF = ZF/Z0;
s=x;
k0=400;
f0 = F0*R0^3/k0;
m = length(x);
initSol = init_bvp(x);
loopsol_2 = zeros(6, length(x), length(a0rng));
figure
% loop over the a0rng vector
h = waitbar(0,sprintf('Calculating... \\alpha_0 = %d/%0.3f', 0, 0, max(a0rng)));

for ii = 1:length(a0rng)
    
     % update the status bar
    waitbar(ii/length(a0rng), h, sprintf('Calculating... \\alpha_0 = %0.3f/%0.3f', a0rng(ii), max(a0rng)))
    
    
    try
    
    % solve for the iith value of a0rng
    [~,Sol] = bvp_tube_2(mesh,a0rng(ii),initSol, C0, R0, Z0, F0, ZF);
    
    % catches errors from endoClathrin
    catch ME
        
        display(ME.message);
        
        loopsol_2 = loopsol_2(:,:,1:ii-1);
        
        a0rng = a0rng(1:ii-1);
        
        break;  % breaks out of the loop
        
    end
    
    % assign iith solution
    loopsol_2(:,:,ii) = Sol;
    
    % set solution as initial guess for next iteration
    initSol = Sol;
    
end
close(h)    % close status bar

display(sprintf('Final solution: a0 = %0.3f', a0rng(end)));

% plot the resultant profile of the membrane
coatArea = [0.5*Z0/R0-a0rng(end) 0.5*Z0/R0+a0rng(end)];
fArea = [0.5*Z0/R0-aF 0.5*Z0/R0+aF];
xLim = [0 2*R0];
yLim = [-1.3*Z0 1.3*Z0];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f, farea = %0.4f, force=%0.4f', a0rng(end), C0, aF, f0);
figure(1)
plot_tube(loopsol_2(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(2)
plot_tube_energy(loopsol_2(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, fArea);