function [loopsol,s] = loop_bvp_cat(R0,Z0,mesh, a0rng, C0) 
x = mesh*Z0/R0;
XP=9;
s=x;
m = length(x);
initSol = init_bvp_cat(x);
loopsol = zeros(6, length(x), length(a0rng));
figure
% loop over the a0rng vector
h = waitbar(0,sprintf('Calculating... \\alpha_0 = %d/%0.3f', 0, 0, max(a0rng)));

for ii = 1:length(a0rng)
    
     % update the status bar
    waitbar(ii/length(a0rng), h, sprintf('Calculating... \\alpha_0 = %0.3f/%0.3f', a0rng(ii), max(a0rng)))
    
    
    try
    
    % solve for the iith value of a0rng
    [~,Sol,p] = bvp_tube_cat(mesh,a0rng(ii),initSol, C0, R0, Z0, 0, 0.01, XP);
    
    % catches errors from endoClathrin
    catch ME
        
        display(ME.message);
        
        loopsol = loopsol(:,:,1:ii-1);
        
        a0rng = a0rng(1:ii-1);
        
        break;  % breaks out of the loop
        
    end
    
    % assign iith solution
    loopsol(:,:,ii) = Sol;
    
    XP = p;
    
    % set solution as initial guess for next iteration
    %initSol = Sol;
    
end
close(h)    % close status bar

display(sprintf('Final solution: a0 = %0.3f', a0rng(end)));

% plot the resultant profile of the membrane
coatArea = [0.5*Z0/R0-a0rng(end) 0.5*Z0/R0+a0rng(end)];
xLim = [0 2*R0];
yLim = [-1.3*Z0 1.3*Z0];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f', a0rng(end), C0);
figure(3)
plot_tube(loopsol(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, 0);
figure(4)
plot_tube_energy(loopsol(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, 0);