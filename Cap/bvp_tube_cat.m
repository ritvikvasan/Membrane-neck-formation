function [s,Sol,f] = bvp_tube_cat(mesh,alpha0,initSol, C0, R0, Z0, F0, ZF, XP)
global k0 lam dc c0 a0 g delta p pw wd P iSol gt r0 z0 f0 aF xp aF_neck
P=0;
r0=R0;
z0=Z0;
x=mesh*Z0/R0;
s = x;
c0 = R0*C0;             % dimensionless preferred curvature
g = 20;
gt = x;
delta = 1;
k0 = 320;
p = P*R0^3/k0;  
pw = 1*R0^3/k0;         % dimensionless pressure from rigid wall
wd = 10/R0; 
a0 = alpha0;
lam = 1/4;
dc = 0;
f0 = F0*r0^3/k0;
aF_neck = 0.05;
aF = ZF/Z0;
xp = XP/r0;

if isempty(initSol)
    initSol = init_bvp_cat(x);        % initial guess
end

beta = F0*r0^3/k0;

iSol = initSol;
options = bvpset('NMax', 100*length(x), 'RelTol', 1e-3);
solinit = bvpinit(x,@mat4init,beta);
sol = bvp4c(@rhs_bvp, @bc_bvp_cat, solinit, options);

f = sol.parameters*k0/r0^3;

Sol = deval(sol,x);

coatArea = [0.5*z0/r0-alpha0 0.5*z0/r0+alpha0];
xLim = [0 2*r0];
yLim = [-2*z0 2*z0];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f', alpha0, C0);
figure(1)
plot_tube(Sol, coatArea, x, R0, plotTitle, xLim, yLim,0);
figure(2)
plot_tube_energy(Sol, coatArea, x, R0, plotTitle, xLim, yLim, 0);
%drawnow
%F= getframe(1);
%writeVideo(writerObj,F);


function Yinit = mat4init(x)
 
 global gt iSol
 
 % returns the vector of values for the initial guess for each mesh point
 Yinit = iSol(:,find(gt==x));