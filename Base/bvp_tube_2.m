function [s,Sol_2,sol] = bvp_tube_2(mesh,alpha0,initSol, C0, R0, Z0, F0, ZF, bpoint)
global k0 lam dc c0 a0 g delta p pw wd P iSol gt r0 z0 f0 aF aa tt
P=0;
r0=R0;
z0=Z0;
aa = bpoint;
x=mesh*aa;
s = x;
c0 = r0*C0;             % dimensionless preferred curvature
g = 20;
gt = x;
delta = 1;
k0 = 320;
p = P*R0^3/k0;  
pw = 1*R0^3/k0;         % dimensionless pressure from rigid wall
wd = 10/R0; 
a0 = alpha0;
aF = ZF/Z0;
lam = 0*R0^2/320;
dc = 0;
f0 = F0*r0^3/k0;
tt = (2 - aa)/aa;

if isempty(initSol)
    initSol = init_bvp(x);        % initial guess
end

iSol = initSol;
options = bvpset('NMax', 100*length(x), 'RelTol', 1e-3);
solinit = bvpinit(x,@mat4init);
sol = bvp4c(@rhs_bvp_2, @bc_bvp_2, solinit, options);
Sol_2 = deval(sol,x);

coatArea = [0.5*z0/r0-alpha0 0.5*z0/r0+alpha0];
fArea = [0.5*z0/r0-aF 0.5*z0/r0+aF];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f, farea = %0.4f, force = %0.4f', alpha0, C0,aF, f0);
%plot(Sol_2(1,:)*R0, Sol_2(2,:)*R0); xlim([0 2*r0]);ylim([-z0*1.3 z0*1.3]);hold on
xLim = [0 2*r0];
yLim = [-z0 z0];
figure(1)
FigHandle = figure(1);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_mean(Sol_2, coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(2)
FigHandle = figure(2);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_psi(Sol_2, coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(3)
FigHandle = figure(3);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_lambda(Sol_2, coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(4)
FigHandle = figure(4);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_K(Sol_2, coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(5)
FigHandle = figure(5);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_energy(Sol_2, coatArea, x, R0, plotTitle, xLim, yLim, fArea);


function Yinit = mat4init(x)
 
 global gt iSol
 
 % returns the vector of values for the initial guess for each mesh point
 Yinit = iSol(:,find(gt==x));