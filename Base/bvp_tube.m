function [s,Sol,f] = bvp_tube(mesh,alpha0,initSol, C0, R0, Z0, F0, ZF,XP, bpoint)
global k0 lam dc c0 a0 g delta p pw wd P iSol gt r0 z0 f0 aF xp aa tt aF_neck
P=0;
r0=R0;
z0=Z0;
xp=XP/r0;
aa = bpoint;
x=mesh*aa;
s = x;
c0 = R0*C0;             % dimensionless preferred curvature
g = 20;
gt = x;
delta = 1;
k0 = 320;
p = P*R0^3/k0;  
pw = 0;%1*R0^3/k0;         % dimensionless pressure from rigid wall
wd = 10/R0; 
a0 = alpha0;
aF = ZF/Z0;
lam = 0*R0^2/320;
dc = 0;
% f0 = F0*r0/k0;
f0 = F0*r0^3/k0;
tt = (2 - aa)/aa;
temp = initSol(8,:)*R0;
sec_x = linspace(aa, 2, 1001);
Z_c = abs(temp(end) - temp(find(x>=aa-aF, 1)));
disp(Z_c)
% Z_c = 1
aF_neck = sec_x(find(temp>=min(temp)+1, 1)) - aa; 
aF_neck = 0.05;
display(aF_neck)
%beta = F0*r0/k0;
beta = F0*r0^3/k0;


if isempty(initSol)
    initSol = init_bvp(x);        % initial guess
end

iSol = initSol;
options = bvpset('NMax', 100*length(x), 'RelTol', 1e-3);
solinit = bvpinit(x,@mat4init, beta);
sol = bvp4c(@rhs_bvp, @bc_bvp, solinit, options);

% extract force
%f = sol.parameters*k0/r0;
f = sol.parameters*k0/r0^3;
%f = sol.parameters*k0/R0^3*pi*rF^2;

Sol = deval(sol,x);

coatArea = [0.5*z0/r0-alpha0 0.5*z0/r0+alpha0];
fArea = [0.5*z0/r0-aF 0.5*z0/r0+aF];
xLim = [0 2*r0];
yLim = [-z0 z0];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f, farea = %0.4f, force = %0.4f', alpha0, C0,aF, f0);

figure(1)
FigHandle = figure(1);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_mean(Sol, coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(2)
FigHandle = figure(2);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_psi(Sol, coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(3)
FigHandle = figure(3);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_lambda(Sol, coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(4)
FigHandle = figure(4);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_K(Sol, coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(5)
FigHandle = figure(5);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_energy(Sol, coatArea, x, R0, plotTitle, xLim, yLim, fArea);

function Yinit = mat4init(x)
 
 global gt iSol
 
 % returns the vector of values for the initial guess for each mesh point
 Yinit = iSol(:,find(gt==x));