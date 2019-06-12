function [loopsol,s,A, XP,q,l, B, hinit,mm] = loop_bvp2(R0,Z0,mesh, a0, C0, F0, ZF,bpoint) 
close all;
p=a0;
r=C0;
x = mesh*bpoint;
aF = ZF/Z0;
s=x;
k0=320;
beta = F0*R0^3/k0;
B = [];
mm = bpoint;

m = length(x);
% initSol = init_bvp(x);
loopsol = zeros(12, length(x), length(F0));

[s, Sol_2] = bvp_tube_2(mesh,a0,[], C0, R0, Z0, 0, 0.1, bpoint);
%XP = [V.p(1:28), 1.5199:-0.01:0.3511];

XP = Sol_2(7, 1)*R0:-0.25:2;
%XP = Sol_2(7,1)*R0;
xp = Sol_2(7,1)*R0;

initSol = Sol_2;
A = zeros(1,length(XP)); 

hinit = Sol_2(2,end);

figure
% loop over the a0rng vector
h = waitbar(0,sprintf('Calculating... \\F_0 = %d/%0.3f', 0, 0, max(F0)));

for ii = 1:length(XP)
    
     % update the status bar
    waitbar(ii/length(XP), h, sprintf('Calculating... \\alpha_0 = %0.3f/%0.3f', XP(ii), max(XP)))
    
    
    try
    
    % solve for the iith value of a0rng
    [~,Sol,f] = bvp_tube(mesh,a0,initSol, C0, R0, Z0, F0, ZF, xp,bpoint);
    
    %Solve for energy integral 
    K = (1/R0^2)* Sol(4,:).^2 - (1/R0 * Sol(4,:) - sin(Sol(3,:))./(R0 * Sol(1,:))).^2;

    E = k0/R0^2*Sol(4,:).^2 - 0.53*k0*K;

    for i = 1:1000
        W1(i) = (k0/R0^2*Sol(4,i).^2 - 0.53*k0*(1/R0^2*Sol(4,i).^2 - (1/R0*Sol(4,i) - sin(Sol(3,i))./(R0*Sol(1,i))).^2))*2*pi*R0*Sol(1,i)*Sol(2,end)/1000;
        W2(i) = (k0/R0^2*Sol(10,i).^2 - 0.53*k0*(1/R0^2*Sol(10,i).^2 - (1/R0*Sol(10,i) - sin(Sol(9,i))./(R0*Sol(7,i))).^2))*2*pi*R0*Sol(7,i)*Sol(2,end)/1000;
    end 

    q(ii) = sum(W1) + sum(W2);
    
    % catches errors from endoClathrin
    catch ME
        
        display(ME.message);
        
        loopsol = loopsol(:,:,1:ii-1);
        
        A = A(1: ii-1);
        
        break;  % breaks out of the loop
        
    end
    
    % assign iith solution
    loopsol(:,:,ii) = Sol;
    
    A(ii) = f;
    
%     fun = @(x) f/(2*aF) * ( 1./(1+exp((x-aF)/0.005)));
    fun = @(x) f * ( 1./(1+exp((x-aF)/0.005)));
    l(ii) = integral(fun,0,1);
    
    if ii >1
        if l(ii)<l(ii-1)
            B = 'instability';
        end
    end
    
    F0 = f;
    % set solution as initial guess for next iteration
    initSol = Sol;
    
    temp = [Sol(2,:) Sol(8,:)];
    x = mesh*bpoint;
    ff = [x ((2-bpoint)/bpoint*x + bpoint)];
    %xx = [x ((2-bpoint)/bpoint*x + bpoint)];
    xx = 1:1:2*length(mesh);
    
 
    bp = xx(find(temp>=hinit,1));
    bpoint = ff(bp);
    mm = [mm bpoint];
   
    
    
    
    xp = Sol(7,1)*R0;
    xp = xp - 0.25;
    %XP = Sol(7,1)*R0 - 1;
    
    %bpoint = 0.3 + ii/10;
    
end
close(h)    % close status bar

display(sprintf('Final solution: F0 = %0.3f', F0(end)));

% plot the resultant profile of the membrane
coatArea = [0.5*Z0/R0-a0 0.5*Z0/R0+a0];
fArea = [0.5*Z0/R0-aF 0.5*Z0/R0+aF];
xLim = [0 2*R0];
yLim = [-Z0 Z0];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f, farea = %0.4f, force=%0.4f', a0, C0, aF, XP(end));
figure(6)
FigHandle = figure(6);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
% plot_tube_mean(loopsol(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(7)
FigHandle = figure(7);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_psi(loopsol(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(8)
FigHandle = figure(8);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_lambda(loopsol(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(9)
FigHandle = figure(9);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_K(loopsol(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(10)
FigHandle = figure(10);
set(FigHandle, 'Position', [0, 1000, 400, 400]);
plot_tube_energy(loopsol(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, fArea);

if length(XP) ~= length(A)
    XP = XP(1: ii-1);
end

figure('Position',[0 1000 300 300])
%subplot(2,1,1)
plot(XP, l)
hold on
scatter(XP,l,24)
xlabel('Radius (nm)','interpreter','Tex')
ylabel('Force (pN)','interpreter','Tex')
title('Force vs Position','interpreter','Tex')
set(gca,'fontsize',16)
axis auto

%{
subplot(2,1,2)
plot(A,XP)
hold on
scatter(A,XP,24)
xlabel('Force','interpreter','latex','fontsize',16)
ylabel('Radius (nm)','interpreter','latex','fontsize',16)
title('Position vs Force','interpreter','latex','fontsize',16)
set(gca,'fontsize',30)
%}
figure('Position',[0 1000 300 300])
scatter(XP, q/4, 24);
hold on
plot(XP, q/4);
ylabel('Energy (k_{B}T)','interpreter','Tex')
xlabel('Neck Radius','interpreter','Tex')
title('Energy vs Neck Radius','interpreter','Tex')
set(gca,'Fontsize',16)
axis auto

%save_figures(p,r,R0,Z0)
%save_matfile(p,r,loopsol,A, XP,q,l, B, R0,Z0)

