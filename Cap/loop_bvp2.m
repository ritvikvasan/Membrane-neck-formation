function [loopsol,s,A,l,q] = loop_bvp2(R0,Z0,mesh, a0, C0, F0, ZF) 
%V = importdata('test.mat');
% writerObj = VideoWriter(['R0=10,Z0=10,a0=0,c0=0,zp=logspace(0,-2,250)'],'MPEG-4');
% writerObj.FrameRate = 10;
% open(writerObj);
x = mesh*Z0/R0;
aF = ZF/Z0;
s=x;
k0=320;
beta = F0*R0^3/k0;

m = length(x);
initSol = init_bvp_cat(x);
loopsol = zeros(6, length(x), length(F0));

%[s, Sol_2] = bvp_tube_2(mesh,a0,initSol, C0, R0, Z0, 0, 0.1);
%XP = [V.p(1:28), 1.5199:-0.01:0.3511];

%XP = R0*logspace(0,-2,100);
XP = R0:-0.5:0.5;

A = zeros(1,length(XP)); 

figure
% loop over the a0rng vector
h = waitbar(0,sprintf('Calculating... \\F_0 = %d/%0.3f', 0, 0, max(F0)));

for ii = 1:length(XP)
    
     % update the status bar
    waitbar(ii/length(XP), h, sprintf('Calculating... \\alpha_0 = %0.3f/%0.3f', XP(ii), max(XP)))
    
    
    try
    
    % solve for the iith value of a0rng
    [~,Sol,f] = bvp_tube_cat(mesh,a0,initSol, C0, R0, Z0, F0, ZF, XP(ii));
    
    %Solve for energy integral 
    K = (1/R0^2)* Sol(4,:).^2 - (1/R0 * Sol(4,:) - sin(Sol(3,:))./(R0 * Sol(1,:))).^2;

    E = k0/R0^2*Sol(4,:).^2 - 0.53*k0*K;

    for i = 1:1001
        W(i) = (k0/R0^2*Sol(4,i).^2 - 0.53*k0*(1/R0^2*Sol(4,i).^2 - (1/R0*Sol(4,i) - sin(Sol(3,i))./(R0*Sol(1,i))).^2))*2*pi*R0*Sol(1,i)*Z0/1000;
    end 

    q(ii) = sum(W);
    
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
    
    fun = @(x) f/(2*aF) * ( 1./(1+exp((x-aF)/0.005)));
    l(ii) = integral(fun,0,1);
    
    F0 = f;
    % set solution as initial guess for next iteration
    %initSol = Sol;
    
end
close(h)    % close status bar
%close(writerObj);

display(sprintf('Final solution: F0 = %0.3f', F0(end)));

% plot the resultant profile of the membrane
coatArea = [0.5*Z0/R0-a0 0.5*Z0/R0+a0];
fArea = [0.5*Z0/R0-aF 0.5*Z0/R0+aF];
xLim = [0 2*R0];
yLim = [-2*Z0 Z0*2];
plotTitle = sprintf('Membrane profile, \\alpha_0 = %0.3f, C_0 = %0.4f, farea = %0.4f, force=%0.4f', a0, C0, aF, XP(end));
figure(3)
plot_tube(loopsol(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, fArea);
figure(4)
plot_tube_energy(loopsol(:,:,end), coatArea, x, R0, plotTitle, xLim, yLim, fArea);

if length(XP) ~= length(A)
    XP = XP(1: ii-1);
end

figure
%subplot(2,1,1)
plot(XP, A)
hold on
scatter(XP,A,24)
xlabel('Radius (nm)','interpreter','latex','fontsize',16)
ylabel('Force','interpreter','latex','fontsize',16)
title('Force vs Position','interpreter','latex','fontsize',16)
set(gca,'fontsize',30)

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
figure
scatter(XP, q, 24);
hold on
plot(XP, q);
ylabel('Energy integral','interpreter','latex','fontsize',16)
xlabel('Neck Radius','interpreter','latex','fontsize',16)
title('Energy integral vs Neck Radius','interpreter','latex','fontsize',16)
set(gca,'fontsize',30)
