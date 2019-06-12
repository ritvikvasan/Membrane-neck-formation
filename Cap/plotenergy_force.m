function [qq] = plotenergy_force(loopsol,Z0, R0,l,q,pos,C0)

k0=320;
mesh =0:0.001:1;
t = mesh;
tt=t;
Sol=loopsol;
a0=0.3;
c0=R0*C0;
g=20;
c = -0.5*c0*(1 - tanh(g*(t - a0)));

fontsize = 14;  
fontsize2 = 6;
lineWidth = 4;
lineWidthf = 3;%lineWidth = 12;
axesWidth = 2;

xLabelOn = 1;
yLabelOn = 1;
xTickLabelOn = [];
yTickLabelOn = [];


coatColor = 'red';
fColor = 'black';

% XP(:) = loopsol(1,1,:);
% 
% fighandle = figure(1);
% hold on
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',fontsize, 'fontweight','bold')
% plot(XP, -l*R0/k0,'Linewidth',lineWidth)
% ylabel('FR/\kappa', 'FontSize',fontsize);
% xlabel('x', 'FontSize',fontsize);
% set(gca,'FontSize',fontsize, 'XMinorTick', 'on', 'YMinorTick', 'on');
% 
% 
% % turns off x tick label
% if xTickLabelOn == 0
%     set(gca, 'XTickLabel', []);
% end
% 
% % turns off y tick label
% if yTickLabelOn == 0
%     set(gca, 'YTickLabel', []);
% end

%% energy

mesh = 0:0.001:1;
alpha = Z0/R0;
t = alpha*mesh;tt=t;

H(:,:) = Sol(4,:,:)/R0;
r(:,:) = Sol(1,:,:)*R0;
y(:,:) = Sol(2,:,:)*R0;
lam(:,:) = Sol(6,:,:)*k0/R0^2;
s(:,:) = Sol(3,:,:);

k2 =-0.53*k0;

% for j = 1:pos
%     for i = 1:1000
%     K(i,j) = H(i, j)^2 - (H(i,j) - sin(s(i,j))/r(i,j))^2;
%     W(i,j) = k0*(H(i,j)- c(i))^2 + k2*K(i,j) ; 
%     %     E(i,j) = W(i,j)*2*pi*r(i,j)*abs((y(i+1,j) - y(i,j))) ;
%     E(i,j) = W(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2 + lam(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2;
%     E2(i,j) = W(i,j)*(tt(i+1) - tt(i))*2*pi*R0^2;
%     end
%     Epos(j) = sum(E(:,j));
%     Epos2(j) = sum(E2(:,j));
% end

for j = 1:pos
  Sol = loopsol(:,:,j);
    for i = 1:1001
        W(i) = (k0/R0^2*(Sol(4,i) - c(i)).^2 - 0.53*k0*(1/R0^2*(Sol(4,i)).^2 - (1/R0*(Sol(4,i)) - sin(Sol(3,i))./(R0*Sol(1,i))).^2))*2*pi*R0*Sol(1,i)*Z0/1000;
    end 
 qq(j) = sum(W);
end


% x(:) = loopsol(1,1,:);
% y = qq*20/(2*pi*20*320);
% st=1;
% x = x(st:end);
% y = y(st:end);
% 
% fontsize = 14;
% 
% fighandle = figure(2);
% hold on
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% set(gca, 'fontsize',fontsize, 'fontweight','bold')
% xlabel('x')
% ylabel('ER/(2\piZ\kappa)')
% h = plot(x,y);
% set(h                          , ...
%   'Color'           , [0 0 0.5]    , ...  % use 0 0 .5, 0 .5 0, .5 0 0, 0 .5 .5,.5 .5 0
%   'LineStyle' , '--' , ... % use --
%   'LineWidth', 3, ... % use 3
%   'Marker', '*',...  % use *, <, o, h, >
%   'MarkerSize', 7, ... % use 4
%   'MarkerIndices',1:1:length(y));
% 
% % general graphics, this will apply to any figure you open
% % (groot is the default figure object).
% set(groot, ...
% 'DefaultFigureColor', 'w', ...
% 'DefaultAxesLineWidth', 1, ...
% 'DefaultAxesXColor', 'k', ...
% 'DefaultAxesYColor', 'k', ...
% 'DefaultAxesFontUnits', 'points', ...
% 'DefaultAxesFontSize', 14, ...
% 'DefaultAxesFontName', 'Helvetica', ...
% 'DefaultLineLineWidth', 1, ...
% 'DefaultTextFontUnits', 'Points', ...
% 'DefaultTextFontSize', 8, ...
% 'DefaultTextFontName', 'Helvetica', ...
% 'DefaultAxesBox', 'off', ...
% 'DefaultAxesTickLength', [0.02 0.025]);
%  
% % set the tickdirs to go out - need this specific order
% set(groot, 'DefaultAxesTickDir', 'out');
% set(groot, 'DefaultAxesTickDirMode', 'manual');