function plot_tube_K(Sol, coatArea, x, R0, plotTitle, xLim, yLim, fArea)

global k0

K = (1/R0^2)* Sol(4,:).^2 - (1/R0 * Sol(4,:) - sin(Sol(3,:))./(R0 * Sol(1,:))).^2;

E = k0/R0^2*Sol(4,:).^2 - 0.53*k0*K;

fontsize = 40;  
fontsize2 = 30;
lineWidth = 2;
lineWidthf = 3;%lineWidth = 12;
axesWidth = 5;

xLabelOn = 1;
yLabelOn = 1;
xTickLabelOn = 1;
yTickLabelOn = 1;

str = {'Increasing force','from 0 to 7000 (steps of 1000) over aF=0.001'};

%memColor = 'black'; 
%memColor = [0.9967    0.7816    0.2007];
%memColor = [0.9139    0.7258    0.3063];
coatColor = 'red';%coatColor = 'blue';
fColor = 'black';
%coatColor = [0    0.4470    0.7410];

pointsize = 10;
scatter(Sol(1,:)*R0, -Sol(2,:)*R0, pointsize, K);%plot(Sol(1,:)*R0, Sol(2,:)*R0,'Color', memColor, 'LineWidth', lineWidth);
hold on
scatter(Sol(1,:)*R0, Sol(2,:)*R0, pointsize, K);
colorbar
h = colorbar;
ylabel('Z (nm)', 'FontSize',fontsize, 'FontName', 'Helvetica');
set(get(h,'title'),'string','Gaussian Curvature (/nm^2)');
xlabel('R (nm)', 'FontSize',fontsize, 'FontName', 'Helvetica');
%annotation('textbox', [.22, .325, .10, .075], 'String', str, 'FontName', 'Arial', 'FontSize', 20, 'Linestyle', 'none')
%annotation('arrow', [.24 0.18], [.5 .5], 'LineWidth', 3, 'HeadWidth', 25)
% sets plot aesthetic properties
set(gca,'FontSize',fontsize, 'FontName', 'Helvetica', 'XMinorTick', 'on', 'YMinorTick', 'on', 'Linewidth', axesWidth);

% turns off x tick label
if xTickLabelOn == 0
    set(gca, 'XTickLabel', []);
end

% turns off y tick label
if yTickLabelOn == 0
    set(gca, 'YTickLabel', []);
end

% sets x limits of plot
if ~isempty(xLim)
    xlim(xLim);
end

% sets y limits of plot
if ~isempty(yLim)
    ylim(yLim);
end
hold on
axis square