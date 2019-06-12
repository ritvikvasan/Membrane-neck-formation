function plotcoatPull_new(pos, loopsol, alpha, act, dash)

coatPullSol = loopsol;

rad(:) = coatPullSol(1,:,pos);
height(:) = coatPullSol(2,:,pos);
mean(:) = coatPullSol(4,:,pos);
lam(:) = coatPullSol(6,:,pos);


    

Sol(:,:) = coatPullSol(:,:,pos);


%memColor = 'black'; 
%memColor = [0.9967    0.7816    0.2007];
memColor = [0.9139    0.7258    0.3063];
coatColor = 'red';
%coatColor = [0    0.4470    0.7410];
actColor = 'blue';

%tt = (2 - alpha)/alpha;
coatArea = [0 act];
%actArea = [0 3.125];
actArea  = [0 0.3];

% set visual properties of the plot
fontsize = 36;  %54;  %70;
lineWidth = 6; %6;  %12;
axesWidth = 5;
R0=1;%R0=20; % 1 for non dimensional
mesh = 0:0.001:1;
t = alpha*mesh;

fighandle = figure(8);
set(fighandle, 'Position', [0, 1000, 300, 300]);
% plot the membrane shape
if dash == 1
    
    plot(Sol(1,:)*R0, Sol(2,:)*R0, -Sol(1,:)*R0, Sol(2,:)*R0,':', 'Color', memColor, 'LineWidth', lineWidth);
    
else
    
    plot(Sol(1,:)*R0, Sol(2,:)*R0, -Sol(1,:)*R0, Sol(2,:)*R0, ':','Color', memColor, 'LineWidth', lineWidth);

end

hold on






% plots the regions of applied force
if ~isempty(actArea)
    % loops over multiple regions
    for ii = 1:size(actArea,1)

        aMin = actArea(ii,1);
        aMax = actArea(ii,2);
        
        if dash == 1
            
            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ':', 'Color', actColor, 'LineWidth', lineWidth);
            
        else

            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0,  ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ':','Color', actColor, 'LineWidth', lineWidth);
        
        end

    end
end
hold on

% plots the isotropic-curvature-generating coat
if ~isempty(coatArea)
    % loops over multiple regions
    for ii = 1:size(coatArea,1)

        aMin = coatArea(ii,1);
        aMax = coatArea(ii,2);
        
        if dash == 1
            
            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ':', 'Color', coatColor, 'LineWidth', lineWidth);
            
        else

            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ':','Color', coatColor, 'LineWidth', lineWidth);
            
        end

    end
end
% xlabel('Radius (nm)')
% ylabel('Height (nm)')
xlabel('x')
ylabel('y')
set(gca, 'fontsize',17, 'Fontweight','bold')
daspect([ 1 1 1])

% fighandle = figure(9);
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% pointsize = 10;
% %plot(Sol(1,:)*R0, Sol(2,:)*R0,'Color', memColor, 'LineWidth', lineWidth);
% 
% scatter(Sol(1,:)*R0, Sol(2,:)*R0, pointsize, Sol(4,:));
% hold on
% %scatter(Sol(7,:)*R0, Sol(8,:)*R0, pointsize, Sol(10,:));%plot(Sol(1,:)*R0, Sol(2,:)*R0,'Color', memColor, 'LineWidth', lineWidth);
% hold on
% colorbar
% h = colorbar;
% ylabel('Z (nm)', 'FontSize',fontsize, 'FontName', 'Helvetica');
% set(get(h,'title'),'string','h');
% xlabel('R (nm)', 'FontSize',fontsize, 'FontName', 'Helvetica');
% %annotation('textbox', [.22, .325, .10, .075], 'String', str, 'FontName', 'Arial', 'FontSize', 20, 'Linestyle', 'none')
% %annotation('arrow', [.24 0.18], [.5 .5], 'LineWidth', 3, 'HeadWidth', 25)
% % sets plot aesthetic properties
% set(gca,'FontSize',fontsize, 'FontName', 'Helvetica', 'XMinorTick', 'on', 'YMinorTick', 'on', 'Linewidth', axesWidth);
% 
% fighandle = figure(10);
% set(fighandle, 'Position', [0, 1000, 300, 300]);
% pointsize = 10;
% %plot(Sol(1,:)*R0, Sol(2,:)*R0,'Color', memColor, 'LineWidth', lineWidth);
% 
% scatter(Sol(1,:)*R0, Sol(2,:)*R0, pointsize, Sol(6,:));
% hold on
% %scatter(Sol(7,:)*R0, Sol(8,:)*R0, pointsize, Sol(12,:));%plot(Sol(1,:)*R0, Sol(2,:)*R0,'Color', memColor, 'LineWidth', lineWidth);
% hold on
% colorbar
% h = colorbar;
% ylabel('Z (nm)', 'FontSize',fontsize, 'FontName', 'Helvetica');
% set(get(h,'title'),'string','h');
% xlabel('R (nm)', 'FontSize',fontsize, 'FontName', 'Helvetica');
% %annotation('textbox', [.22, .325, .10, .075], 'String', str, 'FontName', 'Arial', 'FontSize', 20, 'Linestyle', 'none')
% %annotation('arrow', [.24 0.18], [.5 .5], 'LineWidth', 3, 'HeadWidth', 25)
% % sets plot aesthetic properties
% set(gca,'FontSize',fontsize, 'FontName', 'Helvetica', 'XMinorTick', 'on', 'YMinorTick', 'on', 'Linewidth', axesWidth);

end

