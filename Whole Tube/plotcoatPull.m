function plotcoatPull(pos, loopsol, alpha, act, dash)

coatPullSol = loopsol;

rad(:) = [coatPullSol(1,:,pos) coatPullSol(7,:,pos)];
height(:) = [coatPullSol(2,:,pos) coatPullSol(8,:,pos)];

Sol(:,:) = coatPullSol(:,:,pos);
%memColor = 'black'; 
%memColor = [0.9967    0.7816    0.2007];
memColor = [0.9139    0.7258    0.3063];
coatColor = 'red';
%coatColor = [0    0.4470    0.7410];
actColor = 'red';

tt = (200 - alpha)/alpha;
coatArea = [200/(tt+1)-act 200/(tt + 1)];
actArea = [0 3.125];

% set visual properties of the plot
fontsize = 36;  %54;  %70;
lineWidth = 8; %6;  %12;
axesWidth = 5;
R0=1;%R0=20; % 1 for non dimensional
mesh = (0:1/2000:1).^2;
t = alpha*mesh;

fighandle = figure(8);
set(fighandle, 'Position', [0, 1000, 300, 300]);
% plot the membrane shape
if dash == 1
    
    plot(Sol(1,:)*R0, Sol(2,:)*R0,  -Sol(1,:)*R0, Sol(2,:)*R0, Sol(7,:)*R0, Sol(8,:)*R0,-Sol(7,:)*R0, Sol(8,:)*R0, 'Color', memColor, 'LineWidth', lineWidth);
    
else
    
    plot(Sol(1,:)*R0, Sol(2,:)*R0, -Sol(1,:)*R0, Sol(2,:)*R0, ':',Sol(7,:)*R0, Sol(8,:)*R0, -Sol(7,:)*R0, Sol(8,:)*R0, ':','Color', memColor, 'LineWidth', lineWidth);

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
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, 'Color', coatColor, 'LineWidth', lineWidth);
            
        else

            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ':','Color', coatColor, 'LineWidth', lineWidth);
            
        end

    end
end

% plots the regions of applied force
if ~isempty(actArea)
    % loops over multiple regions
    for ii = 1:size(actArea,1)

        aMin = actArea(ii,1);
        aMax = actArea(ii,2);
        
        if dash == 1
            
            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, 'Color', actColor, 'LineWidth', lineWidth);
            
        else

            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0,  ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0,':', 'Color', actColor, 'LineWidth', lineWidth);
        
        end

    end
end

hold on
% xlabel('Radius (nm)')
% ylabel('Height (nm)')
xlabel('x')
ylabel('y')
set(gca, 'fontsize',14, 'Fontweight','bold')
daspect([ 1 1 1])
xlim([-290/20 290/20])

end

