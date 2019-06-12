x(:) = FvsZp(1,:)/20;
%y = FvsZp(3,:)*20/320;
%y = Epos4*20/(2*pi*20*320);
y = imag(sol);
st=1;
f = length(FvsZp(1,:));
x = x(st:f);
y = y(st:f);

fontsize = 14;

fighandle = figure(1);
hold on
set(fighandle, 'Position', [0, 1000, 300, 300]);
set(gca, 'fontsize',fontsize, 'fontweight','bold')
xlabel('x')
%ylabel('F_{radial}R/\kappa')
%ylabel('ER/(2\piZ\kappa)')
ylabel('Imaginary part of wavenumber p')
h = plot(x,y);
set(h                          , ...
  'Color'           , [0.5 0 0.5]    , ...  % use 0 0 .5, 0 .5 0, .5 0 0, 0 .5 .5,.5 .5 0,.5 0 .5
  'LineStyle' , '--' , ... % use --
  'LineWidth', 3, ... % use 3
  'Marker', 'd',...  % use *, <, o, h, >,d
  'MarkerSize', 7, ... % use 7
  'MarkerIndices',1:25:length(y));

% general graphics, this will apply to any figure you open
% (groot is the default figure object).
set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 1, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 14, ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 8, ...
'DefaultTextFontName', 'Helvetica', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);
 
% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');