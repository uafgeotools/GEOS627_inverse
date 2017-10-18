% Personal start script
home = getenv('HOME');

% ADD Aster library (GEOS 627)
addpath('/usr/local/matlab_toolboxes/aster/cd_5.2/Lib/');

% ADD GEOTOOLS (assumes it is in your home)
%addpath([home '/GEOTOOLS/matlab_util']);
%startup_geotools

% avoid figures splitting two monitors
set(0,'defaultFigurePosition',[1 600 600 600]);
