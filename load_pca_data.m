% get data (modified datafile from from protein.dat)
% X is the n x p predictor (or data) matrix
% n = number of observations (countries)
% p = number of variables (protein sources)

if ~exist('ddir','var'), ddir = './data/'; end
X = load([ddir 'protein_matlab.dat']);      

% columns (variables):
vlabslong = {'RedMeat','WhiteMeat','Eggs','Milk','Fish','Cereals','Starch','Nuts','Fr&Veg'};
vlabs = {'RM','WM','EG','MK','FI','CL','ST','NT','FV'};

% rows (observation points):
dlabslong = {'Albania','Austria','Belgium','Bulgaria','Czechoslovakia',...
           'Denmark','E_Germany','Finland','France','Greece','Hungary',...
           'Ireland','Italy','Netherlands','Norway','Poland','Portugal',...
           'Romania','Spain','Sweden','Switzerland','UK','USSR','W_Germany','Yugoslavia'};
dlabs =  {'ALB','AUS','BEL','BUL','CZK',...
           'DMK','EGR','FIN','FRA','GRE','HUN',...
           'IRE','ITY','NED','NOR','POL','POR',...
           'ROM','SPA','SWE','STZ','UK','USSR','WGR','YUG'};
whos      
