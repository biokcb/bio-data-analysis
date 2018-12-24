%plotting tracking results from split up movies
%rough script, so it doesn't take input, but you can edit it easily to do
%this if you want.
% written 2013 by Kristen Brown

function plot2Dtracks_full

%clear all the variables currently in the matlab workspace
clear all;

%load the matlab file saved after tracking
load('/Users/Kristen/Documents/TamkunLab/5.23.13/E18 rHN not frozen - Nav1.6 + dendra2/Dish 1/Cell 1/Tif - 0.3% -2/Tracking_TW2_SR3_1-2000.mat');

%rename the tracksFinal variable
tracksFinal1=tracksFinal;
clear tracksFinal;

%repeat this set of commands for the number of tracking results you have
%saved. There is probably an actual efficient way of doing this, but I just
%thought it would be ok to do it this way. I just copied + pasted the first
%set of commands and edited the file name by updating the numbers. You only
%need to make sure that it is opening the correct saved Tracking file. And
%that you have all of them opening. Say, if you only had 2 sets of tracks,
%then you could comment out the last 3 sets of commands by putting a % in
%front of each line.

%load the matlab file saved after tracking
load('/Users/Kristen/Documents/TamkunLab/5.23.13/E18 rHN not frozen - Nav1.6 + dendra2/Dish 1/Cell 1/Tif - 0.3% -2/Tracking_TW2_SR3_2001-4000.mat');

%rename the tracksFinal variable
tracksFinal2=tracksFinal;
clear tracksFinal;

%load the matlab file saved after tracking
load('/Users/Kristen/Documents/TamkunLab/5.23.13/E18 rHN not frozen - Nav1.6 + dendra2/Dish 1/Cell 1/Tif - 0.3% -2/Tracking_TW2_SR3_4001-6000.mat');

%rename the tracksFinal variable
tracksFinal3=tracksFinal;
clear tracksFinal;

%load the matlab file saved after tracking
load('/Users/Kristen/Documents/TamkunLab/5.23.13/E18 rHN not frozen - Nav1.6 + dendra2/Dish 1/Cell 1/Tif - 0.3% -2/Tracking_TW2_SR3_6001-8000.mat');

%rename the tracksFinal variable
tracksFinal4=tracksFinal;
clear tracksFinal;

%load the matlab file saved after tracking
load('/Users/Kristen/Documents/TamkunLab/5.23.13/E18 rHN not frozen - Nav1.6 + dendra2/Dish 1/Cell 1/Tif - 0.3% -2/Tracking_TW2_SR3_8001-9538.mat');

%rename the tracksFinal variable
tracksFinal5=tracksFinal;
clear tracksFinal;

%now for the plotting...you don't need to touch this as long as the
%settings for the plots 

hold on %this ensures that the plots stay on the same figure.

%generates the first plot in a new window
plotTracks2D(tracksFinal1, [], 2, [], 0, 1, [], 0, 0, [], 10); 
%each other plot will be on the same figure (changed the 1 to a 0)
plotTracks2D(tracksFinal2, [], 2, [], 0, 0, [], 0, 0, [], 10);
plotTracks2D(tracksFinal3, [], 2, [], 0, 0, [], 0, 0, [], 10);
plotTracks2D(tracksFinal4, [], 2, [], 0, 0, [], 0, 0, [], 10);
plotTracks2D(tracksFinal5, [], 2, [], 0, 0, [], 0, 0, [], 10);

hold off
