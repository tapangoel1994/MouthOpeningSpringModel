%% This code generates Figure 3d-e in the accompanying manuscript

%%Read filenames in the directory containing the matfiles which have the timeseries of the spring network configurations.
%%Calculates the ensamble mean mouth area, timescale and delta_r/r for the given force,forcing rate and strength of non-linearity.
%%Then plots these data as heatmaps.

%% Author: Tapan Goel
%% Date: April 12th 2023

%% Get names of all mat files
files = dir('IndepedentOuterRing8Node4RingMicroAsp5000nN\*i=*.mat');


%% Read each matfile and extract macroscopic properties, flag them and store in Results.
gap = 1;
Results = zeros(length(files)/gap,15); %store the summary statistics in a 

for i = 1:length(Results)
filename = [files(1+(i-1)*gap).folder '\' files(1+(i-1)*gap).name];
Results(i,:) = dataextractor(filename,files(1+(i-1)*gap).name);
end

%% Generate phase diagram
Results = array2table(Results,'VariableNames',{'Force', 'Tau', 'Gamma', 'Iter',...
    'MinArea', 'MaxArea', 'MeanRadius', 'StdRadius', 'MajorRadius',...
    'MinorRadius', 'MajorAxisOrientation', 'OpeningTime', 'flag', 'COMx', 'COMy'});

%% See what the different flags are and how many data points of each type we have
% h = figure;
% scatter3(log10(Results.Force),-log10(Results.Tau),Results.Gamma,25,Results.flag,'filled')
% xlabel('Force (nN)');
% ylabel('Forcing Rate(Hz)');
% zlabel('Non-linearCoefficient');
% title('Control Model')

Results(Results.flag < 0,:) = []; %Remove datasets where a) the max mouth area is > 1.2*(Area of outer ring) OR b) OpeningTime < dt*10 = 0.05s OR if there are less than 50 datapoints in the area time series.
Results.COV = Results.StdRadius./Results.MeanRadius; %define delta_r/r
Results.AspectRatio = Results.MajorRadius./Results.MinorRadius; %define aspect ratio
Results.Angle = abs(atan(tan(Results.MajorAxisOrientation - pi/2)))*180/pi;
%% Take ensamble average
avgresults = grpstats(Results,["Force","Tau","Gamma"]);

%% Generate Heat Maps

%%labels for x and y axes
forcelabels = linspace(1,6,21);
forcelabels(forcelabels ~= floor(forcelabels)) = nan;
ratelabels = linspace(-1,1,7);
ratelabels(ratelabels ~= floor(ratelabels)) = nan;

h = figure;
t = tiledlayout(1,2)


nexttile(1);
g1 = heatmap(avgresults,'Force','Tau','ColorVariable',"mean_Angle");
g1.GridVisible = 'off';
g1.ColorLimits = [0 2];
g1.YDisplayLabels = -ratelabels;
g1.XDisplayLabels = forcelabels;
g1.XLabel = '';
g1.YLabel = '';
g1.Title = 'Angle';


nexttile(2);
g1 = heatmap(avgresults,'Force','Tau','ColorVariable',"mean_AspectRatio");
g1.GridVisible = 'off';
g1.ColorLimits = [0 3];
g1.YDisplayLabels = -ratelabels;
g1.XDisplayLabels = forcelabels;
g1.XLabel = '';
g1.YLabel = '';
g1.Title = 'Aspect Ratio';


xlabel(t,'Force (10^x nN)');
ylabel(t,'Forcing Rate (10^y Hz)')