%% This code generates Figure 2e in the accompanying manuscript

%%Read filenames in the directory containing the matfiles which have the timeseries of the spring network configurations.
%%Calculates the ensamble mean mouth area, timescale and delta_r/r for the given force,forcing rate and strength of non-linearity.
%%Then plots these data as heatmaps.

%% Author: Tapan Goel
%% Date: April 12th 2023

%% Get names of all mat files
files = dir('IndepedentOuterRing8Node4Ring_AzimStrength\*i=*.mat');

%% Read each matfile and extract macroscopic properties, flag them and store in Results.
gap = 1;
Results = zeros(length(files)/gap,16); %store the summary statistics in a 

for i = 1:length(Results)
filename = [files(1+(i-1)*gap).folder '\' files(1+(i-1)*gap).name];
Results(i,:) = dataextractor(filename,files(1+(i-1)*gap).name);
end

%% Generate phase diagram
Results = array2table(Results,'VariableNames',{'Force', 'Tau', 'Gamma','k_azim', 'Iter',...
    'MinArea', 'MaxArea', 'MeanRadius', 'StdRadius', 'MajorRadius',...
    'MinorRadius', 'MajorAxisOrientation', 'OpeningTime', 'flag', 'COMx', 'COMy'});

%% See what the different flags are and how many data points of each type we have
% h = figure;
% scatter3(log10(Results.Force),-log10(Results.Tau),Results.k_azim,25,Results.flag,'filled')
% xlabel('Force (nN)');
% ylabel('Forcing Rate(Hz)');
% zlabel('Non-linearCoefficient');
% title('Control Model')

Results(Results.flag < 0,:) = []; %Remove datasets where a) the max mouth area is > 1.2*(Area of outer ring) OR b) OpeningTime < dt*10 = 0.05s OR if there are less than 50 datapoints in the area time series.
Results.COV = Results.StdRadius./Results.MeanRadius; %define delta_r/r

%% Take ensamble average
avgresults = grpstats(Results,["Force","Tau","Gamma","k_azim"]);

%% Generate Heat Maps
lowk = avgresults.k_azim == 0.5;
midk = avgresults.k_azim == 2;
highk = avgresults.k_azim == 4;

%%labels for x and y axes
forcelabels = linspace(1,6,21);
forcelabels(forcelabels ~= floor(forcelabels)) = nan;
ratelabels = linspace(-1,1,7);
ratelabels(ratelabels ~= floor(ratelabels)) = nan;

h = figure;
t = tiledlayout(3,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% delta_r/r %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(1);
g1 = heatmap(avgresults(lowk,:),'Force','Tau','ColorVariable',"mean_COV");
g1.GridVisible = 'off';
g1.ColorLimits = [0 .15];
g1.YDisplayLabels = -ratelabels;
g1.XDisplayLabels = forcelabels;
g1.XLabel = '';
g1.YLabel = '';
g1.Title = '\kappa _{azim} = 0.5 \kappa_0';
%title('\frac{\delta r}{r}');

nexttile(2);
g1 = heatmap(avgresults(midk,:),'Force','Tau','ColorVariable',"mean_COV");
g1.GridVisible = 'off';
g1.ColorLimits = [0 .15];
g1.YDisplayLabels = -ratelabels;
g1.XDisplayLabels = forcelabels;
g1.XLabel = '';
g1.YLabel = '';
g1.Title = '\kappa _{azim} = 2 \kappa_0';

nexttile(3);
g1 = heatmap(avgresults(highk,:),'Force','Tau','ColorVariable',"mean_COV");
g1.GridVisible = 'off';
g1.ColorLimits = [0 .15];
g1.YDisplayLabels = -ratelabels;
g1.XDisplayLabels = forcelabels;
g1.XLabel = '';
g1.YLabel = '';
g1.Title = '\kappa _{azim} = 4 \kappa_0';

xlabel(t,'Force (10^x nN)');
ylabel(t,'Forcing Rate (10^y Hz)')