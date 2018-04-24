%% LiberateData.m
% This script is intended to convert the .mat files received from Hylke to
% CSVs for analysis in R.
%
% It needs to be run by Matlab from the directory containing the raw data
% (Z:\2.active_projects\Zipper\1.Spatial_data\global\ro_runoff\1original\GlobalDailyStreamflow)

%% Clean up workspace
clear; clc;

%% Step 1: Process the summary of stations
% load .mat file
load('summary_21-Mar-2018.mat')

% get names
out = extractfield(summary_Catchments, 'name')';

% make a big table
T = table(out, summary_aridity, summary_biomes, summary_catchsize, summary_dailyQexists, ...
    summary_drainagedensity, summary_ETPOT_Hargr, summary_forest_gain, summary_forest_loss, ...
    summary_FTC, summary_FW, summary_glaciers, summary_GLWD, summary_HWSD_clay, ...
    summary_HWSD_gravel, summary_HWSD_sand, summary_HWSD_silt, summary_irrig/100, ...
    summary_lat, summary_lon, summary_mean_snow_ground, summary_meanelevation, ...
    summary_meanQobs_entireperiod, summary_meanQobs_studyperiod, summary_meanslope, ...
    summary_meanTa, summary_monthlyQexists, summary_NDVI, summary_nonreference, ...
    summary_P_MSWEP, summary_permafrost, summary_permeability, summary_PPETcorr, ...
    summary_reclength_entireperiod, summary_reclength_studyperiod, summary_res_cap, ...
    summary_res_influence, summary_seasonality_P, summary_seasonality_PET, ...
    summary_snow_fraction, summary_snowfall_fraction, summary_SoilGrids1km_clay, ...
    summary_SoilGrids1km_sand, summary_SoilGrids1km_silt, summary_soilI, ...
    summary_temp_coverage, summary_urban, ...
    NaN(length(summary_urban),1), NaN(length(summary_urban),1), zeros(length(summary_urban),1),...
    'VariableNames', ...
    {'catchment','aridity','biomes','catchsize','dailyQexists', ...
    'drainagedensity','ETPOT_Hargr','forest_gain','forest_loss', ...
    'FTC', 'FW', 'glaciers', 'GLWD', 'HWSD_clay', ...
    'HWSD_gravel','HWSD_sand','HWSD_silt','irrig', ...
    'lat','lon','mean_snow_ground','meanelevation', ...
    'meanQobs_entireperiod','meanQobs_studyperiod','meanslope', ...
    'meanTa','monthlyQexists','NDVI','nonreference', ...
    'P_MSWEP', 'permafrost', 'permeability', 'PPETcorr', ...
    'reclength_entireperiod','reclength_studyperiod','res_cap', ...
    'res_influence', 'seasonality_P', 'seasonality_PET', ...
    'snow_fraction', 'snowfall_fraction', 'SoilGrids1km_clay', ...
    'SoilGrids1km_sand', 'SoilGrids1km_silt', 'soilI', ...
    'temp_coverage', 'urban', ...
    'lat_gage', 'lon_gage', 'natural'});

% ignore first two lines (not actually stations)
T = T(3:end,:);

%% Step 2: Figure out which catchments are natural to save/process output from
% Some parameters which may be useful to filter on:
%  -catchsize [km2] --> Beck et al. (2013, 2016) use <10,000 km2 
%                       Sawicz et al. (2011) use <10,096 km2
%  -nonreference -----> All catchments have nonreference==0; presumably
%                       Hylke removed nonreference gages from dataset
%  -urban [% area] ---> Falcone et al. (2011) GAGES-II uses 5% for
%                       identifying HCDN-2009 reference gages
%  -irrig [% area] ---> 5% for consistency with irrig
%  -res_influence ----> Hylke: >0.1 indicates flows affected by reservoirs
%  -forest_loss [%] --> What cutoff to use? 25% (arbitrary...)
%  -forest_gain [%] --> What cutoff to use? 25% (arbitrary...)
%
% choose what to keep and make some exploratory plots
i_natural = find(T{:,'catchsize'}<10000 & ...
    T{:,'nonreference'}~=1 & ...
    T{:,'urban'}<0.05 & ...
    T{:,'irrig'}<0.05 & ...
    T{:,'res_influence'}<0.1 & ...
    T{:,'forest_loss'}<0.25 & ...
    T{:,'forest_loss'}<0.25);

%i_test = T{i_keep,'res_influence'}>0.1;
%sum(i_test)
%plot(T{i_test,'lon'}, T{i_test,'lat'}, 'o');

%histogram(T{i_natural,'catchsize'});
%plot(T{i_natural,'lon'}, T{i_natural,'lat'}, 'o');

% add column to table
T{i_natural,'natural'} = 1;

%% Step 2: Extract discharge data
% folder containing all the subfolders for individual streams
f_data = 'DATABASE_WORLD_V3';  

% string of datenums
datenums = datenum(1900,1,1):datenum(2039,12,31);

% dates you want to save as output
minyear = 1950;
maxyear = 2017;
datenums_out = datenum(minyear,1,1):datenum(maxyear,12,31);

% find indices of dates you want to save output for
[tf, i_out] = ismember(datenums_out, datenums);
i_out = i_out(i_out ~= 0);

for ind = 1:length(i_natural)
    i = i_natural(ind);

    % clean up workspace
    clearvars -except ind i T datenums datenums_out i_out f_data minyear maxyear i_natural

    % load discharge data
    Qpath = fullfile(f_data, T{i,'catchment'}, 'DISCHARGE.mat');
    load(Qpath{1});
    
    % check if discharge data exist
    if isfield(DISCHARGE, 'Discharge')
        % extract data from dates of interest
        Q = DISCHARGE.Discharge(i_out);
        T{i,'lat_gage'} = DISCHARGE.StationCoords.Lat;
        T{i,'lon_gage'} = DISCHARGE.StationCoords.Lon;

        % load met data
        Ppath = fullfile(f_data, T{i,'catchment'}, 'MSWEP_V220_dly.mat');
        Tpath = fullfile(f_data, T{i,'catchment'}, 'MSMet_V1_Temp_dly.mat');
        Trangepath = fullfile(f_data, T{i,'catchment'}, 'MSMet_V1_Trange_dly.mat');
        Wpath = fullfile(f_data, T{i,'catchment'}, 'MSMet_V1_Wind_dly.mat');

        % all of these met .mat files are a single array named DATA so need to open
        % them one-at-a-time
        load(Ppath{1})
        P = DATA(i_out);

        load(Tpath{1})
        Temp = DATA(i_out);

        load(Trangepath{1})
        Trange = DATA(i_out);

        load(Wpath{1})
        W = DATA(i_out);

        % check all arrays to make sure they are vertical, not horizontal
        if (size(Q,1) < size(Q,2)); Q=Q'; end
        if (size(P,1) < size(P,2)); P=P'; end
        if (size(Temp,1) < size(Temp,2)); Temp=Temp'; end
        if (size(Trange,1) < size(Trange,2)); Trange=Trange'; end
        if (size(W,1) < size(W,2)); W=W'; end
        
        % save as text file
        path_out = fullfile(f_data, T{i,'catchment'}, ...
            strcat('discharge+met_', num2str(minyear), '-', num2str(maxyear), '.csv'));
        fid=fopen(path_out{1},'w');
        fprintf(fid, 'discharge.m3_s,prec.mm_d,Tmean.C,Trange.C,wind.m_s\n');
        fprintf(fid, '%f,%3.3f,%2.2f,%2.2f,%2.2f\n', [Q P Temp Trange W]');
        fclose(fid);
    else
        % if no discharge data: set to 0 to ignore in future analysis
        T{i,'natural'} = 0;
    end
    
    % status update
    disp(strcat(num2str(ind), ' of  ', num2str(length(i_natural)), ' complete'))
    
end

%% Save output table
writetable(T,'summary_21-Mar-2018.csv');