%% Initialize
cd /Users/joshuadimasaka/Desktop/PhD/GitHub/DeepTemporal
addpath('/Users/joshuadimasaka/Desktop/PhD/GitHub/DeepTemporal')

clear
close all
clc

%% Data

[BrgyID, mapR] = readgeoraster("data/EMI/Geospatial References/Barangay Boundaries/brgyBoundary.tif");

% Points
utility_points_layer = ...
    ["school", "market", "worship", "gas", "hotel", ...
     "mall", "bridge", "daycare", "brgyhall", "chemfactory", ...
     "factory", "evacuation", "hospital", "healthcenter", "footbridges", ...
     "multihall", "mrf", "coveredcourt"];
[u, ~] = readgeoraster("data/EMI/Built Environment/utility_points.tif");
for i = 1:length(utility_points_layer)
    geotiffwrite("data/EMI/Built Environment/utility_points_distance_from_"+utility_points_layer(i)+".tif", ...
                 bwdist(u==i).*(BrgyID>0), mapR) %, 'TiffType','bigtiff')
end

% Lines
utility_lines_layer = ...
    ["bikelane", "citybusroute", "jeepneyroute", "uvroute", ...
     "railway", "roadprimary", "roadsecondary", "roadothers"];
[u, ~] = readgeoraster("data/EMI/Built Environment/utility_lines.tif");
for i = 1:length(utility_lines_layer)
    geotiffwrite("data/EMI/Built Environment/utility_lines_distance_from_"+utility_lines_layer(i)+".tif", ...
                 bwdist(u==i).*(BrgyID>0), mapR) %, 'TiffType','bigtiff')
end

% Areas
utility_areas_layer = ...
    ["cemetery", "commercial", "industrial", "informalsettlement", "institutional", ...
     "military", "openspace", "recreational", "reservoir", "residential", ...
      "roads", "socializedhousing", "utility", "vacant", "waterways"];
[u, ~] = readgeoraster("data/EMI/Built Environment/utility_areaZones.tif");
for i = 1:length(utility_areas_layer)
    geotiffwrite("data/EMI/Built Environment/utility_areas_distance_from_"+utility_areas_layer(i)+".tif", ...
                 bwdist(u==i).*(BrgyID>0), mapR) %, 'TiffType','bigtiff')
end

% Disaster Proxy

% Flooding N-year Maps
% <= 0.2: no significant effect
% 0.2-0.5: affect the stability of moving vehicle
% 0.5-1.5: could result in drowning for small children and threaten the
% stability of adults
% 1.5-3.0: building utilities and services no longer functional (water and
% sanitation, electrical) or possibly cutoff from supply
% >3.0: space for human occupancy is lost; high probability of disease
% infection

% We will use 100-year return period because according to QCDRA:
% The selection of the ‘100-year flood’ term used in this report was made because the patterns of
% inundation and damages that can be expected are closer to TS Ondoy and the City Quezon City
% stakeholders can relate to this event

returnperiod = "100";
[u, ~] = readgeoraster("data/EMI/Hazards Vulnerability and Risk/" + ...
    "Baseline Flood Simulation (QC DMP)/" + ...
    "Baseline Flood Simulation - Flow Depth/" + ...
    "floodDepth_baseline_"+returnperiod+"yr.tif");
depthrange = [  0.2 0.5; ...
                0.5 1.5; ...
                1.5 3.0; ...
                3 1000];
for i = 1:size(depthrange,1)
    i, tic
    geotiffwrite("data/EMI/Hazards Vulnerability and Risk/" + ...
                 "Baseline Flood Simulation (QC DMP)/" + ...
                 "Baseline Flood Simulation - Flow Depth/" + ...
                 returnperiod+"yr_flood_proximity_"+ ...
                 "greaterthan_"+string(depthrange(i,1))+"m.tif", ...
                 bwdist(u>depthrange(i,1)).*(BrgyID>0), ...
                 mapR) %, 'TiffType','bigtiff')
    toc
end