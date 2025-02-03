%% Initialize
cd /Users/joshuadimasaka/Desktop/PhD/GitHub/DeepTemporal
addpath('/Users/joshuadimasaka/Desktop/PhD/GitHub/DeepTemporal')

clear
close all
clc

%% Data

% Brgy ID and Label
[BrgyID, mapR] = readgeoraster("data/EMI/Geospatial References/Barangay Boundaries/brgyBoundary.tif");
LabelBrgyID = readtable("data/EMI/Geospatial References/Barangay Boundaries/brgyBoundary.xlsx");

% L4 and L5 - Land Use Class
[L4, ~] = readgeoraster("data/GMMA_Exposure_Compilation_August2013/L4.tif", 'OutputType','uint8');
[L5, ~] = readgeoraster("data/GMMA_Exposure_Compilation_August2013/L5.tif", 'OutputType','uint8');
L4(BrgyID==0) = 0; 
L4 = uint8(L4);
geotiffwrite("data/GMMA_Exposure_Compilation_August2013/L4.tif",L4,mapR)
L5(BrgyID==0) = 0; 
L5 = uint8(L5);
geotiffwrite("data/GMMA_Exposure_Compilation_August2013/L5.tif",L5,mapR)

% Table of Prior Probability Distribution, p
% Category x (Number of Classes = 19)
% Classes are:
classes = [ "W1", "W2", "W3", ...
            "N", "CHB", "URA", "URM", ...
            "RM1", "RM2", "MWS", "CWS", ...
            "C1", "C2", "C4", "PC2", ...
            "S1", "S2", "S3", "S4"]';

% NonResidential
% Note (Consraint): p(MWS) = 0
NR = readtable("data/NonResiPrior.xlsx");

% Residential
% Note (Constraint): only W1, N, CHB, URA, URM, MWS, CWS, C1, S1, S3
R = readtable("data/ResiPrior.xlsx");
validR = ["W1", "N", "CHB", "URA", "URM", "MWS", "CWS", "C1", "S1", "S3"]';

% Generate and export prior maps
rL4 = table2array(NR(NR.NonResi==0, "IndexL4"));
rL5 = table2array(NR(NR.NonResi==0, "IndexL5"));
nrL4 = table2array(NR(NR.NonResi==1, "IndexL4"));
nrL5 = table2array(NR(NR.NonResi==1, "IndexL5"));
irL45 = (   (L4==rL4(1) & L5==rL5(1))   |   ...
            (L4==rL4(2) & L5==rL5(2))   |   ...
            (L4==rL4(3) & L5==rL5(3))   |   ...
            (L4==rL4(4) & L5==rL5(4))   |   ...
            (L4==rL4(5) & L5==rL5(5))   |   ...
            (L4==rL4(6) & L5==rL5(6))   );

for iClass = 2:length(classes)

    tic
    % Initialize
    prior = single(zeros(size(BrgyID)));

    % 2. Check barangay
    if ismember(classes(iClass), validR)
        for iBrgy = 1:max(LabelBrgyID.fid)
            iBrgy
            % iBrgy, tic
            idx =   (BrgyID == iBrgy) & irL45;
            prior(idx) = table2array(R(R.fid==iBrgy, classes(iClass))); % toc
        end
    end

    % NonResidential
    if classes(iClass) ~= "MWS"
        for idx_nr = 1:length(nrL4)
            idx_nr
            % idx_nr, tic
            idx =   ( L4==nrL4(idx_nr) & L5==nrL5(idx_nr) );
            prior(idx) = table2array(NR(NR.IndexL4==nrL4(idx_nr) & NR.IndexL5==nrL5(idx_nr), classes(iClass))); % toc
        end
    end
    toc, iClass

    tic
    geotiffwrite("data/GMMA_Exposure_Compilation_August2013/prior_"+classes(iClass)+".tif", prior, mapR)
    toc

end



