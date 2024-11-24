clear
close all
clc

[ob16,r] = readgeoraster('data/openBLDGtemporal/2016_v3.tif');
[ob17,~] = readgeoraster('data/openBLDGtemporal/2017_v3.tif');
[ob18,~] = readgeoraster('data/openBLDGtemporal/2018_v3.tif');
[ob19,~] = readgeoraster('data/openBLDGtemporal/2019_v3.tif');
[ob20,~] = readgeoraster('data/openBLDGtemporal/2020_v3.tif');
[ob21,~] = readgeoraster('data/openBLDGtemporal/2021_v3.tif');
[ob22,~] = readgeoraster('data/openBLDGtemporal/2022_v3.tif');
[ob23,~] = readgeoraster('data/openBLDGtemporal/2023_v3.tif');

geotiffwrite('data/openBLDGtemporal/diff_2017_2016.tif', ob17-ob16, r);
geotiffwrite('data/openBLDGtemporal/diff_2018_2017.tif', ob18-ob17, r);
geotiffwrite('data/openBLDGtemporal/diff_2019_2018.tif', ob19-ob18, r);
geotiffwrite('data/openBLDGtemporal/diff_2020_2019.tif', ob20-ob19, r);
geotiffwrite('data/openBLDGtemporal/diff_2021_2020.tif', ob21-ob20, r);
geotiffwrite('data/openBLDGtemporal/diff_2022_2021.tif', ob22-ob21, r);
geotiffwrite('data/openBLDGtemporal/diff_2023_2021.tif', ob23-ob22, r);
geotiffwrite('data/openBLDGtemporal/diff_2023_2016.tif', ob23-ob16, r);