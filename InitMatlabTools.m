% Initialize MatlabTools

clear; clc; close all;

global rootDir; 

% change this path for your computer (ispc() check if matlab runs on
% Microsoft or not). Dont know if it's necessary.
if ispc()
    rootDir = 'C:\Users\theomariotte\Documents\01_work\Matlab\MatlabTools\';
else
    rootDir = '/Users/theomariotte/Documents/01_work/Matlab/MatlabTools/'; 
end

% add path to functions to the matlab path
% WARNING : update this when you create a new fomder within MatlabTools
freqpath = [rootDir 'freq/'];
filterpath = [rootDir 'filter/'];
figurepath = [rootDir 'figures/'];

addpath(freqpath);
addpath(filterpath);
addpath(figurepath);

fprintf('Initialization done !')
