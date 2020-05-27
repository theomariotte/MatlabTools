% test figure formating for latex use or other

clear; clc; close all;

path = '/Users/theomariotte/Documents/01_work/CSTB/PFE/MATLAB/Codes/Fig/';
[figname,dapath] = uigetfile('*.fig','Get .fig file...',path);
figname = [dapath '/' figname];
fig2latex(figname);
