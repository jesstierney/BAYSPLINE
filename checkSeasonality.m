function [months] = checkSeasonality(lat,lon)
% Script to check whether a given location has a seasonal alkenone
% signature
%
% months = checkSeasonality(lat,lon)
%
% ----- Inputs -----
% latitude and longitude in decimal degrees. can be scalars or vectors.
%
% ----- Outputs -----
% months of alkenone production in a cell array. if annual it will be 1
% through 12. If Mediterranean, it will be [1 2 3 4 5 11 12]. If N.
% Atlantic it will be [8 9 10]. If N. Pacific it will be [6 7 8].

%ensure column vectors
lat = lat(:);
lon = lon(:);

%calculate number locations
Nloc = length(lat);

% Define polygons
% Mediterranean polygon: Nov-May
poly_m_lat=[36.25; 47.5; 47.5; 30; 30];
poly_m_lon=[-5.5; 3; 45; 45; -5.5];

% North Atlantic polygon: Aug-Oct
poly_a_lat=[48; 70; 70; 62.5; 58.2; 48];
poly_a_lon=[-55; -50; 20; 10; -4.5; -4.5];

% North Pacific polygon: Jun-Aug
poly_p_lat=[45; 70; 70; 52.5; 45];
poly_p_lon=[135; 135; 250; 232; 180];

%calculate lon360 for N Pacific
lon360=lon;
lon360(lon360<0)=360+lon360(lon360<0);

% set default to annual
months = cell(Nloc,1);
[months{:}] = deal(1:12);

% reset to seasonal as needed
indM = inpolygon(lon,lat,poly_m_lon,poly_m_lat);
[months{indM}] = deal([1 2 3 4 5 11 12]);

indA = inpolygon(lon,lat,poly_a_lon,poly_a_lat);
[months{indA}] = deal([8 9 10]);

indP = inpolygon(lon360,lat,poly_p_lon,poly_p_lat);
[months{indP}] = deal([6 7 8]);