%% Main

clc
clear
close all

%%

%Basic parameters%

K = [5000 0.1  2000; %Intrinsic parameters, K, 3*3, 5 degrees of freedom%
    0   4000 1000;
    0   0    1];
    
l = 0.2; %Mirror parameter, (0,1)%


%%
%Calibration

Image_num = 2; % lines image should >= 2
[line_image] = generate_line_image(Image_num,K,l); % generate the simulation lines images

d0 = 0; % distance between mirror contour and sphere center
[mirror_contour] = generate_mirror_contour(d0,K,l); % generate the simulation mirror contour


[est_T_K, est_T_l] = T(line_image,mirror_contour); % use lines images to calibration besed on method T



%%
%Error

rfe_error = (est_T_K(1,1) - K(1,1))/K(1,1); %relative error
fe_error = (est_T_K(2,2) - K(2,2))/K(2,2);
s_error = (est_T_K(1,2) - K(1,2))/K(1,2);
u0_error = (est_T_K(1,3) - K(1,3))/K(1,3);
v0_error = (est_T_K(2,3) - K(2,3))/K(2,3);
l_error = (est_T_l - l)/l;
