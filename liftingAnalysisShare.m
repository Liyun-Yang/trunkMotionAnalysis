%% Symmetric and asymmetric lifting - in Flexion, Bending and Rotation 
close all
clear variables
% Read in IMU data files
% Data from IMU sensors 3:18, 3 Ax Ay Az, 6 Gx Gy Gz, 9 Mx My Mz,...
% 15-18 Qw Qx Qy Qz 19-21 Euler-x -y -z

T4 = csvread('T4.csv',2,0);
C7 = csvread('C7.csv',2,0);
St = csvread('Sternum.csv',2,0);
Sa = csvread('Sacrum.csv',2,0);

% Synchronize IMUs based on the jump signal from ACC
C7Acc = C7(:,3:5);
T4Acc = T4(:,3:5);
StAcc = St(:,3:5);
SaAcc = Sa(:,3:5);

[~, idxc7] = findpeaks(-C7Acc(:,2),'MinPeakHeight',2);
C7 = C7(idxc7(3):end,:);

[~, idxsa] = findpeaks(-SaAcc(:,2),'MinPeakHeight',2);
Sa = Sa(idxsa(4):end,:);

[~, idxst] = findpeaks(-StAcc(:,2),'MinPeakHeight',2);
St = St(idxst(2):end,:);

[~, idxt4] = findpeaks(-T4Acc(:,2),'MinPeakHeight',2);
T4 = T4(idxt4(2):end,:);


%Find shortest length and cut to the same
ln = min([length(Sa(:,1)), length(C7(:,1)), length(St(:,1)), length(T4(:,1))]);

C7 = C7(1:ln,:); 
Sa = Sa(1:ln,:);
T4 = T4(1:ln,:);
St = St(1:ln,:);

C7Acc = C7(:,3:5);
T4Acc = T4(:,3:5);
StAcc = St(:,3:5);
SaAcc = Sa(:,3:5);

figure('Name','ACC-y') %Confirm that the IMUs are synchronized, comparing Acc_y.
plot([C7(:,4),T4(:,4),Sa(:,4),St(:,4)]);
legend('C7','T4','Sa','St')

% Use two fusion algorithm on IMU raw data
[Q_Madgwick, Q_noMag]= signalFusion(C7(:,3:11),1);
C7Q_Mad = Q_Madgwick; C7Q_noMag = Q_noMag; 

[Q_Madgwick, Q_noMag]= signalFusion(T4(:,3:11),1);
T4Q_Mad = Q_Madgwick; T4Q_noMag = Q_noMag; 

[Q_Madgwick, Q_noMag]= signalFusion(St(:,3:11),1);
StQ_Mad = Q_Madgwick; StQ_noMag = Q_noMag;

[Q_Madgwick, Q_noMag]= signalFusion(Sa(:,3:11),1);
SaQ_Mad = Q_Madgwick; SaQ_noMag = Q_noMag;

%Quaternions from built-in IMUs
SaQ_built  = Sa(:,15:18);
C7Q_built  = C7(:,15:18);
T4Q_built  = T4(:,15:18);
StQ_built  = St(:,15:18);

euler_C7 = calculateEuler(C7Q_built,'XYZ');
euler_T4 = calculateEuler(T4Q_built,'XYZ');
euler_St = calculateEuler(StQ_built,'XYZ');

euler_C7_Mad = calculateEuler(C7Q_Mad,'XZY');
euler_T4_Mad = calculateEuler(T4Q_Mad,'XZY');
euler_Sa_Mad = calculateEuler(SaQ_Mad,'XZY');

% Plot to see when can be used as reference frame
figure('Name', 'IMU Euler Angles uncalibrated');

for i = 1:3
    subplot(3,1,i);
hold on;
plot([euler_C7(:,i),euler_T4(:,i),euler_St(:,i)]);
plot([euler_C7_Mad(:,i),euler_T4_Mad(:,i),euler_Sa_Mad(:,i)],'--');
hold off
legend('C7', 'T4','St', 'C7-Mad', 'T4-Mad','Sa-Mad');
xlabel('Time (s)');
ylabel('Angle (deg)');
end

%% calibrate
% This was meant for Madgwick signals due to the long drift of fusion
% algorithm that return to the correct value in the start (the first 10s
% cannot be used).

t1 = 1050; 
t2 = t1 + 10;
% The built-in quaternion reference frame

t1_built  = 1050; 
t2_built  = t1_built + 10; 
   

C7Q_Mad_cal = calibrateQ(C7Q_Mad,t1, t2); 
SaQ_Mad_cal = calibrateQ(SaQ_Mad,t1, t2);
T4Q_Mad_cal = calibrateQ(T4Q_Mad,t1, t2);
StQ_Mad_cal = calibrateQ(StQ_Mad,t1, t2);

% Euler angle from calibrated sensor fusion signal using Madgwick
euler_C7_Mad = calculateEuler(C7Q_Mad_cal,'XZY');
euler_T4_Mad = calculateEuler(T4Q_Mad_cal,'XZY');
euler_St_Mad = calculateEuler(StQ_Mad_cal,'XZY');
euler_Sa_Mad = calculateEuler(SaQ_Mad_cal,'XZY');

% Built in sensor fusion and calibration
C7Q_built_cal = calibrateQ(C7Q_built,t1_built,t2_built); 
SaQ_built_cal = calibrateQ(SaQ_built,t1_built,t2_built);
T4Q_built_cal = calibrateQ(T4Q_built,t1_built,t2_built);
StQ_built_cal = calibrateQ(StQ_built,t1_built,t2_built);

euler_C7built = calculateEuler(C7Q_built_cal,'XYZ'); %% OBS: the order is different for the built in Q! Due to the different Axis!!!
euler_T4built = calculateEuler(T4Q_built_cal,'XYZ');
euler_Stbuilt = calculateEuler(StQ_built_cal,'XYZ');
euler_Sabuilt = calculateEuler(SaQ_built_cal,'XYZ');

% Rotate sternum IMU around Y axis 180 degree 
% OBS: Rotate the sternum IMU, built-in around Z axis [0 0 0 1], 
% Madgwick around Y axis [0 0 1 0]! 
StQ_built_rot = quatmultiply(StQ_built,[0 0 0 1]);
StQ_built_rot_cal = calibrateQ(StQ_built_rot, t1_built, t2_built);
euler_Stbuilt_rot = calculateEuler(StQ_built_rot_cal,'XYZ');

StQ_Mad_rot = quatmultiply(StQ_Mad,[0 0 1 0]);
StQ_Mad_rot_cal = calibrateQ(StQ_Mad_rot, t1, t2);
euler_St_Mad_rot = calculateEuler(StQ_Mad_rot_cal,'XZY');

Fs =25;
time = (1:length(euler_C7_Mad))./Fs;
figure('Name', 'Single IMU Euler Angles from Madgwick quaternion');
for i =1:3
subplot(3,1,i);
hold on;
plot(time, [euler_C7_Mad(:,i),euler_T4_Mad(:,i) ]);
plot(time, euler_Sa_Mad(:,i),'-.','linewidth',2);
plot(time, euler_St_Mad_rot(:,i));
plot(time, -[euler_C7built(:,i),euler_T4built(:,i) ],'--','linewidth',1.5);
plot(time, -euler_Sabuilt(:,i),'--','linewidth',1.5);
plot(time, -euler_Stbuilt_rot(:,i),'--','linewidth',1.5); 
hold off
legend('C7Mad', 'T4Mad', 'SaMad', 'StMadRot', 'C7built', 'T4built','Sabuilt', 'StbuiltRot');
end
xlabel('Time (s)');
ylabel('Angle (deg)');


% Calculate the realtive euler angles between quaternions and calibrate to
% time frame t1:t2
Sa_C7_Mad = relativeQ(SaQ_Mad, C7Q_Mad, 'XZY',t1,t2);
Sa_T4_Mad = relativeQ(SaQ_Mad, T4Q_Mad, 'XZY',t1,t2);
Sa_St_Mad = relativeQ(SaQ_Mad, StQ_Mad_rot_cal, 'XZY',t1,t2); 

Sa_C7_built = relativeQ(SaQ_built, C7Q_built, 'XYZ',t1_built,t2_built);
Sa_T4_built = relativeQ(SaQ_built, T4Q_built, 'XYZ',t1_built,t2_built);
Sa_St_built = relativeQ(SaQ_built, StQ_built_rot_cal, 'XYZ',t1_built,t2_built);

%
time = (1:length(Sa_C7_Mad))./Fs;
figure('Name', 'Relative Euler Angles comparison');
for i = 1:3
subplot(3,1,i);
hold on
plot(time, [Sa_C7_Mad(:,i),Sa_T4_Mad(:,i), Sa_St_Mad(:,i)]);
plot(time, -[Sa_C7_built(:,i),Sa_T4_built(:,i),Sa_St_built(:,i) ], '--');
hold off
legend('Sa-C7-Mad', 'Sa-T4-Mad', 'Sa-St-Mad','Sa-C7-built', 'Sa-T4-built', 'Sa-St-built');
grid on
end


%% Test on functional calibration!
[C7Q_built_fc, R] = trunk_neutral_func(C7Q_built, [t1_built,t2_built], [200,201]);
euler_C7built_fc = calculateEuler(C7Q_built_fc,'XYZ');

figure
for i = 1:3
    subplot(3,1,i)
    hold on
    plot(euler_C7built(:,i));
    plot(euler_C7built_fc(:,i),'--');
    hold off
end
legend('Ori','fc')
