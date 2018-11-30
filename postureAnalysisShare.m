%% Postures - Flexion, Bending and Rotation       
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
C7Acc = C7(:,[3,4,5]);
C7Gyr = C7(:,[6:8]);
T4Acc = T4(:,[3,4,5]);
StAcc = St(:,[3,4,5]);
SaAcc = Sa(:,[3,4,5]);

[~, idxc7] = findpeaks(-C7Acc(:,2),'MinPeakHeight',2);
C7 = C7(idxc7(2):end,:);

[~, idxt4] = findpeaks(-T4Acc(:,2),'MinPeakHeight',2);
T4 = T4(idxt4(2):end,:);

[~, idxst] = findpeaks(-StAcc(:,2),'MinPeakHeight',2);
St = St(idxst(2):end,:);

[~, idxsa] = findpeaks(-SaAcc(:,2),'MinPeakHeight',2);
Sa = Sa(idxsa(2):end,:);


%Find shortest length and cut to the same
ln = min([length(Sa(:,1)), length(C7(:,1)), length(St(:,1)), length(T4(:,1))]);

C7 = C7(1:ln,:); 
Sa = Sa(1:ln,:);
T4 = T4(1:ln,:);
St = St(1:ln,:);

% figure('Name','ACC-y') %Confirm that the IMUs are synchronized, comparing Acc_y.
% plot([C7(:,4),T4(:,4),Sa(:,4),St(:,4)]);
% legend('C7','T4','Sa','St')

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
plot([euler_C7(:,i),euler_T4(:,i)]);
plot([euler_C7_Mad(:,i),euler_T4_Mad(:,i),euler_Sa_Mad(:,i)],'--');
hold off
legend('C7', 'T4', 'C7-Mad', 'T4-Mad','Sa-Mad');
end 

%% calibrate
% This was meant for Madgwick signals due to the long drift of fusion
% algorithm that return to the correct value in the start (the first 10s
% cannot be used).

t1 = 1121; 
t2 = t1 + 10;
C7Q_Mad_cal = calibrateQ(C7Q_Mad,t1, t2); 
SaQ_Mad_cal = calibrateQ(SaQ_Mad,t1, t2);
T4Q_Mad_cal = calibrateQ(T4Q_Mad,t1, t2);
StQ_Mad_cal = calibrateQ(StQ_Mad,t1, t2);


t1_built  = 182; 
t2_built  = t1_built + 10;         
t1_fc = 275;
t2_fc = 333;

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

%%
% Test on functional calibration!!!
[C7Q_built_fc, R] = trunk_neutral_func(C7Q_built, [t1_built,t2_built], [t1_fc,t2_fc]);

C7Q_built_fc_cal = calibrateQ(C7Q_built_fc,t1_built,t2_built); 

euler_C7built_fc = calculateEuler(C7Q_built_fc,'XYZ');
euler_C7built_fc_cal = calculateEuler(C7Q_built_fc_cal,'XYZ');

euler_C7built_fc = quat2eul(C7Q_built_fc,'XYZ')*180/pi;
euler_C7built_fc_cal = quat2eul(C7Q_built_fc_cal,'XYZ')*180/pi;


figure
for i = 1:3
    subplot(3,1,i)
    hold on
    plot(euler_C7built(:,i));
    plot(euler_C7built_fc(:,i),'--');
    plot(euler_C7built_fc_cal(:,i),'.-');
    hold off
end


%%

imuC7.acc = C7Acc;
imuC7.gyr = C7Gyr;

[R] = spineCalibrationModif(imuC7,  [t1_fc,t2_fc], [t1_built,t2_built]);

Rori = quat2rotm(C7Q_Mad);
Rfc = zeros(3,3,length(Rori));

for l = 1:length(Rori)
  Rfc(:,:,l) = Rori(:,:,l) * R';
end


C7Q_fc = rotm2quat(Rfc);
euler_C7built_fc = rotm2eul(Rfc,'XYZ')*180/pi;


C7Q_Mad_cal = calibrateQ(C7Q_fc, t1, t2); 
euler_C7_fc = calculateEuler(C7Q_Mad_cal,'XZY');


figure
for i = 1:3
    subplot(3,1,i)
    hold on
    plot(euler_C7built(:,i));
    plot(euler_C7_fc(:,i),'--');
    hold off
end
