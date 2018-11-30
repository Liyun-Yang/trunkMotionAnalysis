function [ euler ] = calculateEuler( q , string)
%This function calculates the euler angle from quarternion in the order of
% string (e.g.'XZY'). Output is in degree.

matrix = quat2dcm(q);
[r1,r2,r3] = dcm2angle(matrix, string);
euler = [rad2deg(r1) rad2deg(r2) rad2deg(r3)];

end

