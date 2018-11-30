function [ Q_cal ] = calibrateQ( q, t1, t2 )
% This function calculate the quaternion relative to the reference time
% frame t1:t2

 Q_cal = quatmultiply(quatconj(mean(q(t1:t2,:))), q);

end

