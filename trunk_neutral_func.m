function [quat_neutral,offset] = trunk_neutral_func(quat,neu_ind, flex_ind)
% TRUNK_NEUTRAL_FUNC Aligns sensor using functional calibration 
%   Author: Howard Chen, PhD
%   Note: quaternion is in nx4 array, first element is real. 
%
%   Inputs: quat- orientation quaternion
%           neu_ind- indice of neutral position 
%           flex_ind- indice of forward flexion for functional calibration 
%
%  note: both neu_ind and flex_ind can be 1x2 array with
%           start and end indice, of which the function will average it)

%   Output: quat_neutral- quaternion with offset included
%           offset: offset calculated with functional calibration

% use averaging if neu_ind is greater than 1
if length(neu_ind) == 1
    R = quatToDCM(quat(neu_ind,:));

else
    R = quatToDCM(quat(neu_ind(1):neu_ind(2),:));
    R = mean(R,3);       
end
mcnz = R(3,:)'./norm(R(3,:)');

% use averaging if flex_ind is greater than 1
if length(flex_ind) == 1 
    R = quatToDCM(quat(flex_ind,:));

else
    R = quatToDCM(quat(flex_ind(1):flex_ind(2),:));
    R = mean(R,3);  
end
mcnxt = R(3,:)'./norm(R(3,:)');

mcny = cross(mcnxt,mcnz);
mcnx = cross(mcny,mcnz);
offset = [mcnx, mcny, mcnz];

%apply offset
quat_neutral = quatMultiply(quat,dcmToQuat(offset));


end