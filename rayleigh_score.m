function [rayleigh_score, rayleigh_angle, phd]=rayleigh_score...
    (pt,px1,py1,px2,py2,sx1,sy1,sx2,sy2)

dt=median(diff(pt));
%head direction
phd = atan2(py2-py1,px2-px1);

% the head direction of the animal when spike has ocurred
shd = atan2(sy2-sy1,sx2-sx1);

bin_size=deg2rad(3);
radbin = -pi:bin_size:pi;
count = hist(shd,radbin);
time = hist(phd,radbin)*dt;
rate_ang=count./time;

win=hamming(10);
win=win/sum(win);

rate_ang_smoothed=conv([rate_ang rate_ang rate_ang],win,'same');

rate_ang=rate_ang_smoothed(length(count)+1:2*length(count));

x_bin = cos(radbin);
y_bin = sin(radbin);
norm_val = nansum(rate_ang);
x_vec = nansum(x_bin.*rate_ang);
y_vec = nansum(y_bin.*rate_ang);
vec_len = sqrt(x_vec.^2+y_vec.^2);

rayleigh_score = vec_len/norm_val;
rayleigh_angle=(atan2(y_vec,x_vec));

%   figure;bar(radbin,time); title('time spent in each direction');
%   figure;bar(radbin,count);title('number of spikes fired in each direction');
%   figure;bar(radbin,rate_ang);title('number of spikes fired in each direction');
  
end