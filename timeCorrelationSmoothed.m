%lag in ms
function [ pcor ] = timeCorrelationSmoothed( train1,train2,win,lag,sigma)
    train1Ham = train1; train2Ham = train2;
    if win>0
        win=hamming(win);
        train1Ham=conv(train1,win,'same');train2Ham=conv(train2,win,'same'); %win = win/sum(win);
    end
    %  figure, plot(train1Ham), hold on, plot(train1 + 1);
    %correlation
   %pcor = xcorr(train1Ham-mean(train1Ham), train2Ham-mean(train2Ham),lag,'coef');%SAME AS NEXT LINE
    lag=round(lag);
    pcor = xcov(train1Ham, train2Ham,lag,'coef'); % TRY AS MATRIX 'coef'
    %p{i,3} = [pmx (pmxi - ((length(p{i,1})-1)/2) - 1)/1000]; % dont ask
    %smoothing correlation
    if sigma ~= 0
        [b,a] = butter(7,sigma); %LOW PASS 0.15%6th order, fc/fs/2 determined empiracally
        pcor = filtfilt(b,a,pcor);
    end
end

