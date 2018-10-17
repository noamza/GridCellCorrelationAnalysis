% function [t1,t2] = createMsSpikeTrain(s1, offms, s2)
%     [t1,t2] = makeMsSpikeTrain(s1, offms, s2);
% end

function [t1,t2] = createMsSpikeTrain(s1, offms, s2)
t1 = []; t2 = [];
s1=ceil(s1*1000); %from s to ms
if sum(s1<1)>0 %BS THAT I NEED THIS, NEED TO PREPROCESS
    fprintf('negative spikes! %d\n',sum(s1<1));
    s1(s1<1)=1;
end
t1 = zeros(1,max(s1)); t1(s1)=1;
if exist('s2','var')
    s2=ceil(s2*1000);
    mint=min(min(s1),min(s2)); 
    s1 = s1-mint+offms; 
    s2 = s2-mint+offms; 
    maxt=max(max(s1),max(s2))+offms; 
    t1 = zeros(1,maxt); t1(s1)=1;
    t2 = zeros(1,maxt); t2(s2)=1;
elseif exist('offms','var')
    s1 = s1-min(s1)+offms;
    t1 = zeros(1,max(s1)+offms); t1(s1)=1;
end

end