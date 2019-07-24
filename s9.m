
dbstop if error

load('C:\Noam\Data\muscimol\aclls15min','aclls'); %personal
load('.\data\roomdose','room'); %personal

s = {'before';'midall';'after'};tss = {'pre','dur','post'};
vs='rayleigh_score';va='rayleigh_angle';


%room figure
figure(2900);rst=0.3; st='all';
clf; pr = [200 10 800 800]; set(gcf,'position',pr); 
for x=1:3 %session
    t = [aclls.(s{x})]; ts = [t.(vs)]'; tl=arrayfun(@(z) len(z.st),t)';
    z=unique(room);d= arrayfun(@(x) [t(ts>rst & room==x & tl>100).(va)],z,'uni',0);
    for y=1:len(z) %z to d
        subplot(3,len(z),(x-1)*len(z)+y);cla; 
        a=rad2deg(d{y});histogram(a,-180:20:180);xlim([-180 180]); 
        xlabel('r-angle'); axis square;
        title(sprintf('%s rm%d n=%d a=%.0f+-%.0f',tss{x},z(y),len(a),mean(a),std(a)));
    end
end
suptitle([st ' r-angle by room(rm#) rscore> ' n2(rst)]);