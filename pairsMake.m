
%load('C:\Noam\Data\muscimol\cells15nan');
%aclls=cellsn; clear cellsn;
b=[aclls.before]; l =     arrayfun(@(z) len(z.st)>=100,b); 
m=[aclls.midall]; l = l & arrayfun(@(z) len(z.st)>=100,m); 
b=[b.gridscore]; m=[m.gridscore];
%only take cells within threshold and spikes>100
g= b>0.3 & m<0.25; clear b; clear m;
gclls=aclls(g&l); clear g; clear l;

%group by recording date and animal
gs =                 unique(str2double(strcat({gclls.id},{gclls.date}))    )';
group = arrayfun(@(x) gclls(str2double(strcat({gclls.id},{gclls.date}))==x ),gs,'uni',0);
groupis= cellfun (@(x) [x.ind],group,'uni',0)';
%sort by group id;
[~,i]=sort(cellfun(@(x) x(1),groupis));
group=groupis(i)';
clear groupis; clear gs; clear i; clear gccls;

threshRemoveSpikes=-1; 
threshRemoveCell=0.05;
rgs = [];
percentremoved=[];
for gi = 1:len(group)
    gcis=group{gi};
    g = aclls(gcis);
    rcs = [];
    for j = 1:len(g)-1
        for k = j+1:len(g)
            c1 = g(j).before; c2 = g(k).before;
            s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; min_time=min(min(s1),min(s2)); %TIME IN MS
            offset = 0; s1 = s1-min_time+offset; s2 = s2-min_time+offset; max_time=max(max(s1),max(s2))+offset;
            [s1i, s2i] = removeOverlappingSpikes(s1,s2, 1); rs1=s1(s1i); rs2=s2(s2i);
            %fprintf('####b.g%2d i%3dxi%3d: %3d :: %3d\n',ri, g(j).ind ,g(k).ind,length(s1)-length(rs1),length(s2)-length(rs2));
            rem = length(s1)-length(rs1); %assert(abs( (length(s1)-length(rs1)) - (length(s2)-length(rs2)) )<3)
            brem = max(rem/length(s1),rem/length(s2));
            if brem>threshRemoveCell || false
                fprintf('b.g%2d: spikes removed %3dx%3d = %3d, --- c1 %2d%% c2 %2d%% \n',gi, g(j).ind,g(k).ind,rem,...
                    100*round(1-length(rs1)/length(s1),2),100*round(1-length(rs2)/length(s2),2));
            else
                percentremoved(end+1)=brem;
            end
            if brem >  threshRemoveCell;
                %remove cell with higher percentage spikes removed
                rid = gcis(j);
                %s2 got reduced more
                if length(rs2)/length(s2) < length(rs1)/length(s1)
                    rid = gcis(k);
                end
                %fprintf('g.%2d: b c.%3d: %3f%%\n',gids(ri), rid, 100*max(brem,mrem));
                rcs = [rcs rid];
            end
            c1 = g(j).midall; c2 = g(k).midall;
            s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; min_time=min(min(s1),min(s2)); %TIME IN MS
            offset = 0; s1 = s1-min_time+offset; s2 = s2-min_time+offset; max_time=max(max(s1),max(s2))+offset;
            [s1i, s2i] = removeOverlappingSpikes(s1,s2, 1); rs1=s1(s1i); rs2=s2(s2i);
            rem = length(s1)-length(rs1); assert(abs( (length(s1)-length(rs1)) - (length(s2)-length(rs2)) )<3)
            mrem = max( rem/length(s1), rem/length(s2) );
            if mrem>threshRemoveCell || false
                fprintf('m.g%2d: spikes removed %3dx%3d = %3d, --- c1 %2d%% c2 %2d%% \n',gi, g(j).ind,g(k).ind,rem,...
                        100*round(rem/length(s1),2),100*round(rem/length(s2),2));        
                       %100*round(1-length(rs1)/length(s1),2),100*round(1-length(rs2)/length(s2),2));
            else
                percentremoved(end+1)=mrem;
            end
            if mrem > threshRemoveCell
                %remove cell with higher percentage spikes removed
                rid = gcis(j);
                %s2 got reduced more
                if length(rs2)/length(s2) < length(rs1)/length(s1)
                    rid = gcis(k);
                end
                %fprintf('g.%2d: m c.%3d: %3f%%\n',gids(ri), rid, 100*max(brem,mrem));
                rcs = [rcs rid];
            end
        end
    end
    %remove cells
    if ~isempty(rcs)
        rcs = unique(rcs);
        fprintf('$ removing %d cells from g.%d\n',length(rcs),gi);
        [rcs]
        group{gi} = setdiff(gcis, rcs);
    end
end

%remove 1 cell groups
group=group(cellfun(@(x) len(x)>1,group));
%get pairs
pairs=cellfun(@(x) nchoosek(x,2),group,'uni',false);
pairs=vertcat(pairs{:});







