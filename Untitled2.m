
b='';
for i = 1:length(files) %CHECK 36
    a=cells{i};
    if ~isequal(b,strcat(a.id,a.date))
        b=strcat(a.id,a.date);a=a.midall.pt; c(end+1,:)=[min(a),max(a)];
        %fprintf('%d %.1f %.1f \n',i, min(a)/60, max(a)/60);
    end
end
c=[]; b='';
for j = 1:length(cels) %CHECK 36
    i=cels(j);
    a=cells{i};
    if ~isequal(b,strcat(a.id,a.date))
        b=strcat(a.id,a.date);a=a.midall.pt; c(end+1,:)=[min(a),max(a)];
        fprintf('%d %.1f %.1f \n',i, min(a)/60, max(a)/60);
    end
end
    