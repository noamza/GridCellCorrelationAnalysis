function pcor0 = shuffleTimeCorrelations (c1, c2, n, lag, sigma)
    %n=100;
    %lags = 0.2; sigma = 0.006; %lag in s; determined empirically
    %convert to ms and start at 1
    %c2 = cells{71}.before; c1 = cells{75}.before; %BEFORE   
    s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; 
    mint=min(min(s1),min(s2));
    s1 = s1-mint+10; s2 = s2-mint+10; %pad edges
    maxt=max(max(s1),max(s2))+10; %MAKE CORRELATIONS SMALLER? pad ends
    %remove overlapping spikes
    [s1i s2i] = removeOverlappingSpikes(s1,s2, 1); s1=s1(s1i); s2=s2(s2i);
    %create train for c1
    train1 = zeros(1,maxt); train1(s1)=1;
    train2 = zeros(1,maxt); train2(s2)=1;
    %delta to shift train by
    d = floor(length(train1)/n) ;%right?  %amount to shift by each step
    %trains=zeros(n+1,maxt); trains(i,:)=train2shifted; save all
    pcor0 = []; 
    %t2shifted = zeros(n+1,maxt);
    %win=hamming(15);
    for i = 1:n+1 %n+1 because wi = (i-1)*d; to include 0 shift
        %t2shifted(i,:) = conv(train2([(i-1)*d+1:length(train1) 1:(i-1)*d]),win,'same');
    end 
    %pcor=xcov(train1Ham, train2Ham,lag,'coef');
    
    for i = 1:n %n+1 because wi = (i-1)*d; to include 0 shift
        sprintf('done %.2f\%',round(i/(n+1),2));
        wi = (i-1)*d; %first one calculated is 0
        shiftedi = [wi+1:length(train1) 1:wi]; %cyclic
        train2shifted = train2(shiftedi); %trains(i,:)=train2shifted;
        %pcor = timeCorrelationSmoothed(train1,train2shifted,15,lag,sigma);
        pcor = timeCorrelationSmoothed(sparse(train1),sparse(train2shifted),15,lag,sigma);
        pcor0(i) = pcor(round(length(pcor)/2));
        %{
        %smoothing spike train
        win=hamming(15);
        %train1Ham = train1; train2ham = train2;
        train1Ham=conv(train1,win,'same');train2Ham=conv(train2shifted,win,'same'); %win = win/sum(win);
        %  figure, plot(train1Ham), hold on, plot(train1 + 1);
        %correlation
        %pcor=xcov(train1Ham, train2Ham,lag,'coef'); % TRY AS MATRIX
        pcor=xcorr(train1Ham, train2Ham,lag,'coef');
        %p{i,3} = [pmx (pmxi - ((length(p{i,1})-1)/2) - 1)/1000]; % dont ask
        %smoothing correlation
        [b,a] = butter(6,sigma); %LOW PASS 0.15%6th order, fc/fs/2 determined empiracally 
        pcor = filtfilt(b,a,pcor);
        
        pcor0(i) = pcor(round(length(pcor)/2)); %corr at 0
        %count_xcor{i,1}=xcorr(train1,train2, round(1000*lags))'; %how many times spike at same point
        %}
    end
    %disp(pcor0);
    %[mx, imx] = max(abs(pcor0))
    %time_scale=( (1:length(p{i,1})) - ((length(p{i,1})-1)/2) - 1)/1000; %goes from -lag 0 +lag in ms
    %figure;
     %plotmatrix(time_scale',cell2mat(p(:,1))','-');
     %figure
     %plotmatrix((1:maxt)',trains','-');
    %plotmatrix(cell2mat(p(1,2))',cell2mat(p(:,1))')
end