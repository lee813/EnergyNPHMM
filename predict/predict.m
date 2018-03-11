%% HSMM Prediction via negtive binomial duration param

% seq = featureVector';
% pUnique = unique(seq);
% NN = length(seq);
% 
% newSeq = zeros(size(seq));
% for idx = 1:length(pUnique)
%     findIndex = find(seq == pUnique(idx));
%     newSeq(1,findIndex) = idx;
% end
% AA = transitionProb(newSeq);
% 
% pSizeUnique = size(pUnique);
% durParams = zeros(pSizeUnique(2),2);
% 
% 
% % for idx = 1:pSizeUnique(2)
% %     aa = find(newSeq == idx);
% %     A = zeros(size(seq));
% %     A(aa) = 1;
% % 
% %     ne0 = find(A~=0);                                   % Nonzero Elements
% %     ix0 = unique([ne0(1) ne0(diff([0 ne0])>1)]);        % Non-Zero Segment Start Indices
% %     eq0 = find(A==0);                                   % Zero Elements
% %     ix1 = unique([eq0(1) eq0(diff([0 eq0])>1)]);        % Zero Segment Start Indices
% %     ixv = sort([ix0 ix1 length(A)]);     % Consecutive Indices Vector
% %     section = [];
% % 
% %     for k1 = 1:length(ixv)-1
% %         element = A(ixv(k1):ixv(k1+1)-1);
% %         if element(1) == 1
% %            section(length(section) + 1) = length(element);
% %         end
% %     end
% %     parmhat = nbinfit(section);
% %     durParams(idx,:) = parmhat;
% %     
% % end
% 
% for idx = 1:pSizeUnique(2)
%     aa = diff(find(newSeq == idx));
%     parmhat = nbinfit(aa);
%     durParams(idx,:) = parmhat;
% end
% 
% TN = 1000;
% 
% hsmmStates = hsmmPredict(durParams, AA, TN);
% realSeq = pUnique(hsmmStates);
% figure;
% plot(realSeq, 'r*');
% figure;
% plot(featureVector, 'r*');

%% HSMM Prediction via poisson duration param

% seq = featureVector';
% pUnique = unique(seq);
% NN = length(seq);
% 
% newSeq = zeros(size(seq));
% for idx = 1:length(pUnique)
%     findIndex = find(seq == pUnique(idx));
%     newSeq(1,findIndex) = idx;
% end
% AA = transitionProb(newSeq);
% 
% pSizeUnique = size(pUnique);
% durParams = zeros(pSizeUnique(2),1);
% 
% % CALCULATE LENGTH OF DIFF == 1 !!
% 
% for idx = 1:pSizeUnique(2)
%     aa = find(newSeq == idx);
%     A = zeros(size(seq));
%     A(aa) = 1;
% 
%     ne0 = find(A~=0);                                   % Nonzero Elements
%     ix0 = unique([ne0(1) ne0(diff([0 ne0])>1)]);        % Non-Zero Segment Start Indices
%     eq0 = find(A==0);                                   % Zero Elements
%     ix1 = unique([eq0(1) eq0(diff([0 eq0])>1)]);        % Zero Segment Start Indices
%     ixv = sort([ix0 ix1 length(A)]);     % Consecutive Indices Vector
%     section = [];
% 
%     for k1 = 1:length(ixv)-1
%         element = A(ixv(k1):ixv(k1+1)-1);
%         if element(1) == 1
%            section(length(section) + 1) = length(element);
%         end
%     end
%     parmhat =  poissfit(section);
%     durParams(idx) = parmhat;
%     
% end
% 
% TN = 1000;
% 
% hsmmStates = hsmmPredict(durParams, AA, TN);
% realSeq = pUnique(hsmmStates);
% close all;
% figure;
% plot(realSeq, 'r*');
% figure;
% plot(featureVector, 'r*');


%% HMM Prediction 
seq = featureVector';
pUnique = unique(seq);
NN = length(seq);

newSeq = zeros(size(seq));
for idx = 1:length(pUnique)
    findIndex = find(seq == pUnique(idx));
    newSeq(1,findIndex) = idx;
end
AA = transitionProb(newSeq);

N = 603;
predictSeq = zeros(N,1);
for idx = 1:N
    if idx ~= 1
        predictSeq(idx,1) = predictState(AA,predictSeq(idx-1,1));
    else
        predictSeq(1,1) = 1;
    end
end

realSeq = pUnique(predictSeq);
figure;
plot(realSeq, 'r+');
hold;
plot(featureVector, 'b*');
legend('Prediction', 'Ground truth')
%% Prediction accuracy percentage
% close all;
% figure;
uniqP = unique(realSeq);
countP = histc(realSeq, uniqP);

% subplot(1,2,1);
% pie(countP,[1 1 0 0 0 0 0 0 0 0]);

uniqT = unique(featureVector');
countT = histc(featureVector', uniqT);

rmse=sqrt(sum((countT(:)-countP(:)).^2)/numel(countT));
confusionmat(countP,countT)
% subplot(1,2,2);
% pie(countP);
% labels = {'1','2','3','4'};
% legend(labels,'Location','southoutside','Orientation','horizontal');
% hBar = bar([countP;countT], 'stacked');
% xt = get(gca, 'XTick');
% set(gca, 'XTick', xt, 'XTickLabel', {'Prediction' 'Ground Truth'})




%% Generate MVGaussian data from state 

% GeneratedData = zeros(6,N);
% 
% for idx = 1:length(realSeq)
%     mean = Psi200.theta(realSeq(idx)).mu;
%     sigma = inv(Psi200.theta(realSeq(idx)).invSigma);
%     tempGen =  mvnrnd(mean,sigma,1);
%     while ~all(tempGen <1 & tempGen >-1)
%         tempGen =  mvnrnd(mean,sigma,1);
%     end
%     if(all(tempGen <1 & tempGen >-1))
%         GeneratedData(:,idx) = tempGen;
%     end
% end
% 
% % rmse = sqrt(sum(datax-GeneratedData').^2)./datax;
% % rmse = xcorr(datax(:,1),GeneratedData(1,:)');
% 
% close all;
% figure;
% plot(datax);
% figure;
% plot(GeneratedData');
% 
%% ARMIA https://www.mathworks.com/help/ident/ug/forecasting-predator-prey-populations.html
% 
% y = datax(:,1);
% % seq = featureVector';
% pUnique = unique(seq);
% % Psi200.stateSeq(5).z;
% NN = length(seq);
% 
% newSeq = zeros(size(seq));
% for idx = 1:length(pUnique)
%     findIndex = find(seq == pUnique(idx));
%     newSeq(1,findIndex) = idx;
% end
% 
% % Psi200.stateSeq(5).z;
% pSizeUnique = size(pUnique);
% durParam = zeros(2,pSizeUnique(2));
% 
% for idx = 1:pSizeUnique(2)
%     aa = diff([1 find(newSeq == idx) NN]);
%     parmhat = nbinfit(aa);
%     durParam(:,idx) = parmhat;
% end
% 
% rndNB = nbinrnd(0.25,0.005);
% T = length(y);
% 
% Mdl = arima(2,0,6);
% EstMdl = estimate(Mdl,y);
% 
% 
% [yF,yMSE] = forecast(EstMdl,60,'Y0',y);
% upper = yF + 1.96*sqrt(yMSE);
% lower = yF - 1.96*sqrt(yMSE);
% 
% figure
% plot(y,'Color',[.75,.75,.75])
% hold on
% h1 = plot(T+1:T+60,yF,'r','LineWidth',2);
% h2 = plot(T+1:T+60,upper,'k--','LineWidth',1.5);
% plot(T+1:T+60,lower,'k--','LineWidth',1.5)
% xlim([0,T+60])
% title('Forecast and 95% Forecast Interval')
% legend([h1,h2],'Forecast','95% Interval','Location','NorthWest')
% hold off
%% estimate duration parameter
% seq = featureVector';
% pUnique = unique(seq);
% % Psi200.stateSeq(5).z;
% NN = length(seq);
% 
% newSeq = zeros(size(seq));
% for idx = 1:length(pUnique)
%     findIndex = find(seq == pUnique(idx));
%     newSeq(1,findIndex) = idx;
% end
% 
% % Psi200.stateSeq(5).z;
% pSizeUnique = size(pUnique);
% durParam = zeros(2,pSizeUnique(2));
% 
% for idx = 1:pSizeUnique(2)
%     aa = diff([1 find(newSeq == idx) NN]);
%     parmhat = nbinfit(aa);
%     durParam(:,idx) = parmhat;
% end
% 
% rndNB = nbinrnd(0.25,0.005);