uniq = unique(pecan(:,1));
count = histc(pecan(:,1), uniq);
notCompleteIdx = find(count~=722);
ids = uniq(notCompleteIdx)';
for idx = 1:length(ids)
    findIdx = find(pecan(:,1) == ids(idx));
     pecan(findIdx,:) = [];
end
%% ------------------------------------------------- data init---
% data = SeqData();
% % % Loop over 3 different sequences
% % for seqID = 1:2
% %     % Create a (random) matrix that represents all the data for current sequence.
% %     %   has T timesteps, each with D-dimensions of observed measurements
% %     T = 100;
% %     D = 3;
% %     curX = rand(D,T);
% %     % Add the sequence to the collection, and name it "1" or "2" or "3"
% %     data = data.addSeq( curX, num2str(seqID));
% % end
% 
% 
% %6 dishwaser 11 microwave
% house1 = zscore(reddDevices{1,1}([6 11],:),0,2);
% %10 dishwaser 6 microwave 
% house2 = zscore(reddDevices{2,1}([10 6],:),0,2);
% %9 dishwaser 16 microwave
% house3 = zscore(reddDevices{3,1}([9 16],:),0,2);

% %6 dishwaser 11 microwave
% house1 = mapminmax(reddDevices{1,1}([6 11],:));
% %10 dishwaser 6 microwave 
% house2 = mapminmax(reddDevices{2,1}([10 6],:));
% %9 dishwaser 16 microwave
% house3 = mapminmax(reddDevices{3,1}([9 16],:));

% selectedIdx = [93,2242,2974,3456,7536];
% selectPecanData = cell(length(pecanData),1);
% for idx = 1:length(pecanData)
%     tmpData = pecanData{idx,1};
%     c = ismember(tmpData(:,1),selectedIdx);
%     indexes = find(c);
%     selectedData = tmpData(indexes,2);
%     selectPecanData{idx,1} = selectedData;
%     selectedData = tmpData(indexes,3);
%     selectPecanDatamatlab{idx,2} = selectedData;
% end


% uniq = unique(ovenandmic(:,1));
% count = histc(ovenandmic(:,1), uniq);
% 1 air 2 cloth 3 oven
consumationIdx = 3;


% firstData = selectPecanData{consumationIdx,1};
% secondData = selectPecanData{consumationIdx,2};
% data3 = selectPecanData{2,1};
% data4 = selectPecanData{2,2};
% data5 = selectPecanData{1,1};
% % data6 = selectPecanData{1,2};
% 
firstData = pecan(:,2);
secondData = pecan(:,3);
data3 = pecan(:,4);
data4 = pecan(:,5);
data5 = pecan(:,6);
data6 = pecan(:,7);


data = SeqData();

stepDay = 30;
stepLength = stepDay * 24 + 2;
userCount = 34; 
userData = cell(userCount,1);
%phouseTmp= zeros(6,722);

for idx = 0:userCount-1
    fprintf('index %d \n',idx * stepLength + 1 );
    phouseTmp(1,:) = mapminmax(firstData(idx * stepLength + 1:(idx + 1) * stepLength)');
    phouseTmp(2,:) = mapminmax(secondData(idx * stepLength + 1:(idx + 1) * stepLength)');
    phouseTmp(3,:) = mapminmax(data3(idx * stepLength + 1:(idx + 1) * stepLength)');
    phouseTmp(4,:) = mapminmax(data4(idx * stepLength + 1:(idx + 1) * stepLength)');
    phouseTmp(5,:) = mapminmax(data5(idx * stepLength + 1:(idx + 1) * stepLength)');
    phouseTmp(6,:) = mapminmax(data6(idx * stepLength + 1:(idx + 1) * stepLength)');
    % unnormalized
%     phouseTmp(1,:) = firstData(idx * stepLength + 1:(idx + 1) * stepLength)';
%     phouseTmp(2,:) = secondData(idx * stepLength + 1:(idx + 1) * stepLength)';
%     phouseTmp(3,:) = data3(idx * stepLength + 1:(idx + 1) * stepLength)';
%     phouseTmp(4,:) = data4(idx * stepLength + 1:(idx + 1) * stepLength)';
%     phouseTmp(5,:) = data5(idx * stepLength + 1:(idx + 1) * stepLength)';
%     phouseTmp(6,:) = data6(idx * stepLength + 1:(idx + 1) * stepLength)';
    %houseTmp = mapminmax(phouseTmp);
    userData{idx+1,1} = phouseTmp;
    data = data.addSeq(phouseTmp);
end

% 
% 
% figure( 'Units', 'normalized', 'Position', [0.1 0.25 0.75 0.5] );
% for plotIdx = 1:userCount
%     subplot(userCount, 1, plotIdx );
%     plotData( data, plotIdx );
% end


%% -------------------------------------------------   RUN MCMC INFERENCE!
modelP = {'bpM.gamma', 0.1, 'bpM.c', 0.1}; 
algP   = {'Niter', 200, 'HMM.doSampleHypers',0,'BP.doSampleMass',0,'BP.doSampleConc',0}; 
% Start out with just one feature for all objects
initP  = {'F.nTotal', 1}; 
CH = runBPHMM( data, modelP, {1, 1}, algP, initP );
% 
% 
Psi200 = CH.Psi( CH.iters.Psi == 200 );
plotStateSeq( Psi200);
% % 
% % % 
figure( 'Units', 'normalized', 'Position', [0.5 0.5 0.5 0.5] );
plotEmissionParams( Psi200 );
title( 'Theta (@ iter 90)', 'FontSize', 20 );
% 

% alignedPsi200 = alignPsiToTruth_OneToOne( Psi200, data );
% % Estimated feature matrix F
% figure( 'Units', 'normalized', 'Position', [0 0.5 0.5 0.5] );
% plotFeatMat( alignedPsi200 )
% 
% figure;
% plot([-0.990852777426279;-0.973289791504751;-0.996797098594943;-0.769952337627393;-0.996435389822202;-0.001584525839192]);
% 
% figure;
% plot([[-0.993876714691710;-0.977922890807884;-0.979547874401020;-0.980160898376687;-0.384903698038956;-0.434124676572104]]);




