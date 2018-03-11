sizeOfTheta = size(Psi200.theta);
datax = mapminmax(data2965')';
%data2510;
%data2965;
%data7982;
dim = 6;
sizeofDatax = size(datax);
llTable = zeros(sizeofDatax(1,1),sizeOfTheta(1,2));


for idx = 1:sizeOfTheta(1,2)
    mean = Psi200.theta(idx).mu;
    invSigma = Psi200.theta(idx).invSigma;
    ll = 0;
    llTable(:,idx) = multiVNormal(datax,mean',inv(invSigma));
   
end

featureVector = zeros(sizeofDatax(1,1),1);
for dataLengthIdx = 1:sizeofDatax(1,1)
    [A,B] = max(llTable(dataLengthIdx,:));
    featureVector(dataLengthIdx,1) = B;
end


figure;
plot(featureVector, 'r*');




