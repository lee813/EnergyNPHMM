function [ new_state ] = predictState(A, state)
%PREDICTSTATE Summary of this function goes here
%   Detailed explanation goes here

    %ra = randi([1,3],1);
    
    pb = A(state,:);
    
    new_state = find(mnrnd(1,pb) == 1);


end

