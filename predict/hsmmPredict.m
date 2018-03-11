function predictSeq = hsmmPredict(durParams,AA, NN)
        idx = 1;
        predictSeq = zeros(NN,1);
        predictSeq(1) = 1;
        while(idx < NN)
            if idx - 1 == 0
                state = 1;
            else
                state = predictState(AA,predictSeq(idx-1));
            end
            durCount = randomNBDur(durParams(state,:));
            predictSeq(idx:idx+durCount) = state;
            %state = predictState(AA,state);
            idx = idx + durCount + 1;
        end
end


function durCount = randomNBDur(durParam)
    durCount = nbinrnd(durParam(1),durParam(2));
end

function durCount = randomPoiDur(durParam)
    durCount = poissrnd(durParam(1));
%     durCount = round(mean(poissrnd(durParam(1),2,1)));
end

%         
% function [state, dur] =  rle(stateseq)
%     %position that state changed
%     pos, = np.where(np.diff(stateseq) != 0)
%     pos = np.concatenate(([0],pos+1,[len(stateseq)]))
%     %sequence value at changed position(till last one - 1), and the duration
%     return stateseq[pos[:-1]], np.diff(pos)
% end