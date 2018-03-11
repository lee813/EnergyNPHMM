NN = 34;
seqs = zeros(NN, 722);
idx = 1;

while(idx < NN)
 
    seqs(idx,:) = Psi200.stateSeq(idx).z;
    idx = idx + 1;
            
end 

%% ----Calculate similarity


s1y = pdist(seqs,'minkowski');

