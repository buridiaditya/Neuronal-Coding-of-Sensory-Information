function procPSTH = processPSTH(psth, binwidth, overlap)
    indices = 1:(binwidth-overlap):length(psth);
    procPSTH = zeros(1,length(indices));
    
    for i=1:length(indices)
        if indices(i)+binwidth > length(psth)
            procPSTH(i) = sum(psth(indices(i):end));
        else        
            procPSTH(i) = sum(psth(indices(i):1:(indices(i)+binwidth)));
        end
    end
    
end