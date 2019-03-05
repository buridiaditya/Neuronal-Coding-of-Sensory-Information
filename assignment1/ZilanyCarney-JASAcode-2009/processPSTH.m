function procPSTH = processPSTH(psth, binwidth, overlap)
    indices = 1:(binwidth-overlap):length(psth);
    procPSTH = zeros(1,length(indices));
    
    for i=1:length(indices)
        if indices(i)+binwidth-1 > length(psth)
            procPSTH(i) = sum(psth(indices(i):end));
        else        
            procPSTH(i) = sum(psth(indices(i):(indices(i)+binwidth-1)));
        end
    end
    
end