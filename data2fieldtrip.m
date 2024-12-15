function data = data2fieldtrip(y,elecLabel,fsample,t)
    % Converts data into a fieldtrip friendly format.
    data = struct ; 
    data.fsample = fsample ; 
    data.sampleinfo = [1,length(y)] ; 
    data.trial{1} = y ; 
    data.time{1} = linspace(0,20,length(y)) ;
    data.label = elecLabel ; 
end