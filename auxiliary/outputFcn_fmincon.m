function [ stop ] = outputFcn_fmincon( x, optimValues,state, timeStart )
%OUTPUTFCN_FMINCON: This function is to check the time spend in the 
% optimization function

stop = false;

tEnd = toc(timeStart);

if tEnd > 1000
    
    stop = true;
    
end



end

