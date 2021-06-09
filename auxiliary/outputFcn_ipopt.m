function [ continue_ ] = outputFcn_ipopt( x, fvalue,state, timeStart )
%OUTPUTFCN_IPOPT: This function is to check the time spend in the 
% optimization function

continue_ = true;

tEnd = toc(timeStart);

if tEnd > 3600%1000
    
    continue_ = false;
    
end



end

