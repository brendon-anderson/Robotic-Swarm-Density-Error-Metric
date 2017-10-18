function [ stop ] = saveIteration( x,optimValues,~,file,timeLim,firstIter )
%   saveIteration returns fmincon iteration data.
%   Saves the iteration number, objective function value (error), and
%   decision variable vector (robot positions) to a txt file.
    
    stop = false;                   % do not stop the iteration
    iter = optimValues.iteration;   % current iteration number
    iterationGap = 1;               % store data every iterationGap iterations
    
    fileX = '';
    for i = 1:length(x)
        fileX = [fileX,'\t%d'];
    end
    
    if (iter > 0) && (mod(iter,iterationGap) == 0)	% check if iter is a multiple of gap
        fprintf(file,['%d\t%d',fileX,'\r\n'],[iter+firstIter,optimValues.fval,x]);
    end
    
    timeElapsed = toc;
    if timeElapsed > timeLim
        stop = true;
    end

end

