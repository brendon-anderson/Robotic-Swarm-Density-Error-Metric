function [ stop ] = saveIteration( ~,optimValues,~,file )
%saveIteration Returns fmincon iteration data
%   Saves the iteration number, objective function value (error), and
%   decision variable vector (robot positions)
    
    stop = false;               % do not stop the iteration
    iter = optimValues.iteration;   % current iteration number
    iterationGap = 1;         % store data every iterationGap iterations
    
    if mod(iter,iterationGap) == 0	% check if iter is a multiple of gap
        fprintf(file,'%d\t%d\r\n',[iter,optimValues.fval]);
    end

end

