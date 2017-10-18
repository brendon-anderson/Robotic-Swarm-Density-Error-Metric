function parsave(delta,decisionVector,N,error,X,Y,rho_N,xmin,xmax,ymin,ymax,intermediateData,timeElapsed)
    % This function saves the workspace after each robot swarm is
    % optimized.
    fname = ['optimizationWorkspaceN',num2str(N),'.mat'];
    save(fname);
end
