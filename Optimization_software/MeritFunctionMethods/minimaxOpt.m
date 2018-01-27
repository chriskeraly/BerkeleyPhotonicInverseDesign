

% MINIMAX OPTIMIZATION
function[deltaX,deltaFpred] = minimaxOpt(Fmin,dFdx,lb,ub,c,options)

f = @(x)( -(Fmin+dFdx*x) ); % function to minimize

[~,minInd] = min(Fmin);
x0 = dFdx(minInd,:).'; % initial guess = dFdx of worst case

[deltaX,fNewPred] = fminimax(f,x0,[],[],[],[],lb,ub,c,options);

deltaFpred = -fNewPred - Fmin;
deltaX = deltaX.';

end