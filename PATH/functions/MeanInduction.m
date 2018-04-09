function [a,a_p] = MeanInduction(A_axial,A_tangential)

% This function time averages the induction factors.
% Out of range values are set to NaN.
 %% START WHILE ITERATION LOOP HERE   
    
        % TIDY UP TIME !

        
        % FIND POORLY CONVERGED VALUES AND SET TO NAN
        InRange=((abs(A_axial)>0.999)|(A_axial<0.0001));
        A_axial(InRange)=NaN;

        InRange2=((abs(A_tangential)>0.999)|(A_tangential<0.0001));
        A_tangential(InRange2)=NaN;
        
        % CALCULATE THE MEAN IGNORING NaNs
        a=nanmean(A_axial,2);
        a_p=nanmean(A_tangential,2);
        
        % JUST INCASE .......nanmean([NaN,NaN,NaN])= NaN!!!!
        % INTERPOLATE THE BAD NaNs
        
        nanx = isnan(a);
        tx    = 1:numel(a);
        a(nanx) = interp1(tx(~nanx), a(~nanx), tx(nanx),'spline');
 
        nany=  isnan(a_p);
        ty    = 1:numel(a_p);
        a_p(nany) = interp1(ty(~nany), a_p(~nany), ty(nany),'spline');
        
end