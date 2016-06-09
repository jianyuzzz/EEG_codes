function [d dp df dfp] = feature_eval_nParametric(c1,c2)

    % discart all data that exceeds distribution threshold
    m1 = mean(c1);
    fwhm1 = sqrt(2*log(2))*std(c1);
    m2 = mean(c2);
    fwhm2 = sqrt(2*log(2))*std(c2);
    
    c1p = c1((c1>(m1-fwhm1)) & (c1<(m1+fwhm1)));
    c2p = c2((c2>(m2-fwhm2)) & (c2<(m2+fwhm2)));
    
   
    c1min = min(c1);
    c1max = max(c1);
    c2min = min(c2);
    c2max = max(c2);
    
    if (c1max>=c2max && c1min<=c2min) || (c2max>=c1max && c2min<=c1min) 
        ndc1 = 0;
        ndc2 = 0;
    else
        ndc1 = sum((c1>c2max) | (c1<c2min));
        ndc2 = sum((c2>c1max) | (c2<c1min));
    end
    
    d = (ndc1+ndc2)/(length(c1)+length(c2));

    c1min = min(c1p);
    c1max = max(c1p);
    c2min = min(c2p);
    c2max = max(c2p);
    
    if isempty(c1p) || isempty(c2p)
        ndc1 = 0;
        ndc2 = 0;   
    else
    if (c1max>=c2max && c1min<=c2min) || (c2max>=c1max && c2min<=c1min) 
        ndc1 = 0;
        ndc2 = 0;
    else
        ndc1 = sum((c1p>c2max) | (c1p<c2min));
        ndc2 = sum((c2p>c1max) | (c2p<c1min));
    end
    end
    
    dp = (ndc1+ndc2)/(length(c1p)+length(c2p));
    
    
    % Fisher score
    df = (mean(c1)-mean(c2))^2/(var(c1)+var(c2));
    dfp = (mean(c1p)-mean(c2p))^2/(var(c1p)+var(c2p));

end