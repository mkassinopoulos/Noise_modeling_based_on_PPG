function  [CC_all, bounds_all] = estimate_CC(input, output, N_IR)

CC_all = zeros(N_IR,4);
bounds_all = zeros(2,4);
for i = 1 : size(output,2)
    x = output(:,i);
    [val, ~, bounds] = crosscorr(input,x,(N_IR-1)/2);  
    bounds_all(:,i) = bounds;    
    CC_all(:,i) = val;
end
end

