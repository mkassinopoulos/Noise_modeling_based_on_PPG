function IR_all = func_CPM_cos(Ts, memory, M)

t_win= 0 :Ts:memory; nIR = length(t_win); 

IR_all = zeros(nIR,M*2);
for m = 1:M    
    IR_all(:,(m-1)*2+1) = cos(m*2*pi()*t_win/memory)-1;
    IR_all(:,(m-1)*2+2) = sin(m*2*pi()*t_win/memory);
end


%%  -------------------------------