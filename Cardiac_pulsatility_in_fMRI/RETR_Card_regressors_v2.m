function [Regr]=RETR_Card_regressors_v2(time,locsECG,M)
%% ===========================
% RETROICOR
% voxel_f=voxel; locsECG=ECGlocs;

NV=length(time);
Phi=zeros(NV,1);
for i=1:NV    
    t=time(i);    
    [~,minI]=min(abs(locsECG-t));
    
    minOnLeft=t-locsECG(minI)>0;
    if (minI==1 && ~minOnLeft)
        t2=locsECG(minI);
        t1=t2-1;        
    elseif (minI==length(locsECG) && minOnLeft)
        t1=locsECG(minI);
        t2=t1+1;
    elseif minOnLeft       
        t1=locsECG(minI);
        t2=locsECG(minI+1);
    else
        t1=locsECG(minI-1);
        t2=locsECG(minI);
    end        
    Phi(i)=2*pi*(t-t1)/(t2-t1);    
end
    
Regr=zeros(NV,M*2);
for i=1:M
    Regr(:,(i-1)*2+1)=cos(i*Phi);
    Regr(:,i*2)=sin(i*Phi);
end

% Regr=[zeros(1,M*2);diff(Regr)];



