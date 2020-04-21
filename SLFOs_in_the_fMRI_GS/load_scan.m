
function [Ts,Fs,TR, trig,PPGlocs, HR, resp, DVARS, cardiac, GS, FD, RF, BR,movRegr] = load_scan(subject,task,baseDir)


filepath_MRacq=[baseDir,'/Physio/',subject,'_',task,'/phys.mat'];
filepath_movRegr=[baseDir,'/Physio/',subject,'_',task,'/Movement_Regressors_dt.txt'];

load([baseDir,'/Physio/',subject,'_',task,'/Phys_sum.mat']);
load(filepath_MRacq,'Fs','trig','TR','cardiac'); Ts = 1/Fs;

HR = HRV(:);
RF = IHF;
filepath_input=[baseDir,'Atlas/',subject,'_',task,'/'];
load([filepath_input,'TissueBasedRegressors_1199.mat'],'DVARS','WB','FD')
GS = WB.MA;

movRegr=load(filepath_movRegr);  movRegr=[movRegr, movRegr.^2];
% movRegr = zscore(movRegr);

time = 0:Ts:(length(resp)-1)*Ts;
time_10 = 0:0.1:time(end);
indPPGlocs = zeros(size(PPGlocs));
for i = 1 :length(PPGlocs)
    [val loc] = min(abs(PPGlocs(i)-time));
    indPPGlocs(i) = loc;
end
pks = cardiac(indPPGlocs);
uePPG = interp1([0,time(indPPGlocs),time(end)],[pks(1),pks',pks(end)],time_10);
% cardiac = zscore(uePPG); 
cardiac = uePPG;
cardiac = cardiac(:);

resp = interp1(time,resp,time_10);

end




