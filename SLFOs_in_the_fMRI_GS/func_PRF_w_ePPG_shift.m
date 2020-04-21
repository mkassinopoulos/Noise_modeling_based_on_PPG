function [obj_function,CRF_sc,RRF_sc,BRF_sc, r_PRF_sc,yPred, yPred_card, yPred_resp, yPred_PPG]=func_PRF_w_ePPG_shift(P,Ts_10,HR,RF,ePPG, ind_BOLD_10,GS,only_OF,ind_train, ind_valid, filt_b, filt_a,shift)

t1c=P(1) ; d1c=P(2);
t2c=P(3);  d2c=P(4);
t1r=P(5); d1r = P(6);
t2r=P(7); d2r=P(8);
t1p=P(9); d1p = P(10);
t2p=P(11); d2p=P(12);


NV = length(GS);
t_win= 0 :Ts_10:60;

a1= sqrt(t1c)/d1c; a2= sqrt(t1c)*d1c ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardF=IR/max(IR);  
a1= sqrt(t2c)/d2c; a2= sqrt(t2c)*d2c ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardS=IR/max(IR);     

a1= sqrt(t1r)/d1r; a2= sqrt(t1r)*d1r ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respF=IR/max(IR);   
a1= sqrt(t2r)/d2r; a2= sqrt(t2r)*d2r ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respS=IR/max(IR);  

a1= sqrt(t1p)/d1p; a2= sqrt(t1p)*d1p ;
IR = t_win.^a1.*exp(-t_win/a2); IR_ppgF=IR/max(IR);   
a1= sqrt(t2p)/d2p; a2= sqrt(t2p)*d2p ;
IR = t_win.^a1.*exp(-t_win/a2); IR_ppgS=IR/max(IR);  


r_PRF_sc=zeros(3,1);

HR = zscore(HR); RF = zscore(RF);
HR_Fconv=conv(HR,IR_cardF); HR_Fconv_MR=HR_Fconv(ind_BOLD_10);
HR_Sconv=conv(HR,IR_cardS); HR_Sconv_MR=HR_Sconv(ind_BOLD_10);
RF_Fconv=conv(RF,IR_respF); RF_Fconv_MR=RF_Fconv(ind_BOLD_10);
RF_Sconv=conv(RF,IR_respS); RF_Sconv_MR=RF_Sconv(ind_BOLD_10);

PPG_Fconv=conv(ePPG,IR_ppgF); PPG_Fconv_MR=PPG_Fconv(ind_BOLD_10+shift*10);
PPG_Sconv=conv(ePPG,IR_ppgS); PPG_Sconv_MR=PPG_Sconv(ind_BOLD_10+shift*10);

% PPG_MR = ePPG(ind_BOLD_10);

regr = [HR_Fconv_MR,HR_Sconv_MR,RF_Fconv_MR,RF_Sconv_MR, PPG_Fconv_MR, PPG_Sconv_MR];  regr = detrend(regr,'linear');
% regr = [HR_Fconv_MR,HR_Sconv_MR,RF_Fconv_MR,RF_Sconv_MR, PPG_MR ];  regr = detrend(regr,'linear');

regr = filtfilt(filt_b,filt_a,regr ) ;

regr = [regr, ones(NV,1)];

B = regr(ind_train,:)\GS(ind_train);     yPred = regr*B;

obj_function = 1 - corr(yPred(ind_train),GS(ind_train)) ;
    

if only_OF == 0
    CRF_sc = B(1) * IR_cardF + B(2) * IR_cardS;     CRF_sc = CRF_sc/max(abs(CRF_sc));
    RRF_sc = B(3)*IR_respF + B(4)*IR_respS; RRF_sc = RRF_sc/max(abs(RRF_sc));
    BRF_sc = B(5)*IR_ppgF + B(6)*IR_ppgS; BRF_sc = BRF_sc/max(abs(BRF_sc));
        
%     HR_conv = B(1)*HR_Fconv + B(2)*HR_Sconv;    HR_conv = HR_conv(1:length(HR));   HR_conv_MR = HR_conv(ind_BOLD_10);
%     RF_conv = B(3)*RF_Fconv + B(4)*RF_Sconv;     RF_conv = RF_conv(1:length(RF));     RF_conv_MR = RF_conv(ind_BOLD_10);
%     PPG_conv = B(5)*PPG_Fconv + B(6)*PPG_Sconv;     PPG_conv = PPG_conv(1:length(RF));     PPG_conv_MR = PPG_conv(ind_BOLD_10);
        
    r_PRF_sc(1) = corr(yPred(ind_valid),GS(ind_valid));
    yPred_card = regr(:,1:2)*B(1:2);  r_PRF_sc(2) = corr(yPred_card(ind_valid),GS(ind_valid));
    yPred_resp = regr(:,3:4)*B(3:4);  r_PRF_sc(3) = corr(yPred_resp(ind_valid),GS(ind_valid));
    
%     GS_clean = GS - regr(:,1:4)*B(1:4);        
    yPred_PPG = regr(:,5:6)*B(5:6);  r_PRF_sc(4) = corr(yPred_PPG(ind_valid),GS(ind_valid));
    
end




%%   ----------------------------------------------------