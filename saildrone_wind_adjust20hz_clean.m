function saildrone_wind_adjust20hz_clean(sdnum,tc_name,pth)


loca=['lhf_input_' int2str(sdnum) tc_name '_qcwind.txt'];
    
coare_output=coare35vn(loca,pth);
ws_10m_nostab=coare_output(:,29);   %wind speed adjusted to 10 m without stability correction 
ws_10m=coare_output(:,30);          %wind speed adjusted to 10 m height with stability correction (U10N - neutral stability)
cd_10m=coare_output(:,31);
    
save([pth 'ws_10m_sd' int2str(sdnum) tc_name '_cdadj_zo25mps_qcwind'],'ws_10m','ws_10m_nostab');    %adjust wind speed for leveling of Zo, Cd in coare35vn.m

