% Gait_Para_All_excel=[{'filename'},'stride_time','stride_time','stride_time','stride_time','stride_time','stride_time'];
% Gait_Para_All_excel(1,:)=[{'filename'},'stride_time','stride_length','step_length','stride_speed','stance_phase %','swing_phase%'];
% Gait_Para_All_excel(1,:)=[{'filename'},2,3,4,5,6];
% Gait_Para_All_excel=[1,2];
 [file,path]=uigetfile('*.xls','MultiSelect','on');
% [file,path]=uigetfile('*.xlsx')

for KA=1:length(file)
    try
    class(file(KA))
    clear Gait_Para_All_excel
   file_name= cell2mat(file(KA));
   class(file_name)
%    w=xlsread(file_name,'Sheet1');
    try
    w=xlsread(file_name,erase(file_name,".xls"));
    catch 
    w=xlsread(file_name,'Sheet1');    
    end
   A_functtion=function_knee_flexion_Automated_Right(w);
  recycle on % Send to recycle bin instead of permanently deleting.
delete(file_name);
%    file_name=cellstr(file_name)
   class(file_name)
%    file_name=char(mat2cell(file_name,1));
   Gait_Para_All_excel(:,1:2)=A_functtion;
   xlswrite(file_name,Gait_Para_All_excel);
%      Gait_Para_All_excel(KA,:)=[x,y];
     catch
    end
   
end


% w=xlsread(file,'Sheet1')
% [p1,p2,p3,p4,p5,p6]=Gait_function_all_code(w);
% Gait_Para_All_excel=[p1,p2,p3,p4,p5,p6];
% xlswrite('Gait_Para_All_excel.xls',Gait_Para_All_excel);

    
   