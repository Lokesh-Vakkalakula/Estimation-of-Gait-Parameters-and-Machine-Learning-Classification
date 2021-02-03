Gait_Para_All_excel=[{'filename'},'stride_time','stride_time','stride_time','stride_time','stride_time','stride_time','cadence'];
Gait_Para_All_excel(1,:)=[{'filename'},'stride_time','stride_length','step_length','stride_speed','stance_phase %','swing_phase%','cadence'];
% Gait_Para_All_excel(1,:)=[{'filename'},2,3,4,5,6];
% Gait_Para_All_excel=[1,2];
 [file,path]=uigetfile('*.xls','MultiSelect','on');
% [file,path]=uigetfile('*.xls') 
delete('Gait_Para_All_excel.xlsx')
for KA=1:length(file)
   try
    class(file(KA))
   file_name= cell2mat(file(KA));
   class(file_name)
%    w=xlsread(file_name,'sheet1');
%     w=xlsread(file_name,erase(file_name,".xls"));

 try
    w=xlsread(file_name,erase(file_name,".xls"));
    catch
    w=xlsread(file_name,'Sheet1');    
    end
   [p1,p2,p3,p4,p5,p6,p7]=Gait_function_all_code3_Rightside(w);
%    [q1,q2,q3,q4,q5,q6,q7]=Gait_function_all_code3_Leftside(w);
%    p1=(p1+q1)/2;
%    p2=(p2+q2)/2;
%    p3=(p3+q3)/2;
%    p4=(p4+q4)/2;
%    p5=(p5+q5)/2;
%    p6=(p6+q6)/2;
%    p7=(p7+q7)/2;
   file_name=cellstr(file_name)
   class(file_name)
%    file_name=char(mat2cell(file_name,1));
   Gait_Para_All_excel(KA+1,:)=[file_name,p1,p2,p3,p4,p5,p6,p7];
%      Gait_Para_All_excel(KA,:)=[x,y];
   clear p1,p2,p3,p4,p5,p6,p7;
   catch
  end
end


% w=xlsread(file,'Sheet1')
% [p1,p2,p3,p4,p5,p6]=Gait_function_all_code(w);
% Gait_Para_All_excel=[p1,p2,p3,p4,p5,p6];
xlswrite('Gait_Para_All_excel.xls',Gait_Para_All_excel);
Gait_Para_All_excel_modify=1
[ra,rb]=size(Gait_Para_All_excel)
Gait_Para_All_excel_modify(ra,rb)=1
for i=1:1:ra
Gait_Para_All_excel_modify(i,:)=Gait_Para_All_excel(ra-i-1,:);
end
   