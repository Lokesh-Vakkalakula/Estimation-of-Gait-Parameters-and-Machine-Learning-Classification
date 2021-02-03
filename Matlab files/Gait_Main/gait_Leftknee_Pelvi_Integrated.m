%  w=xlsread('Proband_6_CB1.xls','Proband_6_CB1');
  w=xlsread('Proband490_CB1.xls','Sheet1');
  w=sortrows(w,1);
p=w
%w=w(all(w,2),:)
% w=[w(:,1),w(:,47:61)];
w=[w(:,1),w(:,62:76)];
% p=w
[m,n]=size(w)
for i = 2:m-1
    for j=1:n-1
      if w(i,j)==0
        w(i,j)=(w(i+1,j)+w(i+1,j))/2;
      end 
    end
end
j=0
x(1)=0
% for i = 1:length(w(:,1))-10
%     if -(w(i,1)-w(i+10,1))>35
%         w(i:i+10,1)=nan
%     end
% end

% w(x)
w(w==56)= 0;
% % w(w==0)= NaN;
w=w(all(w,2),:);
% w(:,16)=sgolayfilt(w(:,16),2,11)
% w(:,11)=sgolayfilt(w(:,11),2,11)
% w(:,6)=sgolayfilt(w(:,6),2,11)

[m,n]=size(w);
j=0;
x(1)=0;
j=0;
x(1)=0;


% for i = 1:length(w(:,1))-10
%     if -(w(i,1)-w(i+10,1))>15
%         w(i:i+10,:)=[];
%     end
% end
for i = 1:length(w(:,1))-10
    if -(w(i,1)-w(i+10,1))>15
        w(i:i+10,:)=nan;
    end
end
% for i = 1:length(w(:,1))
%     if w(i,1)==nan
%         w(i)=[]
%     end
% end
% w(any(isnan(w), 2), :) = [];


 for i = 2:m-1
    if (abs(w(i,16)-w(i-1,16)))>6
%        if (w(i+1,1)-w(i,1))>3
        w(i,16)=nan;
%        end
    end
end

AngleNum=[((w(:,4)-w(:,9)).*(w(:,14)-w(:,9)))+((w(:,5)-w(:,10)).*(w(:,15)-w(:,10)))+((w(:,6)-w(:,11)).*(w(:,16)-w(:,11)))];
AngleDenom=[(sqrt(((w(:,4)-w(:,9)).^2)+((w(:,5)-w(:,10)).^2)+((w(:,6)-w(:,11)).^2))).*(sqrt(((w(:,14)-w(:,9)).^2)+((w(:,15)-w(:,10)).^2)+((w(:,16)-w(:,11)).^2)))];
% AngleNum=[((w(:,49)-w(:,54)).*(w(:,59)-w(:,54)))+((w(:,50)-w(:,55)).*(w(:,60)-w(:,55)))+((w(:,51)-w(:,56)).*(w(:,61)-w(:,56)))]
% AngleDenom=[(sqrt(((w(:,49)-w(:,54)).^2)+((w(:,50)-w(:,55)).^2)+((w(:,51)-w(:,56)).^2))).*(sqrt(((w(:,59)-w(:,54)).^2)+((w(:,60)-w(:,55)).^2)+((w(:,61)-w(:,56)).^2)))]

% AngleNum=[((w(:,49)-w(:,54)).*(w(:,59)-w(:,54)))+((w(:,50)-w(:,55)).*(w(:,60)-w(:,55)))]
% AngleDenom=[(sqrt(((w(:,49)-w(:,54)).^2)+((w(:,50)-w(:,55)).^2))).*(sqrt(((w(:,59)-w(:,54)).^2)+((w(:,60)-w(:,55)).^2)))]
Angle = [acosd(AngleNum./AngleDenom)];
y=(smooth(180-real(Angle)));
for i = 1:size(y)
     if y(i)>55
         y(i)=50;
     end
     if y(i)<0
         y(i)=50;
     end 
    
end


% for i = 3:size(y)-3
%      if y(i)==y(i+1)
%          if y(i)==(i-1)
%            y(i)= nan
%          end
%      end
%     
% end
% windowSize = 1; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% A = filter(b,a,x)
x = w(:,1);
% figure(1)
% plot(y)
A = ((sgolayfilt(y,2,11)));
%A=y
TF = islocalmin(A);
TY= islocalmax(A);
figure(1);
plot(x,A,x(TF),A(TF),'r*',x(TY),A(TY),'r^')
% figure(2)
% plot(w(:,1),y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%gait calculation
Nan_index=0;
Min_indices=find(TF);
Max_indices=find(TY);
j=1;
Max_indices_sel=1; 
%to calculate the real primary Maxima
for i=1:length(find(TY))
    if A(Max_indices(i))>20 & ~isnan(x(Max_indices(i)))
        Max_indices_sel(j)= Max_indices(i);
        j=j+1;
        disp(Max_indices(i));
    end
end

%to calculate the Nan break
for i=Max_indices_sel(1):length(x)
    if isnan(x(i))& isnan(x(i+1))&isnan(x(i+2))&isnan(x(i+3))&isnan(x(i+4))
        Nan_index=i-1;
        break
    end
        
end
 
%to remove the max selected peaks which are very near by



Min_indices=find(TF);
Max_indices=find(TY);
j=1;
%to calculate the real primary Maxima
% for i=1:length(find(TY))
%     if A(Max_indices(i))>25 
%         Max_indices_sel(j)= Max_indices(i);
%         j=j+1
%     end
% end

%to calculate the real primary minimas. & ~isnan(x(Max_indices(i)))
k=1;
Min_indices_sel=1;
for i=1:length(Min_indices)
    if A(Min_indices(i))<15 & ~isnan(x(Min_indices(i)))
        
        Min_indices_sel(k)= Min_indices(i);
        k=k+1;
    end
end

%to calculate the real primary maximas that are less than the nan break.
k=1
Max_indices_sel_Before_nan=1
for i=1:length(Max_indices_sel)
    if Max_indices_sel(i)<Nan_index
        Max_indices_sel_Before_nan(k)=Max_indices_sel(i);
        k=k+1;
    end

end

%removing the maximaas that do not have minimas on left side
r=1
Max_indices_sel_Before_nan_non_extremas=1
for i=1:length(Max_indices_sel_Before_nan)
    count=0;
    for j=1:length(Min_indices_sel)
    if x(Min_indices_sel(j))<x(Max_indices_sel_Before_nan(i))
        count=count+1;
    end
    end
    if count>=1
        Max_indices_sel_Before_nan_non_extremas(r)=Max_indices_sel_Before_nan(i);
        r=r+1;
    end
        
end
%removing the maximaas that do not have minimas on right side
r=1
Max_indices_sel_Before_nan_non_extremas_all=1
for i=1:length(Max_indices_sel_Before_nan_non_extremas)
    count=0
    for j=1:length(Min_indices_sel)
    if x(Max_indices_sel_Before_nan_non_extremas(i))<x(Min_indices_sel(j)) & x(Min_indices_sel(j))<x(Max_indices_sel_Before_nan_non_extremas(i))+21
        count=count+1;
        disp(x(Min_indices_sel(j)));
    end
    end
    if count>=1
        Max_indices_sel_Before_nan_non_extremas_all(r)=Max_indices_sel_Before_nan_non_extremas(i);
        r=r+1
    end
        
end

%removing the multiple maximas of a single gait 
M=1
M(1)=Max_indices_sel_Before_nan_non_extremas_all(1)
j=2;
for i=2:length(Max_indices_sel_Before_nan_non_extremas_all)
    if (x(Max_indices_sel_Before_nan_non_extremas_all(i))-x(Max_indices_sel_Before_nan_non_extremas_all(i-1)))>10
       M(j)= Max_indices_sel_Before_nan_non_extremas_all(i);
       j=j+1;
    end
end
        
%to remove the gait which do not have primary maximas on left
 r=1
 M_final=1
 for i=1:length(M)
     count=0
     for j=1:length(Max_indices_sel_Before_nan)
     if x(Max_indices_sel_Before_nan(j))< x(M(i))
         count=count+1;
         Gait_Left_Max=Max_indices_sel_Before_nan(j);
     end
     end
     if count>=1
         M_final(r)=M(i);
         r=r+1;
     end
         
 end
 
 
 
 Gait_left_min_index=1
 Gait_left_max_index=1
 Gait_right_min_index=1
 %%%calculating multiple Gaits
 for i=1:length(M_final)
     %calculate Left Primary Maximas
    for j=1:length(Max_indices_sel_Before_nan)
        if(Max_indices_sel_Before_nan(j)==M_final(i))
            if(length(Max_indices_sel_Before_nan)>1)
            Gait_left_max_index(i)=Max_indices_sel_Before_nan(j-1);
            end
        end
    end
        %calculate the left Primary Minima
        for k=1:length(Min_indices_sel)
            if(Min_indices_sel(k)>Gait_left_max_index(i))
                Gait_left_min_index(i)=Min_indices_sel(k);
                break;
            end
        end
        %calculate the left Primary Minima
        for l=1:length(Min_indices_sel)
            if(Min_indices_sel(l)>M_final(i))
                Gait_right_min_index(i)=Min_indices_sel(l);
                break;
            end
        end
        if i==2 
            break
        end
            
 end
 
 
%  for i=1:length(Gait_left_max_index)
%      stride_time(i)=(x(Gait_right_min_index(i))-x(Gait_left_min_index(i)))/25;
% %  stride_length(i)=w(Gait_left_min_index(i),16)-w(Gait_right_min_index(i),16);
%  stride_length(i)=sqrt((w(Gait_left_min_index(i),16)-w(Gait_right_min_index(i),16))^2+(w(Gait_left_min_index(i),15)-w(Gait_right_min_index(i),15))^2+(w(Gait_left_min_index(i),14)-w(Gait_right_min_index(i),14))^2);
%  step_length(i)=stride_length(i)/2;
%  stride_speed(i)=stride_length(i)/stride_time(i);
%  
%  end
%  
% p1= stride_time(1)
% p2= stride_length(1)
% p3=step_length(1)
% p4=stride_speed(1)
%stance pahse and swing phase
count=0;
count_final=1;
index=1;
 Gait_index_start=Gait_left_min_index(1)
 Gait_index_end=Gait_right_min_index(1)
%  for in=Gait_index_start:Gait_index_end
%      count=count+1;
%  if y(in)>30
%      index=in;
%      count_final=count;
%      break;
%  end
%  end
%  stance_phase=0;
%  swing_phase=0;
%  if (Gait_index_end-Gait_index_start)>0
%  stance_phase=(count_final/(Gait_index_end-Gait_index_start))*100
%  swing_phase=100-stance_phase
%  end

 


% &########################################pelvic tilt inclusion

w_pelvic=p;
w_pelvic=[w_pelvic(:,1),w_pelvic(:,7:11),w_pelvic(:,42:51)];
[m,n]=size(w_pelvic);
for i = 2:m-1
    for j=1:n-1
      if w_pelvic(i,j)==0
        w_pelvic(i,j)=(w_pelvic(i+1,j)+w_pelvic(i+1,j))/2;
      end
    end
end
w_pelvic(w_pelvic==56)= 0;
w_pelvic=w_pelvic(all(w_pelvic,2),:);
figure(2);
AngleNum_pelvic=[((w_pelvic(:,2)-w_pelvic(:,7)).*(w_pelvic(:,12)-w_pelvic(:,7)))+((w_pelvic(:,3)-w_pelvic(:,8)).*(w_pelvic(:,13)-w_pelvic(:,8)))];
AngleDenom_pelvic=[(sqrt(((w_pelvic(:,2)-w_pelvic(:,7)).^2)+((w_pelvic(:,3)-w_pelvic(:,8)).^2))).*(sqrt(((w_pelvic(:,12)-w_pelvic(:,7)).^2)+((w_pelvic(:,13)-w_pelvic(:,8)).^2)))];
Angle_pelvic = [acosd(AngleNum_pelvic./AngleDenom_pelvic)];
y_pelvic=90-Angle_pelvic;
for i=1:length(Angle_pelvic)
    if abs(y_pelvic(i))>15
        y_pelvic(i)=NaN;
    end
end
A_pelvic =(movmean(sgolayfilt(y_pelvic,2,11),2));
 figure(2);
plot(w_pelvic(:,1),A_pelvic);
x_pelvic=w_pelvic(:,1);

Gait_left_min_index_pelvic=find(x_pelvic==x(Gait_index_start));
Gait_right_min_index_pelvic=find(x_pelvic==x(Gait_index_end));
% array_gait_pelvic_tilt=A_pelvic(Gait_left_min_index_pelvic+15:Gait_right_min_index_pelvic-5);
% array_frame_gait=x_pelvic(Gait_left_min_index_pelvic:Gait_right_min_index_pelvic);

% [val, idx] = min(array_gait_pelvic_tilt);
% transition_frame=idx+x_pelvic(Gait_left_min_index_pelvic+16)-1;
% transition_frame_index_knee=find(x==transition_frame);
Gait_index_difference=Gait_right_min_index_pelvic-Gait_left_min_index_pelvic;
Gait_index_difference_left=Gait_left_min_index_pelvic+Gait_index_difference*0.55;
Gait_index_difference_right=Gait_left_min_index_pelvic+Gait_index_difference*0.75;
TF_pelvic = islocalmin(A_pelvic);
TF_pelvic_index=find(TF_pelvic);
j=1;
pelvic_index_Gait_Minimas=1;
for i=1:length(TF_pelvic_index)
if Gait_index_difference_left<TF_pelvic_index(i) & Gait_index_difference_right>TF_pelvic_index(i)
    pelvic_index_Gait_Minimas(j)=TF_pelvic_index(i);
    j=j+1
end
end
pelvic_index_Gait_Minimas=min(pelvic_index_Gait_Minimas);
transition_frame=x_pelvic(pelvic_index_Gait_Minimas)
transition_frame_index_knee=find(x==transition_frame);
transition_frame_index_knee_p=find(p(:,1)==transition_frame)
figure(2);
plot(x_pelvic,A_pelvic,x_pelvic(TF_pelvic),A_pelvic(TF_pelvic),'r*');


 
stance_phase=0;
swing_phase=0;
if (Gait_index_end-Gait_index_start)>0
   
    stance_phase=((transition_frame_index_knee-Gait_index_start)/(Gait_index_end-Gait_index_start))*100
    swing_phase=100-stance_phase
end
 if stance_phase>50 & swing_phase>25
p5=stance_phase;
p6=swing_phase;
 end
 p_footlength=sqrt((p(:,124)-p(:,119)).^2+(p(:,125)-p(:,120)).^2+(p(:,126)-p(:,121)).^2);
%  footlength=median(p_footlength);
 footlength=0;
%  footlength=sqrt((p(transition_frame_index_knee_p,124)-p(transition_frame_index_knee_p,119))^2+(p(transition_frame_index_knee_p,125)-p(transition_frame_index_knee_p,120))^2+(p(transition_frame_index_knee_p,126)-p(transition_frame_index_knee_p,121))^2)
 for i=1:length(Gait_left_max_index)
     stride_time(i)=(x(Gait_right_min_index(i))-x(Gait_left_min_index(i)))/25;
%  stride_length(i)=w(Gait_left_min_index(i),16)-w(Gait_right_min_index(i),16);
 stride_length(i)=sqrt((w(Gait_left_min_index(i),16)-w(Gait_right_min_index(i),16))^2+(w(Gait_left_min_index(i),15)-w(Gait_right_min_index(i),15))^2+(w(Gait_left_min_index(i),14)-w(Gait_right_min_index(i),14))^2)+80;
 step_length(i)=stride_length(i)/2;
 stride_speed(i)=stride_length(i)/stride_time(i);
 
 end
 
p1= stride_time(1)
p2= stride_length(1)
p3=step_length(1)
p4=stride_speed(1)/100
p_footlength=sqrt((p(:,124)-p(:,119)).^2+(p(:,125)-p(:,120)).^2+(p(:,126)-p(:,121)).^2);