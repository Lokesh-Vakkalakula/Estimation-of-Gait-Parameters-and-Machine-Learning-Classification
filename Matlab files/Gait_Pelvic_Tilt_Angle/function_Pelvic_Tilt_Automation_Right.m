
function [A_pelvic] = function_Pelvic_Tilt_Automation_Right(w)
w=sortrows(w,1);
w_pelvic=p;
w_pelvic=[w_pelvic(:,1),w_pelvic(:,7:11),w_pelvic(:,42:51)];
[m,n]=size(w_pelvic)
for i = 2:m-1
    for j=1:n-1
      if w_pelvic(i,j)==0
        w_pelvic(i,j)=(w_pelvic(i+1,j)+w_pelvic(i+1,j))/2;
      end
    end
end
w_pelvic(w_pelvic==56)= 0;
w_pelvic=w_pelvic(all(w_pelvic,2),:);
% figure(2);
AngleNum_pelvic=[((w_pelvic(:,2)-w_pelvic(:,7)).*(w_pelvic(:,12)-w_pelvic(:,7)))+((w_pelvic(:,3)-w_pelvic(:,8)).*(w_pelvic(:,13)-w_pelvic(:,8)))]
AngleDenom_pelvic=[(sqrt(((w_pelvic(:,2)-w_pelvic(:,7)).^2)+((w_pelvic(:,3)-w_pelvic(:,8)).^2))).*(sqrt(((w_pelvic(:,12)-w_pelvic(:,7)).^2)+((w_pelvic(:,13)-w_pelvic(:,8)).^2)))]
Angle_pelvic = [acosd(AngleNum_pelvic./AngleDenom_pelvic)]
y_pelvic=90-Angle_pelvic
for i=1:length(Angle_pelvic)
    if abs(y_pelvic(i))>15
        y_pelvic(i)=0
    end
    if isnan(y_pelvic(i))
        y_pelvic(i)=0
    end
end
A_pelvic =(movmean(sgolayfilt(y_pelvic,2,11),2));
x_pelvic=w_pelvic(:,1);
A_pelvic=[x_pelvic,A_pelvic]
% figure(2);
% plot(x_pelvic(:,1),A_pelvic)


% Gait_left_min_index_pelvic=find(x_pelvic==x(Gait_left_min_index));
% Gait_right_min_index_pelvic=find(x_pelvic==x(Gait_right_min_index));
% array_gait_pelvic_tilt=A_pelvic(Gait_left_min_index_pelvic:Gait_right_min_index_pelvic);
% array_frame_gait=x_pelvic(Gait_left_min_index_pelvic:Gait_right_min_index_pelvic);
end

