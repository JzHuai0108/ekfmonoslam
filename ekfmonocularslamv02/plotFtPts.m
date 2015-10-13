chi_095_2 = 5.9915;
% chi_099_2 = 9.2103;
chi_095_3 = 7.8147;
% chi_099_3 = 11.3449;

% plot image stuff
figure(figure_all);
subplot(im_fig);
hold off;
imagesc(tracker.lastImg);
colormap gray;
hold on;
title('Thick red: low innovation inliers. Thin red: high innovation inliers. \newline Magenta: rejected by 1-point RANSAC. Blue: No match found by pyramid KLT.');
half_patch_size_when_initialized = 20;
half_patch_size_when_matching = 6;
features_info=filter.features_info;

for i=1:length(features_info)
    
    if (~isempty(features_info(i).h)&&~isempty(features_info(i).S))
        
        %         imagesc( features_info(i).h(1) - half_patch_size_when_matching,...
        %             features_info(i).h(2) - half_patch_size_when_matching, features_info(i).patch_when_matching);
        
        if features_info(i).low_innovation_inlier
            plotUncertainEllip2D( features_info(i).S,...
                features_info(i).h, chi_095_2, 'r', 4 )
            plot( features_info(i).h(1), features_info(i).h(2),'r+','Markersize',10);
        end
        if features_info(i).high_innovation_inlier
            plotUncertainEllip2D( features_info(i).S,...
                features_info(i).h, chi_095_2, 'r', 2 )
            plot( features_info(i).h(1), features_info(i).h(2),'r+','Markersize',10);
        end
        if (~isempty(features_info(i).z))&&(features_info(i).high_innovation_inlier==0)&&(features_info(i).low_innovation_inlier==0)
            plotUncertainEllip2D( features_info(i).S,...
                features_info(i).h, chi_095_2, 'm', 2 )
            plot( features_info(i).h(1), features_info(i).h(2),'m+','Markersize',10);
        end
        
        if (isempty(features_info(i).z)&&~isempty(features_info(i).h))
            plotUncertainEllip2D( features_info(i).S,...
                features_info(i).h, chi_095_2, 'b', 2 )
            plot( features_info(i).h(1), features_info(i).h(2),'b+','Markersize',10);
        end
        
        if (~isempty(features_info(i).z))
            plot( features_info(i).z(1), features_info(i).z(2),'g+','Markersize',10);
            plot([features_info(i).z(1),features_info(i).h(1)],[features_info(i).z(2),features_info(i).h(2)],'Color','r','LineWidth',2);
        end        
    end
    
end

% plot predicted and measured
% which_are_predicted = predicted_measurements(:,1)>0;
% plot( predicted_measurements(which_are_predicted,1), predicted_measurements(which_are_predicted,2),'r+','Markersize',10);
% which_are_measured = measurements(:,1)>0;
% plot( measurements(which_are_measured,1), measurements(which_are_measured,2),'g+');

axes_handler = get(gcf,'CurrentAxes');
set(axes_handler,'XTick',[],'YTick',[]);

% plot 3D stuff
figure(figure_all);
subplot(near3D_fig);
hold off;

imgctr=step-firstGroupFrameId;
qs02c0=filter.camPose(1:4);
% obtain qc2c0 current camera frame c to start camera frame c0, T c in c0
qs02c=quatmult_v001(filter.camPose(1:4),filter.rvqs0(7:10),0);
qc2c0=quatmult_v001(qs02c0,qs02c,2);
Tc2c0=-quatrot_v000(qc2c0, filter.camPose(5:7),0)+...
    quatrot_v000(qs02c0,filter.rvqs0(1:3),0)+filter.camPose(5:7);
trajectory(:,imgctr) =[Tc2c0; qc2c0];
draw_camera( trajectory(:,imgctr), 'k' );
hold on;
title(sprintf('Camera motion and scene [m] in top view at frmId(1 based)%d.',step));

plot3( trajectory(1, 1:imgctr), trajectory(2, 1:imgctr),...
    trajectory(3, 1:imgctr), 'k', 'LineWidth', 2 );

% for each group frame, compute Cci2c0 and Tcj2c0
grpLeg=length(filter.groupPose);
groupTransform=zeros(7, length(filter.groupPose));
qs02c0=filter.camPose(1:4);
for apple=1:grpLeg
    qs02cj=quatmult_v001(filter.camPose(1:4),filter.groupPose(apple).pose(1:4),0);
    groupTransform(1:4,apple)=quatmult_v001(qs02c0,qs02cj,2);
    groupTransform(5:7,apple)=-quatrot_v000(groupTransform(1:4,apple), filter.camPose(5:7),0)+...
        quatrot_v000(qs02c0,filter.groupPose(apple).pose(5:7),0)+...
        filter.camPose(5:7);
end
groupIds=[filter.groupPose.grpId];
covDim=filter.groupFrameSIP-1+length(filter.groupPose)*6+length(filter.features_info);
featSIP=filter.groupFrameSIP+length(filter.groupPose)*6;
%for each feature with lostNo==0, obtain its position in c0 frame
for wasp = 1:length(features_info)
    if (features_info(wasp).lostNo==0)
        % the feature is still being tracked, note some point
        % may be lost. But before a new group comes, it is still in
        % the states
        apple=find(groupIds==filter.features_info(wasp).grpId);
        
        % because grpId contains gap, this instruction can be sped up by making use of bisection search
        rho=filter.features_info(wasp).invDepth;
        XYZ=quatrot_v000(groupTransform(1:4,apple), [filter.features_info(wasp).xyn; 1]/rho,0)+...
            groupTransform(5:7,apple);
        if rho-3*sqrt(filter.p_k_k(featSIP+wasp-1,featSIP+wasp-1))<0
            if ( rho>0)
                plot3(XYZ(1),XYZ(2),XYZ(3),'k^','Markersize',10)
            end
        else
            if ( rho>0)
                plot3(XYZ(1),XYZ(2),XYZ(3),'r+','Markersize',10)
            end
        end
    end
end

axes_handler = get(gcf,'CurrentAxes');
axis([-100 150 -5 6 -200 210]);
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
hold off
view([0, 0]);% only x and z