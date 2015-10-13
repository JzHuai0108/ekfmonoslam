classdef KLT_tracker < handle
    properties (Hidden)
        maxFrames = 60; % giving the table row number
        maxPointsInTable=60; % maximum points in the table, giving the table col number
        minPointsInStates=30; % minimum number of tracked points in
        % the features_info, i.e.,in the states, before creating a new group of
        % features by triangulation from the feature table
        minFramesToTriang=3;
        half_patch_size_when_initialized=12;
        half_patch_size_when_matching=12;
        excluded_band; % the border area that is allowed to have a feature point
    end
    % The following properties can be set only by class methods
    properties (SetAccess = private)
        featureTable; % row corresponds to frames, columns to points
        firstFrmId; % the frame id w.r.t the video, 1 based
        lastFrmPos; % the frm pos w.r.t the feature table, 1 based
        nextFeatId=1; % the next feature id to be used
        nextGroupNum=1; % the next group number/ID to be used
        outFilePath; % the full path of the file to store the tracking table
        fileHandle; % the file pointer to the output file
        initRho=0.1; % unit 1/m
        lastImg; % the last gray scale single channel image
        videoName; % the name of the video file or image sequence folder
        seqType; % image sequence type, 0 for image sequence and 1 for video
        downScale; % down sample the input image, so that each edge is
        % reduced to 1/(1<<downScale) of its original length
        myVid;        
    end
    % Define an event
    events
        InsufficientPoints
    end
    methods
        function tracker = KLT_tracker(filename,  sequencePath, imgseqtype, downScale)
            tracker.firstFrmId=-1;
            tracker.lastFrmPos=-1;
            tracker.featureTable = zeros(tracker.maxFrames,9+6*tracker.maxPointsInTable);
            tracker.outFilePath=filename;
            tracker.fileHandle=fopen(filename,'W'); % feature table
            tracker.videoName=sequencePath;
            tracker.seqType=imgseqtype;
            tracker.myVid= VideoReader(sequencePath);
            numFrames= tracker.myVid.NumberOfFrames;
            tracker.downScale=downScale;
            tracker.excluded_band=tracker.half_patch_size_when_initialized+1; % the border area that is allowed to have a feature point
        end
        % In the current implementation, we directly provides this function
        % with features_info list and output it, which can be modified for
        % speed. pointGroup is also a list with element similar to features_info
        % input: frmId of gray image, nextImg, 1 based video frame number,
        % camera intrinsics and distortions, cam,
        % predPoints, the predicted point attributes list by a filter
        % qrs0, q s0 2s, r s in s0, s0 is a earth fixed
        % frame, not necessarily e-frame or the IMU start frame
        % output: the measurement of each predicted points, trackedPoints,
        % the new point group with attributes if a new group is required
        function  [trackedPoints, pointGroup]=TrackandDetectFeaturePoints(tracker, predPoints, grpPoses, qrs0, camPose, cam, frmId, timestamp )
            budgetPredPts=0; % how many points have predictions
            featureRow=zeros(1, 9+6*tracker.maxPointsInTable);
            featureRow(1,1)=frmId;
            featureRow(1,2)=timestamp;
            featureRow(1,3:9)=qrs0';
            % predict position of points in the feature table with positive
            % tracked frames No. and not in the states, i.e., features_info
            if(tracker.lastFrmPos~=-1) % not the first image
                % fill in the filter predicted positions
                for petal=1: length(predPoints)
                    if(~isempty(predPoints(petal).h))
                        featureRow(1, predPoints(petal).colFeatTable+(0:4))=...
                            [predPoints(petal).featId, predPoints(petal).h', predPoints(petal).invDepth, 1];
                        budgetPredPts=budgetPredPts+1;
                    end
                end
                % predict the points in the feature table not in the states
                if(tracker.lastFrmPos>1) % at least two frames available
                    % c(j+1) incoming camera frame, c(j) last camera frame
                    % c(j+1) to c(j-1)
                    Pcomp(:,:,1)=getRTcj2cjmk(qrs0, tracker.featureTable(tracker.lastFrmPos-1, 3:9)', camPose);
                    % c(j+1) to c(j)
                    Pcomp(:,:,2)=getRTcj2cjmk(qrs0, tracker.featureTable(tracker.lastFrmPos, 3:9)',camPose);
                end
                % c(j) to c(j+1)
                Pcj2cjp1=getRTcj2cjmk(tracker.featureTable(tracker.lastFrmPos, 3:9)',qrs0,camPose);
                rooster=1;
                while((rooster<tracker.maxPointsInTable+1)&&...
                        tracker.featureTable(tracker.lastFrmPos, 9+(rooster-1)*6+1)>0)
                    % all the feature points filled in the table even lost
                    % and not replaced still pass on feature index
                    trackBeacon=tracker.featureTable(tracker.lastFrmPos, 9+rooster*6);
                    if(tracker.featureTable(tracker.lastFrmPos, 9+rooster*6-1)==0 &&trackBeacon>0)
                        % was not in the states and was not lost
                        % predict the coordinates of the point that has been tracked
                        % at least twice
                        if(trackBeacon>1)
                            % homogeneous normalized image cooridnates without undistortion
                            xyni=  cam.K\[tracker.featureTable(tracker.lastFrmPos+(-1:0), 9+(rooster-1)*6+(2:3)), [1;1]]';
                            Xcjp1 = get_X_from_xP_lin(xyni,Pcomp); % the homogeneous coordinates of the point in c(j+1) frame
                            % predict the points that has been tracked only once,
                            % just detected
                        elseif(trackBeacon==1)
                            % approximate point coordinates in the cj frame
                            Xcj=(cam.K\[tracker.featureTable(tracker.lastFrmPos, 9+(rooster-1)*6+(2:3)),1]')/tracker.initRho;
                            Xcjp1=Pcj2cjp1*[Xcj;1];
                        end
                        uvd=hi_inverse_depth_v002(Xcjp1(1:3), cam); % projection with distortion
                        if(~isempty(uvd))
                            featureRow(1, 9+(rooster-1)*6+(1:5))=...
                                [tracker.featureTable(tracker.lastFrmPos,9+(rooster-1)*6+1), uvd', 1/Xcjp1(3), 0];
                            budgetPredPts=budgetPredPts+1;
                        end
                    end
                    rooster=rooster+1;
                end
            end
            % get the image of index frmId 1 based
            switch tracker.seqType
                case 0
                    nextImg = takeImage_v002( sequencePath, frmId, 'JPG' );
                case 1
                    nextImg=takeimagefromvideo(tracker.myVid, frmId, tracker.downScale);
                otherwise
            end
            assert(~isempty(nextImg));
            % tracking all the points given their previous image coordinates
            % and predicted image coordinates
            trackedPtsNo=0;
            if(budgetPredPts>0)
                ftMap=zeros(budgetPredPts,1);  % records the col number of the predicted points' index
                prevPts=cell(1,budgetPredPts);
                initFlow=cell(1,budgetPredPts);
                attendant=0; % iterator for predicted points
                rooster=1; % iterator for points in the feature table
                while((rooster<tracker.maxPointsInTable+1)&&attendant<budgetPredPts)
                    if(featureRow(1,9+(rooster-1)*6+1)>0) % have been predicted
                        attendant=attendant+1;
                        prevPts{attendant}=tracker.featureTable(tracker.lastFrmPos, 9+(rooster-1)*6+(2:3))-1;
                        % -1 because opencv image coord starts from 0,0
                        initFlow{attendant}=featureRow(1,9+(rooster-1)*6+(2:3))-1;
                        ftMap(attendant)=9+(rooster-1)*6+1;
                    end
                    rooster=rooster+1;
                end
                [prevPts, status, ~] = cv.calcOpticalFlowPyrLK(tracker.lastImg, nextImg, prevPts, ...
                    'InitialFlow',initFlow, 'WinSize', [21, 21], ...
                    'MaxLevel', 3,'Criteria', struct('type', 'Count+EPS', 'maxCount', 20, 'epsilon', 0.03));
                trackedPoints=predPoints;
                for attendant=1: budgetPredPts
                    if(status(attendant)) % i.e., ==1                   
                        % update feature Row that is to be added to the Table
                        featureRow(1,ftMap(attendant)+(1:2))=prevPts{attendant}+1; % because matlab image coord starts from 1,1
                        assert(tracker.featureTable(tracker.lastFrmPos, ftMap(attendant)+5)>0);
                        featureRow(1,ftMap(attendant)+5)=tracker.featureTable(tracker.lastFrmPos, ftMap(attendant)+5)+1; % the tracked No increment
                        trackedPtsNo=trackedPtsNo+1;
                        % else featureRow(1,ftMap(attendant)+5)=0; % lost points will
                        % likely be repleted in point detection phase
                    end
                end
                % update trackedPoints of output based on featureRow
                if(~isempty(grpPoses))
                    grpIds=[grpPoses.grpId];
                end
                correlation_threshold = 0.80;
                for atom=1: length(predPoints)
                    if(featureRow(1,trackedPoints(atom).colFeatTable+5)>tracker.minFramesToTriang&&...
                            featureRow(1,trackedPoints(atom).colFeatTable)==trackedPoints(atom).featId)
                        % this feature is tracked in the incoming frame,
                        % but this point may be a replacing one to a
                        % defunct point in the state, so make sure it is
                        % tracked at least minFramesToTriang+1
                        assert(trackedPoints(atom).lostNo==0&&(~isempty(trackedPoints(atom).h)));
                        % affine consistency check before updating the
                        % features_info
%                         assert(~isempty(trackedPoints(atom).silhou));
                        % Warp patches according to predicted motion,
                        % requiring groupPose, predPatch has size (2k+1)*(2j+1)
%                         center=round(featureRow(1,trackedPoints(atom).colFeatTable+(1:2)));
%                         grpFrmPose=grpPoses(find(grpIds==trackedPoints(atom).grpId,1)).pose;
%                         predPatch=pred_patch_fc( cam,trackedPoints(atom).xyn,...
%                             center', grpFrmPose, qrs0, camPose, trackedPoints(atom).invDepth, trackedPoints(atom).silhou);
%                         
%                         % match predicted apperance and the original one using normalized
%                         % cross-correlation
%                         if(~isempty(predPatch))
%                             half_patch_size_when_matching=floor(size(predPatch)/2);
%                             nextPatch= nextImg(...
%                                 center(2)+(-half_patch_size_when_matching(1):half_patch_size_when_matching(1)),...
%                                 center(1)+(-half_patch_size_when_matching(2):half_patch_size_when_matching(2)));
%                             % although pred_patch_fc is correct,
%                             sometimes it produces zero zones cause low
%                             correlation, as a result of the scale problem,
% we need to be careful about
%                             correlation_matrix = corrcoef(predPatch(:), double(nextPatch(:)));
%                             if( correlation_matrix(1,2) >= correlation_threshold)
                                trackedPoints(atom).z=featureRow(1,trackedPoints(atom).colFeatTable+(1:2))';
                                trackedPoints(atom).trackedNo = trackedPoints(atom).trackedNo + 1;
%                                 continue;
%                             end
%                         end
%                         trackedPoints(atom).lostNo=1;
%                         featureRow(1,trackedPoints(atom).colFeatTable+5)=0; % so that this point may be replaced later
%                         trackedPtsNo=trackedPtsNo-1;
                    else % the feature point is either just lost having h or lost for a while without h
                        %featureRow(1,trackedPoints(atom).colFeatTable+5)
                        %is not necessarily 0 as explained above
                        assert(isempty(trackedPoints(atom).z));
                        if(~isempty(trackedPoints(atom).h)) 
                            % the lostNo for points without prediction are taken care in predictPoints()
                            assert(trackedPoints(atom).lostNo==0);  
                            trackedPoints(atom).lostNo=1;                       
                        % for just lost features in the states, predict its
                        % inverse depth, use for KF update, to be implemented later
%                             trackedPoints(atom).invDepth=get_X_from_xP_lin();
                        end
                    end
                end
            else
                trackedPoints=[];
            end
            % Detect new Points if needed and fill in the featureRow, i.e.,
            % add more points to the feature table if lose too many
            if(trackedPtsNo<tracker.maxPointsInTable*2/3)
                %create mask, avoid detecting already tracked features, 8UC1
                mask=uint8(ones(size(nextImg)));
                mask(1:tracker.excluded_band,:)=0;
                mask(end-tracker.excluded_band+1:end,:)=0;
                mask(:,1:tracker.excluded_band)=0;
                mask(:,end-tracker.excluded_band+1:end)=0;
                
                initializing_box_size = [60,40];
                initializing_box_semisize = initializing_box_size/2;
                % mask predicted features by a box
                if(budgetPredPts>0)
                    for jack=1:budgetPredPts
                        if(status(jack)) % has been tracked
                            uv_pred=uint32(prevPts{jack})+1; % make values of uv_pred be positive integers, +1 becasue prevPts in Opencv starting from 0,0
                            mask(max(1,uv_pred(2)-initializing_box_semisize(2)):min(uv_pred(2)+initializing_box_semisize(2),size(nextImg,1)),...
                                max(1,uv_pred(1)-initializing_box_semisize(1)):min(uv_pred(1)+initializing_box_semisize(1), size(nextImg,2)))=0;
                        end
                    end
                end
                % the argument order and vacancy does not matter in matlab in constrast to
                % opencv C++
                % 4.The remaining corners are sorted by the quality measure in the descending order.
                nextPts = cv.goodFeaturesToTrack(nextImg,'MaxCorners',200, ...
                    'QualityLevel', 0.02,'MinDistance', 20,'Mask', mask,...
                    'BlockSize', 3, 'UseHarrisDetector', 0, 'K',0.04);
                nextPts=cv.cornerSubPix(nextImg, nextPts);
                newPtsNo=size(nextPts,2);
                jerk=1; % iterator for points in feature table
                router=0; % iterator for detected points
                while (jerk<tracker.maxPointsInTable+1)
                    if(featureRow(1, 9+jerk*6)<1) % (i.e., ==0)
                        assert(featureRow(1, 9+jerk*6)==0);
                        % the feature is either not predicted or just lost, we
                        % replace it with new points or pass on its values
                        if(router<newPtsNo)
                            % Opencv image starts from (0,0) upper left corner and matlab from (1,1)
                            router=router+1;
                            featureRow(1, 9+(jerk-1)*6+(1:6))=[tracker.nextFeatId,nextPts{router}+1, tracker.initRho, 0, 1];
                            tracker.nextFeatId=tracker.nextFeatId+1;
                        elseif(tracker.lastFrmPos~=-1 && (tracker.featureTable(tracker.lastFrmPos,9+(jerk-1)*6+1)>0)) % defunct ancestors pass on vlaues
                            trackBeacon=tracker.featureTable(tracker.lastFrmPos,9+jerk*6);
                            if(trackBeacon>0) % either just not predicted successfully or tracking lost
                                featureRow(1, 9+(jerk-1)*6+(1:6))=...
                                    [tracker.featureTable(tracker.lastFrmPos,9+(jerk-1)*6+(1:5)),-1];
                            else % has been lost at least in the last frame
                                assert(trackBeacon<0);
                                featureRow(1, 9+(jerk-1)*6+(1:6))=...
                                    [tracker.featureTable(tracker.lastFrmPos,9+(jerk-1)*6+(1:5)),trackBeacon-1];
                            end
                        end
                    end
                    jerk=jerk+1;
                end
            else % update the feature table for features that are either
                % not predicted or just lost, we pass on its values and
                % reset the trackedNo
                jerk=1; % iterator for points in feature table
                while (jerk<tracker.maxPointsInTable+1)
                    if(featureRow(1, 9+jerk*6)<1) % (i.e., ==0)
                        assert(featureRow(1, 9+jerk*6)==0);
                        % the feature is either not predicted or just lost, we pass on its values
                        if(tracker.lastFrmPos~=-1 && (tracker.featureTable(tracker.lastFrmPos,9+(jerk-1)*6+1)>0)) % defunct ancestors pass on vlaues
                         	trackBeacon=tracker.featureTable(tracker.lastFrmPos,9+jerk*6);
                            if(trackBeacon>0) % either just not predicted successfully or tracking lost
                                featureRow(1, 9+(jerk-1)*6+(1:6))=...
                                    [tracker.featureTable(tracker.lastFrmPos,9+(jerk-1)*6+(1:5)),-1];
                            else % has been lost at least in the last frame
                                assert(trackBeacon<0);
                                featureRow(1, 9+(jerk-1)*6+(1:6))=...
                                    [tracker.featureTable(tracker.lastFrmPos,9+(jerk-1)*6+(1:5)),trackBeacon-1];
                            end
                        end
                    end
                    jerk=jerk+1;
                end
            end
            % up to now, featureRow has the latest information of each
            % feature point, (1) the feature point is tracked
            % successfully, (2) the point is lost and replaced by new
            % points, (3) the point inherits its parent's properties.
            % use tracked points in featureRow to initialize a pointGroup
            trackedPinS=0; % No. of successfully tracked points in states
            jerk=1; % iterator for points in feature table
            while ((jerk<tracker.maxPointsInTable+1) && (featureRow(1, 9+(jerk-1)*6+1)>0)) % all the points are in the front of each row of the table
                if(featureRow(1, 9+(jerk-1)*6+5)==1 && featureRow(1, 9+jerk*6)>tracker.minFramesToTriang)
                    % was in states and really being tracked consecutively
                    trackedPinS=trackedPinS+1;
                end
                jerk=jerk+1;
            end
            pointGroup=[];
            if(trackedPinS<tracker.minPointsInStates && tracker.lastFrmPos>tracker.minFramesToTriang-2)
                dummys=struct('featId', 0,'grpId', 0,'initFrmNo', 0,'colFeatTable',0,...
                    'xyn', [-1;-1], 'invDepth', 0,'h', [], 'z', [], 'H', [], 'S', [],...
                    'low_innovation_inlier', 0, 'high_innovation_inlier', 0, ...
                    'lostNo', 0, 'trackedNo', 0, 'silhou', []);
                pointGroup=repmat(dummys, 1,floor(tracker.minPointsInStates*2/3));
                % the projection matrix for the last few frames w.r.t the
                % lastest/incoming camera frame, cj, e.g., [Rcj2cj-k, Tcj2cj-k]
                Pcomp=zeros(3,4,tracker.minFramesToTriang);
                Pcomp(:,:, 1)=[eye(3), zeros(3,1)]; % qcj2cj and Tcj in cj
                for twig=2:tracker.minFramesToTriang
                    Pcomp(:,:, twig)=getRTcj2cjmk(featureRow(1, 3:9)', tracker.featureTable(tracker.lastFrmPos-twig+2, 3:9)', camPose);
                end
                xyns=ones(3,tracker.minFramesToTriang);                
                jerk=1; % iterator for points in feature table's one row
                howManySoFar=0; % how many points has been triangulated?
                % all the established feature points are in front of a row
                while ((jerk<tracker.maxPointsInTable+1) && (featureRow(1, 9+(jerk-1)*6+1)>0))
                    if(featureRow(1, 9+(jerk-1)*6+5)==0 && featureRow(1, 9+jerk*6)>=tracker.minFramesToTriang)
                        % not in states and has been tracked long enough
%                         center=round(tracker.featureTable(tracker.lastFrmPos, 9+6*(jerk-1)+(2:3)));                  
%                         if(~(center(1)+half_patch_size_when_initialized>cam.nCols||...
%                                 center(1)-half_patch_size_when_initialized<1||...
%                                 center(2)+half_patch_size_when_initialized>cam.nRows||...
%                                 center(2)-half_patch_size_when_initialized<1))
                            howManySoFar=howManySoFar+1;
                            % It is sort of contradictory that I initialize
                            % a feature with its last appearance and set
                            % its group frame as the present frame, but I
                            % believe it won't make much difference in
                            % results
%                             pointGroup(howManySoFar).silhou= tracker.lastImg(...
%                                 center(2)+(-half_patch_size_when_initialized:half_patch_size_when_initialized),...
%                                 center(1)+(-half_patch_size_when_initialized:half_patch_size_when_initialized));                           
                            xyns(1:2, 1)=normalize(featureRow(1, 9+6*(jerk-1)+(2:3))',cam.fc,cam.cc,cam.kc,cam.alpha);
                            for twig=2:tracker.minFramesToTriang
                                xyns(1:2, twig)=normalize(...
                                    tracker.featureTable(tracker.lastFrmPos-twig+2, 9+6*(jerk-1)+(2:3))',...
                                    cam.fc,cam.cc,cam.kc,cam.alpha);
                            end
                            Xcj=get_X_from_xP_lin(xyns, Pcomp);
                            pointGroup(howManySoFar).featId=featureRow(1, 9+(jerk-1)*6+1);
                            pointGroup(howManySoFar).grpId=tracker.nextGroupNum;
                            pointGroup(howManySoFar).initFrmNo=frmId;
                            pointGroup(howManySoFar).colFeatTable=9+(jerk-1)*6+1;
                            pointGroup(howManySoFar).xyn=Xcj(1:2)/Xcj(3);
                            pointGroup(howManySoFar).invDepth=1/Xcj(3);
                            pointGroup(howManySoFar).h=[];
                            pointGroup(howManySoFar).z=[];
                            pointGroup(howManySoFar).H=[];
                            pointGroup(howManySoFar).S=[];
                            pointGroup(howManySoFar).low_innovation_inlier=0;
                            pointGroup(howManySoFar).high_innovation_inlier=0;
                            pointGroup(howManySoFar).lostNo=0;
                            pointGroup(howManySoFar).trackedNo=1;
                            featureRow(1, 9+(jerk-1)*6+5)=1; % this feature is in the states
%                         end
                    end
                    jerk=jerk+1;
                end
                pointGroup=pointGroup(1:howManySoFar);
            end
            % maintain the tracker
            if(tracker.lastFrmPos~=-1) % not the first image
                if (tracker.lastFrmPos<tracker.maxFrames)
                    tracker.lastFrmPos=tracker.lastFrmPos+1;
                    tracker.featureTable(tracker.lastFrmPos, :)= featureRow;
                else
                    SaveFeatureTableRow(tracker.featureTable(1,:), tracker.fileHandle,...
                        tracker.maxPointsInTable); % write one line of feature table into a file
                    tracker.featureTable=[tracker.featureTable(2:end, :); featureRow];
                end
            else
                tracker.featureTable(1, :)= featureRow;
                tracker.firstFrmId=frmId;
                tracker.lastFrmPos=1;
            end
            tracker.lastImg=nextImg;
            if(~isempty(pointGroup))
                tracker.nextGroupNum=tracker.nextGroupNum+1;
            end
        end
        % frmId, 1 based video frame Number, keep tracker.featureTable
        % posted after filter update
        function PostStates(tracker, ftDepthCol, rvqs0, frmId)
            rowNo=min(frmId-tracker.firstFrmId+1, tracker.maxFrames);
            assert(rowNo==tracker.lastFrmPos);
            tracker.featureTable(rowNo, 3:9)=rvqs0([7:10, 1:3])';
            % post inverse depth for each feature
            for tractor=1:size(ftDepthCol,2)
                tracker.featureTable(rowNo, ftDepthCol(1,tractor)+3)=ftDepthCol(2,tractor);
            end
        end
        % convert the e frame coordination of rs in e and qe2s (not qs2e,
        % for extensibility of TrackandDetectFeaturePoints() with filter
        % EKF_INS_GPS, so be careful about its input rvqs0), to s0 frame
        % formulation, rs in s0, q s0 2s.
        % input: rqs02e, rs0 in e and q s0 2e.
        function FeatureTable_e2s0(tracker, rqs02e)
            terminal=min(tracker.minFramesToTriang, tracker.lastFrmPos);
            for clerk=1:terminal
                tracker.featureTable(tracker.lastFrmPos-clerk+1, 3:6)=...
                    quatmult_v001(tracker.featureTable(tracker.lastFrmPos-clerk+1, 3:6),rqs02e(4:7), 0)';
                tracker.featureTable(tracker.lastFrmPos-clerk+1, 7:9)=...
                    quatrot_v000(rqs02e(4:7), tracker.featureTable(tracker.lastFrmPos-clerk+1, 7:9)'-rqs02e(1:3), 1)';
            end
        end
    end % methods
end % classdef
% input: (1)qTs0nsj, qs02sj and Tsj2s0, (2) qTs0nsjmk, qs02sj-k and Tsj-k2s0,
% (3) qTs2c, qs2c, Ts2c, s0 is a earth fixed frame, s(j) is the sensor
% frame at epoch j, c(j) is the associated camera frame
% output: Pcj2cjmk, Rcj2cj-k and Tcj2cj-k, i.e., Tcj in c(j-k)
function Pcj2cjmk=getRTcj2cjmk(qTs0nsj, qTs0nsjmk, qTs2c)
% the camera's pose c(j-k) and c(j) in s0 frame
qs02cjmk=quatmult_v001(qTs2c(1:4),qTs0nsjmk(1:4), 0);
qs02cj=quatmult_v001(qTs2c(1:4),qTs0nsj(1:4), 0);
% the relation between the frame c(j) and the incomiing
% frame c(j-k)
qcj2cjmk=quatmult_v001(qs02cjmk, qs02cj,2);
Tcj2cjmk=qTs2c(5:7)+quatrot_v000(qs02cjmk,qTs0nsj(5:7)-qTs0nsjmk(5:7), 0)-quatrot_v000(qcj2cjmk, qTs2c(5:7),0);
Pcj2cjmk=[quat2dcm_v000(qcj2cjmk),Tcj2cjmk];
end
% write one line of feature table into a file
function SaveFeatureTableRow(featureRow,outFP, maxPointsInTable)
fprintf(outFP, '%d\t%15.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%14.6f\t%14.6f\t%14.6f\n', featureRow(1:9));
jerk=1;
while((jerk<maxPointsInTable+1)&&featureRow(1,9+(jerk-1)*6+1)>0 ) % the feature points are all in front of the table
    fprintf(outFP, '\t%d\t%10.5f\t%10.5f\t%10.8f\t%d\t%d\n', featureRow(1,9+(jerk-1)*6+(1:6)));
    jerk=jerk+1;
end
end
% given two camera poses and its appearance in the first camera frame where the
% world coordinate system sits, predict the patch appearance in the second
% camera, the predicted patch. We denote the first camera frame as c1 that
% is the world frame in this function domain, and the second camera pose, c2 
% input: cam, camera intrinsic and distortion parameters, xync1, the 
% normalized image coordinates of the point in c1 frame, 2x1, uvdc2, 2x1, the 
% distorted image coordinates of the point in c2 frame, poses1, qs0 2 s1, 
% T s1 in s0, s0 is a earth fixed sensor frame, poses2, q s0 2 s2 and T s2
% in s0, qTs2c, the sensor frame and the camera frame relation, qs2c and T
% s in c, rho is the inverse depth of the point in c1 frame, patchc1, the
% appearance of the point in the c1 frame
% the homography transform between two point patches are referred to the
% Ch.13 of Hartley and Zisserman multiple view geometry
function patch_pred=pred_patch_fc(cam,xync1,uvdc2,poses1,poses2, qTs2c, rho, patchc1)
half_patch_size_when_initialized=12;
half_patch_size_when_matching=half_patch_size_when_initialized;
halfW_pred=half_patch_size_when_matching;
halfW_fea=half_patch_size_when_initialized;
if((uvdc2(1)>halfW_pred) && (uvdc2(1)<cam.nCols-halfW_pred+1)&&...
        (uvdc2(2)>halfW_pred) && (uvdc2(2)<cam.nRows-halfW_pred+1))
    % note if a minus sign is added before the two 1's in the 
    % following two lines as done in Civera's this function, 
    % This sign change seems to improve the performance a little bit, but
    % it brings about a bug with the line
    % uv_p_pred_patch=rotate_with_dist_fc_c2c1_v001(); 
    % which sometimes gives huge predicted image coordinates and the 
    % following meshgrid cannot guarantee the correct grid size with large numbers
    Pc12c2=getRTcj2cjmk(poses1, poses2, qTs2c);  
    uvdc1=project_points2([xync1;1],zeros(3,1),zeros(3,1),cam.fc,cam.cc,cam.kc,cam.alpha); 
    n1 = [xync1;1];
    n2 = cam.K\[uvdc2;1];
    n2 = Pc12c2(1:3,1:3)'*n2;
    n1 = n1/norm(n1);
    n2 = n2/norm(n2);
    n = n1+n2;
    n = n/norm(n);
    % we put the c1 as the fundamental frame for now 
    d = -n'*[xync1;1]/rho; 
    % the predicted uv image coordinates in c2 frame
    uv_p_pred_patch=rotate_with_dist_fc_c2c1_v001(cam,uvdc1',Pc12c2(1:3,1:3),Pc12c2(1:3,4), n,d);
    % if a bug occurs with the following two lines, it means the geometry
    % of the camera w.r.t the environment is mistaken. Double check it!
    [u_pred,v_pred]=meshgrid(uv_p_pred_patch(1)-halfW_pred:uv_p_pred_patch(1)+halfW_pred,uv_p_pred_patch(2)-halfW_pred:uv_p_pred_patch(2)+halfW_pred);
    uv_pred=[reshape(u_pred,(halfW_pred*2+1)^2,1),reshape(v_pred,(halfW_pred*2+1)^2,1)];
    % the reprojected image coordinates in c1 frame
    uv_pred_imak_dist=rotate_with_dist_fc_c1c2_v001(cam,uv_pred,Pc12c2(1:3,1:3),Pc12c2(1:3,4), n,d);
    % align the predicted patch in c2 to the intensity patch source in c1 frame
    uv_pred_imak_dist(:,1)=uv_pred_imak_dist(:,1)-(uvdc1(1)-halfW_fea-1);% + 0.5*ones(size(uv_pred_imak_dist,1),1);
    uv_pred_imak_dist(:,2)=uv_pred_imak_dist(:,2)-(uvdc1(2)-halfW_fea-1);% + 0.5*ones(size(uv_pred_imak_dist,1),1);
    u_pred_imak_dist=reshape(uv_pred_imak_dist(:,1),halfW_pred*2+1,halfW_pred*2+1);
    v_pred_imak_dist=reshape(uv_pred_imak_dist(:,2),halfW_pred*2+1,halfW_pred*2+1);
    % interpolation to determine the intensity of the patch in c2 according
    % to its c1 corresnpondences
    [u_fea,v_fea]=meshgrid(1:size(patchc1,1),1:size(patchc1,2));
    patch_pred=interp2(u_fea,v_fea,double(patchc1),u_pred_imak_dist,v_pred_imak_dist);    
else
    patch_pred=[];
end
end
