function matches=associate2(first_list, second_list, offset, max_difference)
% first_list first column are timestamps, so are second_list, both
% first_list and second_list are of same columns, offset is to be applied
% to second_list when doing association, 
% matches first col the index in the first_list and second col, index in
% second_list
% author Jianzhu Huai 2014

    firstArray=[first_list(:,1), ones(size(first_list,1),1), (1:size(first_list,1))'; inf, 1, -1];
    secondArray=[second_list(:,1)+offset, ones(size(second_list,1),1)*2, (1:size(second_list,1))'; inf, 2, -1];
    resultArray=merge_sorted_arrays(firstArray, secondArray);
    secondUsed=zeros(size(second_list,1),1);% whether second data was matched
    matches=zeros(size(firstArray,1),2);
    matchCount=0;
    for i=1:size(resultArray, 1)
        if(resultArray(i, 2)==1) %from first_list
            % check both sides
            leftdiff=inf;
            rightdiff=inf;
            
            if(i>1 && resultArray(i-1,2)==2 && secondUsed(resultArray(i-1,3))==0)
                   leftdiff= abs(firstArray(resultArray(i,3), 1)-secondArray(resultArray(i-1,3), 1));           
            end
            if(i<size(resultArray,1)&& resultArray(i+1,2)==2 && secondUsed(resultArray(i+1,3))==0)
                   rightdiff= abs(firstArray(resultArray(i,3), 1)-secondArray(resultArray(i+1,3), 1));           
            end
            if(leftdiff<max_difference || rightdiff< max_difference)
                if(leftdiff<rightdiff)
                    assert(leftdiff<max_difference);
                    matchCount=matchCount+1;
                    matches(matchCount, :)=[resultArray(i,3), resultArray(i-1,3)];
                    secondUsed(resultArray(i-1,3))= resultArray(i,3);
                else
                    matchCount=matchCount+1;
                    matches(matchCount, :)=[resultArray(i,3), resultArray(i+1, 3)];
                    assert(rightdiff<max_difference);
                    secondUsed(resultArray(i+1,3))= resultArray(i,3);
                end                
            end     
        end
    end
    matches=matches(1:matchCount,:);
    
    function A=merge_sorted_arrays(L, R)
        % L mxn, R pxn, ascendingly sorted arrays,
        % n fields, m-1 and p-1 observations, last row is [inf, ...] as sentinel
        % see introduction to algorithms 2nd edition p29     
        r=size(L,1)+size(R,1)-2;       
        iota=1; 
        j=1;
        A=zeros(r,size(L,2));
        for k=1:r
            if(L(iota,1)<=R(j,1))
                A(k,:)=L(iota,:);
                iota=iota+1;
            else
                A(k, :)=R(j,:);
                j=j+1;
            end
        end
    end
    

end