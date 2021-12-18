function matches=associate(first_list, second_list,offset,max_difference)
%     Associate two dictionaries of (stamp,data). As the time stamps never match exactly, we aim
%     to find the closest match for every input tuple.
%     Input:
%     first_list -- first dictionary of (stamp,data) tuples
%     second_list -- second dictionary of (stamp,data) tuples
%     offset -- time offset between both dictionaries (e.g., to model the delay between the sensors)
%     max_difference -- search radius for candidate generation
%
%     Output:
%     matches -- list of matched tuples ((index1),(index2))

first_keys = first_list(:,1); % timestamps
second_keys = second_list(:,1); % timestamps
matchCount=0;
potential_matches=zeros(length(first_keys)+length(second_keys), 3);   % curDiff, a, b
for i=1:length(first_keys)
    for j=1:length(second_keys)
        curDiff=abs(first_keys(i)-(second_keys(j)+offset));
        if(curDiff<max_difference)
            matchCount=matchCount+1;
            potential_matches(matchCount,:)=[curDiff, i, j];            
        end
    end
end
potential_matches=potential_matches(1:matchCount,:);
sort(potential_matches);
matchCount=0;
matches = zeros(length(first_keys), 2);
first_keys_flag=zeros(length(first_keys)); % whether the keys has been matched and inserted
second_keys_flag=zeros(length(second_keys));
for i=1:size(potential_matches,1)
    if (first_keys_flag(potential_matches(i, 2))==0 && second_keys_flag(potential_matches(i, 3))==0)
        first_keys_flag(potential_matches(i, 2))=1;
        second_keys_flag(potential_matches(i, 3))=1;
        matchCount=matchCount+1;
        matches(matchCount,:)=potential_matches(i, 2:3);
    end
end
matches=matches(1:matchCount, :);
sort(matches);
end

