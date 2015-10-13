
% get the content of q java linkedlist, assume that each element of this
% queue is a column vector of the same size
function array=listcontent(list)
stork= list.size();
dimen=length(list.getFirst());
array=zeros(dimen,stork);
for iota=0:stork-1
    array(:,iota+1)=list.get(iota);
end

% function testqueue()
% import java.util.LinkedList
% q = LinkedList();
% 
% q.isempty()
% for i=1:10000
%     q.addLast(rand(10,1));
% end
% q.getFirst()
% 
% q.size()
% arr=listcontent(q);
% arr(:,1)
% q.removeFirst();
% q.getFirst()
% arr(:,2)

