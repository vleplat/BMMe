function X=remove0col(X,m,n,r)
 X1=[];
 %  remove 0 column
   for i=1:n
       if max(X(:,i))>0 
       X1=[X1 X(:,i)];    
       end
   end
   % remove 0 row
   X=[];
   for i=1:size(X1,1)
       if max(X1(i,:))>0
           X=[X; X1(i,:)];
       end
   end
    [m,n]=size(X);
   if r>min(m,n) 
       error('decrease r');
   end 
end