function a=combin(n,k,p)
%p=[0 0];
%n=4;
%k=2;
p1=floor(log2(nchoosek(n,k)));
z=bi2de(p,'left-msb');
ckk=0:n-1;
m1=0;
m2=0;
l=1;
b=zeros(1:n);
for i=1:k
    b=ckk;
    len=length(b);
   for j=1:len;
       ck=max(b);
       if(ck<k-m1) 
           a(l)=ck+1;
           l=l+1;
           ckk= ckk(ckk~=ck);
           m1=m1+1;
           break;
       else
           if(nchoosek(ck,k-m1)<=z-m2)
               m2=m2+nchoosek(ck,k-m1);
               a(l)=ck+1;
               l=l+1;
               ckk= ckk(ckk~=ck);
               m1=m1+1;
               break;
           else
               b=b(b~=ck);
           end
       end
   end
end
a=fliplr(a);

   
        
