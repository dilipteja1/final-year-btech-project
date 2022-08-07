
function p1=invcomb(k,p)
%n=4;
%k=2;
%p=[1 4];
p=fliplr(p);
z=0;
l=k;
p=p-1;
for i=1:k
    if(p(i)<l) l=l-1;break;
    else
    z=z+nchoosek(p(i),l);
    l=l-1;
    end
end
p1=de2bi(z,8,'left-msb');
