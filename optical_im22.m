clear all;
addpath('C:\Users\Guest\Downloads\');

warning('on','MATLAB:RandStream:ActivatingLegacyGenerators') 
warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState')
clear all
NgType=1;
if NgType==1, nt='CP'; elseif NgType==2, nt='ZP'; end
chType='AWGN'; Target_neb=100; 
Nbps=2; M=2^Nbps;
%Nfft=16; 
N=64;
Nfft=(N+1)*2;
g=16;
n=4;
k=2;
p1=floor(log2(nchoosek(n,k)));
p2=k*log2(M);
Ng=0;
Nsym=Nfft+Ng;   
EbN0=0:1:30; % EbN0
N_iter=1000; % Number of iterations for each EbN0
Nframe=1; % Number of symbols per frame
sigPow=0; 
file2=['OFDM_BER2_' chType '_' nt '_' 'GL' num2str(Ng) '.dat'];
fid=fopen(file2, 'w+');

const=qammod([0:M-1],M);
scale=modnorm(const,'avpow',1);
s1=[0:M-1];
s=qammod(s1,M,0,'gray')*scale;
%norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)]; % BPSK 4-QAM 16-QAM
for i=1:length(EbN0)
randn('state',0); rand('state',0); 
Neb=0; Ntb=0; 
for m=1:N_iter
% Tx
X=randi([0 1],1,Nframe*g*(p1+p2));
if NgType~=2, x_GI=zeros(1,Nframe*Nsym);
elseif NgType==2, x_GI= zeros(1,Nframe*Nsym+Ng);
end
kk1=g*(p1+p2);
kk2=1:kk1;kk3=1:Nsym;
for q=1:Nframe
    B2=zeros(1,N);
    XX1=X(kk2);
    kk4=1:p1+p2;
    kk5=1:n;
    for m=1:g
        B1=zeros(1,n);
        XX2=XX1(kk4);
        X1=XX2(1:p1);
        X2=XX2(p1+1:end);
        %X3=lookup(1,X1);
        X3=combin(n,k,X1);
        A1=reshape(X2,[Nbps,k]);
        A2=transpose(A1);
        A3=num2str(A2);
        A4=bin2dec(A3);
        A5=transpose(A4);
        Xmod=qammod(A5,M,0,'gray')*scale;
        r=1;
        for l=1:n
            for j=1:length(X3);
                if(X3(j)==l)
                    B1(l)=Xmod(r);
                    r=r+1;
                end
            end
        end
        kk4=kk4+(p1+p2);
        B2(kk5)=B1;
        kk5=kk5+n;
    end
    X2mod=[0,B2,0,fliplr(conj(B2))];
      
x=ifft(X2mod);
x_GI(kk3)= guard_interval(Ng,Nfft,NgType,x);
kk2=kk2+kk1;kk3=kk3+Nsym;
end

bdc=12;
clip=sqrt((10.^(bdc/10))-1);
x1_GI=sum(x_GI.^2);
x11_GI=x1_GI/length(x_GI);
b=clip*sqrt(x11_GI);

x2_GI=x_GI+b;
for l=1:length(x2_GI)
if x2_GI(l)<0,
    x2_GI(l)=0;
end
end

y= x2_GI; % No channel
sigPow=sum(y.^2)/length(y);

snr = EbN0(i)+10*log10(Nbps);
noise_mag = sqrt((10.^(-snr/10))*sigPow/2);
y_GI = y + noise_mag*(randn(size(y)));

kk1=(NgType==2)*Ng+[1:Nsym];
kk2=1:Nfft;
kk3=1:N;
kk4=2:Nfft/2;kk7=1:g*(p1+p2);
Y=zeros;
Xmod_r=zeros;
N0=(10.^(-snr/10))*sigPow*2;
for q=1:Nframe
Y(kk2)= fft(remove_GI(Ng,Nsym,NgType,y_GI(kk1)));
Xmod_r(kk3) = Y(kk4);
kk5=1:n;kk6=1:(p1+p2);
for m=1:g
    r1=1;
    r2=1;
    a1=zeros(1,k);
    a2=zeros(1,k);
    Xmod2_r=Xmod_r(kk5);
    for j=1:n
        f=0;
        for t=1:M
            f=f+exp((-1/N0)*(abs(Xmod2_r(j)-s(t))^2));
        end
        llr(j)=log(k)-log(n-k)+((abs(Xmod2_r(j))^2)/N0)+log(f);
       
    end
     for j=1:n
          if(llr(j)>0) a1(r1)=j; a2(r2)=Xmod2_r(j);r1=r1+1;r2=r2+1;
          end
     end
     a1=a1(1:k);
     a2=a2(1:k);
     %pp1=lookup(0,a1);
     pp1=invcomb(k,a1);
     pp1=pp1(length(pp1)-p1+1:end);
     pp2=qamdemod(a2/scale,M,0,'gray');
     pp2bin  = dec2bin(pp2,Nbps)-'0' ;
     pp2bin=transpose(pp2bin);
     pp2bin=reshape(pp2bin,[1,p2]);
     Xmod3_r(kk6)=[pp1 pp2bin];

kk5=kk5+n;kk6=kk6+(p1+p2);
end
Xmod4_r(kk7)=Xmod3_r;
kk1=kk1+Nsym; kk2=kk2+Nfft; kk3=kk3+N; 
kk4=kk4+Nfft;kk7=kk7+g*(p1+p2);
end

Neb=Neb+sum(sum(de2bi(Xmod4_r,Nbps)~=de2bi(X,Nbps)));
Ntb=Ntb+g*(p1+p2)*Nframe;
if Neb>Target_neb, break; end
end
%if i==0, sigPow= sigPow/Nsym/Nframe/N_iter;
%else
Ber = Neb/Ntb;
fprintf('EbN0=%3d[dB], BER=%4d/%8d =%11.3e\n', EbN0(i),Neb,Ntb,Ber)
fprintf(fid, '%d\t%11.3e\n', EbN0(i), Ber);
if Ber<1e-8, break; end
end

if (fid~=0), fclose(fid); end
%plot_ber_op(file2);
file1=optical2();
%[file1,paprSymboldB1]=optical2();
plot_ber_op2(file1,file2);




