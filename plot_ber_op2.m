function plot_ber_op2(file1,file2)
a1= load(file1); semilogy(a1(:,1),a1(:,2),'-.r*'); hold on
a2= load(file2); semilogy(a2(:,1),a2(:,2),'--mo'),grid on
legend('optical ofdm','optical ofdm with index modulation');
xlabel('EbN0[dB]'), ylabel('BER'); axis([a1(1,1) a1(end,1) 1e-4 1])
