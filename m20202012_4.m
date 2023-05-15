%%% Project 3:
% Name: Aibek Bekbergen, ID: 20202012
close all; clear all;  clc;

%% Problem 1.1
load m20202012_4.mat

H1 = freqz(B1,A1);
% figure; plot(abs(h1)); % this is to obtain linear scale, not necessary
figure; semilogy(abs(H1), 'b'); hold on % dB scale

H2 = freqz(B2,A2);
% figure; plot(abs(h2)); % this is to obtain linear scale, not necessary
semilogy(abs(H2), 'y'); hold on % dB scale

H3 = freqz(B3,A3);
% figure; plot(abs(h3)); % this is to obtain linear scale, not necessary
semilogy(abs(H3), 'r'); hold on % dB scale

H4 = freqz(B4,A4);
% figure; plot(abs(h4)); % this is to obtain linear scale, not necessary
semilogy(abs(H4), 'g'); hold on % dB scale

H5 = freqz(B5,A5);
% figure; plot(abs(h5)); % this is to obtain linear scale, not necessary
semilogy(abs(H5), 'c'); hold on % dB scale

H6 = freqz(B6,A6);
% figure; plot(abs(h6)); % linear scale
semilogy(abs(H6), 'm'); title('Magnitude response'); hold off 
ylim([0.01 10]);

%% Problem 1.2
clear all; close all; clc;

load drum.mat
load m20202012_4.mat

y1 = filter(B1, A1, a);
y2 = filter(B2, A2, a);
y3 = filter(B3, A3, a);
y4 = filter(B4, A4, a);
y5 = filter(B5, A5, a);
y6 = filter(B6, A6, a);
y = y1+y2+y3+y4+y5+y6;

sound(y, Fs);
y = y(100000 : 100000+35536);
figure; plot(y); title("drum.mat");
figure; freqz(y,1,35536); title("Magnitude Response");

%%  Problem 1.3
clear all; close all; clc;

load drum.mat
load m20202012_4.mat

H1 = freqz(B1*2,A1);
% figure; plot(abs(H1)); % this is to obtain linear scale, not necessary
figure; semilogy(abs(H1), 'b'); hold on % dB scale

H2 = freqz(B2,A2);
% figure; plot(abs(H2)); % this is to obtain linear scale, not necessary
semilogy(abs(H2), 'y'); hold on % dB scale

H3 = freqz(B3*0.5,A3);
% figure; plot(abs(H3)); % this is to obtain linear scale, not necessary
semilogy(abs(H3), 'r'); hold on % dB scale

H4 = freqz(B4*0.5,A4);
% figure; plot(abs(H4)); % this is to obtain linear scale, not necessary
semilogy(abs(H4), 'b'); hold on % dB scale

H5 = freqz(B5,A5);
% figure; plot(abs(H5)); % this is to obtain linear scale, not necessary
semilogy(abs(H5), 'm'); hold on % dB scale

H6 = freqz(B6*2,A6);
% figure; plot(abs(H6)); % this is to obtain linear scale, not necessary
semilogy(abs(H6), 'c'); title('Magnitude Response'); hold off % dB scale
ylim([0.01 10]);

%% Problem 1.4
y1 = filter(B1*10, A1, a);
y2 = filter(B2/2, A2, a);
y3 = filter(B3/5, A3, a);
y4 = filter(B4*6, A4, a);
y5 = filter(B5/14, A5, a);
y6 = filter(B6*9, A6, a);

y = y1+y2+y3+y4+y5+y6;
sound(y, Fs);

y = y(100000 : 100000+35536);
figure; plot(y); title("Input with Gain");
figure; freqz(y,1,35536); title("Magnitude Response");


%% Problem 2.1

clear all; close all; clc;

h = [1,1,1,1]/4; % our filter
x = [ones(1,32),zeros(1,32),zeros(1,32),ones(1,32)]; % input signal

% create zeros of length 32 * 4 = 128
q = zeros(1, 128);

% we create boundaries  conditions by adding q to the left 
% and right of our input signal
w = [q,x,q];

y = zeros(1, 160);

% now convolve our signal and filter
for n = 1 : 160
    sum = 0;
    for k = 1 : 4
        sum = sum + h(k)*w(n-k+128+1);
    end
    y(n) = sum;
end
figure; stem(y); title("convolution by for loop")
%% Problem 2.2
Hw = fft(h, 160);
Xw = fft(x, 160);

Yw = Hw.*Xw;
y = real(ifft(Yw));
figure; stem(y); title("convolution by fft/ifft");

%% Problem 2.3
w2  = [x, x, x]; % since input is now periodic

y3 = zeros(1, 128);

for n = 1 : 128
    sum = 0;
    for k = 1 : 4
        sum = sum + h(k)*w2(n-k+128+1);
    end
    y3(n) = sum;
end
figure; stem(y3(1:128)); title("convolution by using for loop");

%% Problem 2.4
Hw2 = fft(h, 128);
Xw2 = fft(w2, 128);
result = ifft(Hw2.*Xw2);
figure; stem(result); title("convolution by fft/ifft");
%% Problem 3.1

clear all; close all; clc; 

L = 1025;
M = 512;

P = M*3 + L;
w = hann(L);
w1 = zeros(P, 1);
w2 = zeros(P, 1);
w3 = zeros(P, 1);
w4 = zeros(P, 1);

w1((1:L)+3*M) = w;
w2((1:L)+2*M) = w;
w3((1:L)+M) = w;
w4((1:L)) = w;

x = w1 + w2 + w3 + w4;
t = 1:P;

figure; plot(t, w1, t, w2, t, w3, t, w4, t, x)

V = axis;
V(4) = 1.5;
axis(V)
%% Problem 3.2
clear all; close all; clc;

load flute.mat
figure; plot(a)

N = length(a);
L = 1025;
P = 512;
w = hann(L);

K = length(1:P:N)-2;
S = zeros(2048, K);

for p=1:K
    x = a((1:L) + (p-1)*P);
    x = w.*x;
    X = fft(x, 2048);
    S(:, p) = X;
end

figure; mesh(abs(S(1:1024,:))); view([0,0,1]);

%% Problem 3.3:

y = zeros(N, 1);

for p = 1:K
    X = S(:, p);
    x = real(ifft(X));
    y((1:L) + (p-1)*P) = y((1:L) + (p-1)*P) + x(1:L);
end

figure; plot(y)

figure; plot(y - a);

%% Problem 3.4
clear all; close all; clc;

load flute.mat
figure; plot(a);

N = length(a);
x = (1:N)/Fs;
[s,w,t] = spectrogram(a, hann(1025), 512, 2048, Fs, 'yaxis');

% % to find the sample of the first note with maximum peak
% for i = 1:1:50
%     time = i;
%     S = s(:, time);
%     figure(2); stem(w, abs(S));
%     title(i)
%     pause(1);
% end
% % the 25th sample has the maximum peak
% 
% % to find the sample of the fifth note with maximum peak
% for i = 210:1:250
%     time = i;
%     S = s(:, time);
%     figure(2); stem(w, abs(S));
%     title(i)
%     pause(1);
% end
% % the 210th sample has the maximum peak

figure(1); plot(x, a)
figure(2); plot(t)

time1 = 25;
time2 = 210;

S1 = s(:, time1);
S2 = s(:, time2);

figure(3); spectrogram(a, hann(1025), 512, 2048, Fs, 'yaxis');
figure(4); plot(w, abs(S1)); title("First Note")
figure(5); plot(w, abs(S2)); title("Fifth Note")
figure(6); plot(w, abs(S1), w, abs(S2))

t(time1)
t(time2)
% sound(a, Fs);

% first note is B4
% fifth note is F#5/Gb5