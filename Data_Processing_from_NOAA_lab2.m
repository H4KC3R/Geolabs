%% Считываем сигнал
clear;

[y, fs] = audioread("content/19000020.wav");
t = (1 : length(y))' / fs;
f = (1 : length(y))' / length(y) * fs;
plot(t,y); 
xlabel('Время,с'); 
ylabel('Амлитуда');
title('АМ модулированный сигнал')

%% Первичное изображение

l = length(y);
n = floor(fs/2);
m = floor(l/n);
B = zeros(m,n);
for i = 1:l
    B(floor(i/n)+1,mod(i,n)+1) = y(i);
end
figure("name", "Первичное изображение");
imshow(B);

%% Находим несущую частоту

ffty = fft(y);

figure();
subplot(2, 2, 1);
plot(f, abs(ffty));
title("Исходный спектр");
xlim([0, fs / 2]);
ylim([0, 10000]);
xlabel("Частота, Гц");
ylabel("Амплитуда");

[M, fc] = max(abs(ffty(1 : floor((length(y) / 2)))));
fc = fc / length(ffty) * fs;
disp("fc = " + fc);

%% Преобразуем сигнал

y = y .* sign(y);

%% Фильтруем нижние частоты

ffty = fft(y);
ff = floor(fc / 2 / fs * length(ffty));
ffty(ff : length(ffty) - ff) = 0;

subplot(2, 2, 2);
plot(f, abs(ffty));
title("Новый спектр");
xlim([0, fs / 2]);
ylim([0, 10000]);
xlabel("Частота, Гц");
ylabel("Амплитуда");

%% Обратное преобразование и нормировка

y = abs(ifft(ffty));
y = y / max(y);
plot(t,y); 
xlabel('Seconds'); 
ylabel('Amplitude');
%% Находим синхроимпульсы

syncI = round((0 : 14)' / 2080 * fs);
syncF = [0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0; 1; 0];
syncT = zeros(floor(t(length(t))) * 2, 1);
h = 0;
%an = (sign(y - 0.5) + 1) / 2;
an = ((y > 0.7) - (y < 0.3) + 1) / 2;
for sync = 1 : (length(an) - max(syncI))
    if syncF(:) == an(syncI + sync)
        if h > 0
            h = h + round((sync - syncT(h)) / fs * 2);
        else
            h = 1;
        end
        syncT(h) = sync;
    end
end
syncT = syncT(1 : h);

syncPlot = zeros(h, 1);
j = 0;
for i = 1 : h
    if syncT(i) > 0
        j = j + 1;
        syncPlot(j) = syncT(i);
    end
end
syncPlot = syncPlot(1 : j);

subplot(2, 2, 3);
plot(t(syncPlot));
title("Найденные синхроимпульсы");
xlabel("Номер импульса");
ylabel("Время от начала записи, с");

%% Восстанавливаем ненайденные синхроимпульсы

last = 1;
for i = 2 : h
    if syncT(i) > 0 && i - last > 1
        for j = (last + 1) : (i - 1)
            syncT(j) = round(syncT(last) + (syncT(i) - syncT(last)) * (j - last) / (i - last));
        end
        last = i;
    end
end

subplot(2, 2, 4);
plot(t(syncT));
title("Найденные и восстановленные синхроимпульсы");
xlabel("Номер импульса");
ylabel("Время от начала записи, с");

%% Считываем данные из записи

w = 2080;
im = zeros(h, w);
for j = 1 : (h - 1)
    for i = 0 : (w - 1)
        im(j, i + 1) = y(round(syncT(j) + (syncT(j + 1) - syncT(j)) * i / w));
    end
end

%% Телеметрия А

telemetry = w / 2 - 10;
[telemetryMax, telemetryInd] = max(im(:, telemetry));
telemetryMin = im(telemetryInd + 8, telemetry);
telemetryK = 0.87 / (telemetryMax - telemetryMin);
telemetryB = - telemetryK * telemetryMin;
im(:, 1 : (w / 2)) = telemetryK * im(:, 1 : (w / 2)) + telemetryB;

%% Телеметрия B

telemetry = w - 10;
[telemetryMax, telemetryInd] = max(im(:, telemetry));
telemetryMin = im(telemetryInd + 8, telemetry);
telemetryK = 0.87 / (telemetryMax - telemetryMin);
telemetryB = - telemetryK * telemetryMin;
im(:, (w / 2) : w) = telemetryK * im(:, (w / 2) : w) + telemetryB;

%% Вывод изображения

figure("name", "result");
imshow(im);

%% Сохраняем результат

imwrite(im, "result.bmp");

%%
photo = imread("pic.bmp");
figure('name', 'Initial picture');
imshow(photo);
[h, w] = size(photo);

%%

imTelemetry = photo(:, 2040 : 2070);
figure('name', 'Telemetry');
imshow(imTelemetry);

%%

im = photo(300 : 800, 1127 : 2035);
figure('name', 'Picture');
imshow(im);

%%

telemetry = zeros(1, h);
for i = 1 : h
    telemetry(i) = sum(imTelemetry(i, :)) / 40;
end

figure();
plot(telemetry);

%%

ind = 731;

for i = 1 : 16
    telemetry(i) = sum(telemetry((ind + (i - 1) * 8) : (ind + (i - 1) * 8 + 7))) / 8;
end
telemetry = telemetry(1 : 17);
telemetry(17) = telemetry(1);

figure();
stairs(telemetry);

%%

sz = size(im);

d0 = [276.6067; 276.6119; 276.6311; 276.6168];
d1 = [0.051111; 0.05101; 0.051033; 0.051058];
d2 = [1.405783e-6; 1.496037e-6; 1.496990e-6; 1.493110e-6];

CPRT = telemetry(2 : 5)';
disp('CPRT:');
disp(CPRT);
TPRT = d0 + d1 .* CPRT + d2 .* (CPRT .* CPRT);
disp('TPRT:');
disp(TPRT);

TBB = sum(TPRT) / 4;
disp('TBB:');
disp(TBB);

A = 0.53959;
B = 0.998534;
VC = 928.9;

TBB_star = A + B * TBB;
disp('TBB_star:');
disp(TBB_star);

C1 = 1.1910427e-5;
C2 = 1.4387752;

NBB = C1 * VC * VC * VC / (exp(C2 * VC / TBB_star) - 1);
disp('NBB:');
disp(NBB);

NS = - 5.49;

CSsz = size(photo(451 : 568, 1081 : 1125));
CS = sum(sum(photo(451 : 568, 1081 : 1125))) / CSsz(1) / CSsz(2);
CBB = telemetry(7);

NLIN = NS + (NBB - NS) * (CS - double(im(:, :))) / (CS - CBB);

b0 = 5.7;
b1 = -0.11187;
b2 = 0.00054668;
NCOR = b0 + b1 * NLIN + b2 * (NLIN .* NLIN);

NE = NLIN + NCOR;

TE_star = C2 .* VC ./ log(1 + C1 * VC * VC * VC ./ NE);

TE = (TE_star - A) / B - 273.15;

%%

figure();
imagesc(abs(TE), [0 40]);
colorbar
