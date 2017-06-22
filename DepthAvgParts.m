function DepthAvgParts(fname, pname)
% Hourly depth bin averaged Pb for particles

cd([pname,'/Output']);

load(fname)

PbpredNDR = Pbpot.*PinhNDR*3600;
PbpredEr1 = Pbpot.*squeeze(PinhEr(:,:,1))*3600;
PbpredEr2 = Pbpot.*squeeze(PinhEr(:,:,2))*3600;

% Assess trajectories. 
zmax = max(max(zp)); %max depth of particles

if ceil(zmax) == round(zmax)
    zbin = 0:0.5:ceil(zmax);
else
    zbin = 0:0.5:round(zmax)+0.5;
end

hr = ceil(time/3600); hr(1) = 1; hr = repmat(hr, [1, size(zp,2)]);
Pbhrz = nan*ones(24, length(zbin)); Pbhrzu = nan*ones(24, length(zbin)); PbhrzNDR = nan*ones(24, length(zbin));
PbhrzEr = nan*ones(24, length(zbin), 2);
for i = 1:24;
    for j = 1:length(zbin)-1
    ix = find(hr == i & zp>=zbin(j) & zp<zbin(j+1));
    Pbhrz(i,j) = mean(Pbpred(ix));
    PbhrzNDR(i,j) = mean(PbpredNDR(ix));
    Pbhrzu(i,j) = mean(Pbpot(ix)*3600);  %/s --> /hr
    PbhrzEr(i,j, 1) = mean(PbpredEr1(ix));
    PbhrzEr(i,j, 2) = mean(PbpredEr2(ix));
    end
end

figure; set(gcf, 'position', [73 553 1482 570]); xl = [-0.2 2]; yl = [-35 0.2];
subplot(161); hold on; plot(Pbhrzu(1:4,:), -zbin, 'k--'); plot(Pbhrz(1:4,:), -zbin);  plot(PbhrzNDR(1:4,:), -zbin, ':');  ylim(yl); xlim(xl); title('0-4 hr'); ylabel('Depth (m)');
subplot(162); hold on; plot(Pbhrzu(5:8,:), -zbin, 'k--'); plot(Pbhrz(5:8,:), -zbin); plot(PbhrzNDR(5:8,:), -zbin, ':');   ylim(yl); xlim(xl); title('5-8 hr');
subplot(163); hold on; plot(Pbhrzu(9:12,:), -zbin, 'k--'); plot(Pbhrz(9:12,:), -zbin); plot(PbhrzNDR(9:12,:), -zbin, ':');  ylim(yl); xlim(xl); title('9-12 hr'); xlabel('     Pbpred');
subplot(164); hold on; plot(Pbhrzu(13:16,:), -zbin, 'k--'); plot(Pbhrz(13:16,:), -zbin); plot(PbhrzNDR(13:16,:), -zbin, ':'); ylim(yl); xlim(xl); title('13-16 hr'); xlabel('(gC/gChl/hr)');
subplot(165); hold on; plot(Pbhrzu(17:20,:), -zbin, 'k--'); plot(Pbhrz(17:20,:), -zbin); plot(PbhrzNDR(17:20,:), -zbin, ':'); ylim(yl); xlim(xl); title('17-20 hr');
subplot(166); hold on; plot(Pbhrzu(21:24,:), -zbin, 'k--'); plot(Pbhrz(21:24,:), -zbin); plot(PbhrzNDR(21:24,:), -zbin, ':'); ylim(yl); xlim(xl); title('21-24 hr');

figure; set(gcf, 'position', [73 553 1482 570]); xl = [-0.2 2]; yl = [-35 0.2];
subplot(141); hold on; plot(squeeze(PbhrzEr(7,:,1)), -zbin, 'g--'); plot(Pbhrz(7,:), -zbin); plot(squeeze(PbhrzEr(7,:,2)), -zbin, 'r--');   ylim(yl); xlim(xl); title('7 hr'); ylabel('Depth (m)');
subplot(142); hold on; plot(squeeze(PbhrzEr(11,:,1)), -zbin, 'g--'); plot(Pbhrz(11,:), -zbin); plot(squeeze(PbhrzEr(11,:,2)), -zbin, 'r--');  ylim(yl); xlim(xl); title('11 hr'); xlabel('     Pbpred');
subplot(143); hold on; plot(squeeze(PbhrzEr(13,:,1)), -zbin, 'g--'); plot(Pbhrz(13,:), -zbin); plot(squeeze(PbhrzEr(13,:,2)), -zbin, 'r--'); ylim(yl); xlim(xl); title('13 hr'); xlabel('(gC/gChl/hr)');
subplot(144); hold on; plot(squeeze(PbhrzEr(16,:,1)), -zbin, 'g--'); plot(Pbhrz(16,:), -zbin); plot(squeeze(PbhrzEr(16,:,2)), -zbin, 'r--'); ylim(yl); xlim(xl); title('16 hr');
legend('Rmax+50%', num2str(Rmax), 'Rmax-50%') 


hour = 1:24;
fnameout = [fname(1:end-4),'_hrzavg'];

save(fnameout, 'Pbhrz', 'Pbhrzu', 'PbhrzNDR','zbin', 'hour', 'PbhrzEr')

cd(pname)

% load Pb_Tmod
% figure; set(gcf, 'position', [73 553 1482 570]); xl = [-0.2 2]; yl = [-35 0.2];
% subplot(161); hold on; plot(Pbi_hr(1:4,:), -z, '--'); plot(Pbhrz(1:4,:), -zbin);   ylim(yl); xlim(xl); title('0-4 hr'); ylabel('Depth (m)');
% subplot(162); hold on; plot(Pbi_hr(5:8,:), -z, '--'); plot(Pbhrz(5:8,:), -zbin);   ylim(yl); xlim(xl); title('5-8 hr');
% subplot(163); hold on; plot(Pbi_hr(9:12,:), -z, '--'); plot(Pbhrz(9:12,:), -zbin); ylim(yl); xlim(xl); title('9-12 hr'); xlabel('     Pbpred');
% subplot(164); hold on; plot(Pbi_hr(13:16,:), -z, '--'); plot(Pbhrz(13:16,:), -zbin); ylim(yl); xlim(xl); title('13-16 hr'); xlabel('(gC/gChl/hr)');
% subplot(165); hold on; plot(Pbi_hr(17:20,:), -z, '--'); plot(Pbhrz(17:20,:), -zbin); ylim(yl); xlim(xl); title('17-20 hr');
% subplot(166); hold on; plot(Pbi_hr(21:24,:), -z, '--'); plot(Pbhrz(21:24,:), -zbin); ylim(yl); xlim(xl); title('21-24 hr');




% examine particles stuck below LC 
% for j=1:size(zp,2)
%     check(j) = any(zp(:,j) <=zeu); %also tried 25m and 20m. see notes. 
% end
% stuck = find(check == 0);
% 
% figure;
% subplot(211); plot(time/3600, -zp(:,stuck)); ylabel('Depth');
% subplot(212); plot(time/3600, Pbpred(:,stuck)); ylabel('Pb pred');




