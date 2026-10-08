function []=SolMovie(FolderName,style)
% Movie - solutions along one branch

currentFolder = pwd;

if style==2
    
    figPath = fullfile(currentFolder,'BD3_psigma_case4.fig');  % name of the file of the diag we want to use
    h1 = openfig(figPath,'reuse'); % open figure
    ax1 = gca; % get handle to axes of figure
    fig1 = get(ax1,'children'); 

    figure(3); 
    s1 = subplot(1,2,1); 
    box on
    xlabel('\sigma/d_R')
    ylabel('||R||_{L^1}')
    copyobj(fig1,s1); 
    axis square
    close(h1)
    subplot(1,2,2)
    box on
    xlabel('x')
    ylabel('R,T')
    axis square  
end

cd(FolderName)
ptFileList=dir('pt*.mat');

[~,inx]=sort({ptFileList.date});
ptFileList = ptFileList(inx);
[nn,~]=size(ptFileList);
cd(currentFolder)

% FIX
maxVal = 0;
for i=1:nn
    a=strrep(ptFileList(i).name,'.mat','');
    cd(FolderName); S=load(a); p=S.p; cd(currentFolder);
    maxVal = max(maxVal, max(p.u(1:2*p.np)));
end
yUpper = maxVal*1.1;

VideoName=FolderName;
v=VideoWriter(VideoName,'MPEG-4');
v.FrameRate=2;
open(v);

nRepeat = 2;
nInterp = 3;

allU1 = cell(nn,1); allU2 = cell(nn,1); allDp = zeros(nn,1);

for i=1:nn
    a=ptFileList(i).name;
    a=strrep(a,'.mat','');
    cd(FolderName);
    S=load(a);
    p=S.p;
    cd(currentFolder)
    % u1=p.u(1: p.np);
    % u2=p.u(p.np+1:2*p.np);
    % dp=p.u(2*p.np+10);
    allU1{i}=p.u(1:p.np);
    allU2{i}=p.u(p.np+1:2*p.np);
    allDp(i)=p.u(2*p.np+10);
 end
[po,~ ,~]= getpte(p); %mesh
 
for i=1:nn-1
    u1a=allU1{i};   u1b=allU1{i+1};
    u2a=allU2{i};   u2b=allU2{i+1};
    dpa=allDp(i);   dpb=allDp(i+1);

    for t=linspace(0,1,nInterp)
        u1 = (1-t)*u1a + t*u1b;   % linear interpol. 
        u2 = (1-t)*u2a + t*u2b;
        dp = (1-t)*dpa + t*dpb;

    switch style
        case 1 % solution
            %subplot(1,2,1)
            figure(10)
            hold on
            plot1=plot(po,u1,'color','k','LineWidth',1);
            plot2=plot(po,u2,'color','b','LineWidth',1);
            plot3=plot(po,po*0+4.5,':','color','b','LineWidth',1);
            axis([po(1),po(end),0,10])
            set(gcf,'color','w')
            drawnow;
            frame=getframe(gcf);
            for r=1:nRepeat, writeVideo(v,frame); end
            delete(plot1)
            delete(plot2)
            delete(plot3)
        case 2 % bif diagram on the left
            figure(3) 
            subplot(1,2,1)
            hold on
            % L1norm_R=p.Om*sum(p.mat.M(1:p.np,1:p.np)*abs(u(1:p.np)));
            L1norm_R = p.Om*trapz(po,abs(u1));  % u->u1, uniform grid
            plot1=plot(dp,L1norm_R,'ok','MarkerFaceColor','y');
            subplot(1,2,2)
            hold on
            plot2=plot(po,u1,'color',[31/244,97/244,14/244],'LineWidth',2);
            plot3=plot(po,u2,'color','k','LineWidth',2);
            axis([po(1),po(end),0,yUpper])  
            set(gcf,'color','w')
            drawnow;
            frame=getframe(gcf);
            for r=1:nRepeat, writeVideo(v,frame); end
    
            %
            subplot(1,2,1)
            delete(plot1)
            subplot(1,2,2)
            delete(plot2)
            delete(plot3)
    end

    end
end

close(v);
cd(currentFolder)
end

