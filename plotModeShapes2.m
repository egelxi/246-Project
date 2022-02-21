function plotModeShapes2(Phin,H,figureNumber,saveOrNot)
% Plots mode shapes
[n,~]=size(Phin);
topHeight=max(H);
figure(figureNumber);

for counter=1:n
    subplot(1,n,counter);
    
    height = H;
    modeShape=zeros(n+1,1);
    
    % Normalizing the mode vectors
    modeShape(2:n+1,1)=Phin(:,counter);
    modeShape=modeShape/max(abs(modeShape));
    % Plot
    plot(modeShape,height,'r');
    hold on;
    plot(modeShape,height,'ob');
    title(['Mode ',num2str(counter)]);
    xlim([-1 1])
    ylim([0 topHeight]);
    set(gca,'ytick', H);
    set(gca,'xtick',-1:0.5:1);
    grid on;
end

if saveOrNot==1
    saveas(figure(figureNumber),'Mode shapes.fig');
    close(figureNumber);
end