
function out = tesa_medfilt(data,win,filtOrd)

out = [];
for a = 1:size(win,2);
    if mod(filtOrd,2) == 0
        filtWin = data(:,win(1,a)-30/2 : win(1,a)+filtOrd/2-1);
    else
        filtWin = data(:,win(1,a)-(filtOrd-1)/2 : win(1,a)+(filtOrd-1)/2);
    end
    out(:,a) = median(filtWin')';
end