function [theta_scores] = get_rot_similarity(ref,img,thetas)

% based on code found at /matlabcentral/newsreader/view_thread/254141
% probably not the right way to do this
% see: http://www.peterkovesi.com/matlabfns/

P1 = radon(ref,0);
P2 = radon(img, thetas);

for tx = 1:length(thetas)
    b = P2(:,tx);
    this_corr = abs(ifft(PT.*conj(b)));
    theta_scores(tx) = max(this_corr(:));   
end

end
