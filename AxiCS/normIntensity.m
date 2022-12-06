clear all;

aaa = [];

%% 
for ii = 1:size(aaa,2)
   aaa(:,ii) = aaa(:,ii)./max(aaa(:,ii)) ;
end