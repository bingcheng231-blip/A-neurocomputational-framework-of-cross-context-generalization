%function it4=compit4
function it4=wd_compit4_1

colormap(gray(256));
%orient tall
it=jet;
it4=jet;
for t=1:170
        it4(t,3)=it(t+32,2);
end
for l=129:202
        it4(l,3)=(202-l)/73;
end
for t=44:160
        it4(t,1)=it(t-26,2);
end
for l=225:256
   it4(l,1)=(256-l)/96+65/96;
end
for l=1:128
        it4(l,1)=l/128;
end

for l=25:100
        it4(l,2)=(l-25)/75;
end
for l=144:249
        it4(l,2)=(249-l)/105;
end



