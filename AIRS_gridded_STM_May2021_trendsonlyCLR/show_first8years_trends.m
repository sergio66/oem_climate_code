a01 = load('iType_13_iQAX_4_convert_sergio_clearskygrid_obsonly_Q01.mat');
a02 = load('iType_13_iQAX_4_convert_sergio_clearskygrid_obsonly_Q02.mat');
a03 = load('iType_13_iQAX_4_convert_sergio_clearskygrid_obsonly_Q03.mat');

i1231 = find(a01.h.vchan >= 1231,1);
i2616 = find(a01.h.vchan >= 2616,1);

strow.a01.desc(:,:,1) = squeeze(a01.b_desc(:,:,i1231));
strow.a01.desc(:,:,2) = squeeze(a01.b_desc(:,:,i2616));
strow.a01.asc(:,:,1)  = squeeze(a01.b_asc(:,:,i1231));
strow.a01.asc(:,:,2)  = squeeze(a01.b_asc(:,:,i2616));

strow.a02.desc(:,:,1) = squeeze(a02.b_desc(:,:,i1231));
strow.a02.desc(:,:,2) = squeeze(a02.b_desc(:,:,i2616));
strow.a02.asc(:,:,1)  = squeeze(a02.b_asc(:,:,i1231));
strow.a02.asc(:,:,2)  = squeeze(a02.b_asc(:,:,i2616));

strow.a03.desc(:,:,1) = squeeze(a03.b_desc(:,:,i1231));
strow.a03.desc(:,:,2) = squeeze(a03.b_desc(:,:,i2616));
strow.a03.asc(:,:,1)  = squeeze(a03.b_asc(:,:,i1231));
strow.a03.asc(:,:,2)  = squeeze(a03.b_asc(:,:,i2616));

junk = squeeze(strow.a01.desc(:,:,1)-strow.a01.desc(:,:,2));
aslmap(1,rlat65,rlon73,smoothn((junk'),1), [-90 +90],[-180 +180]);  title('AVGobs DESC 1231-2616');  colormap(jett); caxis([-1 1]/4); colormap(usa)
junk = squeeze(strow.a01.asc(:,:,1)-strow.a01.asc(:,:,2));
aslmap(2,rlat65,rlon73,smoothn((junk'),1), [-90 +90],[-180 +180]);  title('AVGobs ASC 1231-2616');  colormap(jett); caxis([-1 1]/4); colormap(usa)

junk = squeeze(strow.a02.desc(:,:,1)-strow.a02.desc(:,:,2));
aslmap(1,rlat65,rlon73,smoothn((junk'),1), [-90 +90],[-180 +180]);  title('Q0.03 DESC 1231-2616');  colormap(jett); caxis([-1 1]/2); colormap(usa)
junk = squeeze(strow.a02.asc(:,:,1)-strow.a02.asc(:,:,2));
aslmap(2,rlat65,rlon73,smoothn((junk'),1), [-90 +90],[-180 +180]);  title('Q0.03 ASC 1231-2616');  colormap(jett); caxis([-1 1]/2); colormap(usa)

junk = squeeze(strow.a03.desc(:,:,1)-strow.a03.desc(:,:,2));
aslmap(1,rlat65,rlon73,smoothn((junk'),1), [-90 +90],[-180 +180]);  title('Q0.97 DESC 1231-2616');  colormap(jett); caxis([-1 1]/10); colormap(usa)
junk = squeeze(strow.a03.asc(:,:,1)-strow.a03.asc(:,:,2));
aslmap(2,rlat65,rlon73,smoothn((junk'),1), [-90 +90],[-180 +180]);  title('Q0.97 ASC 1231-2616');  colormap(jett); caxis([-1 1]/10); colormap(usa)

save strow_first8years_trends.mat strow
