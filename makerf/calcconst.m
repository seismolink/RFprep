function sel_data = calcconst(sel_data)

sel_data.dt = 1.0/sel_data.drate;
sel_data.dpts = round(sel_data.tdur*sel_data.drate);
sel_data.dpts2 = round(sel_data.tdur2*sel_data.drate);
sel_data.tprior=sel_data.tdur/100*15;
sel_data.nmt = round(sel_data.tprior*sel_data.drate);
sel_data.npre = round(sel_data.tpre*sel_data.drate);
sel_data.npost = sel_data.dpts2 - sel_data.npre;
sel_data.dln2 = log(sel_data.tdur*sel_data.drate)/log(2.);
sel_data.dpad = round(2.^(ceil(sel_data.dln2)));
sel_data.tpad = sel_data.dpad*sel_data.dt;
sel_data.df = 1.0/sel_data.tpad;
sel_data.nf = round(sel_data.fmax/sel_data.df + 1);

end