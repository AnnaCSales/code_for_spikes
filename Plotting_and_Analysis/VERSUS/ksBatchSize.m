ks=get(gcf, 'UserData')

%change batch size. This must be an integer*32  + ntbuff. Usually the
%number is 2048.  Increase batch size to address memory problems. 

ks.ops.NT=(2300*32 )  + ks.ops.ntbuff;

ks.ops.throw_out_channels=0;

ks.ops.minFR=0.0002; %want to have clusters that are sparsely firing
ks.ops.AUCsplit=0.4; %allow KS to split more clusters. (0.9 default)
% Threshold on the area under the curve (AUC) criterion for performing a 
% split in the final step. If the AUC of the split is higher than this,
% that split is considered good. 

% Thresholds on spike detection used during the optimization Th(1) or during 
% the final pass Th(2). These thresholds are applied to the template projections,
% not to the voltage. Typically, Th(1) is high enough that the algorithm only
% picks up sortable units, while Th(2) is low enough that it can pick all of 
% the spikes of these units. Default [10, 4] 

ks.ops.Th=[10, 4];


% The individual spike amplitudes are biased towards the mean of the cluster
% by this factor; 50 is a lot, 0 is no bias.
% at 50 - 4 multi unit, 1 clear, very few noise clusters.
ks.ops.lam = 10 %(10 default).