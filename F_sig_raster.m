% F_sig_raster() - Plots 2D (electrode x time/frequency) raster diagram of 
%                variables that show significant effects according to F-tests. 
%                
% Usage:
%  >> [img, h_ax] = F_sig_raster(GND_GRP_specGND_or_fname, test_id, varargin);
%
% Required Inputs:
%  GND_GRP_specGND_or_fname  - A GND/GRP/specGND structure variable or the filename of 
%                      such a variable that has been saved to disk.   
%                      To create a GND variable from EEGLAB *.set files use 
%                      sets2GND.m.  To create a GRP structure use GNDs2GRP.m. 
%                      See Mass Univariate ERP Toolbox documentation for 
%                      detailed information about the format of GND and GRP 
%                      variables. If you specifiy a filename 
%                      be sure to include the file's path, unless the file is
%                      in the current working directory.              
%  test_id           - [integer] The index # of the F-test results 
%                      stored in the GND/GRP/specGND variable that you wish to visualize.  
%                      To see what test results are available, look at 
%                      the "F_tests" field of your variable (e.g., GND.F_tests)
%
% Optional Inputs:
%  effect           - The effect within the test to plot. Required for
%                     factorial ANOVA.
%  x_ticks          - [vector] The times/frequencies you want marked with ticks
%                     on the x-axis.  Note, because the EEG has been sampled 
%                     at discrete time points, not all times between the 
%                     minimum and maximum analyzed times/frequencies can be
%                     used.  When a tick is requested that is not available, 
%                     the closest possible value is used. If not
%                     specified, x_ticks are automatically chosen based
%                     on the time/frequency range covered by the diagram. 
%                     This option has no effect if the F-test was
%                     executed based on mean voltages/power within specified 
%                     time windows/frequency bands.
%  fig_id           - [integer] The index # of the MATLAB figure in which
%                     the diagram will be produced.  Useful for overwriting
%                     old figures. {default: lowest unused index}
%  use_color        - ['n','rb','rgb'] If 'n', raster will be black, white,
%                     and grey.  If 'rb', raster will be red (significant) 
%                     and white(nonsignificant). If 'rgb', raster will use a 
%                     temperature scale to indicate graded values of 
%                     significance. {default: 'n'} 
%  plot_vert_lines  - [integer] If non-zero, vertical lines separating all
%                     time points included in the F-test will be 
%                     shown. If 0, the only vertical lines drawn will be 
%                     those that separate times/frequencies included in a test 
%                     from those not included. This option has no effect if 
%                     the F-test was executed based on mean voltage/power 
%                     within specified time windows. {default: 1}
%  lr_sym           - [integer] If non-zero, left hemisphere electrodes are
%                     plot from most anterior to most posterior and right 
%                     hemisphere electrodes are plot from most posterior to 
%                     most anterior.  This may make it easier to visualize 
%                     the degree of effect left-right symmetry.  If 0, all
%                     electrodes are plot from most anterior to most
%                     posterior.
%  scale_limits     - [cmin cmax] The minimum and maximum value of the
%                     colorscale when 'use_color' option is set to 'rgb'. 
%                     If 'use_color' is NOT set to 'rgb', this option has 
%                     no effect. If 'scale_limits' is not specified the 
%                     limits are +/- the maximum absolute value being shown.
%                     This option is useful for making multiple sig_raster
%                     plots on the same scale.
%  fontsize         - [integer] Fontsize of tick mark text. Other text
%                     (e.g., axis labels) are a function of this. {default: 12}
%  mask_opacity     - [0<=scalar<=1] The degree to which the grey masking
%                     of statistically insignificant varialbles will be
%                     transparent. 0=fully transparent. 1=fully opaque.
%                     {default: 1}
%  verblevel        - An integer specifiying the amount of information you want
%                     this function to provide about what it is doing during runtime.
%                     Options are:
%                      0 - quiet, only show errors, warnings, and EEGLAB reports
%                      1 - stuff anyone should probably know
%                      2 - stuff you should know the first time you start working
%                          with a data set {default value}
%                      3 - stuff that might help you debug (show all
%                          reports)
%
% Outputs:
%  img              - Two-dimensional matrix (channel x time/frequency) 
%                     representing the visualized raster.  This could be
%                     useful for making your own personalized raster plot
%                     (e.g., >> imagesc(img);).  Note that the matrix
%                     includes rows for any blank rows separating groups of
%                     channels (e.g., left electrodes from midline
%                     electrodes).
%  h_ax             - h_ax(1) is the handle of the raster axis. h_ax(2) is
%                     the handle of the colorscale (if used).
%         
%
% Global Variable:
%   VERBLEVEL - Mass Univariate ERP Toolbox level of verbosity (i.e., tells functions
%               how much to report about what they're doing during
%               runtime) set by the optional function argument 'verblevel'
%
% Notes:
% -The printed/exported figure will have the same dimensions as the figure
% on the display.  Thus you can undo figure clutter by re-sizing it.
%
% -Some image viewers (e.g., Apple's Preview) may have trouble viewing 
% postscript raster images since they downsample the image to make the file 
% smaller, which blurs the lines.  Ghostview  and Adobe Illustrator should
% be able to view the files just fine.  If viewing the raster as a
% postscript file is producing problems, you can save the figure in jpg 
% format.
%
% -Assignment of electrode to left, midline, or right grouping of
% electrodes is based on the GND.chanlocs(x).theta coordinate.  Anterior to
% posterior organization of electrodes is based on GND.chanlocs(x).radius.
%
%VERSION DATE: 7 December 2017
%AUTHOR: Eric Fields (modified from sig_raster.m by David Groppe)
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
%This function is based on code from the Mass Univariate Toolbox by David Groppe: 
%Copyright (c) 2015, David Groppe: https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox/blob/master/LICENSE

function [img, h_ax] = F_sig_raster(GND_GRP_specGND_or_fname,test_id,varargin)

p=inputParser;
p.addRequired('GND_GRP_specGND_or_fname',@(x) isstruct(x) || ischar(x)); 
p.addRequired('test_id',@(x) isnumeric(x) && (length(x)==1));
p.addParameter('effect', [], @ischar);
p.addParameter('x_ticks',[],@isnumeric);
p.addParameter('fig_id',[],@(x) isnumeric(x) && (length(x)==1));
p.addParameter('use_color','n',@(x) strcmpi(x,'n') || strcmpi(x,'rb') || strcmpi(x,'rgb'));
p.addParameter('units','F',@(x) strcmpi(x,'F'));
p.addParameter('plot_vert_lines',1,@(x) isnumeric(x) || ischar(x));
p.addParameter('lr_sym',0,@(x) isnumeric(x) || ischar(x));
p.addParameter('verblevel',[],@(x) isnumeric(x) && (length(x)==1));
p.addParameter('scale_limits',[],@(x) isnumeric(x) && (length(x)==2));
p.addParameter('unused_color',[1 1 1]*.7,@(x) isnumeric(x) && (length(x)==3)); % this option isn't all that useful actually
p.addParameter('fontsize',12,@(x) isnumeric(x) && (length(x)==1));
p.addParameter('mask_opacity',1,@(x) isnumeric(x) && (length(x)==1)); 

p.parse(GND_GRP_specGND_or_fname,test_id,varargin{:});

% Manage VERBLEVEL
if isempty(p.Results.verblevel),
    VERBLEVEL=2; %not global, just local
else
    VERBLEVEL=p.Results.verblevel;
end

use_color=p.Results.use_color;
plot_vert_lines=str2bool(p.Results.plot_vert_lines);
lr_sym=str2bool(p.Results.lr_sym);

%Load GND, GRP, or specGND struct
if ischar(GND_GRP_specGND_or_fname),
    fprintf('Loading GND, GRP, or specGND struct variable from file %s.\n',GND_GRP_specGND_or_fname);
    load(GND_GRP_specGND_or_fname,'-MAT');
    if ~exist('GND','var') && ~exist('GRP','var') && ~exist('specGND','var')
        error('File %s does not contain a GND, GRP, or specGND variable.',GND_GRP_specGND_or_fname);
    end
    if exist('GRP','var'),
        GND=GRP; %for simplicity GRP variables are re-named GND
        clear GRP;
    elseif exist('specGND','var'),
        GND=specGND; %for simplicity specGND variables are re-named GND
        clear specGND;
    end
    VerbReport(sprintf('Experiment: %s',GND.exp_desc),2,VERBLEVEL);
else
    GND=GND_GRP_specGND_or_fname; %for simplicity GRP and specGND variables are re-named GND
    clear GND_GRP_specGND_or_fname;
end

n_test=length(GND.F_tests);
if test_id>n_test,
   error('There are only %d F-tests stored with these data, but you requested Test #%d',n_test,test_id); 
end

%Extract results for plotting
results = GND.F_tests(test_id);
if isempty(p.Results.effect)
    if isstruct(results.F_obs)
        error('You must specify which effect in the factorial ANOVA you want to plot.');
    else
        null_test   = results.null_test;
        Fval        = results.F_obs;
        pval        = results.adj_pval;
        effect_name = results.factors{1};
        estimated_alpha = results.estimated_alpha;
    end
else
    if ~isfield(results.F_obs, p.Results.effect)
        error('''%s'' does not name an effect in GND.F_tests(%d).', p.Results.effect, test_id); 
    else
        null_test   = results.null_test.(p.Results.effect);
        Fval        = results.F_obs.(p.Results.effect);
        pval        = results.adj_pval.(p.Results.effect);
        effect_name = p.Results.effect;
        if isnan(results.n_perm)
            estimated_alpha = results.estimated_alpha;
        else
            estimated_alpha = results.estimated_alpha.(p.Results.effect);
        end
    end
end

%FDR or permutation test correction for multiple comparisons
if isnan(results.n_perm)
    fdr_crct=1;
    if VERBLEVEL,
       fprintf('Correcting for multiple comparisons via FDR procedure: %s\n', ...
               results.mult_comp_method); 
    end
else
    fdr_crct=0;
    if VERBLEVEL,
       fprintf('Correcting for multiple comparisons via permutation test.\n');  
    end
end

if isfield(results,'freq_band')
%     %create temporary fields to make frequency domain plotting
%     % easily compatible with time domain plot
%     results.time_wind=results.freq_band;
%     results.mean_wind=results.mean_band;
%     results.used_tpt_ids=results.used_freq_ids;
%     GND.grands_t=GND.grands_pow_dB_t;
%     freq_step=GND.freqs(2)-GND.freqs(1);
%     freq_ord=orderofmag(freq_step);
%     ord_pow=log10(freq_ord); %useful for formatting x-tick labels
%     if ord_pow>=0,
%         n_dig_past_dot=0;
%     else
%         n_dig_past_dot=-ord_pow;
%     end
%     GND.time_pts=rnd_orderofmag(GND.freqs);
%     freq_domain=1;
else
    freq_domain=0;
end

%Clean up code by copying some GND.F_tests variables to their own
%variable
use_chans = results.used_chan_ids;
n_wind    = size(results.time_wind, 1);
if strcmpi(results.mean_wind, 'yes'),
    use_tpts=zeros(1,length(GND.time_pts)); %preallocate mem
    mean_pval = pval;
    mean_Fval = Fval;
    pval=repmat(use_tpts,length(results.used_chan_ids),1);
    Fval=pval;
    for w=1:n_wind,
        ids=results.used_tpt_ids{w};
        use_tpts(ids)= ...
            use_tpts(ids)+1;
        if fdr_crct,
            pval(:,ids)=repmat(1-null_test(:,w),1,length(ids));
        else
            pval(:,ids)=repmat(mean_pval(:,w),1,length(ids));
        end
        if strcmpi(p.Results.units,'F')
            Fval(:,ids)=repmat(mean_Fval(:,w),1,length(ids));
        else
%             temp_mn=mean(GND.grands(results.used_chan_ids, ...
%                 GND.F_tests(3).used_tpt_ids{w},use_bin),2);
%             Fval(:,ids)=repmat(temp_mn,1,length(ids));
        end
    end
    if max(use_tpts)>1,
        if freq_domain,
            error('Tested frequency bands overlap. F_sig_raster.m can''t handle this currently.');
        else
            error('Tested time windows overlap. F_sig_raster.m can''t handle this currently.');
        end
    end
    ids=find(use_tpts);
    use_tpts=ids;
    pval=pval(:,ids);
    Fval=Fval(:,ids);
else
    use_tpts = results.used_tpt_ids;
    if fdr_crct,
        pval = 1 - null_test;
    end
end

% Figure out electrode ordering 
n_use_chans = length(use_chans);

% Standardize theta to -180<theta<=180 degrees
for c=1:n_use_chans,
    if ~isfield(GND.chanlocs(c),'theta'),
        watchit(sprintf('Channel #%d (Label: %s) does not have a theta coordinate. To produce raster, I''m assigning it temporary coordinate of: theta=0.',c,GND.chanlocs(use_chans(c)).labels));
        theta=0;
    else
        theta=mod(GND.chanlocs(use_chans(c)).theta,360);
        if isempty(theta),
            watchit(sprintf('Channel #%d (Label: %s) does not have a theta coordinate. To produce raster, I''m assigning it temporary coordinate of: theta=0.',c,GND.chanlocs(use_chans(c)).labels));
            theta=0;
        end
    end
    if theta>180,
        theta=theta-360;
    elseif theta<=-180,
        theta=360+theta;
    end
    GND.chanlocs(use_chans(c)).theta=theta;
    
    if ~isfield(GND.chanlocs(c),'radius'),
        GND.chanlocs(use_chans(c)).radius=2;
        watchit(sprintf('Channel #%d (Label: %s) does not have a radius coordinate. To produce raster, I''m assigning it temporary coordinate of: radius=2.',c,GND.chanlocs(use_chans(c)).labels));
    else
        if isempty(GND.chanlocs(use_chans(c)).radius),
            GND.chanlocs(use_chans(c)).radius=2;
            watchit(sprintf('Channel #%d (Label: %s) does not have a radius coordinate. To produce raster, I''m assigning it temporary coordinate of: radius=2.',c,GND.chanlocs(use_chans(c)).labels));
        end
    end
end

% Get midline electrodes from front to back
midline=[];
midline_dist=[]; %distance from MiCe
for c=1:n_use_chans,
    theta=GND.chanlocs(use_chans(c)).theta;
    radius=GND.chanlocs(use_chans(c)).radius;

    if theta==0 || radius==0,
        midline=[midline use_chans(c)];
        midline_dist=[midline_dist radius];
    elseif theta==180,
        midline=[midline use_chans(c)];
        midline_dist=[midline_dist -radius];
    end
end
[midline_dist, id]=sort(midline_dist,2,'descend');
midline=midline(id);
n_midline=length(midline);

% Get left electrodes from front to back
left=[];
left_dist=[];
for c=1:n_use_chans,
    %standardize theta to a range between +/-180
    theta=GND.chanlocs(use_chans(c)).theta;
    
    if (theta<0) && (theta>-180)
        if GND.chanlocs(use_chans(c)).radius,
            %make sure electrode isn't midline or right
            left=[left use_chans(c)];
            ang=90-abs(theta);
            dist=GND.chanlocs(use_chans(c)).radius*sin(ang*2*pi/360);
            left_dist=[left_dist dist];
        end
    end
end
[left_dist, id]=sort(left_dist,2,'descend');
left=left(id);
n_left=length(left);

% Get right electrodes from front to back
right=[];
right_dist=[];
for c=1:n_use_chans,
    %standardize theta to a range between +/-180
    theta=GND.chanlocs(use_chans(c)).theta;
    
    if (theta>0) && (theta<180)
        if GND.chanlocs(use_chans(c)).radius,
            %make sure electrode isn't midline or right
            right=[right use_chans(c)];
            ang=90-abs(theta);
            dist=GND.chanlocs(use_chans(c)).radius*sin(ang*2*pi/360);
            right_dist=[right_dist dist];
        end
    end
end
% Make right electrodes back to front if requested
if lr_sym,
   right_dist=-right_dist; 
end
[right_dist id]=sort(right_dist,2,'descend');
right=right(id);
n_right=length(right);

% Compute time/frequency tick mark locations and labels
show_tpts=use_tpts(1):use_tpts(end); %show_tpts differs from use_tpts if multiple non-contiguous time windows were used
n_show_tpts=length(show_tpts);
if strcmpi(results.mean_wind,'yes'),
    % Test was based on mean voltage in time window(s)/frequency band(s)
    xtick=NaN;
    n_xtick=0;
    for a=1:n_wind,
        n_xtick=n_xtick+1;
        mn_tpt=mean(results.used_tpt_ids{a});
        xtick(n_xtick)=mn_tpt-show_tpts(1)+1;
        if freq_domain,
            lab_form=['%.' int2str(n_dig_past_dot) 'f-%.' int2str(n_dig_past_dot)  'f'];
            xtick_lab{n_xtick}=sprintf(lab_form,results.time_wind(a,1), ...
                results.time_wind(a,2));
        else
            xtick_lab{n_xtick}=sprintf('%d-%d',results.time_wind(a,1), ...
                results.time_wind(a,2));
        end
        if a~=n_wind,
            skipped_wind=[(results.used_tpt_ids{a}(end)+1): ...
                (results.used_tpt_ids{a+1}(1)-1)];
            if ~isempty(skipped_wind),
                n_xtick=n_xtick+1;
                mn_tpt=mean(skipped_wind);
                xtick(n_xtick)=mn_tpt-show_tpts(1)+1;
                if freq_domain,
                    if length(skipped_wind)>1,
                        lab_form=['%.' int2str(n_dig_past_dot) 'f-%.' int2str(n_dig_past_dot)  'f'];
                        xtick_lab{n_xtick}=sprintf(lab_form,GND.time_pts(skipped_wind(1)), ...
                            GND.time_pts(skipped_wind(end)));
                    else
                        lab_form=['%.' int2str(n_dig_past_dot) 'f'];
                        xtick_lab{n_xtick}=sprintf(lab_form,GND.time_pts(skipped_wind(1)));
                    end
                else
                    if length(skipped_wind)>1,
                        xtick_lab{n_xtick}=sprintf('%d-%d',GND.time_pts(skipped_wind(1)), ...
                            GND.time_pts(skipped_wind(end)));
                    else
                        xtick_lab{n_xtick}=sprintf('%d',GND.time_pts(skipped_wind(1)), ...
                            GND.time_pts(skipped_wind(end)));
                    end
                end
            end
        end
    end
else
    % Test was based on each time point in time window(s)/frequency band(s)
    if isempty(p.Results.x_ticks),
        tm_range=GND.time_pts(use_tpts(end))-GND.time_pts(use_tpts(1));
        omag=orderofmag(tm_range);
        tm_step=round(omag*GND.srate/1000);
        xtick=1:tm_step:n_show_tpts;
        xtick_lab=cell(1,length(xtick));
        for a=1:length(xtick),
            xtick_lab{a}=num2str(GND.time_pts(show_tpts(xtick(a))));
        end
    else
        tk_ct=0;
        xtick=zeros(1,length(p.Results.x_ticks));
        for t=p.Results.x_ticks,
            tk_ct=tk_ct+1;
            xtick(tk_ct)=find_tpt(GND.time_pts(show_tpts),t);
        end
        xtick=unique(xtick); %get rid of any redundant ticks
        xtick_lab=cell(1,length(xtick));
        for a=1:length(xtick),
            xtick_lab{a}=num2str(GND.time_pts(show_tpts(xtick(a))));
        end
    end
end

if VERBLEVEL,
    if fdr_crct,
        fprintf('q level of critical F-scores: %f.\n',results.desired_alphaORq);
    else
        fprintf('Estimated alpha level of critical F-scores: %f.\n',estimated_alpha);
    end
end

%Initialize variables
n_blank=(n_right>0)+(n_left>0)+(n_midline>0)-1;
%if ~strcmpi(use_color,'n'),
if strcmpi(use_color,'rb'),
    img=-.25*ones(n_use_chans+n_blank,n_show_tpts);%makes skipped lines grey (non-sig rectangles are white)
else
    img=zeros(n_use_chans+n_blank,n_show_tpts); %makes skipped lines grey
end
if strcmpi(use_color,'rgb')
    mask=img;
end
chan_lab=cell(1,n_use_chans);
dat.chan_lab=cell(1,n_use_chans+n_blank); %channel labels including empty cells for skipped lines
%dat_chan_lab is needed to make sig_raster plot interactive

%Get IDs of time points/frequencies used in perm test (some shown time points/frequencies may not
%have been included
used_tpt_ids=zeros(1,n_show_tpts);
for a=1:n_show_tpts,
   if ismember(show_tpts(a),use_tpts),
      used_tpt_ids(a)=1; 
   end
end
used_tpt_ids=find(used_tpt_ids);

%Construct raster for left electrodes
ct=0;
skipped=[];
for a=left,
    ct=ct+1;
    elec_id=find(use_chans==a);
    if fdr_crct,
        if strcmpi(use_color,'rgb')
            mask(ct,used_tpt_ids)=(pval(elec_id,:)<=results.desired_alphaORq);
            img(ct,used_tpt_ids)=Fval(elec_id,:);
        else
            img(ct,used_tpt_ids)=(pval(elec_id,:)<=results.desired_alphaORq).*sign(Fval(elec_id,:));
        end
    else
        if strcmpi(use_color,'rgb')
            mask(ct,used_tpt_ids)=(pval(elec_id,:)<results.desired_alphaORq);
            img(ct,used_tpt_ids)=Fval(elec_id,:);
        else
            img(ct,used_tpt_ids)=(pval(elec_id,:)<results.desired_alphaORq).*sign(Fval(elec_id,:));
        end
    end
    chan_lab{ct}=GND.chanlocs(use_chans(elec_id)).labels;
    dat.chan_lab{ct}=GND.chanlocs(use_chans(elec_id)).labels;
end

if ~isempty(midline),
    if ct>0,
        %skip line to separate left and midline electrodes
        ct=ct+1;
        skipped=ct;
    end
    for a=midline,
        ct=ct+1;
        elec_id=find(use_chans==a);
        if fdr_crct,
            if strcmpi(use_color,'rgb')
                mask(ct,used_tpt_ids)=(pval(elec_id,:)<=results.desired_alphaORq);
                img(ct,used_tpt_ids)=Fval(elec_id,:);
            else
                img(ct,used_tpt_ids)=(pval(elec_id,:)<=results.desired_alphaORq).*sign(Fval(elec_id,:));
            end
        else
            if strcmpi(use_color,'rgb')
                mask(ct,used_tpt_ids)=(pval(elec_id,:)<results.desired_alphaORq);
                img(ct,used_tpt_ids)=Fval(elec_id,:);
            else
                img(ct,used_tpt_ids)=(pval(elec_id,:)<results.desired_alphaORq).*sign(Fval(elec_id,:));
            end
        end
        chan_lab{ct-length(skipped)}=GND.chanlocs(use_chans(elec_id)).labels;
        dat.chan_lab{ct}=GND.chanlocs(use_chans(elec_id)).labels;
    end
end

if ~isempty(right),
    if ct>0,
        %skip line to separate right from left or midline electrodes
        ct=ct+1;
        skipped=[skipped ct];
    end
    for a=right,
        ct=ct+1;
        elec_id=find(use_chans==a);
        if fdr_crct,
            if strcmpi(use_color,'rgb')
                mask(ct,used_tpt_ids)=(pval(elec_id,:)<=results.desired_alphaORq);
                img(ct,used_tpt_ids)=Fval(elec_id,:);
            else
                img(ct,used_tpt_ids)=(pval(elec_id,:)<=results.desired_alphaORq).*sign(Fval(elec_id,:));
            end
        else
            if strcmpi(use_color,'rgb')
                mask(ct,used_tpt_ids)=(pval(elec_id,:)<results.desired_alphaORq);
                img(ct,used_tpt_ids)=Fval(elec_id,:);
            else
                img(ct,used_tpt_ids)=(pval(elec_id,:)<results.desired_alphaORq).*sign(Fval(elec_id,:));
            end
        end
        chan_lab{ct-length(skipped)}=GND.chanlocs(use_chans(elec_id)).labels;
        dat.chan_lab{ct}=GND.chanlocs(use_chans(elec_id)).labels;
    end
end


if isempty(p.Results.fig_id)
    fig_h=figure;
else
    fig_h=figure(p.Results.fig_id); clf;
end

set(fig_h, 'Name', effect_name, 'paperpositionmode','auto');

%setting paperpositionmode to 'auto' means that if the figure is manually
%resized, the printed version of the figure will reflect the whatever the
%shown size was (at the time of printing)

if strcmpi(use_color,'rb'),
    %skipped rows/columns=grey, nonsig=white, +sig=red, -sig=blue
    h_img=imagesc(img,[-1 1]);colormap([0 0 1; .5 .5 .5; 1 1 1; 1 0 0]);
elseif strcmpi(use_color,'rgb')
    if verLessThan('matlab','8.4')
        cmap=colormap('jet');
    else
        cmap=colormap('winter');
    end
    colormap(cmap);
    abs_mx=max(max(abs(img)));
    if isempty(p.Results.scale_limits)
        %masked_img=img.*mask;
        %masked_img(masked_img==0)=1e-10; %This is necessary to make masked 
        % values slightly greater than zero, since cmap has an even # of 
        % entries and we need cmap(33,:) to map to 0. 
        %h_img=imagesc(masked_img,[-1 1]*abs_mx);
        h_img=imagesc(img,[0 1]*abs_mx);
    else
        %h_img=imagesc(img.*mask,p.Results.scale_limits);
        h_img=imagesc(img,p.Results.scale_limits);
    end
else
    h_img=imagesc(img,[-1 1]);colormap([0 0 0; .6 .6 .6; 1 1 1]);
end
set(gca,'fontsize',p.Results.fontsize);
h_ax(1)=gca;
if freq_domain,
    bdfcn=['Cp = get(gca,''CurrentPoint''); ' ...
        'Xp=Cp(2,1);', ...
        'Yp=Cp(2,2);', ...
        'dat=get(gcbo,''userdata'');', ...
        'ht=text(Xp,Yp,[num2str(dat.times(round(Xp))) '' Hz, '' dat.chan_lab{round(Yp)}]);' ...
        'set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''buttondownfcn'',''delete(gcbo);'');'];
else
    bdfcn=['Cp = get(gca,''CurrentPoint''); ' ...
        'Xp=Cp(2,1);', ...
        'Yp=Cp(2,2);', ...
        'dat=get(gcbo,''userdata'');', ...
        'ht=text(Xp,Yp,[int2str(dat.times(round(Xp))) '' ms, '' dat.chan_lab{round(Yp)}]);' ...
        'set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''buttondownfcn'',''delete(gcbo);'');'];
end
dat.times=GND.time_pts(show_tpts); %dat.chan_lab has been set above
set(h_img,'buttondownfcn',bdfcn,'userdata',dat);
set(gca,'xtick',xtick,'xticklabel',xtick_lab,'tickdir','out');

set(gca,'ytick',setdiff(1:(n_use_chans+n_blank),skipped),'yticklabel',chan_lab,'box','off');

if freq_domain,
    h=xlabel('Hz');
else
    h=xlabel('Time (ms)');
end
set(h,'fontsize',p.Results.fontsize+2,'fontweight','bold');

h=ylabel('Electrode');
set(h,'fontsize',p.Results.fontsize+2,'fontweight','bold');



h=title(effect_name);

set(h,'fontsize',p.Results.fontsize+3,'fontweight','bold');

hold on;
v=axis;
%rightmost vertical line to cover up weird white lip that peeks out from
%under imagesc
h=plot([1 1]*v(2),v(3:4));
set(h,'color',[1 1 1]*.6);


%% Mask non-significant squares if using rgb
if strcmpi(use_color,'rgb'),
    [n_x, n_y]=size(img);
    % Fill in blank rows that divide different groups of electrodes
    for x=skipped,
        hp=fill([0 n_y+.5 n_y+.5 0],[-0.5 -0.5 0.5 0.5]+x,[1 1 1]*.78);
        set(hp,'linestyle','none');
    end
    for x=1:n_x,
        if ~ismember(x,skipped),
            for y=1:n_y,
                if mask(x,y)==0,
                    %plot(y,x,'m*');
                    hp=fill([-0.5 0.5 0.5 -0.5]+y,[-0.5 -0.5 0.5 0.5]+x,[1 1 1]*.78);
                    set(hp,'facealpha',p.Results.mask_opacity,'linestyle','none');
                end
            end
        end
    end
end

%horizontal lines
show_times=GND.time_pts(show_tpts);
for a=1:n_wind,
    start_tpt=find_tpt(results.time_wind(a,1),show_times);
    stop_tpt=find_tpt(results.time_wind(a,2),show_times);
    for c=0:(n_use_chans+n_blank-1),
        h=plot([start_tpt-.5 stop_tpt+.5],[1 1]*c+.5);
        if ismember(c+1,skipped) || ismember(c,skipped),
            %make lines near boundaries different??  currently not used
            set(h,'color',[0 0 0],'buttondownfcn',bdfcn,'userdata',dat);
        else
            set(h,'color',[0 0 0],'buttondownfcn',bdfcn,'userdata',dat);
        end
    end
end

%vertical lines
if plot_vert_lines && ~strcmpi(results.mean_wind,'yes'),
    for t=1:(n_show_tpts),
        if ismember(show_tpts(t),use_tpts),
            if ismember(t+1,xtick) || ismember(t,xtick),
                %draw vertical lines at time tick marks differently??
                %currently option isn't used
                c=[0 0 0];
            else
                c=[0 0 0];
            end
            if isempty(skipped),
                h=plot([1 1]*t+.5,v(3:4));
                set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                if t==1 || ~ismember(show_tpts(t-1),use_tpts),
                    h=plot([1 1]*t-.5,v(3:4));
                    set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                end
            else
                h=plot([1 1]*t+.5,[v(3) skipped(1)-.5]);
                set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                if length(skipped)==2,
                    h=plot([1 1]*t+.5,[skipped(1)+.5 skipped(2)-.5]);
                    set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                end
                h=plot([1 1]*t+.5,[skipped(end)+.5 v(4)]);
                set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                
                if t==1 || ~ismember(show_tpts(t-1),use_tpts),
                    h=plot([1 1]*t-.5,[v(3) skipped(1)-.5]);
                    set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                    if length(skipped)==2,
                        h=plot([1 1]*t-.5,[skipped(1)+.5 skipped(2)-.5]);
                        set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                    end
                    h=plot([1 1]*t-.5,[skipped(end)+.5 v(4)]);
                    set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                end
            end
        end
    end
else
    c=[0 0 0];
    for w=1:n_wind,
        tt(1)=find_tpt(show_times,results.time_wind(w,1));
        tt(2)=find_tpt(show_times,results.time_wind(w,2));
        ct=-1;
        for t=tt,
            ct=ct+1;
            if isempty(skipped),
                if ct,
                    t_plt=[1 1]*t+.5;
                else
                     t_plt=[1 1]*t-.5;
                end
                h=plot(t_plt,v(3:4));
                set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
            else
                if ct,
                    t_plt=[1 1]*t+.5;
                else
                    t_plt=[1 1]*t-.5;
                end
                
                h=plot(t_plt,[v(3) skipped(1)-.5]);
                set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                if length(skipped)==2,
                    h=plot(t_plt,[skipped(1)+.5 skipped(2)-.5]);
                    set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
                end
                h=plot(t_plt,[skipped(end)+.5 v(4)]);
                set(h,'color',c,'buttondownfcn',bdfcn,'userdata',dat);
            end
        end
    end
end


    
%% Draw colorbar is using rgb
if strcmpi(use_color,'rgb'),
    cbar_ax=cbarDG();
    set(gca,'xtick',[]);
    if verLessThan('matlab','8.4')
        axes(cbar_ax);
    end
    v=axis;
    axis(v);
    rngX=v(2)-v(1);
    rngY=v(4)-v(3);
    if strcmpi(p.Results.units,'F')
        ht=text(v(2)+rngX*.8,mean(v(3:4)),'F');
        if verLessThan('matlab','8.4')
            set(ht,'fontsize',p.Results.fontsize+2,'rotation',0,'fontweight','bold', ...
                'horizontalalignment','left','verticalalignment','middle');
        else
            set(ht,'fontsize',p.Results.fontsize,'rotation',0, ...
                'horizontalalignment','left','verticalalignment','middle');
        end
    else
%         ht=text(v(2)+rngX*.8,rngY*.005,'\muV');
%         if verLessThan('matlab','8.4')
%             set(ht,'fontsize',p.Results.fontsize+2,'rotation',0,'fontweight','bold', ...
%                 'horizontalalignment','left','verticalalignment','middle');
%         else
%             set(ht,'fontsize',p.Results.fontsize,'rotation',0, ...
%                 'horizontalalignment','left','verticalalignment','middle');
%         end
    end
    h_ax(2)=cbar_ax;
end


%
%% %%%%%%%%%%%%%%%%%%%%% function orderofmag() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ord=orderofmag(val)
%function ord=orderofmag(val)
%
% Returns the order of magnitude of the value of 'val' in multiples of 10
% (e.g., 10^-1, 10^0, 10^1, 10^2, etc ...)
% used for computing erpimage trial axis tick labels as an alternative for
% plotting sorting variable

val=abs(val);
if val>=1
    ord=1;
    val=floor(val/10);
    while val>=1,
        ord=ord*10;
        val=floor(val/10);
    end
    return;
else
    ord=1/10;
    val=val*10;
    while val<1,
        ord=ord/10;
        val=val*10;
    end
    return;
end


%
%% %%%%%%%%%%%%%%%%%%%%% function str2bool() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function bool=str2bool(str)
%function bool=str2bool(str)

if ischar(str),
    if strcmpi(str,'yes') || strcmpi(str,'y')
        bool=1;
    else
        bool=0;
    end
else
   bool=str; 
end
