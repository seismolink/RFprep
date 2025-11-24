function hf=fill_trace(a,b,c,d,e,f,max_trace_ratio)

% []=fill_trace(tx,trace,is_mean,is_neg,couleur,couleur_contour,max_trace_ratio)
%
%
% version 2.1 : Janvier 2001, modifie Mars 2003, par A. KAVIANI
%
% repr�sentation d'une trace sismique en coloriant
% la partie positive du signal.
% tx : vecteur base de temps
%      !! INUTILE si trace est une structure SAC !!!
% trace : vecteur des donn�es (structure SAC ou vecteur)
% is_mean : 'm' si on enl�ve la moyenne (optionnel / d�faut : 'm')
%	    [] si on ne fait rien
%           val si on rempli les valeurs au dessus de "val"
% is_neg : 1 si on veut colorier la partie n�gative (optionnel / d�faut : 0)
% couleur : couleur du remplissage (optionnel / d�faut : noire)
% couleur_contour : couleur du contour (optionnel / d�faut : idem couleur)
%
% RQ : des points sont rajout�s sur la courbe pour mieux tenir compte des
% passages � 0
% max_trace_ratio: To adjust the fill level from zero base. defualt 100.
% Initialisations
hf=[];
%return
if isstruct(a)
    trace=a.trace;
    beg_time=0;
    if isfield(a,'bt');beg_time=a.bt;end
    tx=[0:1:length(trace)-1]*a.tau+beg_time;
    if nargin<2,is_mean='m';else is_mean=b;end
    if nargin<3,is_neg=0;else is_neg=c;end
    if nargin<4,couleur='k';else couleur=d;end
    if nargin<5,couleur_contour=couleur;else couleur_contour=e;end
    if nargin<6, max_trace_ratio=100;end
else
    tx=a;trace=b;
    if nargin<3,is_mean='m';else is_mean=c;end
    if nargin<4,is_neg=0;else is_neg=d;end
    if nargin<5,couleur='k';else couleur=e;end
    if nargin<6,couleur_contour=couleur;else couleur_contour=f;end
    if nargin<7, max_trace_ratio=100;end
end

% Tests sur la taille des vecteurs
if size(tx,1)~=size(trace,1)
    trace=trace';
end
if size(tx,1)~=size(trace,1) || size(tx,2)~=size(trace,2)
    disp('vector x and y of different size...')
    return
end

% Moyenne
if isempty(is_mean)
    mtrace=0;
else
    if ischar(is_mean) && strcmp(is_mean,'m')==1
        mtrace=mean(trace);
        trace=trace-mean(trace);
    end
    if ~ischar(is_mean)
        if is_mean<min(trace) || is_mean>max(trace)
            %            disp(['is_mean=' num2str(is_mean)]);
            %     disp(['max(trace)=' num2str(max(trace))]);
            %    disp(['min(trace)=' num2str(min(trace))]);
            %	disp('The parameters of filling are out of border of the trace');
            %			return
            mtrace=mean(trace);
            trace=trace-mean(trace);
        else
            mtrace=is_mean;
            trace=trace-is_mean;
        end
    end
end




% on remplace les 0 par des petites valeurs
ii0=find(abs(trace)<max(abs(trace))/1e10);trace(ii0)=max(abs(trace))/1e10;

% selection des valeurs positives (ou negatives)
if is_neg==1
    trace=-trace;
end

%if max(abs(trace))
% recherche des zeros pour une plus belle repr�sentation ....

%descendaning (downward) parts
t1=(trace>0).*trace;
t2=[trace(1) trace];
t2=t2(1:length(trace));
t2=(t2<0).*t2;
i0=find(abs(t1+t2)<1e-15);
if ~isempty(i0) && min(i0)>1
    %disp([min(i0) max(i0) length(t1) length(tx)])
    T=[tx(i0-1);tx(i0)];
    Y=[trace(i0-1);trace(i0)];
    nTd=abs((T(1,:)-T(2,:))).*(Y(1,:)./(Y(1,:)-Y(2,:)))+T(1,:);
else
    nTd=[];
end

%ascendaning (upward) parts
t1=(trace<0).*trace;
t2=[trace(1) trace];
t2=t2(1:length(trace));
t2=(t2>0).*t2;
i0=find(abs(t1+t2)<1e-15);
if ~isempty(i0) && min(i0)>1
    T=[tx(i0-1);tx(i0)];
    Y=[trace(i0-1);trace(i0)];
    nTa=abs((T(1,:)-T(2,:))).*(Y(1,:)./(Y(1,:)-Y(2,:)))+T(1,:);
else
    nTa=[];
end


% rajout des points sur la trace
tx=[tx,nTd,nTa];
trace=[trace zeros(size(nTd)) zeros(size(nTa))];
[tx,it]=sort(tx);
trace2=trace(it);

% on ne retient que les parties positives pour le remplissage
trace3=(trace2>=0).*trace2;
% inversion de la trace si dessin des negatifs

if is_neg==1
    trace3=-trace3;
    trace2=-trace2;
end

% Fermeture du polygone
tx2=[tx(1) tx tx(length(tx))];
trace3=[0 trace3 0];

%hf=fill(tx2,trace3+mtrace,'g');
%set(hf,'EdgeColor',[1 1 1]);
%set(hf,'EdgeColor','none');

% Dessin du polygone rempli
%figure;
max_trace=max(abs(trace));
criter=(max_trace/max_trace_ratio);
criter2=criter*1.5;
%disp([max_trace max_trace_ratio criter criter2]);
Indx=find(abs(trace3)>criter);
if ~isempty(Indx)
    Indx2=Indx(2:end)-Indx(1:end-1);
    indx=find(Indx2>1);
    if ~isempty(indx)
        Indx3=Indx(indx);
        Indx4=Indx3+Indx2(indx);
        
        tr=trace3(Indx(1)-1:Indx3(1)+1);
        tr(1)=0.0;tr(end)=0.0;
        tx1=tx2(Indx(1)-1:Indx3(1)+1);
        %        disp([ max(abs(tr)) criter2]);
        if max(abs(tr))>criter2
            hf=fill(tx1,tr+mtrace,couleur);
            %        set(hf,'EdgeColor',couleur);
            set(hf,'EdgeColor','none');
            hold on;
        end
        for i=2:length(Indx3)
            tr=trace3(Indx4(i-1)-1:Indx3(i)+1);
            tr(1)=0.0;tr(end)=0.0;
            tx1=tx2(Indx4(i-1)-1:Indx3(i)+1);
            if max(abs(tr))>criter2
                hf=fill(tx1,tr+mtrace,couleur);
                set(hf,'EdgeColor','none');
                %        set(hf,'EdgeColor',couleur);
            end
        end
        tr=trace3(Indx4(end)-1:Indx(end)+1);
        tr(1)=0.0;tr(end)=0.0;
        tx1=tx2(Indx4(end)-1:Indx(end)+1);
        if max(abs(tr))>criter2
            hf=fill(tx1,tr+mtrace,couleur);
            set(hf,'EdgeColor','none');
            %         set(hf,'EdgeColor',couleur);
        end
    else
        tr=trace3(Indx(1)-1:Indx(end)+1);
        tr(1)=0.0;tr(end)=0.0;
        tx1=tx2(Indx(1)-1:Indx(end)+1);
        if max(abs(tr))>criter2
            hf=fill(tx1,tr+mtrace,couleur);
            set(hf,'EdgeColor','none');
            %        set(hf,'EdgeColor',couleur);
            hold on;
        end
    end
end

%get(hf)
%pause;
%xlim([0 15]);

%hi=couleur_contour;
% dessin de la courbe
%hold on
%hp=plot(tx,trace2+mtrace);
%set(hp,'Color',couleur_contour)
% hold off
