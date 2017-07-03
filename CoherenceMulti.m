%Version 1.00 
%**************************************************************************
%***************************** CoherenceMulti GUI ******************************
%**************************************************************************
%---------------------------Credits----------------------------------------
%Wavelet Transform: Dmytro Iatsenko
%----------------------------Documentation---------------------------------
%Comnig Soon



function varargout = CoherenceMulti(varargin)
% COHERENCEMULTI MATLAB code for CoherenceMulti.fig
%      COHERENCEMULTI, by itself, creates a new COHERENCEMULTI or raises the existing
%      singleton*.
%
%      H = COHERENCEMULTI returns the handle to a new COHERENCEMULTI or the handle to
%      the existing singleton*.
%
%      COHERENCEMULTI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COHERENCEMULTI.M with the given input arguments.
%
%      COHERENCEMULTI('Property','Value',...) creates a new COHERENCEMULTI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CoherenceMulti_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CoherenceMulti_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help CoherenceMulti

% Last Modified by GUIDE v2.5 03-Jul-2017 17:59:05
%*************************************************************************%
%                BEGIN initialization code - DO NOT EDIT                  %
%                ----------------------------------------                 %
%*************************************************************************%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CoherenceMulti_OpeningFcn, ...
                   'gui_OutputFcn',  @CoherenceMulti_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%*************************************************************************%
%                END initialization code - DO NOT EDIT                    %
%*************************************************************************%


function CoherenceMulti_OpeningFcn(hObject, eventdata, handles, varargin)

movegui('center') 
axes(handles.logo)
matlabImage = imread('physicslogo.png');
image(matlabImage)
axis off
axis image
h = findall(0,'Type','uicontrol');
set(h,'FontUnits','normalized');
handles.calc_type = 1;
guidata(hObject,handles);
drawnow;

handles.output = hObject;
guidata(hObject, handles);
function varargout = CoherenceMulti_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
function plot_type_CreateFcn(hObject, eventdata, handles)
function wavlet_transform_CreateFcn(hObject, eventdata, handles)
function surrogate_type_Callback(hObject, eventdata, handles)
function surrogate_type_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function surrogate_count_Callback(hObject, eventdata, handles)
function surrogate_count_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function surrogate_analysis_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function signal_length_Callback(hObject, eventdata, handles)
function signal_length_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function orientation_Callback(hObject, eventdata, handles)
function orientation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sampling_freq_Callback(hObject, eventdata, handles)
function sampling_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function max_freq_Callback(hObject, eventdata, handles)
function status_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Please Import Signal');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function max_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function min_freq_Callback(hObject, eventdata, handles)
function min_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function wavelet_type_Callback(hObject, eventdata, handles)
function wavelet_type_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function central_freq_Callback(hObject, eventdata, handles)
function central_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function time_series_1_ButtonDownFcn(hObject, eventdata, handles)
function preprocess_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function cutedges_Callback(hObject, eventdata, handles)
function cutedges_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sampling_rate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function length_Callback(hObject, eventdata, handles)
function length_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function intervals_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xlim_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ylim_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function detrend_signal_popup_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function display_type_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function surrogate_percentile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function signal_list_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%--------------------------------------------------Unused Callbacks--------
function status_Callback(hObject, eventdata, handles, msg)
set(handles.status,'String',msg);
drawnow;
%--------------------------------------------------------------------------

function sampling_rate_Callback(hObject, eventdata, handles)
%Replots after changin sampling rate
   % display_selected(hObject, eventdata, handles);

function intervals_Callback(hObject, eventdata, handles)
%Marking lines on the graphs    
    intervals = csv_to_mvar(get(handles.intervals,'String'));    
    
    if(size(intervals)>0)
        zval = 1;
        child_handles = allchild(handles.wt_pane);
        for i = 1:size(child_handles,1)
            
            if(strcmp(get(child_handles(i),'Type'),'axes'))
                set(child_handles(i),'Ytick',intervals);
                hold(child_handles(i),'on');
                warning('off');
                
                for j = 1:size(intervals,2)
                    xl = get(child_handles(i),'xlim');
                    x = [xl(1) xl(2)];        
                    z = ones(1,size(x,2));
                    z = z.*zval;
                    y = intervals(j)*ones(1,size(x,2));
                    plot3(child_handles(i),x,y,z,'--k');
                end
                
                warning('on');
                hold(child_handles(i),'off');
            end
            
        end
    end
    
    
% function preprocess_Callback(hObject, eventdata, handles)
% %Detrending Part Visualisation
%     cla(handles.plot_pp,'reset');
%     data = guidata(hObject);
%     sig = data.sig; 
%     time_axis = data.time_axis;
%     L = size(sig,2);
%     fs = str2double(get(handles.sampling_freq,'String'));
%     fmax = str2double(get(handles.max_freq,'String'));
%     fmin = str2double(get(handles.min_freq,'String'));
%     
%     contents = cellstr(get(handles.detrend_signal_popup,'String'));
%     i = contents{get(handles.detrend_signal_popup,'Value')};
%     i = str2double(i);
%     
%     %Detrending
%     cur_sig = sig(i,:);
%     cur_sig = cur_sig(:);
%     X=(1:length(cur_sig))'/fs; XM=ones(length(X),4); 
% 
%     for pn=1:3 
%         CX=X.^pn; 
%         XM(:,pn+1)=(CX-mean(CX))/std(CX); 
%     end
% 
%     w=warning('off','all'); 
%     new_signal=cur_sig-XM*(pinv(XM)*cur_sig); 
%     warning(w);
% 
%     %Filtering
%     fx=fft(new_signal,L); % Fourier transform of a signal
% 
%     Nq=ceil((L+1)/2); 
%     ff=[(0:Nq-1),-fliplr(1:L-Nq)]*fs/L; 
%     ff=ff(:); % frequencies in Fourier transform
% 
%     fx(abs(ff)<=max([fmin,fs/L]) | abs(ff)>=fmax)=0; % filter signal in a chosen frequency domain
%     new_signal=ifft(fx);
%     %Plotting
%     
%     plot(handles.plot_pp,time_axis,cur_sig);
%     hold(handles.plot_pp,'on');
%     plot(handles.plot_pp,time_axis,new_signal,'-r');
%     
%     
%     %legend(handles.plot_pp,'Original','Pre-Processed','Location','Best');
%     xlim(handles.plot_pp,[0,size(sig,2)./fs]);
    
%-------------------------------------------------------------------------    

function wavlet_transform_Callback(hObject, eventdata, handles)
%Does the wavelet transform 

% Get user input from GUI
    status_Callback(hObject, eventdata, handles, 'Calculating Wavelet Transform...');
    fs = str2double(get(handles.sampling_freq,'String'));
    fmax = str2double(get(handles.max_freq,'String'));
    fmin = str2double(get(handles.min_freq,'String'));
    fc =  str2double(get(handles.central_freq,'String'));
    surrogate_count = floor(str2double(get(handles.surrogate_count,'String')));
    
    if (surrogate_count == 1)
        errordlg('Number of surrogates must be greater than 1','Parameter Error');
        set(handles.status,'String','Enter Valid Parameters before continuing');
        drawnow;
        return;
    end
    
    if isnan(fs)
      errordlg('Sampling frequency must be specified','Parameter Error');
      set(handles.status,'String','Enter Valid Parameters before continuing');
      drawnow;
      return;
    end
    
    if ~isfield(handles,'sig')
      errordlg('Signal not found','Signal Error');
      set(handles.status,'String','Enter Valid Parameters before continuing');
      drawnow;
      return;
    end
    
    items = get(handles.wavelet_type,'String');
    index_selected = get(handles.wavelet_type,'Value');
    wavelet_type_selected = items{index_selected};
    
    items = get(handles.preprocess,'String');
    index_selected = get(handles.preprocess,'Value');
    preprocess_selected = items{index_selected};
    
    items = get(handles.cutedges,'String');
    index_selected = get(handles.cutedges,'Value');
    cutedges_selected = items{index_selected};
    
    
    sig = handles.sig;    % do not optimise as the signal is being cut in this case
    
    xl = csv_to_mvar(get(handles.xlim,'String'));
    xl = xl.*fs;
    xl(2) = min(xl(2),size(sig,2));
    xl(1) = max(xl(1),1);
    xl = xl./fs;
    time_axis = xl(1):1/fs:xl(2);
    
    if length(time_axis)>=2000
        screensize = max(get(groot,'Screensize'));
        under_sample = floor(size(sig,2)/screensize);%TODO improve reliability with screens
    else 
        under_sample = 1;
    end

    handles.time_axis_us = time_axis(1:under_sample:end);
    n = size(handles.sig,1)/2 ;
    
%Taking only selected part of the signal
    xl = get(handles.xlim,'String');
    xl = csv_to_mvar(xl);
    xl = xl.*fs;
    xl(2) = min(xl(2),size(handles.sig,2));
    xl(1) = max(xl(1),1);
    sig = sig(:,xl(1):xl(2));
    
    set(handles.status,'String','Calculating Wavelet Transform...');
    
    %Calculating wavelet transform
    for p = 1:n
        status_Callback(hObject, eventdata, handles, sprintf('Calculating Wavelet Tranform of Signal %d/%d',p,n));
        if(isnan(fmax)&& isnan(fmin))
            if(isnan(fc))                              
                    [wt_1,handles.freqarr]=wt(sig(p,:),fs,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 

                    [wt_2,handles.freqarr]=wt(sig(p+n,:),fs,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
            else
                    [wt_1,handles.freqarr]=wt(sig(p,:),fs,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 

                    [wt_2,handles.freqarr]=wt(sig(p+n,:),fs,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);      
            end
        elseif(isnan(fmax))
            if(isnan(fc))
                    [wt_1,handles.freqarr]=wt(sig(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 

                    [wt_2,handles.freqarr]=wt(sig(p+n,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
            else
                    [wt_1,handles.freqarr]=wt(sig(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 

                    [wt_2,handles.freqarr]=wt(sig(p+n,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
            end
        elseif(isnan(fmin))
            if(isnan(fc))
                    [wt_1,handles.freqarr]=wt(sig(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 

                    [wt_2,handles.freqarr]=wt(sig(p+n,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
            else
                    [wt_1,handles.freqarr]=wt(sig(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 

                    [wt_2,handles.freqarr]=wt(sig(p+n,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
            end
        else
            if(isnan(fc))
                    [wt_1,handles.freqarr]=wt(sig(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 

                    [wt_2,handles.freqarr]=wt(sig(p+n,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
            else
                    [wt_1,handles.freqarr]=wt(sig(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 

                    [wt_2,handles.freqarr]=wt(sig(p+n,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                        'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
            end
        end
        
        handles.TPC{p,1} = tlphcoh(wt_1,wt_2,handles.freqarr,fs);
        handles.time_avg_wpc{p,1} = nanmean(handles.TPC{p,1}.');
        handles.TPC{p,1} = handles.TPC{p,1}(:,1:under_sample:end);
    end   
    
    %---------------
    
    %Surrogate Calculation
     
    items = get(handles.surrogate_type,'String');
    index_selected = get(handles.surrogate_type,'Value');
    surrogate_type = items{index_selected};
    if (surrogate_count > 1) && (floor(surrogate_count) == surrogate_count)
        
        handles.surrogates = cell(size(handles.sig,1),1);
        for p = 1:size(handles.sig,1)/2
            handles.surrogates{p,1} = surrogate(sig(p,:),surrogate_count,surrogate_type);        
        end    
        handles.TPC_surr_avg_arr = cell(surrogate_count,size(handles.sig,1)/2);

        for idx = 1:size(handles.sig,1)/2
            status_Callback(hObject, eventdata, handles, sprintf('Calculating Wavelet Tranform of Signal %d/%d',idx,size(handles.sig,1)/2));
            if(isnan(fmax)&& isnan(fmin))
                if(isnan(fc))                              
                        [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
                else
                        [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);     
                end
            elseif(isnan(fmax))
                if(isnan(fc))
                        [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
                else
                        [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
                end
            elseif(isnan(fmin))
                if(isnan(fc))
                        [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
                else
                        [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
                end
            else
                if(isnan(fc))
                        [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
                else
                        [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);
                end
            end
            for p = 1:surrogate_count
                status_Callback(hObject, eventdata, handles, ...
                    sprintf('Calculating Wavelet Transform for surrogate:%d/%d for signal pair %d',p,surrogate_count,idx));

                if(isnan(fmax)&& isnan(fmin))
                    if(isnan(fc))                  
                        [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);     
                    else
                        [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);    
                    end
                elseif(isnan(fmax))
                    if(isnan(fc))
                        [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);  
                    else
                        [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);  
                    end
                elseif(isnan(fmin))
                    if(isnan(fc))
                        [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);    
                    else
                        [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);  
                    end
                else
                    if(isnan(fc))
                        [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);
                    else
                        [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
                            'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);  
                    end
                end

                TPC_surrogate = tlphcoh(wt_1,WT_surrogate,handles.freqarr,fs);
                handles.TPC_surr_avg_arr{p,idx} = nanmean(TPC_surrogate.'); 
            end
        end

        status_Callback(hObject, eventdata, handles, 'Finished calculating surrogates');
    %-------------------------------------------------------------------------

        %handles.TPC_surr_avg_arr = cell2mat(TPC_surr_avg_arr);

        surrogate_analysis = get(handles.surrogate_analysis,'Value');        
        handles.TPC_surr_avg_max = cell(size(handles.sig,1)/2,1);
        for i = 1:size(handles.sig,1)/2
            if(surrogate_analysis == 2)
                t = cell2mat(handles.TPC_surr_avg_arr);
                t = t(:,length(handles.freqarr)*(i-1)+1:length(handles.freqarr)*(i));        
                surrogate_percentile = str2double(get(handles.surrogate_percentile,'String'));
                handles.TPC_surr_avg_max{i,1} = prctile(t,surrogate_percentile);

            elseif(surrogate_analysis == 1)    
                t = cell2mat(handles.TPC_surr_avg_arr);
                t = t(:,length(handles.freqarr)*(i-1)+1:length(handles.freqarr)*(i));
                handles.TPC_surr_avg_max{i,1} = max(t);

            end
        end     
    end
    guidata(hObject,handles);
    xyplot_Callback(hObject, eventdata, handles);
    intervals_Callback(hObject, eventdata, handles)
    guidata(hObject,handles);
    set(handles.intervals,'Enable','on');
    drawnow;

function xyplot_Callback(hObject, eventdata, handles)
%Plots all figures
    signal_selected = get(handles.signal_list,'Value');  
    if any(signal_selected == size(handles.sig,1)/2+1) && isfield(handles,'freqarr')     
        cla(handles.plot3d,'reset');
        cla(handles.plot_pow,'reset');
        cla(handles.cum_avg,'reset');
        set(handles.plot3d,'visible','off');
        set(handles.plot_pow,'visible','off');   
        set(handles.cum_avg,'visible','on');
        hold(handles.cum_avg,'on');
        if size(handles.sig,1)/2 > 1
            plot(handles.cum_avg, handles.freqarr, mean(cell2mat(handles.time_avg_wpc)),'--','Linewidth',3);
            plot(handles.cum_avg, handles.freqarr, median(cell2mat(handles.time_avg_wpc)),'-','Linewidth',3);
        else
            plot(handles.cum_avg, handles.freqarr, cell2mat(handles.time_avg_wpc),'--','Linewidth',3);
            plot(handles.cum_avg, handles.freqarr, cell2mat(handles.time_avg_wpc),'-','Linewidth',3);
        end
        ylabel(handles.cum_avg,'Average Coherence');
        xlabel(handles.cum_avg,'Frequency (Hz)');

        legend(handles.cum_avg,'mean','median')
        for i = 1:size(signal_selected,2)            
            if(signal_selected(i) <= size(handles.sig,1)/2)                                
                plot(handles.cum_avg, handles.freqarr, handles.time_avg_wpc{signal_selected(i),1});                     
                ylabel(handles.cum_avg,'Average Coherence');
                xlabel(handles.cum_avg,'Frequency (Hz)');
                [M,I] = max(handles.time_avg_wpc{signal_selected(i),1});
                text(handles.cum_avg,handles.freqarr(I),M,num2str(signal_selected(i)));  
            end
        end
        
        set(handles.cum_avg,'xscale','log');     
        idx_first = find(sum(~isnan(handles.time_avg_wpc{1,1}),1) > 0, 1 ,'first');
        idx_last = find(sum(~isnan(handles.time_avg_wpc{1,1}),1) > 0, 1 , 'last');   
        xlim(handles.cum_avg,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
        
    elseif isfield(handles,'freqarr') 
        cla(handles.cum_avg,'reset');
        cla(handles.plot3d,'reset');
        cla(handles.plot_pow,'reset');
        set(handles.cum_avg,'visible','off');
        set(handles.plot3d,'visible','on');
        set(handles.plot_pow,'visible','on');
        
        handles.peak_value = max(handles.TPC{signal_selected,1}(:))+.1;
        pcolor(handles.plot3d, handles.time_axis_us , handles.freqarr, handles.TPC{signal_selected,1});                
        plot(handles.plot_pow, handles.time_avg_wpc{signal_selected,1}, handles.freqarr,'LineWidth',2);     
        hold(handles.plot_pow,'on');
        if isfield(handles,'TPC_surr_avg_max')
            plot(handles.plot_pow ,handles.TPC_surr_avg_max{signal_selected,1} , handles.freqarr,'LineWidth',2);
            plot(handles.plot_pow,handles.TPC_surr_avg_max{signal_selected,1} , handles.freqarr,'LineWidth',2);
            color_positive_breach(handles.plot_pow,handles.freqarr,handles.time_avg_wpc{signal_selected,1}, handles.TPC_surr_avg_max{signal_selected,1},'color','red','flipped');
        end
        hold(handles.plot_pow,'off');
        xlabel(handles.plot_pow,'Average Coherence');       
        
        c = colorbar(handles.plot3d,'Location','east');
        set(c, 'position',[0.73 .12 .015 .85],'Linewidth',0.2);
        set(c, 'fontsize',8,'units','normalized');
        shading(handles.plot3d,'interp');       
        set(handles.plot3d,'yscale','log');
        set(handles.plot_pow,'yscale','log');        
        set(handles.plot_pow,'yticklabel',[]);
        set(handles.plot3d,'ylim',[min(handles.freqarr) max(handles.freqarr)]);%making the axes tight
        set(handles.plot3d,'xlim',[handles.time_axis_us(1) handles.time_axis_us(end)]);%making the axes tight
        xlabel(handles.plot3d,'Time (s)');
        ylabel(handles.plot3d,'Frequency (Hz)');    
        ylabel(handles.plot_pow,'Frequency (Hz)');    
        ylim(handles.plot_pow,[min(handles.freqarr) max(handles.freqarr)]);
        set(handles.status,'String','Done Plotting');
    end
    
    set(handles.plot3d,'Fontunits','normalized');
    set(handles.plot_pow,'Fontunits','normalized');
    set(handles.cum_avg,'Fontunits','normalized');
    guidata(hObject,handles);

%---------------------------Surrogate Analysis-----------------
function surrogate_analysis_Callback(hObject, eventdata, handles)
%To enable and disable the Percentile box

surrogate_analysis = get(handles.surrogate_analysis,'Value');  
if(surrogate_analysis == 2)
    set(handles.surrogate_percentile,'Enable','on');
    surrogate_percentile = str2double(get(handles.surrogate_percentile,'String'));
    handles.TPC_surr_avg_max = prctile(handles.TPC_surr_avg_arr,surrogate_percentile);
    display('percentile');
elseif(surrogate_analysis == 1)
    set(handles.surrogate_percentile,'Enable','off');
    handles.TPC_surr_avg_max = max(handles.TPC_surr_avg_arr);
end

guidata(hObject,handles);
subtract_surrogates_Callback(hObject, eventdata, handles)
guidata(hObject,handles);

function surrogate_percentile_Callback(hObject, eventdata, handles)
surrogate_analysis = get(handles.surrogate_analysis,'Value');  
if(surrogate_analysis == 2)
    set(handles.surrogate_percentile,'Enable','on');
    surrogate_percentile = str2double(get(handles.surrogate_percentile,'String'));
    handles.TPC_surr_avg_max = prctile(handles.TPC_surr_avg_arr,surrogate_percentile);
    display('percentile');
elseif(surrogate_analysis == 1)
    set(handles.surrogate_percentile,'Enable','off');
    handles.TPC_surr_avg_max = max(handles.TPC_surr_avg_arr);
end

guidata(hObject,handles);
subtract_surrogates_Callback(hObject, eventdata, handles)
guidata(hObject,handles);

function subtract_surrogates_Callback(hObject, eventdata, handles)
display_selected = get(handles.display_type,'Value');
if display_selected == 3 
    toggle = get(handles.subtract_surrogates,'Value');
    if toggle == get(handles.subtract_surrogates,'Max')
        cla(handles.plot_pow,'reset');
        corrected_coherence = handles.time_avg_wpc - handles.TPC_surr_avg_max;
        corrected_coherence = subplus(corrected_coherence);
        plot(handles.plot_pow ,corrected_coherence, handles.freqarr,'LineWidth',2);
        set(handles.plot_pow,'yscale','log','yticklabel',[]);     
        ylim(handles.plot_pow,[min(handles.freqarr) max(handles.freqarr)]);
        xlabel(handles.plot_pow,{'Average Coherence','(Surrogate Subtracted)'});
        %legend(handles.plot_pow,'Surrogate Subtracted','Location','Best');
    else
        cla(handles.plot_pow,'reset');
        hold(handles.plot_pow,'on');
        plot(handles.plot_pow ,handles.time_avg_wpc, handles.freqarr,'LineWidth',2);
        if(size(handles.TPC_surr_avg_max)>0)
            plot(handles.plot_pow ,handles.TPC_surr_avg_max , handles.freqarr,'LineWidth',2);
        end
        hold(handles.plot_pow,'off');     
        set(handles.plot_pow,'yscale','log','yticklabel',[]);     
        ylim(handles.plot_pow,[min(handles.freqarr) max(handles.freqarr)]);
        color_positive_breach(handles.plot_pow,handles.freqarr,handles.time_avg_wpc, handles.TPC_surr_avg_max,'color','red','flipped');
        set(handles.status,'String','Done Plotting');
        xlabel(handles.plot_pow,'Average Coherence');   
        %legend(handles.plot_pow,'Original Signal','Surrogate','Location','Best');
    end
    
end    
%---------------------------------------Surrogate Analysis------

% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
%Loading data

% --------------------------------------------------------------------
function csv_read_Callback(hObject, eventdata, handles)
%Read csv file
    linkaxes([handles.time_series_1 handles.time_series_2],'x');
    set(handles.status,'String','Importing Signal...');
    fs = str2double(get(handles.sampling_freq,'String'));     
    if isnan(fs)
      errordlg('Sampling frequency must be specified before importing','Parameter Error');
      return;
    end
    
    sig = read_from_csv();
    
    % Construct a questdlg with three options
    choice = questdlg('Select Orientation of Data set?', ...
        'Data Import','Column wise','Row wise','default');
    switch choice
        case 'Column wise'
            sig = sig';
            %msgbox(sprintf('Number of data sets %d',size(sig,1)),'Import Complete');
        case 'Row wise'
            %msgbox(sprintf('Number of data sets %d',size(sig,1)),'Import Complete');
        case 'default'
            errordlg('Data set orientation must be specified')
            return;
    end
    list = cell(size(sig,1)/2+1,1);
    list{1,1} = 'Signal Pair 1';
    
    for i = 2:size(sig,1)/2
        list{i,1} = sprintf('Signal Pair %d',i);
    end
    set(handles.signal_list,'String',list);
    list{i+1,1} = sprintf('Average Plot(All)');
    set(handles.signal_list,'String',list); 
    
    handles.sig = sig;   
    time = 1:size(sig,2);
    time = time./fs;
    handles.time_axis = time;
    
    plot(handles.time_series_1,time,sig(1,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_1,[0,size(sig,2)./fs]);
    plot(handles.time_series_2,time,sig(1+size(sig,1)/2,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_2,[0,size(sig,2)./fs]);
    
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    guidata(hObject,handles);  
    %cla(handles.plot_pp,'reset');
    %preprocess_Callback(hObject, eventdata, handles);%plots the detrended curve
    xlabel(handles.time_series_2,'Time (s)');
    set(handles.status,'String','Select Data And Continue With Wavelet Transform');
    set(handles.signal_length,'String',strcat(num2str(size(sig,2)/fs/60),' minutes'));

% --------------------------------------------------------------------
function mat_read_Callback(hObject, eventdata, handles)
%Read mat file    
    linkaxes([handles.time_series_1 handles.time_series_2],'x');
    set(handles.status,'String','Importing Signal...');
    fs = str2double(get(handles.sampling_freq,'String'));   
    if isnan(fs)
      errordlg('Sampling frequency must be specified before importing','Parameter Error');
      return;
    end
    
    sig = read_from_mat();
    sig = struct2cell(sig);
    sig = cell2mat(sig);
    
    % Construct a questdlg with three options
    choice = questdlg('Select Orientation of Data set?', ...
        'Data Import','Column wise','Row wise','default');
    switch choice
        case 'Column wise'
            sig = sig';
            %msgbox(sprintf('Number of data sets %d',size(sig,1)),'Import Complete');
        case 'Row wise'
            %msgbox(sprintf('Number of data sets %d',size(sig,1)),'Import Complete');
        case 'default'
            errordlg('Data set orientation must be specified')
            return;
    end
    list = cell(size(sig,1)/2+1,1);
    list{1,1} = 'Signal Pair 1';
    i = 2;
    for i = 2:size(sig,1)/2
        list{i,1} = sprintf('Signal Pair %d',i);
    end
    set(handles.signal_list,'String',list);
    list{size(sig,1)/2+1,1} = sprintf('Average Plot(All)');
    set(handles.signal_list,'String',list); 
    
    handles.sig = sig;   
    time = 1:size(sig,2);
    time = time./fs;
    handles.time_axis = time;
    
    plot(handles.time_series_1,time,sig(1,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_1,[0,size(sig,2)./fs]);
    plot(handles.time_series_2,time,sig(1+size(sig,1)/2,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_2,[0,size(sig,2)./fs]);
    
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    guidata(hObject,handles);  
    %cla(handles.plot_pp,'reset');
    %preprocess_Callback(hObject, eventdata, handles);%plots the detrended curve
    xlabel(handles.time_series_2,'Time (s)');
    set(handles.status,'String','Select Data And Continue With Wavelet Transform');
    set(handles.signal_length,'String',strcat(num2str(size(sig,2)/fs/60),' minutes'));

    
%---------------------------Limits-----------------------------
function xlim_Callback(hObject, eventdata, handles)
%When the values of xlim are changed the graphs are updated
    xl = csv_to_mvar(get(handles.xlim,'String'));
    xlim(handles.time_series_1,xl);
    xlim(handles.time_series_2,xl);
    xlim(handles.plot_pp,xl);
    t = xl(2) - xl(1);
    set(handles.length,'String',t);

function ylim_Callback(hObject, eventdata, handles)
%When the values of ylim are changed the graphs are updated  
    yl = csv_to_mvar(get(handles.ylim,'String'));
    ylim(handles.time_series_1,yl);
    ylim(handles.time_series_2,yl);


%---------------------------Updating Value of limits Limits-----------------------------
function refresh_limits_Callback(hObject, eventdata, handles)
%Calcualtes limits of the plot    
    
    x = get(handles.time_series_1,'xlim');
%    xlim(handles.plot_pp,x);
    t = x(2) - x(1);
    x = strcat(num2str(x(1)),' , ',num2str(x(2)));    
    
    y = get(handles.time_series_1,'ylim');
    y = strcat(num2str(y(1)),' , ',num2str(y(2)));
    
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);
    
% ---------------------------Zoom Updating--------------------------
function zoom_in_OffCallback(hObject, eventdata, handles)
%Refreshes the limit values right after the tool is deselected
    x = get(handles.time_series_1,'xlim');
%    xlim(handles.plot_pp,x);
    t = x(2) - x(1);
    x = strcat(num2str(x(1)),' , ',num2str(x(2)));    
    
    y = get(handles.time_series_1,'ylim');
    y = strcat(num2str(y(1)),' , ',num2str(y(2)));
    
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);

% -----------------------------Zoom Updating--------------------------
function zoom_out_OffCallback(hObject, eventdata, handles)
%Refreshes the limit values right after the tool is deselected
    x = get(handles.time_series_1,'xlim');
%    xlim(handles.plot_pp,x);
    t = x(2) - x(1);
    x = strcat(num2str(x(1)),' , ',num2str(x(2)));    
    
    y = get(handles.time_series_1,'ylim');
    y = strcat(num2str(y(1)),' , ',num2str(y(2)));
    
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);
    
function plot_type_SelectionChangeFcn(hObject, eventdata, handles)
%deciding which plot
    switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
        case 'power'
            plot_type = 1;
        case 'amp'
            plot_type = 2;
    end

    data = guidata(hObject);
    data.plot_type = plot_type;
    guidata(hObject,data); 

% ----------------------------------------Saving Files---------------
function save_Callback(hObject, eventdata, handles)
%Honestly you're just here because I don't know how to get rid of you

function save_3dplot_Callback(hObject, eventdata, handles)
%Saves the 3d plot
    Fig = figure;
    copyobj(handles.plot3d, Fig);
    Fig = tightfig(Fig);

function save_power_plot_Callback(hObject, eventdata, handles)
%Saves the power plot
    Fig = figure;
    copyobj(handles.plot_pow, Fig);
    Fig = tightfig(Fig);

function save_pow_arr_Callback(hObject, eventdata, handles)
%Saves the avg power array
    [FileName,PathName] = uiputfile
    save_location = strcat(PathName,FileName)
    data = guidata(hObject);
    pow_arr = data.pow_arr;
    save(save_location,'pow_arr');

function save_wt_Callback(hObject, eventdata, handles)
%Saves the wavelet transform
    [FileName,PathName] = uiputfile
    save_location = strcat(PathName,FileName)
    data = guidata(hObject);
    WT = data.WT;
    freqarr = data.freqarr;
    save(save_location,'freqarr','-v7.3');
    save(save_location,'WT','-v7.3');%Sometimes the compression is faulty

function detrend_signal_popup_Callback(hObject, eventdata, handles)
%Detrends the signal plots the chosen one
    cla(handles.plot_pp,'reset');
    %preprocess_Callback(hObject, eventdata, handles);

function signal_list_Callback(hObject, eventdata, handles)
%Selecting signal and calling other necessary functions
    signal_selected = get(handles.signal_list, 'Value');
    
    if any(signal_selected == size(handles.sig,1)/2+1)
        set(handles.signal_list,'Max',size(handles.sig,1)/2);
    else
        if size(signal_selected,2) == 1
            set(handles.signal_list,'Max',1);
        else
            set(handles.signal_list, 'Value', 1);
            set(handles.signal_list,'Max',1);
            drawnow;
            xyplot_Callback(hObject, eventdata, handles);
        end
    end
    
    if any(signal_selected ~= size(handles.sig,1)/2+1) && length(signal_selected) == 1
        
        plot(handles.time_series_1, handles.time_axis, handles.sig(signal_selected,:));%Plotting the time_series part after calculation of appropriate limits
        xl = csv_to_mvar(get(handles.xlim, 'String'));
        xlim(handles.time_series_1, xl);
        plot(handles.time_series_2, handles.time_axis, handles.sig(signal_selected+size(handles.sig,1)/2,:));%Plotting the time_series part after calculation of appropriate limits
        xlim(handles.time_series_2, xl);        
        xlabel(handles.time_series_2, 'Time (s)');
        refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
        %cla(handles.plot_pp, 'reset');
        %preprocess_Callback(hObject, eventdata, handles);%plots the detrended curve
        %xlabel(handles.time_series, 'Time (s)');
        set(handles.status, 'String', 'Select Data And Continue With Wavelet Transform');
        if isfield(handles,'TPC')
            xyplot_Callback(hObject, eventdata, handles);
        end
        intervals_Callback(hObject, eventdata, handles)
        
    elseif any(signal_selected == size(handles.sig,1)/2+1)
        xyplot_Callback(hObject, eventdata, handles);
        intervals_Callback(hObject, eventdata, handles)
    end
