function varargout = PCA_example_GUI_2D(varargin)
% PCA_EXAMPLE_GUI_2D M-file for PCA_example_GUI_2D.fig
%      PCA_EXAMPLE_GUI_2D, by itself, creates a new PCA_EXAMPLE_GUI_2D or raises the existing
%      singleton*.
%
%      H = PCA_EXAMPLE_GUI_2D returns the handle to a new PCA_EXAMPLE_GUI_2D or the handle to
%      the existing singleton*.
%
%      PCA_EXAMPLE_GUI_2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCA_EXAMPLE_GUI_2D.M with the given input arguments.
%
%      PCA_EXAMPLE_GUI_2D('Property','Value',...) creates a new PCA_EXAMPLE_GUI_2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PCA_example_GUI_2D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PCA_example_GUI_2D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PCA_example_GUI_2D

% Last Modified by GUIDE v2.5 29-Oct-2009 08:49:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PCA_example_GUI_2D_OpeningFcn, ...
                   'gui_OutputFcn',  @PCA_example_GUI_2D_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before PCA_example_GUI_2D is made visible.
function PCA_example_GUI_2D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PCA_example_GUI_2D (see VARARGIN)

% Choose default command line output for PCA_example_GUI_2D
handles.output = hObject;

set(hObject,'toolbar','figure');
% Compact layout: inputs across the top, NMF and PCA side by side below
compact_layout_2d(hObject,handles);
% Show the default synergies and data when the GUI opens
PCA_example_GUI_2D('pushbutton2_Callback',handles.pushbutton2,[],handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PCA_example_GUI_2D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PCA_example_GUI_2D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%***********************************************
% Various initializations & declarations       *
%***********************************************

% Declare Global objects and variables
global W C C1_wt C2_wt nmus DATA


%***********************************************
% Defining and Plotting the muscle synergies   *
%***********************************************

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global W nmus

% DEFINE THE SYNERGIES
% Read them from the W boxes next to the bar plots
% (chapter example: W=[0.5 1;1 0.5], i.e. M2=2*M1 and M2=0.5*M1)
W=read_W(handles);

[nmus nsyn]=size(W);

%PLOT THE SYNERGIES
axes(handles.axes1)
b= bar(W(:,1));set(b,'FaceColor','b','EdgeColor',[1 1 1]);
title('W')
guidata(hObject, handles); %updates the handles
axes(handles.axes2)
b= bar(W(:,2));set(b,'FaceColor','b','EdgeColor',[1 1 1]);
clear b
sync_W_axes(handles); % same y-scale on all W bar plots
guidata(hObject, handles); %updates the handles



%***********************************************
% Defining the synergy coefficients            *
%***********************************************
replot_data(handles); % update the data plot right away


function NumPoints_Callback(hObject, eventdata, handles)
% hObject    handle to NumPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumPoints as text
%        str2double(get(hObject,'String')) returns contents of NumPoints as a double

numpts=str2double(get(hObject,'String')); %returns contents of NumPoints as a double
global C Cb
C=rand(2,numpts); %random numbers for C (data set 1)
Cb=rand(2,numpts); %random numbers for data set 2
replot_data(handles); % update the data plot right away


% --- Executes during object creation, after setting all properties.
function NumPoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global C1_wt

%obtains the slider value from the slider component
sliderValue = get(handles.slider1,'Value');
%puts the slider value into the edit text component
set(handles.slider1_editText,'String', num2str(sliderValue));
% Update handles structure
guidata(hObject, handles);
C1_wt=sliderValue;
replot_data(handles); % update the data plot right away


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global C2_wt

%obtains the slider value from the slider component
sliderValue2 = get(handles.slider2,'Value');
%puts the slider value into the edit text component
set(handles.slider2_editText,'String', num2str(sliderValue2));
% Update handles structure
guidata(hObject, handles);

C2_wt=sliderValue2;
replot_data(handles); % update the data plot right away


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function slider1_editText_Callback(hObject, eventdata, handles)
% hObject    handle to slider1_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slider1_editText as text
%        str2double(get(hObject,'String')) returns contents of slider1_editText as a double

global C1_wt

%get the string for the editText component
sliderValue = get(handles.slider1_editText,'String');
%convert from string to number if possible, otherwise returns empty
sliderValue = str2num(sliderValue);
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(sliderValue) || sliderValue < 0 || sliderValue > 100)
    set(handles.slider1,'Value',0);
    set(handles.slider1_editText,'String','0');
else
    set(handles.slider1,'Value',sliderValue);
end
C1_wt=sliderValue;
replot_data(handles); % update the data plot right away


% --- Executes during object creation, after setting all properties.
function slider1_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slider2_editText_Callback(hObject, eventdata, handles)
% hObject    handle to slider2_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slider2_editText as text
%        str2double(get(hObject,'String')) returns contents of slider2_editText as a double

global C2_wt

%get the string for the editText component
sliderValue2 = get(handles.slider2_editText,'String');
%convert from string to number if possible, otherwise returns empty
sliderValue2 = str2num(sliderValue2);
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(sliderValue2) || sliderValue2 < 0 || sliderValue2 > 100)
    set(handles.slider2,'Value',0);
    set(handles.slider2_editText,'String','0');
else
    set(handles.slider2,'Value',sliderValue2);
end
C2_wt=sliderValue2;
replot_data(handles); % update the data plot right away


% --- Executes during object creation, after setting all properties.
function slider2_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%***********************************************
% Plotting the made-up dataset                 *
%***********************************************

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global W C Cb C1_wt C2_wt DATA nmus DATA1

% Synergies from the W boxes, weights from the sliders
W=read_W(handles); nmus=size(W,1);
C1_wt=get(handles.slider1,'Value');
C2_wt=get(handles.slider2,'Value');
% Random samples: C for data set 1, Cb for data set 2 (new ones if the number of points changed)
n=str2double(get(handles.NumPoints,'String'));
if isnan(n) || n<2, n=100; set(handles.NumPoints,'String','100'); end
n=round(n);
if isempty(C) || size(C,2)~=n, C=rand(2,n); end
if isempty(Cb) || size(Cb,2)~=n, Cb=rand(2,n); end

% Data set 1
DATA1=W*[C(1,:)*C1_wt; C(2,:)*C2_wt];
DATA=DATA1;
ax=handles.axes3;
cla(ax); hold(ax,'off')
plot(ax,DATA1(1,:),DATA1(2,:),'.','Color',[0 0.45 0.74]);
% Data set 2 (if ticked): its own random sample and its own C1/C2 weights
if get(handles.checkbox3,'Value')==1
    s2=findall(handles.figure1,'Tag','set2_C1'); s3=findall(handles.figure1,'Tag','set2_C2');
    DATA2=W*[Cb(1,:)*get(s2,'Value'); Cb(2,:)*get(s3,'Value')];
    DATA=[DATA1 DATA2]; %combine both data sets for synergy extraction
    hold(ax,'on')
    plot(ax,DATA2(1,:),DATA2(2,:),'.','Color',[0.85 0.33 0.1]);
    hold(ax,'off')
end
set(ax,'XLim',[0 1],'YLim',[0 1],'PlotBoxAspectRatio',[1 1 1]); % fixed, square axes for the data plot
guidata(hObject, handles); %updates the handles



%***********************************************
% Extract and Plot NNMF synergies              *
%***********************************************

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DATA

%Extract muscle synergies using NNMF:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
V=DATA;
r=2;

V = V.*(V>0); % Any potential negative entrie in data matrix will be set to zero

test=sum(V,2); % Any potential muscle channnel with only zeros is not included in the iteration 
index=find(test~=0);
ind=find(test==0);
Vnew_m=V(index,:);

test_cond=sum(V,1); % Any potential condition with only zeros is not included in the iteration 
index_cond=find(test_cond~=0);
ind_cond=find(test_cond==0);
Vnew=Vnew_m(:,index_cond);

% Scale the input data to have unit variance %%%%%%%%%
stdev = std(Vnew'); %scale the data to have unit variance of this data set
Vnew = diag(1./stdev)*Vnew;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions. 
[n,m]=size(Vnew);
H=rand(r,m);
W=rand(n,r);

% Finds zeros in the initial matrix and removes them:
h0=find(H==0);
w0=find(W==0);
while isempty(h0)==0
    H=rand(r,m);
    h0=find(H==0);
end
while isempty(w0)==0
    W=rand(n,r);
    w0=find(W==0);
end

% Sum of squared errors between original and reconstruction:
err=sum(sum( (Vnew-W*H).*(Vnew-W*H) )); 

MAX_IT=100000; 

% Error goal - the "err" quantity, as defined, is the squared error.  If we
% want a 1% mse, then, we want .01*prod(size(V))=.01*n*m.
% Update...  For normed data, the max err is n x m
ERR_GOAL=.0001*(n*m); 

% Update...  For normed data, the max err is n x m
err_save=[]; 

while err>ERR_GOAL
    
    H_fac=W'*Vnew;
    
    H_fac=H_fac./(W'*W*H);
    H=H.*H_fac;
    
    W_fac=Vnew*H';
    W_fac=W_fac./(W*H*H');
    
    W=W.*W_fac;

    err=sum(sum( (Vnew-W*H).*(Vnew-W*H) ));
    err_save=[err_save;err];

    if length(err_save)>550
        err_change = (err_save(end-500) - err_save(end))/err_save(end-500);
        if abs(err_change)<.01
            disp('NMF.m: Error has changed less than 1% over the last 500 iterations. Exiting.')
            break;
        end
    end
    
    if length(err_save)>MAX_IT
        disp('NMF.m: Maximum iterations exceeded. Exiting.')
        break;
    end
end

% Re-scale the original data and the synergies; add in zero rows; calculate 
% final error.

%undo the unit variance scaling so synergies are back out of unit variance
%space and in the same scaling as the input data was
Vnew = diag(stdev)*Vnew;
W = diag(stdev)*W;


% Synergy vectors normalization  %

m=max(W);% vector with max activation values 
for i=1:r
    H(i,:)=H(i,:)*m(i);
    W(:,i)=W(:,i)/m(i);
end

% Set to zero the columns or rows that were not included in the iteration
[n_o,m_o]=size(V);
Hnew=[];
Wnew=[];
for l=1:length(ind_cond)
    if ind_cond(l)==1
        Hnew=[zeros(r,1) H];
        H=Hnew;
    elseif ind_cond(l)==m_o
        Hnew=[H zeros(r,1)];
        H=Hnew;
    else 
        for k=1:m_o
            if ind_cond(l)==k
                Hnew=[H(:,1:k-1) zeros(r,1) H(:,k:end)];
                H=Hnew; break
            else
                Hnew=H;
            end
        end
    end
end
for l=1:length(ind)
    if ind(l)==1
        Wnew=[zeros(1,r); W];
        W=Wnew;
    elseif ind(l)==n_o
        Wnew=[W; zeros(1,r)];
        W=Wnew;
    else 
        for k=1:n_o
            if ind(l)==k
                Wnew=[W(1:k-1,:); zeros(1,r); W(k:end,:)];
                W=Wnew; break
            else
                Wnew=W;
            end
        end
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Wnnmf=W;
Hnnmf=H;
err1=err;


%Determine the percentage of variability accounted for by each muscle
%synergy:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
W=Wnnmf; C=Hnnmf;

[nmuscles ncond]=size(DATA);
[nsyn ndum]=size(C);

%Calculate reconstructed values
ReconData=W*C;

Wstr=mat2cell(W,[nmuscles],[ones(1,nsyn)]); % <1x7>
Cstr=mat2cell(C,[ones(1,nsyn)],[ncond]);
Recon_parts=cellfun(@(N1,N2) N1*N2,Wstr,Cstr','UniformOutput',0);

% Calculate contribution of each synergy to overall variability
synconstr=cellfun(@(D1) 100*sum(sum(D1))/sum(sum(ReconData)),Recon_parts,'UniformOutput',0);
syncontr=cell2mat(synconstr);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%Plot the extracted muscle synergies and their activation coefficients
%(weightings)
axes(handles.axes4)
b= bar(Wnnmf(:,1));set(b,'FaceColor','b','EdgeColor',[1 1 1]);
title('W')
guidata(hObject, handles); %updates the handles
axes(handles.axes5)
b= bar(Wnnmf(:,2));set(b,'FaceColor','b','EdgeColor',[1 1 1]);
sync_W_axes(handles); % same y-scale on all W bar plots
clear b
guidata(hObject, handles); %updates the handles

%Print the VAF by each synergy to the GUI
set(handles.VAF_NMF1,'String', num2str(round(syncontr(1))));
set(handles.VAF_NMF2,'String', num2str(round(syncontr(2))));

[nmus nsyn]=size(Wnnmf);
axes(handles.axes6)
cla
hold on;
first=Wnnmf(:,1);
second=Wnnmf(:,2);
if nmus==2;
    plot_data_sets(DATA); % data set 1 blue, data set 2 orange (as in the data plot)
    hold on;
    plot([0 first(1)],[0 first(2)],'Color','b');
    plot([0 second(1)],[0 second(2)],'Color','g');
elseif nmus==3;
    plot3(DATA(1,:),DATA(2,:),DATA(3,:),'.');
    hold on;
    plot3([0 first(1)],[0 first(2)],[0 first(3)],'Color','b');
    plot3([0 second(1)],[0 second(2)],[0 second(3)],'Color','g');
end
set(handles.axes6,'XLim',[-1 1],'YLim',[-1 1],'PlotBoxAspectRatio',[1 1 1]); % fixed, square axes
guidata(hObject, handles); %updates the handles




%***********************************************
% Extract and Plot PCA synergies              *
%***********************************************

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DATA nmus

% princomp was removed from MATLAB; pca is its replacement (same outputs)
[Wpca,SCORE,pcvars] = pca(DATA','Economy',true);
%[residuals,reconstructed] = pcares(x,ndim)
VAFs_cum=cumsum(pcvars./sum(pcvars) * 100);

axes(handles.axes8)
b= bar(Wpca(:,1));set(b,'FaceColor','b','EdgeColor',[1 1 1]);
title('W')
guidata(hObject, handles); %updates the handles
[nmus nsyn]=size(Wpca);
if nsyn>1
    axes(handles.axes9)
    b= bar(Wpca(:,2));set(b,'FaceColor','b','EdgeColor',[1 1 1]);
    clear b
    guidata(hObject, handles); %updates the handles
end
sync_W_axes(handles); % same y-scale on all W bar plots

%Print the VAF by each synergy to the GUI
set(handles.VAF_PCA1,'String', num2str(round(VAFs_cum(1))));
set(handles.VAF_PCA2,'String', num2str(round(VAFs_cum(2)-VAFs_cum(1))));

axes(handles.axes7)
cla
hold on;
first=Wpca(:,1);
second=Wpca(:,2);
if nmus==2;
    plot_data_sets(DATA); % data set 1 blue, data set 2 orange (as in the data plot)
    hold on;
    %quiver([0;0],[0;0],[first],[second]);
    plot([0 first(1)],[0 first(2)],'Color','b');
    plot([0 second(1)],[0 second(2)],'Color','g');
elseif nmus==3;
    plot3(DATA(1,:),DATA(2,:),DATA(3,:),'.');
    hold on;
    %quiver3([0;0;0],[0;0;0],[first],[second]);
    plot3([0 first(1)],[0 first(2)],[0 first(3)],'Color','b');
    plot3([0 second(1)],[0 second(2)],[0 second(3)],'Color','g');
end
set(handles.axes7,'XLim',[-1 1],'YLim',[-1 1],'PlotBoxAspectRatio',[1 1 1]); % fixed, square axes
guidata(hObject, handles); %updates the handles




%***********************************************
% Reset all fields and plots                   *
%***********************************************

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i=1:5
    eval(['axes(handles.axes',num2str(i),')']);
    cla
end
axes(handles.axes6);cla;
axes(handles.axes7);cla;
axes(handles.axes8);cla;
axes(handles.axes9);cla;

set([handles.VAF_NMF1 handles.VAF_NMF2 handles.VAF_PCA1 handles.VAF_PCA2],'String',''); % clear old VAFs
set(handles.NumPoints,'String','100');
global C Cb
C=rand(2,100); Cb=rand(2,100);
reset_W(handles);
reset_weights(handles);
PCA_example_GUI_2D('pushbutton2_Callback',handles.pushbutton2,[],handles); % plots synergies and data


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3



function VAF_NMF1_Callback(hObject, eventdata, handles)
% hObject    handle to VAF_NMF1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VAF_NMF1 as text
%        str2double(get(hObject,'String')) returns contents of VAF_NMF1 as a double


% --- Executes during object creation, after setting all properties.
function VAF_NMF1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VAF_NMF1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VAF_NMF2_Callback(hObject, eventdata, handles)
% hObject    handle to VAF_NMF2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VAF_NMF2 as text
%        str2double(get(hObject,'String')) returns contents of VAF_NMF2 as a double


% --- Executes during object creation, after setting all properties.
function VAF_NMF2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VAF_NMF2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VAF_PCA1_Callback(hObject, eventdata, handles)
% hObject    handle to VAF_PCA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VAF_PCA1 as text
%        str2double(get(hObject,'String')) returns contents of VAF_PCA1 as a double


% --- Executes during object creation, after setting all properties.
function VAF_PCA1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VAF_PCA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VAF_PCA2_Callback(hObject, eventdata, handles)
% hObject    handle to VAF_PCA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VAF_PCA2 as text
%        str2double(get(hObject,'String')) returns contents of VAF_PCA2 as a double


% --- Executes during object creation, after setting all properties.
function VAF_PCA2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VAF_PCA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savefig.
function savefig_Callback(hObject, eventdata, handles)
% hObject    handle to savefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fn = 'testing';
print( gcf, '-depsc2', fn )


% --- Rearrange the controls so the window has less blank space.
% Units are characters. Top band: synergies | data sets and weights | data.
% Bottom band: NMF (left) and PCA (right).
function compact_layout_2d(fig,h)
% Lay out only once: running PCA_example_GUI_2D again reuses the open window
if isappdata(fig,'compactLayout'), return; end
setappdata(fig,'compactLayout',true);
W=100; H=35;
pos=get(fig,'Position');
set(fig,'Units','characters','Position',[pos(1) pos(2)+pos(4)-H W H]);
P=@(o,p) set(o,'Units','characters','Position',p);
L=@(o) set(o,'HorizontalAlignment','left');
yt=33.1;                                   % row of section titles
% ---- Defined synergies, with editable weights
P(h.text5,[3 yt 26 1.2]);
wW=15;                                     % width of every W bar plot (defined, NMF, PCA)
P(h.axes1,[6 26.6 wW 4.6]);
P(h.axes2,[6 20.2 wW 4.6]);
P(h.pushbutton2,[6 16.5 25.5 2.1]);
Wdef=[0.5 1;1 0.5];                        % chapter example
ys=[26.6 20.2];
for k=1:2
    for m=1:2
        y=ys(k)+1.9*(2-m);
        uicontrol(fig,'Style','text','Units','characters','Position',[21.5 y+0.2 4 1.2], ...
            'String',sprintf('M%d',m),'HorizontalAlignment','right');
        uicontrol(fig,'Style','edit','Units','characters','Position',[25.5 y 6 1.6], ...
            'String',num2str(Wdef(m,k)),'Tag',sprintf('W%d%d',m,k), ...
            'TooltipString',sprintf('Weight of muscle %d in synergy %d (0 or more)',m,k), ...
            'Callback',@(o,e) PCA_example_GUI_2D('pushbutton2_Callback',o,e,guidata(o)));
    end
end
% ---- Data sets: number of points, resample, and C1/C2 weights for each set
P(h.text10,[33 yt 31 1.2]); L(h.text10); set(h.text10,'String','Data sets: C1 and C2 weights');
P(h.text1,[33 30.5 8 1.2]); L(h.text1); set(h.text1,'String','Points');
P(h.NumPoints,[41 30.2 9 1.8]);
set(h.NumPoints,'String','100');           % default number of points
uicontrol(fig,'Style','pushbutton','Units','characters','Position',[51 30.1 12.5 2], ...
    'String','Resample','TooltipString','Draw new random C values for both data sets', ...
    'Callback',@(o,e) resample_C(guidata(o)));
set([h.text6 h.text7 h.text8 h.text9],'Visible','off');   % 0/1 slider labels: the boxes show the value
uicontrol(fig,'Style','text','Units','characters','Position',[33 28.3 30 1.3], ...
    'String','Data set 1 (blue)','FontWeight','bold','HorizontalAlignment','left');
P(h.text2,[33 26.4 4 1.3]); P(h.slider1,[37 26.4 17.5 1.5]); P(h.slider1_editText,[55.5 26.2 8 1.8]);
P(h.text3,[33 24.3 4 1.3]); P(h.slider2,[37 24.3 17.5 1.5]); P(h.slider2_editText,[55.5 24.1 8 1.8]);
P(h.checkbox3,[33 22.0 30 1.6]); set(h.checkbox3,'String','Data set 2 (orange)','FontWeight','bold', ...
    'TooltipString','Add a second data set with its own random sample and C1/C2 weights', ...
    'Callback',@(o,e) set2_checkbox(guidata(o)));
for k=1:2
    y=20.0-2.1*(k-1);
    uicontrol(fig,'Style','text','Units','characters','Position',[33 y 4 1.3],'String',sprintf('C%d',k));
    uicontrol(fig,'Style','slider','Units','characters','Position',[37 y 17.5 1.5],'Min',0,'Max',1, ...
        'Tag',sprintf('set2_C%d',k),'Callback',@(o,e) sync_set2(o,'slider'));
    uicontrol(fig,'Style','edit','Units','characters','Position',[55.5 y-0.2 8 1.8], ...
        'Tag',sprintf('set2_C%d_edit',k),'Callback',@(o,e) sync_set2(o,'edit'));
end
% ---- Data plot, Plot Data, Reset and Save Figure
P(h.text4,[70 yt 8 1.2]); L(h.text4);
P(h.pushbutton5,[78 32.8 9 1.9]);
P(h.savefig,[88 32.8 10 1.9]);
P(h.axes3,[71 17.6 26 14.3]);
set(h.pushbutton1,'Visible','off');      % the data replot automatically
% Live updates while a slider is being dragged
for sl=[h.slider1 h.slider2 findall(fig,'Tag','set2_C1') findall(fig,'Tag','set2_C2')]
    addlistener(sl,'ContinuousValueChange',@(o,e) feval(get(o,'Callback'),o,[]));
end
% ---- NMF (left) and PCA (right): button, two synergy bar plots with VAF, scatter
set([h.text11 h.text12],'Visible','off');  % the buttons already name each block
blk={2, h.pushbutton3, h.axes4, h.axes5, h.text16, h.VAF_NMF1, h.text17, h.VAF_NMF2, h.axes6; ...
     51, h.pushbutton4, h.axes8, h.axes9, h.text18, h.VAF_PCA1, h.text19, h.VAF_PCA2, h.axes7};
for k=1:2
    x=blk{k,1};
    P(blk{k,2},[x 13.2 46 2]);
    P(blk{k,3},[x+4 8.3 wW 3.4]);
    P(blk{k,4},[x+4 3.0 wW 3.4]);
    P(blk{k,5},[x+21.5 10.3 5 1.2]); P(blk{k,6},[x+21.5 8.3 7 2]);
    P(blk{k,7},[x+21.5 5.0 5 1.2]);  P(blk{k,8},[x+21.5 3.0 7 2]);
    P(blk{k,9},[x+31 2.6 16.5 9.8]);       % square scatter (see PlotBoxAspectRatio)
end
h.figure1=fig;
reset_weights(h);


% --- Default C1/C2 weights: set 1 mostly C1, set 2 mostly C2; set 2 off.
function reset_weights(handles)
f=handles.figure1;
set(handles.slider1,'Value',1);    set(handles.slider1_editText,'String','1');
set(handles.slider2,'Value',0.15); set(handles.slider2_editText,'String','0.15');
set(findall(f,'Tag','set2_C1'),'Value',0.15); set(findall(f,'Tag','set2_C1_edit'),'String','0.15');
set(findall(f,'Tag','set2_C2'),'Value',1);    set(findall(f,'Tag','set2_C2_edit'),'String','1');
set(handles.checkbox3,'Value',0);
toggle_set2(handles);


% --- Turn the data set 2 controls on or off with its checkbox.
function toggle_set2(handles)
f=handles.figure1;
if get(handles.checkbox3,'Value')==1, en='on'; else, en='off'; end
set([findall(f,'Tag','set2_C1') findall(f,'Tag','set2_C2') ...
     findall(f,'Tag','set2_C1_edit') findall(f,'Tag','set2_C2_edit')],'Enable',en);


% --- Keep each data set 2 slider and its number box in step (values 0 to 1).
function sync_set2(o,src)
f=ancestor(o,'figure'); t=get(o,'Tag');
if strcmp(src,'slider')
    set(findall(f,'Tag',[t '_edit']),'String',num2str(round(get(o,'Value'),2)));
else
    v=str2double(get(o,'String'));
    if isnan(v), v=0; end
    v=min(max(v,0),1);
    set(o,'String',num2str(v));
    set(findall(f,'Tag',strrep(t,'_edit','')),'Value',v);
end
replot_data(guidata(f)); % update the data plot right away


% --- Draw new random C values for both data sets and replot the data.
function resample_C(handles)
global C Cb
n=str2double(get(handles.NumPoints,'String'));
if isnan(n) || n<2, n=100; set(handles.NumPoints,'String','100'); end
C=rand(2,round(n)); Cb=rand(2,round(n));
PCA_example_GUI_2D('pushbutton1_Callback',handles.pushbutton1,[],handles);


% --- Read the defined synergies W (muscles x synergies) from the W boxes.
% Blank, non-numeric or negative entries are set to 0.
function W=read_W(handles)
W=zeros(2,2);
for k=1:2
    for m=1:2
        e=findall(handles.figure1,'Tag',sprintf('W%d%d',m,k));
        if isempty(e), W=[0.5 1;1 0.5]; return; end   % boxes missing: use the chapter example
        e=e(1);
        v=str2double(get(e,'String'));
        if isnan(v) || v<0, v=0; set(e,'String','0'); end
        W(m,k)=v;
    end
end


% --- Put the chapter example synergies back in the W boxes.
function reset_W(handles)
Wdef=[0.5 1;1 0.5];
for k=1:2
    for m=1:2
        set(findall(handles.figure1,'Tag',sprintf('W%d%d',m,k)),'String',num2str(Wdef(m,k)));
    end
end


% --- Plot the data in the NMF/PCA panels with the same colors as the data plot:
% data set 1 blue, data set 2 orange.
function plot_data_sets(DATA)
global DATA1
n1=size(DATA1,2);
if isempty(DATA1) || n1>=size(DATA,2), n1=size(DATA,2); end
plot(DATA(1,1:n1),DATA(2,1:n1),'.','Color',[0 0.45 0.74]);
if n1<size(DATA,2)
    hold on
    plot(DATA(1,n1+1:end),DATA(2,n1+1:end),'.','Color',[0.85 0.33 0.1]);
end


% --- Consistent W bar plot scaling: defined and NMF synergies 0 to 1;
% PCA synergies 0 to 1, or -1 to 1 when a PC has a negative weight.
function sync_W_axes(handles)
nonneg=[handles.axes1 handles.axes2 handles.axes4 handles.axes5];
pcs=[handles.axes8 handles.axes9];
set([nonneg pcs],'XLim',[0.4 2.6]);
set(nonneg,'YLim',[0 1],'YTick',[0 0.5 1]);
v=[];
for a=pcs
    b=findall(a,'Type','bar');
    for k=1:numel(b), v=[v b(k).YData(:)']; end %#ok<AGROW>
end
if any(v<0)
    set(pcs,'YLim',[-1 1],'YTick',[-1 0 1]);
else
    set(pcs,'YLim',[0 1],'YTick',[0 0.5 1]);
end


% --- Replot the data (Plot Data) with the current synergies, weights and points.
function replot_data(handles)
PCA_example_GUI_2D('pushbutton1_Callback',handles.pushbutton1,[],handles);


% --- Data set 2 checkbox: enable its controls and replot.
function set2_checkbox(handles)
toggle_set2(handles);
replot_data(handles);
