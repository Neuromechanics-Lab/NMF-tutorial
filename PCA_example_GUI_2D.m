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

    W=[0.5 1;1 0.5]; %example for chapter: M2=2*M1, M2=0.5*M1
    %W=[0.3 1;1 0.3]; %example for chapter: M2=2*M1, M2=0.5*M1

[nmus nsyn]=size(W);

%PLOT THE SYNERGIES
axes(handles.axes1)
b= bar(W(:,1));set(b,'FaceColor','b','EdgeColor',[1 1 1]);
title('W')
guidata(hObject, handles); %updates the handles
axes(handles.axes2)
b= bar(W(:,2));set(b,'FaceColor','b','EdgeColor',[1 1 1]);
clear b
guidata(hObject, handles); %updates the handles



%***********************************************
% Defining the synergy coefficients            *
%***********************************************

function NumPoints_Callback(hObject, eventdata, handles)
% hObject    handle to NumPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumPoints as text
%        str2double(get(hObject,'String')) returns contents of NumPoints as a double

numpts=str2double(get(hObject,'String')); %returns contents of NumPoints as a double
global C
C=rand(2,numpts); %random numbers for C


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

global W C C1_wt C2_wt DATA nmus DATA1
axes(handles.axes3)
Cnew=[C(1,:)*C1_wt;C(2,:)*C2_wt]; %The user can set different weightings for C1 and C2
dat=W*Cnew;

dimcheck2=get(handles.checkbox3,'Value');
if dimcheck2==1; %If the box labeled "Two Data Sets" is checked
    DATA2=dat;
    clear dat
    DATA=[DATA1 DATA2]; %combine the previous data and the current data
else
    DATA1=dat; %Store the data in DATA1 if only one dataset is used, and DATA1 is never used further if only one dataset is used
    DATA=DATA1;
end

if nmus==2;plot(DATA(1,:),DATA(2,:),'.');
    %axis([0 1 0 1]);set(gca,'Clipping','off');
elseif nmus==3;plot3(DATA(1,:),DATA(2,:),DATA(3,:),'.');end

% DATA1=W*C1;
% DATA2=W*C2;
% plot(DATA1(1,:),DATA1(2,:),'.');
% hold on; plot(DATA2(1,:),DATA2(2,:),'.b');
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
    plot(DATA(1,:),DATA(2,:),'.');
    hold on;
    plot([0 first(1)],[0 first(2)],'Color','b');
    plot([0 second(1)],[0 second(2)],'Color','g');
elseif nmus==3;
    plot3(DATA(1,:),DATA(2,:),DATA(3,:),'.');
    hold on;
    plot3([0 first(1)],[0 first(2)],[0 first(3)],'Color','b');
    plot3([0 second(1)],[0 second(2)],[0 second(3)],'Color','g');
end
guidata(hObject, handles); %updates the handles

figure(1)
if nmus==2;
    plot(DATA(1,:),DATA(2,:),'.');
    hold on;
    plot([0 first(1)],[0 first(2)],'Color','b');
    plot([0 second(1)],[0 second(2)],'Color','g');
elseif nmus==3;
    plot3(DATA(1,:),DATA(2,:),DATA(3,:),'.');
    hold on;
    plot3([0 first(1)],[0 first(2)],[0 first(3)],'Color','b');
    plot3([0 second(1)],[0 second(2)],[0 second(3)],'Color','g');
end



%***********************************************
% Extract and Plot PCA synergies              *
%***********************************************

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DATA nmus

[Wpca,SCORE,pcvars] = princomp(DATA','econ');
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

%Print the VAF by each synergy to the GUI
set(handles.VAF_PCA1,'String', num2str(round(VAFs_cum(1))));
set(handles.VAF_PCA2,'String', num2str(round(VAFs_cum(2)-VAFs_cum(1))));

axes(handles.axes7)
cla
hold on;
first=Wpca(:,1);
second=Wpca(:,2);
if nmus==2;
    plot(DATA(1,:),DATA(2,:),'.');
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
guidata(hObject, handles); %updates the handles

figure(2)
if nmus==2;
    plot(DATA(1,:),DATA(2,:),'.');
    hold on;
    plot([0 first(1)],[0 first(2)],'Color','b');
    plot([0 second(1)],[0 second(2)],'Color','g');
elseif nmus==3;
    plot3(DATA(1,:),DATA(2,:),DATA(3,:),'.');
    hold on;
    plot3([0 first(1)],[0 first(2)],[0 first(3)],'Color','b');
    plot3([0 second(1)],[0 second(2)],[0 second(3)],'Color','g');
end



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
axes(handles.axes7);cla;
axes(handles.axes8);cla;
axes(handles.axes9);cla;

set(handles.slider1,'Value',0);
set(handles.slider1_editText,'String','0');
set(handles.slider2,'Value',0);
set(handles.slider2_editText,'String','0');
set(handles.NumPoints,'String','0');


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


