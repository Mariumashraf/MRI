function varargout = FINAL(varargin)
% FINAL MATLAB code for FINAL.fig
%      FINAL, by itself, creates a new FINAL or raises the existing
%      singleton*.
%
%      H = FINAL returns the handle to a new FINAL or the handle to
%      the existing singleton*.
%
%      FINAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINAL.M with the given input arguments.
%
%      FINAL('Property','Value',...) creates a new FINAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FINAL_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FINAL_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FINAL

% Last Modified by GUIDE v2.5 11-May-2018 12:00:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FINAL_OpeningFcn, ...
                   'gui_OutputFcn',  @FINAL_OutputFcn, ...
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


% --- Executes just before FINAL is made visible.
function FINAL_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FINAL (see VARARGIN)

% Choose default command line output for FINAL
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FINAL wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FINAL_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function Rx=xrot(phi)

Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
function Ry=yrot(phi)

Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];
function Rz=zrot(phi)

Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];
function [Afp,Bfp]=freeprecess(T,T1,T2,df)
%
%	Function simulates free precession and decay
%	over a time interval T, given relaxation times T1 and T2
%	and off-resonance df.  Times in ms, off-resonance in Hz.

phi = 2*pi*df*T/1000;	% Resonant precession,radians aroound z(around B0.
E1 = exp(-T/T1);	
E2 = exp(-T/T2);

Afp = [E2 0 0;0 E2 0;0 0 E1]*zrot(phi);
Bfp = [0 0 1-E1]';
% --- Executes on button press in show_decay.
function show_decay_Callback(hObject, eventdata, handles)
% hObject    handle to show_decay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dT = 1;		% 1ms delta-time.
T = 1000;	% total duration
N = ceil(T/dT)+1; % number of time steps.
df = 10;	% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
% ===== Get the Propagation Matrix ======

[A,B] = freeprecess(dT,T1,T2,df);


% ===== Simulate the Decay ======

M = zeros(3,N);	% Keep track of magnetization at all time points.
M(:,1)=[1;0;0];	% Starting magnetization.

for k=2:N
	M(:,k) = A*M(:,k-1)+B;
end;


% ===== Plot the Results ======

time = [0:N-1]*dT;
figure
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;


% --- Executes on button press in showRFpulse.
function showRFpulse_Callback(hObject, eventdata, handles)
% hObject    handle to showRFpulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
dT = 1;
TR = 500;	% ms.
flip = pi/3;	% radians.
Ntr = round(TR/dT);
Nex = 10;	% 20 excitations.

M = [0;0;1];
Rflip = yrot(flip);
[A1,B1] = freeprecess(dT,T1,T2,df);


M(1,Nex*Ntr)=0;	%	Allocate to record all M's.
		% 	Not necessary, but makes program faster.

Mcount=1;
for n=1:Nex
	M(:,Mcount) = Rflip*M(:,Mcount);	

	for k=1:Ntr
		Mcount=Mcount+1;
		M(:,Mcount)=A1*M(:,Mcount-1)+B1;
	end;
end;

time = [0:Mcount-1]*dT;
figure
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;


% --- Executes on button press in non_homogenous.
function non_homogenous_Callback(hObject, eventdata, handles)
% hObject    handle to non_homogenous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a= 0:1:10;
m = 10*exp(-a/100)+200;
axes(handles.axes1);
plot(a,m)


% --- Executes on button press in plotimageinK.
function plotimageinK_Callback(hObject, eventdata, handles)
% hObject    handle to plotimageinK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [filename pathname] = uigetfile({'*.png';'*.tif'},'File Selector');
% image = strcat(pathname, filename);
% [rows, cols] = size(image);
% image1 = imread(image);
load mri; % Load Matlab's default set of MR images
image1=D
 KS=zeros(size(image1)); % Allocate memory for k-space
 for i=1:size(image1,4) % Compute the Fourier transform of each image
    KS(:,:,:,i)=fftshift(fft2((image1(:,:,:,i))));
 end
% %plot image
% axes(handles.axes2)
% imagesc(squeeze(image1(:,:,:, 1))),colormap gray
%plot image in kspace
axes(handles.axes5)
imagesc(squeeze(abs(KS(:,:,:, 1)))); caxis([0 10000]),colormap gray


                     % Display image
% FF = fft2(image1);                   % Take FFT
% 
% IFF = ifft2(FF);                 % take IFFT
% 
% FINAL_IM = uint8(real(IFF));      % Take real part and convert back to UINT8
% 
% axes(handles.axes4)
% imshow(FINAL_IM),colormap gray





% Image1_FFT=fft2(image1);
% 
% %figure:
% figure
% imshow(abs(fftshift(Image1_FFT)),[10 100]), colormap gray
%  figure
% %figure:
% imshow(angle(fftshift(Image1_FFT)),[-pi pi]), colormap gray





% --- Executes on button press in loadimagee.
function loadimagee_Callback(hObject, eventdata, handles)
% hObject    handle to loadimagee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
IM=imread('D:\TASKMRI\mri.tif');  % Read in a image

[filename pathname] = uigetfile({'*.tif';'*.png'},'File Selector');
image = strcat(pathname, filename);
[rows, cols] = size(image);
image1 = imread(image);
load mri; % Load Matlab's default set of MR images
image1=D
IM=imread('D:\TASKMRI\mri.tif');
axes(handles.axes2)
imagesc(squeeze(IM(:,:,:, 1))),colormap gray

Image1_FFT=fft2(IM);

axes(handles.axes3)
imshow(abs(fftshift(Image1_FFT)),[10 100]), colormap gray
axes(handles.axes4)
imshow(angle(fftshift(Image1_FFT)),[-pi pi]), colormap gray





function P= gradiant(P,Gx,gero,hObject, eventdata, handles)
 [m,n]=size(P);
 handles.phase=360/m;
m=m/2;
n=n/2;
for y=0:m-1
    for x=0:n-1
w=(gero*10^6)*(3-Gx*x+handles.phase*y);
Ry=cos(w);
P(n-x,m+y)=Ry*P(n-x,m+y);
    end
end

for y=0:m-1
    for x=0:n-1
w=(gero*10^6)*(3-Gx*x-handles.phase*y);
Ry=cos(w);
P(n-x,m-y)=Ry*P(n-x,m-y);
    end
end
for y=0:m-1
    for x=0:n-1
w=(gero*10^6)*(3+Gx*x+handles.phase*y);
Ry=cos(w);
P(n+x,m+y)=Ry*P(n+x,m+y);
    end
end
for y=0:m-1
    for x=0:n-1
w=(gero*10^6)*(3+Gx*x-handles.phase*y);
Ry=cos(w);
P(n+x,m-y)=Ry*P(n+x,m-y);
    end
end


% --- Executes on button press in plotfromks.
function plotfromks_Callback(hObject, eventdata, handles)
% hObject    handle to plotfromks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

A = phantom('Modified Shepp-Logan',256);
P = phantom('Modified Shepp-Logan',256);
Gx=.3;
handles.p=gradiant(P,Gx,42.6,hObject, eventdata, handles);
axes(handles.axes6)
imagesc(handles.p)

imagesc(log(abs(fft(handles.p)))) , colormap gray ;
guidata(hObject, handles);


% --- Executes on button press in reconstructed_image.
function reconstructed_image_Callback(hObject, eventdata, handles)
% hObject    handle to reconstructed_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes7)
img = abs(ifft(fft(handles.p)));
imagesc(img) , colormap gray;

h = HeatMap(img);
imagesc(h)
guidata(hObject, handles);


% --- Executes on button press in fat_oil.
function fat_oil_Callback(hObject, eventdata, handles)
% hObject    handle to fat_oil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 A = phantom('Modified Shepp-Logan',256);
 Gx=.7;
 [m,n]=size(A);
 
oil=gradiant(A,Gx,10.71,hObject, eventdata, handles);
  oil=oil(:,1:n/2);
water=gradiant(A,Gx,42.6,hObject, eventdata, handles);
  water=water(:,n/2:n);
  A=cat(2,oil,water);
%imagesc(log(abs(fft(A)))) , colormap gray ;
img = abs(ifft(fft(A)));
axes(handles.axes8)
imagesc(img) , colormap gray;

hmo = HeatMap(img);
imagesc(hmo)

  guidata(hObject, handles);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rotate.
function rotate_Callback(hObject, eventdata, handles)
% hObject    handle to rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
phi = str2num(get(handles.edit1,'string'))
t= 0:1:100;
m0 = 100*exp(-t/200);
m = 1-m0;
y = zeros(size(m0));
m0matrix = [y(:) y(:) m0(:)];

Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];
Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];
rotationaxis = get(handles.edit2,'string')
if rotationaxis == 'x'
    rotation = m0matrix*Rx;
elseif rotationaxis == 'y'
    rotation = m0matrix*Ry;
else rotationaxis == 'z'
    rotation = m0matrix*Rz;
end
disp(rotationaxis(:))
axes(handles.axes9)
plot3(m0matrix(:,1),m0matrix(:,2),m0matrix(:,3));
axes(handles.axes10)
comet3(rotation(:,1),rotation(:,2),rotation(:,3));
axes(handles.axes11)
comet(t,m)
axes(handles.axes12)
comet(t,m0)
