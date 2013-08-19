
close all;

open('costs_orig.fig');
h = gcf; %current figure handle
axesObjs = get(h, 'Children'); %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type'); %type of low-level graphics object
xdata_orig = get(dataObjs, 'XData'); %data from low-level grahics objects
ydata_orig = get(dataObjs, 'YData');
zdata_orig = get(dataObjs, 'ZData');

close(h);

open('costs_vs.fig');
h = gcf; %current figure handle
axesObjs = get(h, 'Children'); %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type'); %type of low-level graphics object
xdata_vs = get(dataObjs, 'XData'); %data from low-level grahics objects
ydata_vs = get(dataObjs, 'YData');
zdata_vs = get(dataObjs, 'ZData');

close(h);

open('costs_vt.fig');
h = gcf; %current figure handle
axesObjs = get(h, 'Children'); %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type'); %type of low-level graphics object
xdata_vt = get(dataObjs, 'XData'); %data from low-level grahics objects
ydata_vt = get(dataObjs, 'YData');
zdata_vt = get(dataObjs, 'ZData');

close(h);



for i=1:length(xdata_vs)
    
   xdata_vs_new(i) = xdata_vs(i);
   ydata_vs_new(i) = ydata_vs(i);
   if zdata_vs(i) ~= 0
       zdata_vs_new(i) = zdata_vs(i);
   else
       zdata_vs_new(i) = Inf;
   end
end


for i=1:length(xdata_vt)
   xdata_vt_new(i) = xdata_vt(i);
   ydata_vt_new(i) = ydata_vt(i);
   if zdata_vt(i) ~= 0
       zdata_vt_new(i) = zdata_vt(i);
   else
       zdata_vt_new(i) = Inf;
   end
end


open('costs_orig.fig');
h = figure('name','State Costs of the System - Filtered (KEEP_VALID_STATES)'); scatter3(xdata_vs_new, ydata_vs_new, zdata_vs_new, 5, zdata_vs_new); colormap(jet); xlabel('State Variable x(1)'); ylabel('State Variable x(2)'); zlabel('State Cost');
figure('name','State Costs of the System - Filtered (KEEP_VALID_TRANSITIONS)'); scatter3(xdata_vt_new, ydata_vt_new, zdata_vt_new, 5, zdata_vt_new); colormap(jet); xlabel('State Variable x(1)'); ylabel('State Variable x(2)'); zlabel('State Cost');

