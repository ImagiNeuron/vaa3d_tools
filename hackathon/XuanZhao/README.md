## neuronQC
### ��Ҫ����
���ڼ����Զ��ؽ���Ԫ��һЩӲ��ָ�����
### ���÷���
```
vaa3d_path /x dll_path /f neuronQC_batch /i swc_dir csv_path /p short_length_thres node_length_thres loop_three_bifircation
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·����Ҳ�ɸ�Ϊ�����neuronQC��<br>
swc_dir: ����swc���ļ���swc���ж�����ò������Ϊ��������<br>
csv_path: csv�ļ�·�������ڼ�¼���<br>
short_length_thres���̷�֧�ĳ�����ֵ���ò�����Ĭ��Ϊ10<br>
node_length_thres��swc����node֮��ĳ�����ֵ���ò�����Ĭ��Ϊ15<br>
loop_three_bifircation���Ƿ���loop�����ֲ棬Ĭ�ϼ��<br>

## imageProcess
### ��Ҫ����
һЩ���ܻ����õ�һЩ��������
### ���������÷���
#### enhanceImage
��ͼ��������ļ򵥴���
```
f(i) = ((i/255)^2/3)*255
```
����
```
vaa3d_path /x dll_path /f enhanceImage /i img_path
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·����Ҳ�ɸ�Ϊ�����,����.pro�ļ��鿴,����Ϊ$$qtLibrary(name),����nameΪ�������<br>
img_path: ͼ��·��<br>
#### get_2d_image2
��ͼ��Ͷ�Ӧswc������һ����һ��mip<br>

����
```
vaa3d_path /x dll_path /f get_2d_image2 /i swc_path1 swc_path2 ... /p img_path
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
swc_pathi: swc·��,�ɶ��<br>
img_path��ͼ��·��<br>
#### getSWCL0image
��ȡswc��Ӧ��ͼ��ͼ��̫����ȡ����<br>
����
```
vaa3d_path /x dll_path /f getSWCL0image /i swc_path /p brain_path times /o out_dir 
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
swc_path: swc·��<br>
brain_path��terafly��ʽȫ��·��<br>
times����һ��ֱ��ʣ�1��2��4��8...֮�������<br>
out_dir: �����ļ���·��<br>
#### convertTeraflyDataTov3draw
��terafly��ʽ���ݵ�ĳһ��ֱ���ƴ�ӳ�һ��ͼ<br>
����
```
vaa3d_path /x dll_path /f convertTeraflyDataTov3draw /i brain_path /p resolution /o out_path 
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
brain_path��terafly��ʽȫ��ĳһ��ֱ��ʵ�·��<br>
resolution����һ��ֱ��ʣ�1��2��4��8...֮�����,1Ϊ��߷ֱ��ʣ�2Ϊ�θ�...��<br>
out_path: ����ͼ���·��<br>
#### bilateralfilter
��ͼ����һ��˫���˲�����<br>
����
```
vaa3d_path /x dll_path /f bilateralfilter /i img_path /p is_normal space_sigma_xy space_sigma_z color_sigma 
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
img_path��ͼ��·��<br>
is_normal���Ƿ��Ƚ�ͼ��ת����0-255��Ĭ����<br>
space_sigma_xy: ��xy����Ĵ�С,Ĭ��Ϊ2<br>
space_sigma_z: ��xy����Ĵ�С��Ĭ��Ϊ1<br>
color_sigma: color sigma��С��Ĭ��Ϊ35<br>
#### changeContrast
��ͼ����һ�����Ա任<br>
����
```
vaa3d_path /x dll_path /f changeContrast /i img_path /p pencert_down pencert_up
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
img_path��ͼ��·��<br>
pencert_down��ת��Ϊ0����㣬�ٷ���<br>
pencert_up: ת��Ϊ255����㣬�ٷ���<br>
#### v3draw2tif
��v3draw��ʽ��ͼ������תΪtif��ʽ<br>
����
```
vaa3d_path /x dll_path /f v3draw2tif /i img_dir /o out_dir
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
img_dir��Ҫת��ͼ����ļ���<br>
out_dir������ļ���<br>

## dynamicApp2
### ��Ҫ����
����app2��һЩ�Ķ�����
### ���������÷���
#### dynamicApp2
��һ��ͼ����ڽ���׷�٣���ÿ��׷�ٵĳ��ȼ���һ�����ƣ����������з���׷�ٽ�����֤��Ȼ���ټ���׷�٣�ֱ��׷�ٵ�tip��<br>
����
```
vaa3d_path /x dll_path /f dynamicApp2 /i img_path /p ... length_thres ...
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
img_path: ͼ��·��<br>
length_thres: �������Ʋ���<br>
����������app2һ��<br>
#### ultratracerAxonTerafly
��terafly���ݸ�ʽ������ͻ׷��<br>
����
```
vaa3d_path /x dll_path /f ultratracerAxonTerafly /i swc_path brain_path tmp_image_dir
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
swc_path: ��ʼswc·������ʼswcΪһ����ͻ���źţ��ú������ڸ��źŵĻ����ϼ�������׷��<br>
brain_path: teraflyȫ��·��<br>
tmp_image_dir���м������������򲻴�<br>

## consensus
### ��Ҫ����
consensus��һЩ����
### ���������÷���
#### consensus
��һ��ͼ����ڸı������ø��ַ�ʽ����׷�٣���׷�ٵ�ͼ������ںϣ�Ȼ�����ںϵ�ͼ�����׷�٣����õ�׷�ٽ��
```
vaa3d_path /x dll_path /f consensus /i img_path marker_path /p kmeans_th
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
img_path: ͼ��·��<br>
marker_path: soma marker·��<br>
kmeans_th: �Ƿ���kmeans�ķ������õ�һ����ֵ��Ĭ�ϲ���<br>

## app2WithPreinfo
### ��Ҫ����
��app2һЩԤ��Ϣ��һЩ����
### ���������÷���
#### app2WithPreinfo
��һ��ǰ������һ��������app2�ṩһЩ��ֵ�ϵ���Ϣ
```
vaa3d_path /x dll_path /f app2WithPreinfo /i input_dir brain_path /o out_dir /p ratio th resolution_times image_flag lower upper is_multi app2_length contrast_th contrast_ratio
```
vaa3d_path: vaa3dӦ��·��<br>
dll_path: dll·��<br>
input_dir: �����Ŀ¼����Ŀ¼���һ��swc�ļ���һ��apo�ļ���swc�ļ�����һ��������һ��ǰ����apoΪsomaλ��<br>
brain_path: terafly���ݸ�ʽȫ��·��<br>
out_dir: ���Ŀ¼<br>
ratio: ratio������ֵ���㣬th = fmean*ratio + bmean*��1-ratio����fmeanΪǰ����bmeanΪ����<br>
th: ���ratioΪ0����bstd/��bstd+fstd���Զ�����ratio����ratioΪ-4����bmean/��bmean+fmean���Զ�����ratio��ratio��Ϊ-1����app2���������ֵ�ķ�ʽ��ratio��Ϊ-2���Լ�����th��ratioΪ-3����ֵ��Ϊbmean+th<br>
resolution_times�����ĸ��ֱ��ʵ�ͼ���ܣ�1Ϊ��ߣ�2Ϊ�θߣ�4Ϊ...<br>
image_flag: ��ͼ����Ԥ����Ĳ�ͬ����<br>
lower����image_flagΪ2��3ʱ����������ӳ��任����Ϊ3ʱ����һ��(0-bmean)�͵ڶ��Σ�bmean-fmean���������Ա任����Ϊ2ʱ��ֻ�еڶ��������Ա任<br>
upper��lower��Ӧbmean��upper��Ӧfmean<br>
is_multi���Ƿ���multi app2<br>
app2_length��app2��֦�ĳ�����ֵ<br>
contrast_th���Աȶ���ֵ<br>
contrast_ratio���Աȶ�ratio<br>

## Retrace
���Retace�ļ������README

## SNAP
���SNAP�ļ������README

## TypeLength
### ��Ҫ����
ͳ��swc��type�ĳ����Լ�������Ϣ
### ���÷���
�˵�����

## judgeBranch
### ��Ҫ����
ͳ��һЩ��֧������ʹ�����ɭ����ѵ�����ж�swc����֧������
### ���÷���
�˵�����