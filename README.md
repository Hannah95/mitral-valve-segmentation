# Mitral Valve Segmentation

## run in main.m:

For running code for mitral valve **detection** and **highlighting**:  <br> <br>

`rnmf(method, #videos, varargin)` <br> 
where  `method` is the respective method to run, that has to be listed in `rnmf.m`. <br>
`#videos` is the number of videos the method as to be applied on.<br>
Example: `rnmf('robustNMF',10)`<br>
If desired, the standard variable values can be replaced by others:<br>
Example: `rnmf('robustNMF',10,sparsity,[0,0.1])`<br>

 <br> <br>
For running code for mitral valve **segmentation**: <br> <br>

`segment(method, rnmf method, ids, varargin)` <br> 
where  `method` is the respective method to run, that has to be listed in `segment.m`. <br>
`rnmf method` should include the name of the method on whose results (that are maked with IDs) the segmentation should run. `ids` identify the input data for the segmentation.<br>
Example: `rnmf('segmentCVPlus','robustNMF',[1,2,3,7])`<br>
If desired, the standard variable values can be replaced by others:<br>
Example: `rnmf('segmentCVPlus','robustNMF',[1,2,3,7],'rank',5)`<br>


 <br> <br>

For running code for mitral valve **segmentation in one step**: <br> <br>

`dual(method, #videos, varargin)`<br>
where  `method` is the respective method to run, that has to be listed in `dual.m`. <br>
`#videos` is the number of videos the method as to be applied on.<br>
Example: `dual('robustNMF_GL',50)`<br>
If desired, the standard variable values can be replaced by others:<br>
Example: `dual('robustNMF_GL',50,'epsilon',[0.1,0.9,1.5])`<br>


