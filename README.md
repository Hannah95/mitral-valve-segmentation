# Mitral Valve Segmentation

## run in main.m:

For running code for mitral valve **detection** and **highlighting**: 

`rnmf(method, #videos, varargin)`

where  `method` is the respective method to run, that has to be listed in `rnmf.m`. 
`#videos` is the number of videos the method as to be applied on.
Example: `rnmf('robustNMF',10)`

If desired, the standard variable values can be replaced by others:
Example: `rnmf('robustNMF',10,sparsity,[0,0.1])`


For running code for mitral valve **segmentation**: 

`segment(method, rnmf method, ids, varargin)`

where  `method` is the respective method to run, that has to be listed in `segment.m`. 
`rnmf method` should include the name of the method on whose results (that are maked with IDs) the segmentation should run. `ids` identify the input data for the segmentation.
Example: `rnmf('segmentCVPlus','robustNMF',[1,2,3,7])`

If desired, the standard variable values can be replaced by others:
Example: `rnmf('segmentCVPlus','robustNMF',[1,2,3,7],'rank',5)`




For running code for mitral valve **segmentation in one step**: 

`dual(method, #videos, varargin)`

where  `method` is the respective method to run, that has to be listed in `dual.m`. 
`#videos` is the number of videos the method as to be applied on.
Example: `dual('robustNMF_GL',50)`

If desired, the standard variable values can be replaced by others:
Example: `dual('robustNMF_GL',50,'epsilon',[0.1,0.9,1.5])`


