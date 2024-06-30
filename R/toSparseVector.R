toSparseVector=function(vv){
  # type conversion as.logical() necessary depending on classs of vv
  sparseVector(vv[vv!=0],which(as.logical(vv!=0)),length=length(vv))
}