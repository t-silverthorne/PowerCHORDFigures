getPcLinCsts=function(Nmeas,Nfine,drts){
  # construct linear constraint matrix
  zero_vec        = rep(0,Nfine+1)
  lin_csts        = list()
  rhs_list        = list()
  sense_list      = list()
  lin_csts[[1]]   = c(rep(1,Nfine),0)
  rhs_list[[1]]   = Nmeas
  sense_list[[1]] = '=' 
  
  # add support constraints if drts finite 
  if (drts<Inf){
    lin_cstr_vec          = rep(0,Nfine+1)
    rt_inds               = seq(1,Nfine,drts)
    lin_cstr_vec[rt_inds] = 1
    lin_csts[[length(lin_csts)+1]] = lin_cstr_vec  
    rhs_list[[2]]              = sum(lin_cstr_vec)
    sense_list[[2]]            = '>'
  }
  A =lin_csts %>% unlist() %>% matrix(byrow=T,ncol=Nfine+1)
  return(list(A          = A,
              sense_list = sense_list,
              rhs_list   = rhs_list))
}