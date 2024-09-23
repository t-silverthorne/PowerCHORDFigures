BEGIN {
  min_len = ARGV[1];
  max_len = ARGV[2];
  delete ARGV[1];
  delete ARGV[2];
}

function my_filt(line,min_len,max_len){
  if (length(line) > min_len && length(line) < max_len) {
    return 1; 
  } else {
    return 0;
  }
}

{
  if(my_filt($0,min_len,max_len)){
    print $0
  }
}
