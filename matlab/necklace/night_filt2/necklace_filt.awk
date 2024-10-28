BEGIN { # read in args and then remove so they are not treated as files
  width   = ARGV[1];
  density = ARGV[2];
  nfine   = ARGV[3];
  delete ARGV[1];
  delete ARGV[2];
  delete ARGV[3];
  count=0;
}

# define custom filter
function contains_spec_window(line_in,width,density, csw_flag, ss , ss_init, line){
  csw_flag=0;
  line=line_in;
  line = line substr(line,1,width-1);
  for (ii=1;ii<=length(line_in);ii++){
    ss=substr(line,ii,width);
    ss_init=substr(line,ii,width);
    gsub(/[^1]/, "",ss); # modifies ss
    if(length(ss)==density){
      csw_flag=1;
    }
  }
  return(csw_flag);
}

# parse file
{
  if(contains_spec_window($0,width,density) && length($0)==nfine){
    #print $0
    count++;
    filename = "output_" int((count-1)/1000000) ".txt"
    print $0 > filename
  }
}
