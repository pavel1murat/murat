#!/usr/bin/ruby
#------------------------------------------------------------------------------
# shift X values in X,Y table (two columns, multiple rows) by half-bin in X
# assume constant bin
# output_fn = #{input_fn}+'.new'
#------------------------------------------------------------------------------

def correct(fn)

  fin = File.open(fn);
  text = fin.readlines;
  fin.close

  fout = File.open(fn+'.new','w');

  x0 = text[0].split()[0].to_f
  x1 = text[1].split()[0].to_f

  bin = x1-x0;

  puts "bin = #{bin}"


  for line in text 
    x = line.split()[0].to_f+bin/2;
    y = line.split()[1]
    fout.printf("%8.3f  %s\n",x,y);
  end

  fout.close
end


if __FILE__ == $0 then
  puts "p1 = #{ARGV}"
  fn = ARGV[0];
  correct(fn);
end
