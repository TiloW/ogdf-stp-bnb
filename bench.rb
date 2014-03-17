program = ARGV[0]

print "building optimized program .."
%x(make clean)
%x(make #{program} DEBUG=0)
puts "done!"

totalTime = 0.0
count = 0
good = true
File.open('bench.txt').each do |line|
  tmp = line.split /\s+/
  instance = tmp[0]
  optCost = tmp[-1]

  if File.file? instance
    count += 1 
    print "starting program for #{instance} .."
    out = %x(./#{program} #{instance})
    tmp = /after (?<seconds>\d+(\.\d+)?) seconds/.match(out)
    time = tmp ? tmp[:seconds].to_f : 0
    totalTime = totalTime + time
  
    tmp = /solution costs (?<cost>\d+(\.\d+)?)./.match(out)
    cost = tmp ? tmp[:cost] : Infinity
  
    # TODO: Assert result is optimal
    if cost != optCost
      puts "found bad solution: #{cost}"
      puts "Optimal solution is #{optCost}"
      puts "ALGORITHM IS FLAWED."
      puts "Halting execution."
      good = false
      break
    else
      puts "optimal solution found: #{cost}"
    end
  else
    puts "skipping #{instance}"
  end
end

if good
  puts ""
  puts "#{count} instances processed."
  puts "Total computation time was #{totalTime} seconds."
end

%x(make clean)

