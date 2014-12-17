#!/bin/csh -f

foreach i (`ls vrp_data/data/Taillard`)
 valgrind --tool=callgrind ./bin/vrp_init -f ./vrp_data/data/Taillard/$i -m 1 >> ./test_sols.out.tmp
 set outf=$i"_SW.out"
 mv callgrind.out.* $outf
end
