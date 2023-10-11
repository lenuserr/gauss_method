CFLAGS = -O3  -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

a.out:	main.o input.cpp solution.cpp output.cpp r.cpp
		g++ $^ $(CFLAGS) -o $@

main.o: main.cpp inc.h
		g++ -c $(CFLAGS) main.cpp

clean:
		rm -f *.out *.o
		
