# ptmt: main.o ptmt.o
# 	g++ -std=c++14 main.o ptmt.o -o ptmt

# main.o: main.cpp
# 	g++ -c -std=c++14 main.cpp

# ptmt.o: ptmt.cpp
# 	g++ -c -std=c++14 ptmt.cpp

# clean:
# 	rm *.o ptmt


ifeq ($(OS),Windows_NT)
    IS_WINDOWS := yes
    RM := del *.o ptmt.exe 
else
    IS_WINDOWS := no
    RM := rm -f *.o ptmt 
endif


CXXFLAGS = -std=c++14 -fopenmp

ptmt: main.o partition.o ptmt.o 
	g++ $(CXXFLAGS) main.o partition.o ptmt.o -o PTMT

main.o: main.cpp
	g++ $(CXXFLAGS) -c main.cpp

partition.o: partition.cpp
	g++ $(CXXFLAGS) -c partition.cpp

ptmt.o: ptmt.cpp
	g++ $(CXXFLAGS) -c ptmt.cpp

clean:
	$(RM)