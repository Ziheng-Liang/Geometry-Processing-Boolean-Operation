CC = g++
INCLUDE = -I ~/eigen3/ -I ~/libigl/include/
CFLAGS = -std=c++11 $(INCLUDE)
LDFLAGS = -lpthread

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)
 
main: main.cpp
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)
	

clean: 
	rm ./main 
	rm *.o
