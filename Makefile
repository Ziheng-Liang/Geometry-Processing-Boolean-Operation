CC = g++
EIGENPATH = ~/lib/
LIBIGLPATH = ~/lib/libigl/include/
BOOSTPATH = ~/lib/


INCLUDE = -I $(EIGENPATH) -I $(LIBIGLPATH) -I $(BOOSTPATH)
CFLAGS = -std=c++11 $(INCLUDE)
LDFLAGS = -lpthread

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)
 
main: main.cpp
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)
	

clean: 
	rm ./main 
	rm *.o
