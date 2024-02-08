CFILE = main.cpp
exec:$(CFILE)      
	g++ -o exec $(CFILE) -lm -lgsl -lgslcblas -Wall
clean:
	rm -rf exec
