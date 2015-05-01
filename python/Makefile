

all:
	clang++ -O0 -g -I ../../.. -lpython2.7 -lboost_python -lboost_system -shared -fPIC -I /usr/include/python2.7/ bindings.cpp -o genetics.so

fast:
	clang++ -O3 -I ../../.. -lpython2.7 -lboost_python -lboost_system -shared -fPIC -I /usr/include/python2.7/ bindings.cpp -o genetics.so

