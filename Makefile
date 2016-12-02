default: espresso

espresso: main.cpp
	g++ -std=c++14 -o espresso *.cpp -I.
