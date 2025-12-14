
all:
	g++ -O3 -DNDEBUG *.cpp -o ./build/build.exe

debug:
	g++ -g -DDEBUG_MODE -O0 *.cpp -o ./build/build.exe

run: all
	./build/build

clean:
	rm -rf build
