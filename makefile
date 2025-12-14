
all:
	g++ -O3 -DNDEBUG *.cpp -o ./build/build.exe

run: all
	./build/build

clean:
	rm -rf build
