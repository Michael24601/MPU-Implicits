
all:
	g++ *.cpp -o ./build/build

run: all
	./build

clean:
	rm -rf build
