
all:
	g++ *.cpp -o build

run: all
	./build

clean:
	rm -rf build
