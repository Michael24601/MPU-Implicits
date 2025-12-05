
all:
	g++ main.cpp -o build/app

run: all
	./build/app

clean:
	rm -rf build
