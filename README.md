# HTStream
A toolset for high throughput sequence analysis using a streaming approach facilitated by Linux pipes.


## build
```
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make all 
```

## debug build
```
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make all
```

## clang build

note if you change the exports on an existing build dir, you must delete the CMakeCache.txt

```
mkdir clang-debug
cd clang-debug
export CXX=clang++
export CC=clang
cmake -DCMAKE_BUILD_TYPE=Debug ..
make all
```

## verbose 
```
make VERBOSE=1

## test

```
make testSD
make test
```
